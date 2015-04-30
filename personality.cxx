/*
 *  Created by Preeti Malakar
 *  Argonne National Laboratory
 *
 *  Determine routing order, bridge node information etc.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <mpix.h>

#include "personality.h"
#include "route.h"

/* 
 * Size of midplane on BGQ
 */
//const int MidplaneSize = 512; 

/* 
 * Structure containing BGQ parameters 
 */
MPIX_Hardware_t hw; 

/*
 * Ranks per node, Core ID (0...15) 
 */ 
int ppn, coreID, nodeID;

/*
 * Size of each dimension 
 */ 
int dimSize[MPIX_TORUS_MAX_DIMS];	// torus dimension size

/*
 * Whether dimension is torus or not	//TODO bool 
 */ 
int isTorus[MPIX_TORUS_MAX_DIMS];	// torus wrap = 0/1 

/*
 *  2D array for each rank, contains bridge node info
 *  [0] - Rank of the bridge node, [1] - Distance from the bridge node 
 */ 
int bridgeNodeInfo[2];						

/*
 *  Routing order : varies based on the partition, node count etc. 
 */
int *routingOrder;

/*
 * Initialize system parameters
 * - get ranks per node
 * - get routing order
 */
void initSystemParameters(int myrank) {

  int i;

	MPIX_Hardware(&hw);
	ppn = hw.ppn;
  coreID = hw.coreID;	
  nodeID = myrank/ppn;

	routingOrder = new int[MPIX_TORUS_MAX_DIMS];
	getRoutingOrder(routingOrder);

	for (i=0; i<MPIX_TORUS_MAX_DIMS ; i++) {
		isTorus[i] = hw.isTorus[i];
		dimSize[i] = hw.Size[i];
	}

	//Rank of bridge node and distance to IO node
	bridgeNodeInfo[0] = MPIX_IO_link_id (); 
	bridgeNodeInfo[1] = MPIX_IO_distance (); 
}

void getPersonality (int myrank) {
		
	//local variables
	int i;

	initSystemParameters(myrank);

#ifdef DEBUG
	printf("Torus dimensions = (%u,%u,%u,%u,%u) Routing order = (%d,%d,%d,%d,%d)\n", hw.Size[0], hw.Size[1], hw.Size[2], hw.Size[3], hw.Size[4], routingOrder[0], routingOrder[1], routingOrder[2], routingOrder[3], routingOrder[4]);
	printf("Rank: %d Node: %d Torus coords = (%u,%u,%u,%u,%u) distance to ION: %d link ID: %d\n", myrank, myrank/ppn, hw.Coords[0], hw.Coords[1], hw.Coords[2], hw.Coords[3], hw.Coords[4], bridgeNodeInfo[1], bridgeNodeInfo[0]);

	if (myrank ==0)
		printf("Torus wraps? %u,%u,%u,%u,%u\n", hw.isTorus[0], hw.isTorus[1], hw.isTorus[2], hw.isTorus[3], hw.isTorus[4]);
#endif

		/*
		 * Find coordinates of bridge node
		 */
		int bridgeCoords[6];
		MPIX_Rank2torus (bridgeNodeInfo[0], bridgeCoords);

		/*
		 * Initialize intermediate nodes in original path to the bridge node
     */
		int intmdt_coords[6];
		for (int dim=0; dim < MPIX_TORUS_MAX_DIMS; dim++) {
			intmdt_coords[dim] = hw.Coords[dim];
		}   
		intmdt_coords[5] = 0;

		int hopnum = 0;
		int hopDiff, intmdt_rank, child, parent;
		child = myrank;
		for (int dim=0; dim<MPIX_TORUS_MAX_DIMS; dim++) {

			int dimID = routingOrder[dim];
			hopDiff = abs(bridgeCoords[dimID] - hw.Coords[dimID]);

			if (hw.isTorus[dimID] == 1 && (hopDiff*2 > hw.Size[dimID])) {
		//			printf("hopDiff earlier rank %d size[%d] %d hops %d\n", myrank, dimID, hw.Size[dimID], hopDiff);
					hopDiff = hw.Size[dimID] - hopDiff ;//+ 1;		//torus dimension
		//			printf("hopDiff later rank %d hops %d current hopnum=%d\n", myrank, hopDiff, hopnum);
			}
#ifdef DEBUG
			printf("SD %d to %d Difference in dim %d = %d\n", myrank, bridgeNodeInfo[0], dimID, hopDiff);
#endif

			for(int diff=0; diff<hopDiff ;diff++) {
				int unitHop = 1;
				if (hw.isTorus[dimID] == 0) {
					if(bridgeCoords[dimID] < hw.Coords[dimID]) intmdt_coords[dimID] -= unitHop;  
					else intmdt_coords[dimID] += unitHop;
				}
				else {		// torus
					if (abs(bridgeCoords[dimID] - hw.Coords[dimID])*2 > hw.Size[dimID]) {
						printf("check > %d bridgecoords[%d]=%d hw.Coords[%d]=%d hw.Size[%d]=%d\n", myrank, dimID, bridgeCoords[dimID], dimID, hw.Coords[dimID], dimID, hw.Size[dimID]);

						if (bridgeCoords[dimID] > hw.Coords[dimID]) 
										intmdt_coords[dimID] = ((intmdt_coords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
						else 
										intmdt_coords[dimID] = (intmdt_coords[dimID] + unitHop) % hw.Size[dimID];
					}
					else if (abs(bridgeCoords[dimID] - hw.Coords[dimID])*2 < hw.Size[dimID]) {
						printf("check < %d bridgecoords[%d]=%d hw.Coords[%d]=%d hw.Size[%d]=%d\n", myrank, dimID, bridgeCoords[dimID], dimID, hw.Coords[dimID], dimID, hw.Size[dimID]);
						if (bridgeCoords[dimID] < hw.Coords[dimID]) 
										intmdt_coords[dimID] = ((intmdt_coords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
						else 
										intmdt_coords[dimID] = (intmdt_coords[dimID] + unitHop) % hw.Size[dimID];
					}
					else {
									//if source coord is even, plus direction
						if (hw.Coords[dimID]%2 == 0)	// see phil's email: Aug 22, 2014
										intmdt_coords[dimID] = (intmdt_coords[dimID] + unitHop) % hw.Size[dimID]; 		//even source coord: traverse in plus direction  
						else 
										intmdt_coords[dimID] = ((intmdt_coords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
					}
				}

				++hopnum;

				//get the rank
				MPIX_Torus2rank (intmdt_coords, &intmdt_rank);
				parent = intmdt_rank;
#ifdef DEBUG
				printf ("Enroute %d (%d %d %d %d %d) to %d (%d %d %d %d %d) Hop %d: in dimension %d: to rank %d (%d %d %d %d %d)\n", \
				myrank, hw.Coords[0], hw.Coords[1], hw.Coords[2], hw.Coords[3], hw.Coords[4], \
				bridgeNodeInfo[0], bridgeCoords[0], bridgeCoords[1], bridgeCoords[2], bridgeCoords[3], bridgeCoords[4], \
				hopnum, dimID, intmdt_rank, intmdt_coords[0], intmdt_coords[1], intmdt_coords[2], intmdt_coords[3], intmdt_coords[4]);

				printf ("Route %d to %d Hop %d\n", myrank, intmdt_rank, hopnum);
				printf (" %d->%d;\n", child, parent);
#endif
				child = parent;
			}
		}   

#ifdef DEBUG
		//Check if everyone was routed to their bridge nodes?
		if (parent != bridgeNodeInfo[0] && bridgeNodeInfo[1] > 1)  printf("Rank %d lost, did not reach %d, instead captivated by %d :D\n", myrank, bridgeNodeInfo[0], parent);
		if (hopnum != bridgeNodeInfo[1]-1)  printf("Rank %d differs in number hops %d!=%d !!!!\n", myrank, bridgeNodeInfo[1], hopnum);
		printf ("%d reaches %d in %d hops\n", myrank, bridgeNodeInfo[0], hopnum);
#endif

		return;

}

