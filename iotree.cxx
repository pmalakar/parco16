/*
 *
 * Author: Preeti Malakar
 * Affiliation: Argonne National Laboratory
 * Email: pmalakar@anl.gov
 *
 */

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
//#include <mpio.h>
#include <mpix.h>
#include <queue>
#include <unistd.h>
#include <inttypes.h>

#include <omp.h>

#include <hwi/include/bqc/A2_inlines.h>		//time base 
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/memory.h>

#include "neighbour.h"
#include "personality.h"
#include "iotree.h"

///projects/Performance/preeti/utils
#include "mem.h"
#include "hwcnt.h"

int rootps;
int numBridgeNodes;
int numBridgeNodesAll;

int numMPInodes, size, bncommsize;

int *bridgeNodeAll; 					//[MidplaneSize*2]				//2 integers per rank
bool *visited, *processed;		//[MidplaneSize][MidplaneSize];
int *newBridgeNode;						//[MidplaneSize]
uint8_t **revisit;
uint8_t **depthInfo; 					//[numBridgeNodes][MidplaneSize];

int *bridgeRanks; 						//[numBridgeNodes];
uint8_t bridgeNodeCurrIdx;

int lb, ub;
int collector;

MPI_Comm MPI_COMM_core, MPI_COMM_NODE;
MPI_Comm MPI_COMM_MIDPLANE;
MPI_Comm COMM_BRIDGE_NODES, COMM_BRIDGE_NODES_core;

#ifdef CETUS
int BAG = 256;
#else
int BAG = 64;
#endif

//Average load per BN
float *avgWeight;			//[numBridgeNodes];
int maxWeight;
uint64_t memAvail;

int *shuffledNodes;
double **shuffledNodesData;
double *dataPerNode;

queue <Node *> nodeList;
queue <Node *> *rootNodeList;	//[numBridgeNodes];

Node **bridgeNodeRootList;

int SKIP = 5;
int MAXTIMES = 10;
int MAXTIMESD = 5;

int coalesced;	//0=false, 1=true
int blocking;	//0=nonblocking, 1=blocking
int type;	//0=optimized independent, 1=independent, 2=collective

int count; 

double tstart, tend;
double Tmax, Tmin, Tmax_test2;

/*
 *  Independent MPI-IO
 *
 *  writeFile 
 *  	- writes to IO node or file system
 *
 * 		all - 0 (optimized, only bridge nodes perform the write)
 * 		all - 1 (default, independent write from all node)
 * 		all - 2 (default, collective write from all nodes)
 */

int writeFile(dataBlock *datum, int count, int all) {

	double start=0.0, end=0.0, test1=0.0, test2=0.0;
	int nbytes=0, totalBytes=0, i;

	MPI_Request sendreq, nodesendreq, noderecvreq[ppn-1];
	MPI_Status sendst, nodesendst, noderecvst[ppn-1];

	int result;

	// Optimized independent I/O
	if (all == 0) {

#ifdef DEBUG
	 	if (coreID == 0) 
			printf("%d:%d: called writeFile: %d %d %d %d\n", myrank, coreID, newBridgeNode[myrank], bridgeRanks[newBridgeNode[myrank]], bridgeNodeInfo[0], bridgeNodeInfo[1]);
#endif

		//If I am not a bridge node 
		if (bridgeNodeInfo[1] > 1) {
			//If I have been assigned a new bridge node 
			  int index = (nodeID*ppn) % midplane ; 
		  	if(newBridgeNode[index] != -1) {	
					int	myBridgeRank = bridgeRanks[newBridgeNode[index]] + coreID; 
#ifdef DEBUG
	printf("%d will send to %d\n", myrank, myBridgeRank);		
#endif

					if (coalesced == 1) {

						if (coreID == 0) 
							if (dataPerNode == NULL) printf("allocation error at %d\n", myrank);

						assert(datum->getAlphaBuffer());

//simple MPI_Gather at core 0
/*
						result = MPI_Gather (datum->getAlphaBuffer(), count, MPI_DOUBLE, dataPerNode, count, MPI_DOUBLE, 0, MPI_COMM_NODE);
						if (result != MPI_SUCCESS) 
								prnerror (result, "coalesced nonBN node MPI_Gather error:");
*/

/* replacing gather by p2p */
/*
						if (coreID > 0) {
							//result = MPI_Isend (datum->getAlphaBuffer(), count, MPI_DOUBLE, 0, coreID, MPI_COMM_NODE, &nodesendreq);	
							result = MPI_Send (datum->getAlphaBuffer(), count, MPI_DOUBLE, 0, coreID, MPI_COMM_NODE);	
							if (result != MPI_SUCCESS) 
								prnerror (result, "coalesced nonBN node MPI_Isend error:");
							//MPI_Wait (&nodesendreq, &nodesendst);
						}

						else if (coreID == 0) {
							dataPerNode[0] = *(datum->getAlphaBuffer());
							for (int i=1; i<ppn; i++)
								MPI_Irecv (dataPerNode+count*i, count, MPI_DOUBLE, i, i, MPI_COMM_NODE, &noderecvreq[i-1]);
							MPI_Waitall (ppn-1, noderecvreq, noderecvst);
						}
*/
/* end - replacing gather by p2p */

// trying 2 collectors instead of just core 0

						//Send data to collectors
						if (coreID != 0 && coreID != ppn/2) {
							if (coreID < ppn/2) collector = 0;
							else collector = ppn/2;
							//printf ("%d (%d, %d) sends %d bytes to core %d on node %d\n", myrank, nodeID, coreID, count, collector, nodeID);
							result = MPI_Send (datum->getAlphaBuffer(), count, MPI_DOUBLE, collector, coreID, MPI_COMM_NODE);	
							//result = MPI_Isend (datum->getAlphaBuffer(), count, MPI_DOUBLE, collector, coreID, MPI_COMM_NODE, &nodesendreq);	
							if (result != MPI_SUCCESS) 
								prnerror (result, "coalesced nonBN node MPI_Isend error:");
							//MPI_Wait (&nodesendreq, &nodesendst);
						}

						//Receive data at collectors from senders in the node
						else if (coreID == 0) {
							dataPerNode[0] = *(datum->getAlphaBuffer());
							for (int i=1; i<ppn/2; i++) {
								MPI_Irecv (dataPerNode+count*i, count, MPI_DOUBLE, i, i, MPI_COMM_NODE, &noderecvreq[i-1]);
								//printf ("%d (%d, %d) waiting for %d\n", myrank, nodeID, coreID, i);
							}

							MPI_Status st;	
							result = MPI_Irecv (dataPerNode+count*ppn/2, count*ppn/2, MPI_DOUBLE, ppn/2, ppn/2, MPI_COMM_NODE, &noderecvreq[ppn/2-1]);
							if (result != MPI_SUCCESS) 
								prnerror (result, "coalesced MPI_Recv error:");
							MPI_Waitall (ppn/2, noderecvreq, noderecvst);	// cores 1 - ppn/2-1
							//MPI_Waitall (ppn/2-1, noderecvreq, noderecvst);	// cores 1 - ppn/2-1
						}
						else if (coreID == ppn/2) {
							dataPerNode[0] = *(datum->getAlphaBuffer());
							for (int i=ppn/2+1; i<ppn; i++) {
								int j = i - ppn/2;
								MPI_Irecv (dataPerNode+count*j, count, MPI_DOUBLE, i, i, MPI_COMM_NODE, &noderecvreq[j-1]);
								//printf ("%d (%d, %d) waiting for %d\n", myrank, nodeID, coreID, i);
							}
							MPI_Waitall (ppn/2-1, noderecvreq, noderecvst);	// cores 1 - ppn/2-1
							result = MPI_Send (dataPerNode, count*ppn/2, MPI_DOUBLE, 0, coreID, MPI_COMM_NODE);	
							if (result != MPI_SUCCESS) 
								prnerror (result, "coalesced MPI_Send error:");
						}

/*
						if (coreID == ppn/2) { 
							result = MPI_Send (dataPerNode, count*ppn/2, MPI_DOUBLE, 0, coreID, MPI_COMM_NODE);	
							if (result != MPI_SUCCESS) 
								prnerror (result, "coalesced MPI_Send error:");
						}
						else if (coreID == 0) { 
							MPI_Status st;	
							result = MPI_Recv (dataPerNode+count*ppn/2, count*ppn/2, MPI_DOUBLE, ppn/2, ppn/2, MPI_COMM_NODE, &st);
							if (result != MPI_SUCCESS) 
								prnerror (result, "coalesced MPI_Recv error:");
						}
*/
// end - trying 2 collectors instead of just core 0

						//only core 0 sends
						if (coreID == 0) {
							result = MPI_Isend (dataPerNode, count*ppn, MPI_DOUBLE, myBridgeRank, myBridgeRank, MPI_COMM_WORLD, &sendreq);	
							if (result != MPI_SUCCESS) 
								prnerror (result, "coalesced nonBN node MPI_Isend error:");
							MPI_Wait (&sendreq, &sendst);
						}
					}
					//no coalescing
					else {
						MPI_Isend (datum->getAlphaBuffer(), count, MPI_DOUBLE, myBridgeRank, myBridgeRank, MPI_COMM_WORLD, &sendreq);	
						MPI_Wait (&sendreq, &sendst);
					}
#ifdef DEBUG
			  	printf("%d sent to %d\n", myrank, bridgeRanks[newBridgeNode[myrank]]);
#endif
		  	}

			//If I have not been assigned a new bridge node, write 
		  	else {
#ifdef DEBUG
				if (coreID == 0) printf("%d will write %d doubles\n", myrank, count);
#endif
				if (blocking == 1) { 
					result = MPI_File_write_at (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				}
				else {
					result = MPI_File_iwrite_at (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &sendreq);
					MPI_Wait (&sendreq, &sendst);
				}	
				if (result != MPI_SUCCESS) 
					prnerror (result, "nonBN MPI_File_write_at Error:");
				
				MPI_Get_elements( &status, MPI_CHAR, &nbytes );
				totalBytes += nbytes;
#ifdef DEBUG
			  	printf("%d wrote directly, BN was %d at dist %d\n", myrank, bridgeNodeInfo[0], bridgeNodeInfo[1]);
#endif
		  	}
		}
		//If I am a bridge node, receive data from the new senders
		else if (bridgeNodeInfo[1] == 1) {

				//int arrayLength = myWeight;	//*ppn;
#ifdef DEBUG
			  printf("%d am a BN: myWeight = %d\n", myrank, myWeight);
#endif

				MPI_Request req[myWeight], wrequest[myWeight+1];
				MPI_Status stat, wstatus[myWeight+1];

				//shuffledNodesData = new double *[myWeight];
				assert(shuffledNodesData);

				if (coalesced == 1) {
					if (coreID == 0) { 
#ifdef DEBUG
#endif
				//		for (int i=0; i<myWeight; i++) shuffledNodesData[i] = new double[count*ppn];
				//		printf ("%d: allocated %d * %d bytes\n", myrank, myWeight, count*ppn);
					}
				}
				else {
#ifdef DEBUG
					if (myrank == bridgeRanks[0])
						printf ("%d: about to allocate %d * %d bytes\n", myrank, myWeight, count);
#endif
//					for (int i=0; i<myWeight; i++) shuffledNodesData[i] = new double[count];
//					printf ("uncoalesced %d: allocated %d * %d bytes\n", myrank, myWeight, count);
				}

				// Write out my data first, before waiting for senders' data 
				if (blocking == 1) {
					assert(fileHandle);
					assert(datum->getAlphaBuffer());
					result = MPI_File_write_at (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
					if (result != MPI_SUCCESS) 
						prnerror (result, "BN own MPI_File_write_at Error:");
				}
				else {
					assert(fileHandle);
					assert(datum->getAlphaBuffer());
#ifdef DEBUG
					if (myrank == bridgeRanks[0])
						printf("BN %d will write %d doubles req %d\n", myrank, count, myWeight);
#endif
					result = MPI_File_iwrite_at (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &wrequest[myWeight]);
				}

				int idx;
				// Post the nonblocking receives for my senders
				if (coalesced == 1) {
					if (coreID == 0) {
 
				 	for (int i=0; i<myWeight ; i++) { 
						assert(shuffledNodesData);
						assert(shuffledNodesData[i]);
						assert(shuffledNodes);
#ifdef DEBUG
						printf("\n%d: myWeight = %d shuffledNodes[%d] = %d\n\n", myrank, myWeight, i, shuffledNodes[i]);
#endif
						MPI_Irecv (shuffledNodesData[i], count*ppn, MPI_DOUBLE, shuffledNodes[i], myrank, MPI_COMM_WORLD, &req[i]);				
				 	}

//#pragma omp parallel for
					for (int i=0; i<myWeight ; i++) {

						MPI_Waitany (myWeight, req, &idx, &stat);
#ifdef DEBUG
						printf ("BN %d received data from %d\n", myrank, shuffledNodes[idx]);
#endif
						if (blocking == 1) {
							result = MPI_File_write_at (fileHandle, (MPI_Offset)shuffledNodes[idx]*count*sizeof(double), shuffledNodesData[idx], count*ppn, MPI_DOUBLE, &status);
							if (result != MPI_SUCCESS) 
								prnerror (result, "BN for nonBN MPI_File_write_at Error:");
						}
						else {
							result = MPI_File_iwrite_at (fileHandle, (MPI_Offset)shuffledNodes[idx]*count*sizeof(double), shuffledNodesData[idx], count*ppn, MPI_DOUBLE, &wrequest[idx]);
						}
					}
					if (blocking == 0) 
						MPI_Waitall (myWeight+1, wrequest, wstatus);
					}
				}
				//non-coalesced
				else {
					for (int i=0; i<myWeight ; i++) {
						assert(shuffledNodesData);
						assert(shuffledNodesData[i]);
						assert(shuffledNodes);
						MPI_Irecv (shuffledNodesData[i], count, MPI_DOUBLE, shuffledNodes[i], myrank, MPI_COMM_WORLD, &req[i]); 
					}
//#pragma omp parallel for
					for (int i=0; i<myWeight ; i++) {
						MPI_Waitany (myWeight, req, &idx, &stat);
#ifdef DEBUG
						printf ("BN %d received data from %d\n", myrank, shuffledNodes[idx]);
#endif
						if (blocking == 1) { 
							result = MPI_File_write_at (fileHandle, (MPI_Offset)shuffledNodes[idx]*count*sizeof(double), shuffledNodesData[idx], count, MPI_DOUBLE, &status);
							if (result != MPI_SUCCESS) 
								prnerror (result, "BN for nonblocking MPI_File_write_at Error:");
						}
						else {
							result = MPI_File_iwrite_at (fileHandle, (MPI_Offset)shuffledNodes[idx]*count*sizeof(double), shuffledNodesData[idx], count, MPI_DOUBLE, &wrequest[idx]);
						}
						if (blocking == 1) {
							MPI_Get_elements( &status, MPI_CHAR, &nbytes );
							totalBytes += nbytes;
						}
					}
					if (blocking == 0) 
						MPI_Waitall (myWeight+1, wrequest, wstatus);
				}
			
				//free 
			//	if ((coalesced == 1 && coreID == 0) || coalesced == 0)
				//if (coalesced == 0)
				//if (coalesced == 1 && coreID == 0)
				//	for (int i=0; i<myWeight; i++) free(shuffledNodesData[i]);
				//free(shuffledNodesData);
			}
		}

		// all = 1, default independent write
		else if (all == 1) {
				MPI_Request wrequest;
				start = MPI_Wtime();

				if (blocking == 1) {
					result = MPI_File_write_at (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
					if (result != MPI_SUCCESS) 
						prnerror (result, "all MPI_File_write_at Error:");
				}
				else {
					MPIO_Request req_iw;
					MPI_Status st_iw;
					result = MPI_File_iwrite_at (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &req_iw);
					MPI_Wait(&req_iw, &st_iw);
				}

				end = MPI_Wtime()-start;
				MPI_Get_elements( &status, MPI_CHAR, &nbytes );
				totalBytes += nbytes;
		}

		// all = 2, default collective write
		else if (all == 2) {
			
				start = MPI_Wtime();
				if (blocking == 1) {
					result = MPI_File_write_at_all (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
					if (result != MPI_SUCCESS) 
						prnerror (result, "collective MPI_File_write_at_all Error:");
				}
				else {
#ifdef MPI3
					result = MPIX_File_iwrite_at_all (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &request);
					result = MPI_Wait(&request, &status);
					if (result != MPI_SUCCESS) 
						prnerror (result, "collective MPIX_File_iwrite_at_all Error:");
#else
					MPI_File_set_view (fileHandle, (MPI_Offset)myrank*count*sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
					result = MPI_File_write_all (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
					if (result != MPI_SUCCESS) 
						prnerror (result, "collective MPI_File_write_all Error:");
#endif
				}
				end = MPI_Wtime()-start;
				MPI_Get_elements( &status, MPI_CHAR, &nbytes );
				totalBytes += nbytes;
		}
		return totalBytes;
}


int build5DTorus (int root) {

	double tStart = MPI_Wtime();
	for (int rank = root ; rank < root+midplane ; rank=rank+ppn) {
#ifdef DEBUG
		printf ("%d: findNeighbours of %d\n", myrank, rank);
#endif
		int result = findNeighbours (myrank, rank);
	}
	double tEnd = MPI_Wtime() - tStart;

#ifdef DEBUG
	printf("%d: findNeighbours time taken: %lf\n", myrank, tEnd);	

	for (int i=0; i<midplane ; i++)
		for (int j=0; j<10 ; j++)
			printf ("%d: Nghbr %d of %d -- %d;\n", myrank, j, i, neighbourRanks[i][j]);
#endif

	return 0;

}

//is BN on the default path?
//check if network path from newNode (which is a neighbour of parentNode) to the BN (destination) is a path in the tree of children from BN 
int checkDefaultRoutes (int destination, Node *parentNode, int newNode, int myrank) {

	int coords[6], parentCoords[6], rootCoords[6], parentId, hopDiff=0;
	MPIX_Rank2torus (newNode, coords);
	MPIX_Rank2torus (destination, rootCoords);
	Node *current = parentNode;

	int intmdt_coords[6];
	for (int dim=0; dim < MPIX_TORUS_MAX_DIMS; dim++) {
		intmdt_coords[dim] = coords[dim];
	}   
	intmdt_coords[5] = 0;

	int hopnum = 0, flag = 1;
	int intmdt_rank, parent;

	for (int dim=0; dim<MPIX_TORUS_MAX_DIMS; dim++) {

		int dimID = routingOrder[dim];
		hopDiff = abs(rootCoords[dimID] - coords[dimID]);

		if (hw.isTorus[dimID] == 1 && (hopDiff*2 > hw.Size[dimID])) 
			hopDiff = hw.Size[dimID] - hopDiff ;

		for (int diff=0 ; diff<hopDiff ; diff++) {

			parentId = current->getNodeId();

			int unitHop = 1;
			if (hw.isTorus[dimID] == 0) {
				if(rootCoords[dimID] < coords[dimID]) intmdt_coords[dimID] -= unitHop;  
				else intmdt_coords[dimID] += unitHop;
			}
			else {		// torus
					
					if (abs(rootCoords[dimID] - coords[dimID])*2 > hw.Size[dimID]) {
						//printf("check > %d bridgecoords[%d]=%d hw.Size[%d]=%d\n", myrank, dimID, rootCoords[dimID], dimID, hw.Size[dimID]);

						if (rootCoords[dimID] > coords[dimID]) 
							intmdt_coords[dimID] = ((intmdt_coords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
						else 
							intmdt_coords[dimID] = (intmdt_coords[dimID] + unitHop) % hw.Size[dimID];
					}
					else if (abs(rootCoords[dimID] - coords[dimID])*2 < hw.Size[dimID]) {
						//printf("check < %d bridgecoords[%d]=%d hw.Size[%d]=%d\n", myrank, dimID, rootCoords[dimID], dimID, hw.Size[dimID]);
						if (rootCoords[dimID] < coords[dimID]) 
							intmdt_coords[dimID] = ((intmdt_coords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
						else 
							intmdt_coords[dimID] = (intmdt_coords[dimID] + unitHop) % hw.Size[dimID];
					}
					else {
									//if source coord is even, plus direction
						if (coords[dimID]%2 == 0)	// see phil's email: Aug 22, 2014
							intmdt_coords[dimID] = (intmdt_coords[dimID] + unitHop) % hw.Size[dimID]; 		//even source coord: traverse in plus direction  
						else 
							intmdt_coords[dimID] = ((intmdt_coords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
					}
			}

			++hopnum;

			//get the rank
			MPIX_Torus2rank (intmdt_coords, &parent);
			
			//is this the rank of the current?
			if (parent == parentId) {
#ifdef DEBUG
				printf("%d: match for %d: intmdt_node %d = parent %d, destination= %d\n", myrank, newNode, parentId, parent, destination);
#endif
				current = current->getParent();
				if (current == NULL) {
#ifdef DEBUG
					if(parentId != destination)
						printf("%d: root reached beforehand? error for %d\n", myrank, newNode);
					else 
						printf("%d: Did %d reach %d?\n", myrank, newNode, destination);
#endif
					break;
				}
			}
			else {
				flag = 0;
#ifdef DEBUG
				printf("%d: mismatch for %d: %d != %d destination %d\n", myrank, newNode, parentId, parent, destination);
#endif
				break;	
			} 
		}
	}   

	return flag;

}

void traverse (int index, int level) {

#ifdef DEBUG
	printf ("%d: index=%d level=%d\n", myrank, index, level);
	fflush(stdout);
#endif

	int i, child, nid, rid, nChildren, depth, bn, newDepth;

	Node *node;
	int count=0;
	int lastCount;
	assert(index>=0 && index<numBridgeNodes);
	assert(rootNodeList != NULL);

	lastCount = rootNodeList[index].size();	//get current count of the queue
		
#ifdef DEBUG
	printf ("%d: lastCount %d\n", myrank, lastCount);
	fflush(stdout);
#endif

	while (count < lastCount) {
		node = rootNodeList[index].front();
		rootNodeList[index].pop();
		nChildren = node->getNumChildren();
		nid = node->getNodeId(), rid=node->getRootId();

#ifdef DEBUG
		printf ("%d: nid %d lastCount %d nChildren %d\n", myrank, nid, lastCount, nChildren);
#endif
		for (i=0; i<nChildren; i++) {

			assert(index>=0 && index<numBridgeNodes);
			assert(rootNodeList != NULL);

			rootNodeList[index].push (node->getChild(i));
			child = node->getChildId(i);
			depth = (node->getChild(i))->getDepth();
#ifdef DEBUG
			printf ("%d: child %d depth %d i %d\n", myrank, child, depth, i);
#endif

#ifdef DEBUG
//TODO check seg fault
//			if (myrank == rid) { 
//				printf ("%d: Tally: %d = %d ?\n", myrank, depth, depthInfo[nid][child]);
//				printf("%d level %d : %d [%d] curr %d [%d] orig %d [%d]\n", rid, level, child, i, depth, nid, bridgeNodeAll[child*2+1], bridgeNodeAll[child*2]);
//			}
#endif
			//find the weighted min depth for child if not already processed
			if (!processed[child]) {

				int childIdx = child % midplane;
#ifdef DEBUG
				printf("%d process child %d %d\n", myrank, child, childIdx);
#endif
			 	newBridgeNode[childIdx] = -1;//index;
				newDepth = bridgeNodeAll[child*2+1]; //254; //depth;

				for (bn=0; bn<numBridgeNodes ; bn++) {
				
#ifdef DEBUG
					printf("%d: %d: currentAvg %4.2f currentSum %4.2f avgWeight[%d]=%4.2f\n", myrank, child, currentAvg, currentSum, bn, avgWeight[bn]);
#endif
				 	if (avgWeight[bn] > currentAvg) {
#ifdef DEBUG
						printf("%d: not assigning %d to %d (%d) to balance\n", myrank, child, bridgeRanks[bn], depthInfo[bn][child]);
#endif
						continue;
				 	}
					//TODO check revisit is true or not
#ifdef DEBUG
					printf ("%d: was this node %d marked %d\n", myrank, child, revisit[child][0]);
#endif
					//if (revisit[child][0] == 1)	
				 	if (depthInfo[bn][child] < newDepth && avgWeight[bn] < maxWeight-1) {//TODO fixme : based on mem
						newBridgeNode[childIdx] = bn;
					 	newDepth = depthInfo[bn][child];
#ifdef DEBUG
						printf("%d: May assign %d to %d (%d) current avgWeight[%d]=%4.2f\n", myrank, child, bridgeRanks[bn], newDepth, bn, avgWeight[bn]);
#endif
				 	}
				}

				if (newBridgeNode[childIdx] != -1) {
					//adjust the weights
					++currentSum;
					avgWeight[newBridgeNode[childIdx]] += 1.0;
					currentAvg = currentSum/numBridgeNodes;

					processed[child] = true;

#ifdef DEBUG
					printf("%d: Processed %d cSum %4.2lf cAvg %4.2f wt[ %d ](%d) = %4.2f\n", myrank, child, currentSum, currentAvg, bridgeRanks[newBridgeNode[childIdx]], newBridgeNode[childIdx], avgWeight[newBridgeNode[childIdx]]);
#endif
				}
				else {
						//else assign to someone else
				}
			}
			else {
#ifdef DEBUG
				printf("%d processed already %d\n", myrank, child);
#endif
			}
		}
		count++;
	}

}


bool isParent (Node *child, int nodeId) {

	if (nodeId == (child->getParent())->getNodeId()) {
#ifdef DEBUG
		printf("%d: %d is parent of %d\n", myrank, nodeId, child->getNodeId());
#endif
		return true;
	}
	else 
		return false;

}


//			visit all children, take global decision, fcfs child formation fails to get all nodes

void expandNode (Node *currentNodePtr) {
		
	int currentNode = currentNodePtr->getNodeId();
	int rootid = root->getNodeId();

#ifdef DEBUG	
	printf ("%d: expandNode currentNode=%d rootid=%d\n", myrank, currentNode, rootid);
#endif

	short childNum=-1;
	for (int j=8; j>=0; j--) {		//start with the lowest differing dimension neighbour to "hope 4" default path

		// check all eligible neighbours
		int localNode = neighbourRanks[currentNode][j];

		if (localNode < lb || localNode >= ub) {
#ifdef DEBUG	
			printf("%d: localNode %d %d %d\n", myrank, localNode, lb, ub);
#endif
		 	continue;
		}

		int localNode_ = neighbourRanks[currentNode][j]%midplane ; 

#ifdef DEBUG
		if (myrank == bridgeRanks[bridgeNodeCurrIdx])
			printf ("%d: check for %d neighbour %d of %d\n", myrank, localNode, j, currentNode);
#endif

		if (currentNodePtr != root && isParent(currentNodePtr, localNode) == true) continue;

		//assert(visited);

		if (visited[localNode] == false) {
			//check if network path from localNode (which is a neighbour of currentNodePtr) to the BN (rootid) is a path in the tree of children from BN 
			int success = checkDefaultRoutes(rootid, currentNodePtr, localNode, myrank);
#ifdef DEBUG
			printf ("%d: success %d for check route from %d\n", myrank, success, localNode);
#endif
			if (success) {

#ifdef DEBUG
			printf ("%d: rootid %d currentNode %d localNode %d\n", myrank, rootid, currentNode, localNode);
#endif
				childNum ++;
				Node *childNodePtr = currentNodePtr->addChild(localNode, rootid);
				nodeList.push (childNodePtr);
				int currDepth = childNodePtr->getDepth();
				depthInfo[bridgeNodeCurrIdx][localNode] = currDepth;
#ifdef DEBUG
				if (myrank == bridgeRanks[bridgeNodeCurrIdx]) {
					printf("%d: %d is new child of %d\n", myrank, localNode, currentNode);
					printf("%d DEBUG: %d has been visited as neighbour of %d\n", myrank, localNode, currentNode); 							}
#endif

				//resolve distance metric later
				if(bridgeNodeAll[localNode_*2+1] > currDepth+1) {
					revisit[localNode][0] = 1;// mark - to be visited later
#ifdef DEBUG
					printf("%d: %d: can change the depth for %d from %d to %d\n", myrank, bridgeRanks[bridgeNodeCurrIdx], localNode, bridgeNodeAll[localNode_*2+1], currDepth);
#endif
				}

#ifdef DEBUG
				if (myrank == bridgeRanks[bridgeNodeCurrIdx]) {
					printf("%d: %d %d Tree %d -> %d [label=\"%d\"];\n", \
					 myrank, bridgeNodeAll[localNode_*2+1], bridgeNodeAll[localNode_*2], currentNode, localNode, currDepth);
					if(bridgeNodeAll[localNode_*2+1] > currDepth+1)
					 printf("%d: %d %d ModTree %d -> %d [label=\"%d\"];\n", \ 
					  myrank, bridgeNodeAll[localNode_*2+1], bridgeNodeAll[localNode_*2], currentNode, localNode, currDepth);
				}
#endif
//build map... ??
			}
		}
	}

	if (numNodes >= 0)
		numNodes += childNum;

#ifdef DEBUG
	if (myrank == bridgeRanks[bridgeNodeCurrIdx]) {
		printf("%d: numNodes %d childNum %d BAG %d\n", myrank, numNodes, childNum, BAG);
		if(!nodeList.empty()) printf("%d: empty queue size %d\n", myrank, nodeList.size());
	}
#endif

	if (numNodes > BAG) {
	//	numNodes = 0;
		//no more expansion for this node?
		while(!nodeList.empty())	nodeList.pop();	
#ifdef DEBUG
	if (myrank == bridgeRanks[bridgeNodeCurrIdx]) 
		printf("%d: queue size %d\n", myrank, nodeList.size());
#endif
		return;
	}
	else {
		if(!nodeList.empty()) {
			Node *next = nodeList.front();
			nodeList.pop();
			expandNode (next);
		}
#ifdef DEBUG
		else
			printf("%d: nodeList empty\n", myrank);
#endif
	}

}


void buildTree (int currentNode) {

	static int iter=-1;
	root = new Node (currentNode);
	bridgeNodeRootList[++iter] = root;
#ifdef DEBUG
	getMemStats(myrank, 1);
	printf("%d: expanding root %d\n", myrank, currentNode);
#endif
	expandNode (root);
#ifdef DEBUG
	printf("%d: expanded root %d\n", myrank, currentNode);
#endif

	//delete root;
}

/*
 * formBridgeNodesRoutes(): 
 *
 * There are 8 bridge nodes in BG/Q per midplane.
 *
 * (1) All bridge nodes get to know other bridge nodes in their partition
 *
 * 	- All bridge nodes send their ranks (1 integer) to process 0
 * 	- Process 0 receives all bridge node ranks
 * 	- Process 0 sends bridge node ranks to all bridge nodes [8 integers]
 * 	- All bridge nodes receive ranks of other bridge nodes
 * 	- Isends and Irecvs used in process 0
 *
 * (2) calls buildTree()
 *
 */

void formBridgeNodesRoutes () {

	int i, j, bn, result;
	double tStart, tEnd;

	for (j=0; j<midplane; j++)
			newBridgeNode[j] = -1;

	if (myrank == rootps) {

		MPI_Request request[numBridgeNodes], requestAll[numBridgeNodes];
		MPI_Status status[numBridgeNodes];

#ifdef DEBUG
		printf("I am the rootps %d\n", myrank);
#endif

		tStart = MPI_Wtime();
		for (int iter = 0; iter < numBridgeNodes ; iter ++) {
			result = MPI_Irecv (&bridgeRanks[iter], 1, MPI_INT, MPI_ANY_SOURCE, rootps, MPI_COMM_WORLD, &request[iter]);//to contain messages within midplanes
			if (result != MPI_SUCCESS) 
				prnerror (result, "MPI_Irecv Error: ");
		}

		result = MPI_Waitall(numBridgeNodes, request, status);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Waitall Error: ");
		
	  for (i=0; i<numBridgeNodes ; i++) {
		
#ifdef DEBUG
			printf("%d sends to BN[%d] %d\n", myrank, i, bridgeRanks[i]);
#endif
			//Introduce all bridge nodes to each other
			int tag = rootps + 1;
			result = MPI_Isend (bridgeRanks, numBridgeNodes, MPI_INT, bridgeRanks[i], tag, MPI_COMM_WORLD, &request[i]);
			if (result != MPI_SUCCESS) 
				prnerror (result, "MPI_Isend Error: ");

			//Processes on same node by default have the same bridge node
			tag = rootps + 2;
			result = MPI_Isend (bridgeNodeAll, 2*midplane, MPI_INT, bridgeRanks[i], tag, MPI_COMM_WORLD, &requestAll[i]);
			if (result != MPI_SUCCESS) 
				prnerror (result, "MPI_Isend Error: ");
		}

		result = MPI_Waitall(numBridgeNodes, request, status);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Waitall Error: ");

		result = MPI_Waitall(numBridgeNodes, requestAll, status);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Waitall Error: ");

		tEnd = MPI_Wtime();
#ifdef DEBUG
		printf("%d: send recv overhead in rootps process %6.3f\n", myrank, tEnd-tStart);
#endif

		MPI_Status st;	
		MPI_Recv(newBridgeNode, midplane, MPI_INT, bridgeRanks[0], bridgeRanks[0], MPI_COMM_WORLD, &st);	
	}
	
	//process on Bridge node core 0
	if (bridgeNodeInfo[1] == 1 && coreID == 0) {

#ifdef DEBUG
		printf("%d am the BN on core %d\n", myrank, coreID); 
#endif

		MPI_Request requestSend, requestRecv, requestRecvAll;
		MPI_Status statusSend, statusRecv, statusRecvAll;

		tStart = MPI_Wtime();

		//result = MPI_Isend (&myrank, 1, MPI_INT, rootps, 100, MPI_COMM_WORLD, &requestSend);
		result = MPI_Isend (&myrank, 1, MPI_INT, rootps, rootps, MPI_COMM_WORLD, &requestSend);		//to contain messages within midplanes
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Isend Error: ");

		int tag = rootps + 1;
		result = MPI_Irecv (bridgeRanks, numBridgeNodes, MPI_INT, rootps, tag, MPI_COMM_WORLD, &requestRecv);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Irecv Error: ");

		tag = rootps + 2;
		result = MPI_Irecv (bridgeNodeAll, 2*midplane, MPI_INT, rootps, tag, MPI_COMM_WORLD, &requestRecvAll);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Irecv Error: ");

		result = MPI_Wait(&requestSend, &statusSend);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Waitall Error: ");

		result = MPI_Wait(&requestRecv, &statusRecv);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Waitall Error: ");

		result = MPI_Wait(&requestRecvAll, &statusRecvAll);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Waitall Error: ");

		tEnd = MPI_Wtime();
#ifdef DEBUG
		printf("%d: send recv overhead %6.3f\n", myrank, tEnd-tStart);
#endif

		for (i=0; i<numBridgeNodes ; i++) { 
			if (myrank == bridgeRanks[i]) myBNIdx = i;
#ifdef DEBUG
			printf("%d got it: %d\n", myrank, bridgeRanks[i]);		
#endif
		}

		result = build5DTorus (rootps);

		//	initialize array
		//	all nodes but the bridge nodes are unvisited
		
		for (i=0; i<commsize ; i=i+1) 
			  visited[i] = false, processed[i] = false;

		for (int bn=0; bn<numBridgeNodes ; bn++) { 
		  visited[bridgeRanks[bn]] = true;
		  for (i=0; i<commsize ; i++) {
				revisit[i][0] = 255, revisit[i][1] = 255;
		  	depthInfo[bn][i] = -1;			
#ifdef DEBUG
			  if (myrank == bridgeRanks[bn])
			    printf ("%d: %d Test init depthInfo[%d][%d] = %d\n", myrank, bridgeRanks[bn], bn, i, depthInfo[bn][i]);
#endif
		  }
		}

		// build the topology local to the partition for all bridge nodes
		tStart = MPI_Wtime();
		for (bridgeNodeCurrIdx=0; bridgeNodeCurrIdx<numBridgeNodes ; bridgeNodeCurrIdx++) {

#ifdef DEBUG
			printf("%d: bridge num %d\n", myrank, bridgeNodeCurrIdx);
#endif
			numNodes=0;

#ifdef DEBUG
			printf("%d: building for %d\n", myrank, bridgeRanks[bridgeNodeCurrIdx]);
#endif
			buildTree(bridgeRanks[bridgeNodeCurrIdx]);
		}

		tEnd = MPI_Wtime();
#ifdef DEBUG
		printf("%d: buildTree overhead %6.3f\n", myrank, tEnd-tStart);
#endif

		MPI_Barrier (COMM_BRIDGE_NODES_core);

		while(!(nodeList.empty())) 
			nodeList.pop();

#ifdef DEBUG
		if (!(nodeList.empty())) 
				printf("%d: Strange: i just emptied nodeList\n", myrank);	
		else
				printf("%d:%d nodeList is empty currentAvg=%d\n", myrank, coreID, currentAvg);	
#endif

		//traverse the trees 
		
		currentAvg=currentSum/numBridgeNodes;

		for (bn=0; bn<numBridgeNodes ; bn++)
			avgWeight[bn]=0;

		for (i=0; i<numBridgeNodes; i++) 
			rootNodeList[i].push(bridgeNodeRootList[i]);

		int flag=1, level=0;
		while (flag == 1) {
			level++;
			for (i=0; i<numBridgeNodes; i++) {
				flag=0;
				traverse (i, level);
				if (rootNodeList[i].size() > 0) flag=1;
			}
		}

#ifdef DEBUG
		for (j=0; j<midplane; j++) 
			if (newBridgeNode[j] >= 0 && myrank == bridgeRanks[newBridgeNode[j]])
				printf("%d: %d (%d) is the new BN for %d at distance %d %d\n", myrank, bridgeRanks[newBridgeNode[j]], newBridgeNode[j], j, bridgeNodeAll[j*2+1], depthInfo[newBridgeNode[j]][j]);
#endif

		for (bn=0; bn<numBridgeNodes ; bn++) 
			if (bridgeRanks[bn] == myrank) myWeight = int(avgWeight[bn]);

	  if (myrank == bridgeRanks[0])
			MPI_Send(newBridgeNode, midplane, MPI_INT, rootps, bridgeRanks[0], MPI_COMM_WORLD);	

/*
 *
 * Attempt to resolve only the long distant ones.. does not help much... 
 *
 *//*
		//resolve duplicates
		//tree has lot of duplicates: each node appears avg 6-7 times in the 8 trees
		//approach 1: take those nodes which have nearer bridge nodes than their default bridge nodes 
		//- in such cases, use FCFS to assign these...
		//- make visited = true for these... and then for the rest...
		
		//revisit 3D array: 
		//	dim 0 - 0 or 1 - whether to revisit or not
		//	dim 1 - smallest depth 
		//	dim 2 - bridgenode with the smallest depth 
*/
	}
}


void distributeInfo() {

	int j;

#ifdef DEBUG
		printf ("\n%d:%d:%d: myWeight = %d\n", myrank, coreID, nodeID, myWeight);
#endif

		MPI_Bcast(&myWeight, ppn, MPI_INT, 0, MPI_COMM_NODE);	
		MPI_Bcast(bridgeNodeAll, 2*midplane, MPI_INT, 0, MPI_COMM_NODE);	
#ifdef DEBUG
		printf ("\n%d:%d:%d: real myWeight = %d\n", myrank, coreID, nodeID, myWeight);
#endif

		shuffledNodes = (int *) malloc (myWeight * sizeof(int));
		if (!shuffledNodes) printf("Panic: shuffledNodes not allocated\n");
		for (j=0; j<myWeight ; j++) 
			shuffledNodes[j] = -1;

		int k=-1;
		for (j=0; j<midplane; j=j+ppn) { 
#ifdef DEBUG
			if (newBridgeNode[j] >= 0) 
				printf ("%d: checking: %d %d %d %d %d\n", myrank, nodeID, coreID, j, bridgeRanks[newBridgeNode[j]], bridgeNodeAll[j*2+1]);
#endif

			if (newBridgeNode[j] >= 0 && bridgeRanks[newBridgeNode[j]] == nodeID*ppn && bridgeNodeAll[j*2+1]>1) 
			{
				k++;
				shuffledNodes[k] = lb + j + coreID;
#ifdef DEBUG
				printf("%d:(%d,%d) k=%d newBridgeNode[%d]=%d bn=%d %d\n", myrank, nodeID, coreID, k, j, newBridgeNode[j], bridgeRanks[newBridgeNode[j]], shuffledNodes[k]); 
#endif
			}
		}

		if (k+1 != myWeight) {
		 printf("%d: Error in number of shuffledNodes: %d %d\n", myrank, k, myWeight);
		//abort
		}

		//Allocation	
		//shuffledNodesData = new double *[myWeight];
		shuffledNodesData = (double **) malloc (myWeight * sizeof(double));
		if (shuffledNodesData == NULL)
			printf("\n%d: Error in allocating %ld bytes\n", myrank, myWeight * sizeof (double));

		if (coalesced == 1) {
				if (coreID == 0) { 
					for (int i=0; i<myWeight; i++) {
						//shuffledNodesData[i] = new double[count*ppn];
						shuffledNodesData[i] = (double *) malloc (count * ppn * sizeof(double));
						if (shuffledNodesData[i] == NULL) printf("\n%d: Error in allocating %ld bytes for %d\n", myrank, count * ppn * sizeof (double) ,  i);
					}
#ifdef DEBUG
					printf ("%d: allocated %d * %d bytes\n", myrank, myWeight, count*ppn);
#endif
				}
		}
		else
		{
				for (int i=0; i<myWeight; i++) { 
					//shuffledNodesData[i] = new double[count];
					shuffledNodesData[i] = (double *) malloc (count * sizeof(double));
				}
#ifdef DEBUG
				printf ("uncoalesced %d: allocated %d * %d bytes\n", myrank, myWeight, count);
#endif
		}

		MPI_Barrier (COMM_BRIDGE_NODES);	

}

void initTree(int n) {

	bridgeNodeRootList = (Node **) malloc (n*sizeof(Node *));
#ifdef DEBUG
	if (bridgeNodeRootList == NULL)	printf("\n%d: malloc fail\n", myrank);
#endif

}

int main (int argc, char **argv) {

//		if (argc<5) return 0;
		MPI_Init (&argc, &argv);

		MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
		MPI_Comm_size (MPI_COMM_WORLD, &commsize);

		int exp = atoi (argv[1]);
		fileSize = exp * oneKB;	
		coalesced = atoi(argv[2]);
		blocking = atoi(argv[3]);
		type = atoi(argv[4]);
		count = fileSize;				//weak scaling

		getPersonality(myrank);

		//form inter-communicator - mainly reqd for core 0 processes 
		MPI_Comm_split (MPI_COMM_WORLD, coreID, myrank, &MPI_COMM_core);
		MPI_Comm_size (MPI_COMM_core, &numMPInodes);

		if (numMPInodes <= 1024)
			midplane = MidplaneSize * ppn;
		else
			midplane = commsize / 2;
		
#ifdef CETUS
		if (numMPInodes <= 1024)
			numBridgeNodes = 8;	
		else if (numMPInodes == 2048)
			numBridgeNodes = 16;	
		else if (numMPInodes == 4096)
			numBridgeNodes = 32;		
		else if (numMPInodes == 8192)
			numBridgeNodes = 64;		
#else
		numBridgeNodes = 32; 	//vesta
#endif

#ifdef DEBUG
		if (myrank == 0 || myrank == 1) getMemStats(myrank, 1);
#endif

		bridgeNodeAll = new int [2*midplane];
		newBridgeNode = new int [midplane];

		visited = new bool [commsize];
		processed = new bool [commsize];

		revisit = new uint8_t *[commsize];
		for (int i=0 ; i<commsize ; i++)
			revisit[i] = new uint8_t[2];

		depthInfo = new uint8_t *[numBridgeNodes];
		for (int i=0 ; i<numBridgeNodes ; i++)
			depthInfo[i] = new uint8_t[commsize];

		bridgeRanks = new int [numBridgeNodes];
		avgWeight = new float [numBridgeNodes];

		rootNodeList = new queue<Node *> [numBridgeNodes];

#ifdef DEBUG
		if (myrank == 0 || myrank == 1) getMemStats(myrank, 1);
#endif

		rootps = floor(myrank/(midplane)) * (midplane);
		lb = floor(myrank/midplane) * (midplane);
		ub = lb + midplane; 

		numMidplanes = (commsize/ppn) / midplane;

#ifdef DEBUG
		if (coreID == 0) printf("Logistics: %d:%d:%d: %d %d %d %d\n", myrank, nodeID, coreID, lb, ub, rootps, midplane);
#endif

		double tStart = MPI_Wtime();	//entire execution

		dataBlock *datum = new dataBlock(count);			//initializes alpha array to random double values

		initNeighbours(commsize);
		initTree(numBridgeNodes);

		//form intra-communicator - mainly reqd for processes on a node
		MPI_Comm_split (MPI_COMM_WORLD, nodeID, myrank, &MPI_COMM_NODE);

		//form intra-communicator per midplane 
		MPI_Comm_split (MPI_COMM_WORLD, rootps, myrank, &MPI_COMM_MIDPLANE);

		//form intra-communicator - mainly reqd for bridge nodes
		MPI_Comm_split (MPI_COMM_WORLD, bridgeNodeInfo[1], myrank, &COMM_BRIDGE_NODES);

		//form intra-communicator - mainly reqd for bridge nodes core wise
		MPI_Comm_split (COMM_BRIDGE_NODES, coreID, myrank, &COMM_BRIDGE_NODES_core);

		MPI_Comm_size (COMM_BRIDGE_NODES, &bncommsize);
		//MPI_Comm_size (COMM_BRIDGE_NODES_core, &size);

#ifdef DEBUG
		if (bridgeNodeInfo[1] == 1)
			if (numBridgeNodes * numMidplanes * ppn != bncommsize)
				printf("%d: PANIC: %d * %d != %d\n", myrank, numBridgeNodes, numMidplanes, ppn, bncommsize);
			else
				printf("%d: numMPInodes = %d numBNnodes = %d\n", myrank, numMPInodes, numBridgeNodes);
#endif

		//gather bridgeNodeInfo at the rootps
		//bridgeNodeAll - per root
		double ts = MPI_Wtime();
		MPI_Gather (bridgeNodeInfo, 2, MPI_INT, bridgeNodeAll, 2, MPI_INT, 0, MPI_COMM_MIDPLANE);
		ts = MPI_Wtime() - ts;

#ifdef DEBUG
		if (coreID == 0) printf("%d: mybridgeNodeInfo: %d %d\n", myrank, bridgeNodeInfo[0], bridgeNodeInfo[1]);
		if (myrank == rootps)
			printf("%d: bridgeNodeInfo %d %d %d %d %lf\n", myrank, bridgeNodeAll[2], bridgeNodeAll[3], bridgeNodeAll[6], bridgeNodeAll[7], ts);
		if (myrank == rootps) printf("%d: gather time %lf\n", nodeID, ts);
#endif

#ifdef DEBUG
		if (myrank == 0 || myrank == 1) getMemStats(myrank, 1);
#endif

		double tOStart = MPI_Wtime();

		Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapAvail); 
		//MPI_Reduce(&heapAvail, &memAvail, 1, MPI_UINT64_T, MPI_MIN, 0, MPI_COMM_WORLD);
		//no need of global reduction
		//MPI_Reduce(&heapAvail, &memAvail, 1, MPI_UINT64_T, MPI_MIN, 0, COMM_BRIDGE_NODES);
		//MPI_Bcast(&memAvail, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);	
		MPI_Reduce(&heapAvail, &memAvail, 1, MPI_UINT64_T, MPI_MIN, 0, MPI_COMM_NODE);
		MPI_Bcast(&memAvail, 1, MPI_UINT64_T, 0, MPI_COMM_NODE);	

		if (coalesced == 1)
			maxWeight = memAvail/(2 * count * ppn * sizeof(double));
		else
			maxWeight = 1000;	//high

#ifdef DEBUG
		if (nodeID < 60 && coreID < 2) printf("%d maxWeight = %d heapAvail = %d memAvail = %d\n", myrank, maxWeight, heapAvail, memAvail);
#endif

		if (coreID == 0) formBridgeNodesRoutes ();

		MPI_Bcast(newBridgeNode, midplane, MPI_INT, 0, MPI_COMM_MIDPLANE);	
		MPI_Bcast(bridgeRanks, numBridgeNodes, MPI_INT, 0, MPI_COMM_MIDPLANE);	

		if (bridgeNodeInfo[1] == 1) distributeInfo();

		double tOEnd = MPI_Wtime();

#ifdef DEBUG
		if (coreID == 0) printf("%d: %d MyNewBN %d %d\n", myrank, ppn, newBridgeNode[myrank], bridgeRanks[newBridgeNode[myrank]]);
		printf("%d: overhead %6.3f\n", myrank, tOEnd-tOStart);
#endif

		MPI_Barrier (MPI_COMM_WORLD);

#ifdef STATS
		//bgpminit(0, 0);
#endif

		//Testing BGQ compute nodes to IO nodes performance
		//Write to /dev/null
		//writeFlag = checkION();

		/*
		 * * * * * * * * * * * * Independent MPI-IO to IO nodes from all compute nodes - shared file * * * * * * * * * *
		 */

		// combine data from all ranks in the node
		if (coalesced == 1 && (coreID == 0 || coreID == ppn/2)) {		//FIXME
			//dataPerNode = new double[count*ppn];		//why does this err?
			dataPerNode = (double *) malloc (count * ppn * sizeof(double));
			if (dataPerNode == NULL) printf("%d: allocation error for %ld bytes\n", myrank, count*ppn*sizeof(double));
		}

		double tION[2];

		/* set file open mode */
		mode = MPI_MODE_CREATE | MPI_MODE_RDWR; //WRONLY;

		/* allocate buffer */
		datum->allocElement (1);

		MPI_Barrier (MPI_COMM_WORLD);

		if (type == 1) {

		MPI_File_open (MPI_COMM_WORLD, fileNameION, mode, MPI_INFO_NULL, &fileHandle);
		for (int i=1; i<=SKIP; i++)
			totalBytes[0][1] += writeFile(datum, count, 1);
		tIOStart = MPI_Wtime();
		for (int i=1; i<=MAXTIMES; i++)
			totalBytes[0][1] += writeFile(datum, count, 1);
		tIOEnd = MPI_Wtime();
		tION[1] = (tIOEnd - tIOStart)/MAXTIMES;
		MPI_File_close (&fileHandle);

		}

		else if (type == 0) {
	
		MPI_File_open (MPI_COMM_WORLD, fileNameION, mode, MPI_INFO_NULL, &fileHandle);
		for (int i=1; i<=SKIP; i++)
			totalBytes[0][0] += writeFile(datum, count, 0);
		tIOStart = MPI_Wtime();
		for (int i=1; i<=MAXTIMES; i++)
			totalBytes[0][0] += writeFile(datum, count, 0);
		tIOEnd = MPI_Wtime();
		tION[0] = (tIOEnd - tIOStart)/MAXTIMES;
		MPI_File_close (&fileHandle);

		}

		/*
		 * * * * * * * * * * * * * * * * * * * * * Independent MPI-IO to file system from all compute nodes - shared file * * * * * * * * * * * * * * * * * *
		 */

		double te[3];

		//mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;

		if (type == 1) {

		MPI_File_open (MPI_COMM_WORLD, fileNameFS, mode, MPI_INFO_NULL, &fileHandle);
		for (int i=1; i<=SKIP; i++)
			totalBytes[1][1] += writeFile(datum, count, 1);
		tstart = MPI_Wtime();
		for (int i=1; i<=MAXTIMESD; i++)
			totalBytes[1][1] += writeFile(datum, count, 1);
		tend = (MPI_Wtime() - tstart)/MAXTIMESD;
		MPI_File_close (&fileHandle);

		}

		else if (type == 2) {

		MPI_File_open (MPI_COMM_WORLD, fileNameFSCO, mode, MPI_INFO_NULL, &fileHandle);
		for (int i=1; i<=SKIP; i++)
			totalBytes[1][2] += writeFile(datum, count, 2);
		tstart = MPI_Wtime();
		for (int i=1; i<=MAXTIMESD; i++)
			totalBytes[1][2] += writeFile(datum, count, 2);
		tend = (MPI_Wtime() - tstart)/MAXTIMESD;
		MPI_File_close (&fileHandle);

		}
	
		else if (type == 0) {

		MPI_File_open (MPI_COMM_WORLD, fileNameFSBN, mode, MPI_INFO_NULL, &fileHandle);
		for (int i=1; i<=SKIP; i++)
			totalBytes[1][0] += writeFile(datum, count, 0);
		tstart = MPI_Wtime();
		for (int i=1; i<=MAXTIMESD; i++)
			totalBytes[1][0] += writeFile(datum, count, 0);
		tend = (MPI_Wtime() - tstart)/MAXTIMESD;
		MPI_File_close (&fileHandle);

		}

		/* free buffer */
		datum->freeElement (1);

		//* * * * * * * * * * * * * * * * * * * * * * * * * * * End of IO * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef STATS
		//bgpmfinalize(0, 0);
#endif

		if (coalesced == 1 && (coreID == 0 || coreID == ppn/2)) 		//FIXME
		free(dataPerNode);
	
		if (bridgeNodeInfo[1] == 1) {
		if ((coalesced == 1 && coreID == 0) || coalesced == 0)
			for (int i=0; i<myWeight; i++) free(shuffledNodesData[i]);
		 free(shuffledNodesData);
		}

		double max[5];

//just testing: turn these on
		//MPI_Reduce(&tION[0], &max[0], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		//MPI_Reduce(&tION[1], &max[1], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		//MPI_Reduce(&tend, &max[2], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		MPI_Finalize ();

#ifdef DEBUG
		if (myrank == 0 || myrank == 1) getMemStats(myrank, 1);
#endif

//		if (myrank == 0) {
			//printf ("Times: %d: %d: %d: %d | %d %d | %6.2f | %4.2lf %4.2lf | %4.2lf\n", type, blocking, coalesced, commsize, ppn, omp_get_num_threads(), 8.0*fileSize/1024.0, max[0], max[1], max[2]);
//		}
			printf ("%d: Times: %d: %d: %d: %d | %d %d | %6.2f | %4.2lf %4.2lf | %4.2lf\n", myrank, type, blocking, coalesced, commsize, ppn, omp_get_num_threads(), 8.0*fileSize/1024.0, tION[0], tION[1], tend);

#ifdef STATS
    PrintCounts("NW", hNWSet, myrank);
    PrintCounts("IO", hIOSet, myrank);
#endif
		return 0;

}

