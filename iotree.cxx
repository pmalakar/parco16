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
#include <math.h>
#include <string.h>
#include <mpi.h>
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

int rootps = 0;					// Default root process

int *bridgeNodeAll; 	//[MidplaneSize*2];				//2 integers per rank
bool *visited, *processed;		//[MidplaneSize], processed[MidplaneSize];
int *newBridgeNode;	//[MidplaneSize];
uint8_t **revisit;
uint8_t **depthInfo; //[numBridgeNodes][MidplaneSize];

int *bridgeRanks; //	[numBridgeNodes];
uint8_t bridgeNodeCurrIdx;

MPI_Comm COMM_BRIDGE_NODES, MPI_COMM_core, MPI_COMM_node;
//#ifdef CETUS
const int numBridgeNodes = 8;	//cetus/mira 
//#else
//const int numBridgeNodes = 32; 	//vesta
//#endif

float *avgWeight;			//[numBridgeNodes];

queue <Node *> nodeList;
queue <Node *> *rootNodeList;	//[numBridgeNodes];

Node **bridgeNodeRootList;

int BLOCKING;

double Tmax, Tmin, Tmax_test2;

/*
 *  Independent MPI-IO
 *
 *  writeFile 
 *  	- writes to IO node or file system
 *
 * 		all - 0 (optimized, only bridge nodes perform the write)
 * 		all - 1 (default, all the nodes perform the write)
 *
 */

int writeFile(dataBlock *datum, int count, int all) {

	double start=0.0, end=0.0, test1=0.0, test2=0.0;
	int nbytes=0, totalBytes=0, i;

	MPI_Request sendreq;
	MPI_Status sendst;

	int result;

	// Optimized independent I/O
	if (all == 0) {

		tIOEnd = 0;
		tIOStart = MPI_Wtime();

#ifdef DEBUG
	 	if (coreID == 0) printf("%d:%d: called writeFile: %d %d %d %d\n", myrank, coreID, newBridgeNode[myrank], bridgeRanks[newBridgeNode[myrank]], bridgeNodeInfo[0], bridgeNodeInfo[1]);
#endif

		//If I am not a bridge node 
		if (bridgeNodeInfo[1] > 1) {

			//  printf("%d am not a bridge node, old BN %d, new BN %d\n", myrank, bridgeNodeInfo[1], newBridgeNode[myrank]);
			//If I have been assigned a new bridge node 
		  	if(newBridgeNode[nodeID*ppn] != -1) {	
					int	myBridgeRank = bridgeRanks[newBridgeNode[nodeID*ppn]] + coreID; 
#ifdef DEBUG
			  	printf("%d will send to %d\n", myrank, myBridgeRank);		//bridgeRanks[newBridgeNode[myrank]]);
#endif
					//MPI_Isend (datum->getAlphaBuffer(), count, MPI_DOUBLE, bridgeRanks[newBridgeNode[myrank]], bridgeRanks[newBridgeNode[myrank]], MPI_COMM_WORLD, &sendreq);	
					MPI_Isend (datum->getAlphaBuffer(), count, MPI_DOUBLE, myBridgeRank, myBridgeRank, MPI_COMM_WORLD, &sendreq);	
	//		  	printf("%d sent to %d\n", myrank, bridgeRanks[newBridgeNode[myrank]]);
					MPI_Wait (&sendreq, &sendst);
		  	}

			//If I have not been assigned a new bridge node, write 
		  	else {
#ifdef DEBUG
				printf("%d will write\n", myrank);
#endif
				result = MPI_File_write_at (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				if ( result != MPI_SUCCESS) {
					prnerror (result, "nonBN MPI_File_write_at Error:");
				}
				MPI_Get_elements( &status, MPI_CHAR, &nbytes );
				totalBytes += nbytes;
//			  	printf("%d wrote\n", myrank);
		  	}

		}

		//If I am a bridge node, receive data from the new senders
		else if (bridgeNodeInfo[1] == 1) {

				int arrayLength = myWeight;	//*ppn;
#ifdef DEBUG
			  printf("%d am a BN: myWeight = %d\n", myrank, arrayLength);
#endif

				MPI_Request req[arrayLength], wrequest[myWeight+1];
				MPI_Status stat, wstatus[myWeight+1];

				double shuffledNodesData[arrayLength][count];

				// Write out my data first, before waiting for senders' data 
				result = MPI_File_write_at (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				if ( result != MPI_SUCCESS) {
					prnerror (result, "BN MPI_File_write_at Error:");
				}
				MPI_Get_elements( &status, MPI_CHAR, &nbytes );
				totalBytes += nbytes;

	//TODO file seek //file set view??	
				//result = MPI_File_write_shared (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);

				//MPI_Offset disp = myrank *count*sizeof(double); 
				//MPI_File_seek (fileHandle, disp, MPI_SEEK_SET);
				//result = MPI_File_write(fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				//if (NONBLOCKING == 1)
				//	result = MPI_File_iwrite(fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &wrequest[myWeight]);


				//result = MPI_File_iwrite_at (fileHandle, myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &wrequest[myWeight]);

#pragma omp parallel for
				// Post the nonblocking receives for my senders
				for (i=0; i<arrayLength ; i++) {
					MPI_Irecv (shuffledNodesData[i], count, MPI_DOUBLE, shuffledNodes[i], myrank, MPI_COMM_WORLD, &req[i]); 
				}

        test1 = MPI_Wtime();
#pragma omp parallel for
				for (i=0; i<arrayLength ; i++) {
					int idx;
					MPI_Waitany (myWeight, req, &idx, &stat);
					//result = MPI_File_iwrite_at (fileHandle, shuffledNodes[idx]*count*sizeof(double), shuffledNodesData[idx], count, MPI_DOUBLE, &wrequest[idx]);
					//result = MPI_File_write_shared (fileHandle, shuffledNodesData[idx], count, MPI_DOUBLE, &status);
					
//					disp = shuffledNodes[idx]*count*sizeof(double); 
//					MPI_File_seek (fileHandle, disp, MPI_SEEK_SET);
//					result = MPI_File_write (fileHandle, shuffledNodesData[idx], count, MPI_DOUBLE, &status);
					//if (NONBLOCKING == 1)
					//	result = MPI_File_iwrite (fileHandle, shuffledNodesData[idx], count, MPI_DOUBLE, &wrequest[i]);
					
					result = MPI_File_write_at (fileHandle, (MPI_Offset)shuffledNodes[idx]*count*sizeof(double), shuffledNodesData[idx], count, MPI_DOUBLE, &status);
					if ( result != MPI_SUCCESS) {
						prnerror (result, "BN for nonBN MPI_File_write_at Error:");
					}
					MPI_Get_elements( &status, MPI_CHAR, &nbytes );
					totalBytes += nbytes;
				}

        test2 += MPI_Wtime()-test1;

				if (NONBLOCKING == 1)
					MPI_Waitall (myWeight+1, wrequest, wstatus);

			/*
				int num = 9;
				double *shuffledData = (double *) malloc (num * count * sizeof (double));
				for (i=0; i<num ; i++)
					for (int j=0; j<count ; j++)
						shuffledData[i+j] = shuffledNodesData[i][j];

				start = MPI_Wtime();
				result = MPI_File_write_at (fileHandle, myBNIdx*count*sizeof(double), shuffledData, num*count, MPI_DOUBLE, &status);
				//result = MPI_File_write_all (fileHandle, shuffledData, num*count, MPI_DOUBLE, &status);
			*/	
				/*for (i=0; i<myWeight ; i++) {
					result = MPI_File_write_at (fileHandle, shuffledNodes[i]*count*sizeof(double), shuffledNodesData[i], count, MPI_DOUBLE, &status);
				}*/

				//result = MPI_File_write_at (fileHandle, myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				//result = MPI_File_write_all (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				//double *datumBuffer = new double[1];
				//datumBuffer[0] = 2.1;
				//result = MPI_File_write_all (fileHandle, datumBuffer, 1, MPI_DOUBLE, &status);

				tIOEnd = MPI_Wtime()-tIOStart;
//			  printf("%d %d B %d BN write time %lf test2=%lf\n", myrank, fileSize, level, tIOEnd, test2);
			}

			//	MPI_Reduce(&tIOEnd, &Tmax, 1, MPI_DOUBLE, MPI_MAX, rootps, MPI_COMM_WORLD);
			//	MPI_Reduce(&test2, &Tmax_test2, 1, MPI_DOUBLE, MPI_MAX, rootps, MPI_COMM_WORLD);
			//	MPI_Reduce(&tIOEnd, &Tmin, 1, MPI_DOUBLE, MPI_MIN, rootps, MPI_COMM_WORLD);
			//  if (myrank == rootps) printf("%d %d B %d BN write time %lf min %lf test2=%lf\n", myrank, fileSize, level, Tmax, Tmin, Tmax_test2);
		}
		// all = 1, default independent write
		else {
			 
				MPI_Request wrequest;
				start = MPI_Wtime();
				result = MPI_File_write_at (fileHandle, (MPI_Offset)myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				if ( result != MPI_SUCCESS) {
					prnerror (result, "all MPI_File_write_at Error:");
				}

				MPI_Get_elements( &status, MPI_CHAR, &nbytes );
				totalBytes += nbytes;
				
				//result = MPI_File_write_shared (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				//result = MPI_File_write_all (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				//result = MPI_File_iwrite_at (fileHandle, myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &wrequest);

	//			MPI_Offset disp = myrank *count*sizeof(double); 
	//			MPI_File_seek (fileHandle, disp, MPI_SEEK_SET); 
	//			result = MPI_File_write(fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);

//				result = MPI_File_iwrite(fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &wrequest);
//				MPI_Wait(&wrequest, &status);

				end = MPI_Wtime()-start;
	//			printf("%d %d %d all write time %lf\n", myrank, fileSize, level, end);

		//		MPI_Reduce(&end, &Tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		//		MPI_Reduce(&end, &Tmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		//		if (myrank == 0) printf("%d %d B %d all write time %lf min %lf\n", myrank, fileSize, level, Tmax, Tmin);
		}

		/* Worse performance
		MPI_Offset disp = myrank*count*sizeof(double); 
		int result = MPI_File_set_view (fileHandle, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
		result = MPI_File_write_all (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
		*/

		return totalBytes;

}


int build5DTorus (int root) {

	double tStart = MPI_Wtime();
	//for (int rank = root ; rank < root+MidplaneSize ; rank++) {
	for (int rank = root ; rank < root+MidplaneSize*ppn ; rank=rank+ppn) {
#ifdef DEBUG
		printf ("%d: findNeighbours of %d\n", myrank, rank);
#endif
		int result = findNeighbours (myrank, rank);
		//if (result != 0) return result;
	}
	double tEnd = MPI_Wtime() - tStart;
#ifdef DEBUG
	printf("%d: findNeighbours time taken: %lf\n", myrank, tEnd);	
#endif

/*
#ifdef DEBUG
	for (int i=0; i<MidplaneSize ; i++)
		for (int j=0; j<10 ; j++)
			printf ("%d: Nghbr %d of %d -- %d;\n", myrank, j, i, neighbourRanks[i][j]);
#endif
*/
	return 0;

}


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

		for (int diff=0; diff<hopDiff ;diff++) {

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
				//printf("%d: match for %d: %d equals %d destination %d\n", myrank, newNode, parentId, parent, destination);
				current = current->getParent();
#ifdef DEBUG
				if (current == NULL)
					printf("%d: root reached beforehand? error for %d\n", myrank, newNode);
#endif
			}
			else {
				flag = 0;
				//printf("%d: mismatch for %d: %d != %d destination %d\n", myrank, newNode, parentId, parent, destination);
				break;	
			} 
				
		}
	}   

	return flag;

}

void traverse (int index, int level) {

	int i, child, nid, rid, nChildren, depth, bn, newDepth;

	Node *node;
	int count=0, lastCount = rootNodeList[index].size();	//get current count of the queue
#ifdef DEBUG
	printf ("%d: lastCount %d\n", myrank, lastCount);
#endif

	while (count < lastCount) {
		node = rootNodeList[index].front();
		rootNodeList[index].pop();
		nChildren = node->getNumChildren();
		nid = node->getNodeId(), rid=node->getRootId();

		for (i=0; i<nChildren; i++) {

			rootNodeList[index].push (node->getChild(i));
			child = node->getChildId(i);
			depth = (node->getChild(i))->getDepth();
			if (myrank == rid) { 
#ifdef DEBUG
				printf ("%d: Tally: %d = %d ?\n", depth, depthInfo[nid][child]);
				printf("%d level %d : %d [%d] curr %d [%d] orig %d [%d]\n", rid, level, child, i, depth, nid, bridgeNodeAll[child*2+1], bridgeNodeAll[child*2]);
#endif
			}
			//find the weighted min depth for child if not already processed
			if (!processed[child]) {

			 	newBridgeNode[child] = -1;//index;
				newDepth = 254; //depth;

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
				 	if (depthInfo[bn][child] < newDepth) {// && avgWeight[bn] <= currentAvg) {
						newBridgeNode[child] = bn;
					 	newDepth = depthInfo[bn][child];
#ifdef DEBUG
						printf("%d: May assign %d to %d (%d) current avgWeight[%d]=%4.2f\n", myrank, child, bridgeRanks[bn], newDepth, bn, avgWeight[bn]);
#endif
				 	}
				}

				if (newBridgeNode[child] != -1) {
					//adjust the weights
					++currentSum;
					avgWeight[newBridgeNode[child]] += 1.0;
					currentAvg = currentSum/numBridgeNodes;

					processed[child] = true;
#ifdef DEBUG
					printf("%d: Processed %d cSum %4.2lf cAvg %4.2f wt[ %d ](%d) = %4.2f\n", myrank, child, currentSum, currentAvg, bridgeRanks[newBridgeNode[child]], newBridgeNode[child], avgWeight[newBridgeNode[child]]);
#endif
				}
				else {
						//else assign to someone else
				}
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
//			visited[localNode] = true;		//trying to find how many satisfy default routes

void expandNode (Node *currentNodePtr) {
		
	int currentNode = currentNodePtr->getNodeId();
	int rootid = root->getNodeId();
	short childNum=-1;
	for (int j=8; j>=0; j--) {		//start with the lowest differing dimension neighbour to "hope 4" default path

		// check all eligible neighbours
		int localNode = neighbourRanks[currentNode][j] ;
		if (myrank == bridgeRanks[bridgeNodeCurrIdx])
#ifdef DEBUG
			printf ("%d: check for %d neighbour %d of %d\n", myrank, localNode, j, currentNode);
#endif

		//TODO continue if localNode is parent of currentNode... should not become cycle...
		if (currentNodePtr != root && isParent(currentNodePtr, localNode) == true) continue;

		if (visited[localNode] == false) {
			int success = checkDefaultRoutes(rootid, currentNodePtr, localNode, myrank);
			//printf ("%d: success %d for check route from %d\n", myrank, success, localNode);

			if (success) {

				childNum ++;
				Node *childNodePtr = currentNodePtr->addChild(localNode, rootid);
				nodeList.push (childNodePtr);
				int currDepth = childNodePtr->getDepth();
				depthInfo[bridgeNodeCurrIdx][localNode] = currDepth;
#ifdef DEBUG
				if (myrank == bridgeRanks[bridgeNodeCurrIdx]) {
					printf("%d: %d is new child of %d\n", myrank, localNode, currentNode);
					printf("%d DEBUG: %d has been visited as neighbour of %d\n", myrank, localNode, currentNode); 
				}
#endif



				


				


				//resolve distance metric later
				//if(depthInfo[bridgeNodeCurrIdx][localNode] > currDepth+1) {
				if(bridgeNodeAll[localNode*2+1] > currDepth+1) {
					revisit[localNode][0] = 1;// mark - to be visited later
#ifdef DEBUG
					printf("%d: %d: can change the depth for %d from %d to %d\n", myrank, bridgeRanks[bridgeNodeCurrIdx], localNode, bridgeNodeAll[localNode*2+1], currDepth);
#endif
					//printf("%d: %d: can change the depth for %d from %d to %d\n", myrank, bridgeRanks[bridgeNodeCurrIdx], localNode, depthInfo[bridgeNodeCurrIdx][localNode], currDepth);
				}




#ifdef DEBUG
				if (myrank == bridgeRanks[bridgeNodeCurrIdx]) {
					printf("%d: %d %d Tree %d -> %d [label=\"%d\"];\n", \
					 myrank, bridgeNodeAll[localNode*2+1], bridgeNodeAll[localNode*2], currentNode, localNode, currDepth);
					if(bridgeNodeAll[localNode*2+1] > currDepth+1)
					 printf("%d: %d %d ModTree %d -> %d [label=\"%d\"];\n", \ 
					  myrank, bridgeNodeAll[localNode*2+1], bridgeNodeAll[localNode*2], currentNode, localNode, currDepth);
				}
#endif

//build map... ??

			}
		}

		//if (childNum == 3 || numNodes>63)		break;
		//if (numNodes>63)		break;		//TEST: discover all
	}

	numNodes += childNum;
#ifdef DEBUG
	printf("%d: numNodes %d childNum %d\n", myrank, numNodes, childNum);
	if(!nodeList.empty()) printf("%d: queue size %d\n", myrank, nodeList.size());
#endif

	//if (childNum != 3) printf ("%d: ERROR: childNum not 3 but %d\n", myrank, childNum);
int BAG = 64;
#ifdef CETUS
	BAG = 256;
#endif

	if (numNodes > BAG) {
	//	numNodes = 0;
		//no more expansion for this node?
		while(!nodeList.empty())	nodeList.pop();	
#ifdef DEBUG
		printf("%d: queue size %d\n", myrank, nodeList.size());
#endif
		return;
	}
	else {
		if(!nodeList.empty()) {
			Node *next = nodeList.front();
			nodeList.pop();
			//printf("%d: %d expanding child %d\n", myrank, root->getNodeId(), next->getNodeId());
			expandNode (next);
			//printf("%d: %d expanded child %d\n", myrank, root->getNodeId(), next->getNodeId());
		}
#ifdef DEBUG
		else
			printf("%d: nodeList empty\n");
#endif
	}

}


void buildTree (int currentNode) {

	static int iter=-1;
	root = new Node (currentNode);
	bridgeNodeRootList[++iter] = root;
#ifdef DEBUG
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
 * (1) All bridge nodes get to know other bridge nodes in their MidplaneSize partition
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
//	rootps = floor(myrank/(MidplaneSize*ppn));

//	int numBridgeNodes = 8; 		// TODO Number of bridge nodes?
//	int bridgeRanks[numBridgeNodes];

	//for (j=0; j<MidplaneSize ; j++)
	for (j=0; j<MidplaneSize*ppn ; j++)
			newBridgeNode[j] = -1;

	if (myrank == rootps) {

		MPI_Request request[numBridgeNodes], requestAll[numBridgeNodes];
		MPI_Status status[numBridgeNodes];

#ifdef DEBUG
		printf("I am the rootps %d\n", myrank);
#endif
		tStart = MPI_Wtime();
		for (int iter = 0; iter < numBridgeNodes ; iter ++) {
			result = MPI_Irecv (&bridgeRanks[iter], 1, MPI_INT, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, &request[iter]);
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
			result = MPI_Isend (bridgeRanks, numBridgeNodes, MPI_INT, bridgeRanks[i], 101, MPI_COMM_WORLD, &request[i]);
			if (result != MPI_SUCCESS) 
				prnerror (result, "MPI_Isend Error: ");

			//Send all nodes' bridgeNodeInfo to all bridge nodes//TODO hardcoded 512
			//Processes on same node by default have the same bridge node
			//result = MPI_Isend (bridgeNodeAll, 2*MidplaneSize, MPI_INT, bridgeRanks[i], 102, MPI_COMM_WORLD, &requestAll[i]);
			result = MPI_Isend (bridgeNodeAll, 2*MidplaneSize*ppn, MPI_INT, bridgeRanks[i], 102, MPI_COMM_WORLD, &requestAll[i]);
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
		//MPI_Recv(newBridgeNode, MidplaneSize, MPI_INT, bridgeRanks[0], bridgeRanks[0], MPI_COMM_WORLD, &st);	
		MPI_Recv(newBridgeNode, MidplaneSize*ppn, MPI_INT, bridgeRanks[0], bridgeRanks[0], MPI_COMM_WORLD, &st);	
	}
	
	//process on Bridge node core 0
	if (bridgeNodeInfo[1] == 1 && coreID == 0) {

		MPI_Request requestSend, requestRecv, requestRecvAll;
		MPI_Status statusSend, statusRecv, statusRecvAll;

		tStart = MPI_Wtime();

		result = MPI_Isend (&myrank, 1, MPI_INT, rootps, 100, MPI_COMM_WORLD, &requestSend);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Isend Error: ");

		result = MPI_Irecv (bridgeRanks, numBridgeNodes, MPI_INT, rootps, 101, MPI_COMM_WORLD, &requestRecv);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Irecv Error: ");

		//result = MPI_Irecv (bridgeNodeAll, 2*MidplaneSize, MPI_INT, rootps, 102, MPI_COMM_WORLD, &requestRecvAll);
		result = MPI_Irecv (bridgeNodeAll, 2*MidplaneSize*ppn, MPI_INT, rootps, 102, MPI_COMM_WORLD, &requestRecvAll);
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
			//printf("%d got it: %d\n", myrank, bridgeRanks[i]);		
		}

		result = build5DTorus (rootps);

		//	initialize MidplaneSize array
		//	all nodes but the bridge nodes are unvisited
		
		//for (i=0; i<MidplaneSize ; i++) 
		for (i=0; i<MidplaneSize*ppn ; i=i+1) 
			  visited[i] = false, processed[i] = false;

		for (int bn=0; bn<numBridgeNodes ; bn++) { 
		  visited[bridgeRanks[bn]] = true;
		  //for (i=0; i<MidplaneSize ; i++) {
		  for (i=0; i<MidplaneSize*ppn ; i=i+1) {
				revisit[i][0] = 255, revisit[i][1] = 255;
		  	depthInfo[bn][i] = -1;			
#ifdef DEBUG
			  if (myrank == bridgeRanks[bn])
			    printf ("%d: %d Test init depthInfo[%d][%d] = %d\n", myrank, bridgeRanks[bn], bn, i, depthInfo[bn][i]);
#endif
		  }
		}

		// build the topology local to the MidplaneSize-node partition for all bridge nodes
		tStart = MPI_Wtime();
		for (bridgeNodeCurrIdx=0; bridgeNodeCurrIdx<numBridgeNodes ; bridgeNodeCurrIdx++) {

#ifdef DEBUG
			printf("%d: bridge num %d\n", myrank, bridgeNodeCurrIdx);
#endif
			numNodes=-1;

#ifdef DEBUG
			printf("%d: building for %d\n", myrank, bridgeRanks[bridgeNodeCurrIdx]);
#endif
			buildTree(bridgeRanks[bridgeNodeCurrIdx]);
		}

		tEnd = MPI_Wtime();
#ifdef DEBUG
		printf("%d: buildTree overhead %6.3f\n", myrank, tEnd-tStart);
#endif

		//MPI_Barrier (COMM_BRIDGE_NODES);

		while(!nodeList.empty()) {
			nodeList.pop();
		}

#ifdef DEBUG
		if (!nodeList.empty()) 
				printf("%d: Strange: i just emptied nodeList\n", myrank);	
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
		//for (j=0; j<MidplaneSize ; j++) 
		for (j=0; j<MidplaneSize*ppn ; j++) 
			if (newBridgeNode[j] >= 0)
				printf("%d: %d (%d) is the new BN for %d\n", myrank, bridgeRanks[newBridgeNode[j]], newBridgeNode[j], j);
#endif

		for (bn=0; bn<numBridgeNodes ; bn++) 
			if (bridgeRanks[bn] == myrank) myWeight = int(avgWeight[bn]);

	  if (myrank == bridgeRanks[0])
			MPI_Send(newBridgeNode, MidplaneSize*ppn, MPI_INT, rootps, bridgeRanks[0], MPI_COMM_WORLD);	//TODO subcomm for the midplane
			//MPI_Send(newBridgeNode, MidplaneSize, MPI_INT, rootps, bridgeRanks[0], MPI_COMM_WORLD);	//TODO subcomm for the midplane
	

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
	
		currentAvg=currentSum/numBridgeNodes;

		for (bn=0; bn<numBridgeNodes ; bn++)
			avgWeight[bn]=0;

		for (j=0; j<MidplaneSize ; j++)
			newBridgeNode[j] = -1;

		int winner=-1;
		for (j=0; j<MidplaneSize ; j++) {
			winner = bridgeNodeAll[2*j];			//default BN
		  for (bn=0; bn<numBridgeNodes ; bn++) {
		  	printf("%d: revisit %d %d %d %d %d\n", myrank, j, revisit[j][1], depthInfo[bn][j], revisit[j][0], bridgeRanks[bn]);
		  	if ((revisit[j][1] > 1+depthInfo[bn][j]) && revisit[j][0] != 255 && depthInfo[bn][j] != -1 && bridgeNodeAll[2*j] != bridgeRanks[bn]) { //TODO check if this is child of 9? i.e. isChild() == true 

					if (avgWeight[bn] > currentAvg) {
						printf("%d: not assigning %d (%d) for %d to balance\n", myrank, bridgeRanks[bn], depthInfo[bn][j], j);
						continue;
					}
					printf("%d: BN may be changing for %d: old: %d (%d) new: %d (%d)\n", myrank, j, revisit[j][1], newBridgeNode[j], bridgeRanks[bn], depthInfo[bn][j]);
					revisit[j][1] = depthInfo[bn][j]; newBridgeNode[j] = bridgeRanks[bn]; winner = bn;
					printf("%d: BN may change for %d: originally assigned: %d (%d) new: %d (%d)\n", myrank, j, bridgeNodeAll[2*j], bridgeNodeAll[2*j+1], newBridgeNode[j], depthInfo[bn][j]);
		 		}
		  }

			if (newBridgeNode[j] != -1) {
				++currentSum;
				currentAvg = currentSum/numBridgeNodes;
				avgWeight[winner] += 1;

				printf("%d: Processed node %d, currentSum %4.2lf currentAvg %4.2f weight[%d](%d) = %d\n", myrank, j, currentSum, currentAvg, winner, bridgeRanks[winner], avgWeight[winner]);
			}
		}

		for (bn=0; bn<numBridgeNodes ; bn++) {
			printf("%d: currentAvg %4.2f weight[%d](%d) = %d\n", myrank, currentAvg, bn, bridgeRanks[bn], avgWeight[bn]);
			if (bridgeRanks[bn] == myrank) myWeight = avgWeight[bn];
		}

		int k=-1;
		shuffledNodes = (int *) malloc (myWeight * sizeof(int));
		for (j=0; j<myWeight ; j++) 
			shuffledNodes[j] = -1;

		for (j=0; j<MidplaneSize ; j++) {
		 	//if (revisit[j][0] != 255 && myrank == revisit[j][2]) 
		 	if (revisit[j][0] != 255 && newBridgeNode[j] != -1) 
				printf("%d: %d (%d) is the new BN for %d: Old: %d (%d)\n", myrank, newBridgeNode[j], revisit[j][1], j, bridgeNodeAll[j*2], bridgeNodeAll[j*2+1]);
			if (newBridgeNode[j] == myrank) 
					shuffledNodes[++k] = j;
		}

		if (k+1 != myWeight) printf("%d: Error in number of shuffledNodes: %d %d\n", myrank, k, myWeight);

	
	MPI_Bcast(&newBridgeNode, MidplaneSize, MPI_INT, bridgeRanks[0], MPI_COMM_WORLD);	//TODO subcomm for the midplane
	MPI_Bcast(&revisit, MidplaneSize*2, MPI_BYTE, bridgeRanks[0], MPI_COMM_WORLD);	//TODO subcomm for the midplane

*/

		
	}


}


void distributeInfo() {

	int j;

#ifdef DEBUG
		printf ("\n%d:%d:%d: myWeight = %d\n", myrank, coreID, nodeID, myWeight);
#endif
//
		MPI_Bcast(&myWeight, ppn, MPI_INT, 0, MPI_COMM_node);	
		MPI_Bcast(bridgeNodeAll, 2*MidplaneSize*ppn, MPI_INT, 0, MPI_COMM_node);	
#ifdef DEBUG
		printf ("\n%d:%d:%d: real myWeight = %d\n", myrank, coreID, nodeID, myWeight);
#endif

		shuffledNodes = (int *) malloc (myWeight * sizeof(int));
		for (j=0; j<myWeight ; j++) 
			shuffledNodes[j] = -1;

		//MPI_Barrier (COMM_BRIDGE_NODES);	

		int k=-1;
		//for (j=0; j<MidplaneSize ; j++) { 
		for (j=0; j<MidplaneSize*ppn ; j=j+ppn) { 
			if (newBridgeNode[j] >= 0) 
#ifdef DEBUG
				printf ("%d: checking: %d %d %d %d %d\n", myrank, nodeID, coreID, j, bridgeRanks[newBridgeNode[j]], bridgeNodeAll[j*2+1]);
#endif
			//if (bridgeRanks[newBridgeNode[j]] == myrank && bridgeNodeAll[j*2+1]>1) {
			if (bridgeRanks[newBridgeNode[j]] == nodeID*ppn && bridgeNodeAll[j*2+1]>1) {
				k++;
				shuffledNodes[k] = j+coreID;
#ifdef DEBUG
				printf("%d:(%d,%d) k=%d newBridgeNode[%d]=%d bn=%d %d\n", myrank, nodeID, coreID, k, j, newBridgeNode[j], bridgeRanks[newBridgeNode[j]], shuffledNodes[k]); 
#endif
			}
		}

#ifdef DEBUG
		if (k+1 != myWeight) printf("%d: Error in number of shuffledNodes: %d %d\n", myrank, k, myWeight);
#endif
//

//		for (j=0; j<MidplaneSize*ppn ; j++) 
//			if (newBridgeNode[j] >= 0)
//				printf("%d: %d (%d) is the new BN for %d\n", myrank, bridgeRanks[newBridgeNode[j]], newBridgeNode[j], j);

}

void initTree(int n) {

	bridgeNodeRootList = (Node **) malloc (n*sizeof(Node *));
#ifdef DEBUG
	if (bridgeNodeRootList == NULL)	printf("\n%d: malloc fail\n", myrank);
#endif

}

int main (int argc, char **argv) {

		MPI_Init (&argc, &argv);

		MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
		MPI_Comm_size (MPI_COMM_WORLD, &commsize);

		if (argc > 1){
			int exp = atoi (argv[1]);
			fileSize = exp * oneKB;	
		}

//TODO
//	NONBLOCKING = atoi (argv[2]);

		int i, c;
		opterr = 0;
  	while ((c = getopt (argc, argv, "d:")) != -1)
    		switch (c)
      		{
      			case 'd':
        			//#define DEBUG = 1;
        			break;
      			default:
        			//DEBUG = 0;
        			break;
		}

		double tAStart, tAEnd; 

//TODO if I am on core 0
		getPersonality(myrank);

		bridgeNodeAll = new int [2*MidplaneSize*ppn];
		visited = new bool [MidplaneSize*ppn];
		processed = new bool [MidplaneSize*ppn];
		newBridgeNode = new int [MidplaneSize*ppn];

		revisit = new uint8_t *[MidplaneSize*ppn];
		for (i=0 ; i<MidplaneSize*ppn ; i++)
			revisit[i] = new uint8_t[2];

		depthInfo = new uint8_t *[numBridgeNodes];
		for (i=0 ; i<numBridgeNodes ; i++)
			depthInfo[i] = new uint8_t[MidplaneSize*ppn];

		bridgeRanks = new int [numBridgeNodes];
		avgWeight = new float [numBridgeNodes];

		rootNodeList = new queue<Node *> [numBridgeNodes];

		rootps = floor(myrank/(MidplaneSize*ppn));

		double tStart = MPI_Wtime();	//entire execution

		initNeighbours();
		initTree(numBridgeNodes);

		//form inter-communicator - mainly reqd for core 0 processes 
		MPI_Comm_split (MPI_COMM_WORLD, coreID, myrank, &MPI_COMM_core);

		//form intra-communicator - mainly reqd for processes on a node
		MPI_Comm_split (MPI_COMM_WORLD, nodeID, myrank, &MPI_COMM_node);

		//form intra-communicator - mainly reqd for bridge nodes
		MPI_Comm_split (MPI_COMM_WORLD, bridgeNodeInfo[1], myrank, &COMM_BRIDGE_NODES);

		//gather bridgeNodeInfo at the rootps
		double ts = MPI_Wtime();
		MPI_Gather (bridgeNodeInfo, 2, MPI_INT, bridgeNodeAll, 2, MPI_INT, rootps, MPI_COMM_WORLD);
		ts = MPI_Wtime() - ts;

#ifdef DEBUG
		printf("%d: mybridgeNodeInfo: %d %d\n", myrank, bridgeNodeInfo[0], bridgeNodeInfo[1]);
		if (myrank == rootps)
			printf("%d: bridgeNodeInfo %d %d %d %d %lf\n", myrank, bridgeNodeAll[2], bridgeNodeAll[3], bridgeNodeAll[6], bridgeNodeAll[7], ts);
#endif

		double tOStart = MPI_Wtime();
		if (coreID == 0) formBridgeNodesRoutes ();

		//MPI_Bcast(newBridgeNode, MidplaneSize, MPI_INT, rootps, MPI_COMM_WORLD);	
		MPI_Bcast(newBridgeNode, MidplaneSize*ppn, MPI_INT, rootps, MPI_COMM_WORLD);	
		MPI_Bcast(bridgeRanks, numBridgeNodes, MPI_INT, rootps, MPI_COMM_WORLD);	

		if (bridgeNodeInfo[1] == 1) distributeInfo();

		double tOEnd = MPI_Wtime();

#ifdef DEBUG
		if (coreID == 0) printf("%d: %d MyNewBN %d %d\n", myrank, ppn, newBridgeNode[myrank], bridgeRanks[newBridgeNode[myrank]]);
		//printf("%d: overhead %6.3f\n", myrank, tOEnd-tOStart);
#endif

		MPI_Barrier (MPI_COMM_WORLD);

		int count = fileSize;				//weak scaling
		dataBlock *datum = new dataBlock(count);			//initializes alpha array to random double values

		//Testing BGQ compute nodes to IO nodes performance
		//Write to /dev/null
		//writeFlag = checkION();

		/*
		 * * * * * * * * * * * * Independent MPI-IO to IO nodes from all compute nodes - shared file * * * * * * * * * *
		 */

		double tION[2];

		/* open file */
		mode = MPI_MODE_CREATE | MPI_MODE_RDWR; //WRONLY;

		/* allocate buffer */
		datum->allocElement (1);

		MPI_File_open (MPI_COMM_WORLD, fileNameION, mode, MPI_INFO_NULL, &fileHandle);
		tIOStart = MPI_Wtime();
		totalBytes[0][0] += writeFile(datum, count, 0);			
		totalBytes[0][0] += writeFile(datum, count, 0);			
		totalBytes[0][0] += writeFile(datum, count, 0);			
		totalBytes[0][0] += writeFile(datum, count, 0);	
		tIOEnd = MPI_Wtime();
		tION[0] = tIOEnd - tIOStart;
		MPI_File_close (&fileHandle);

		MPI_File_open (MPI_COMM_WORLD, fileNameION, mode, MPI_INFO_NULL, &fileHandle);
		tIOStart = MPI_Wtime();
		totalBytes[0][1] += writeFile(datum, count, 1);
		totalBytes[0][1] += writeFile(datum, count, 1);
		totalBytes[0][1] += writeFile(datum, count, 1);
		totalBytes[0][1] += writeFile(datum, count, 1);
		tIOEnd = MPI_Wtime();
		tION[1] = tIOEnd - tIOStart;
		MPI_File_close (&fileHandle);

		/*
		 * * * * * * * * * * * * * * * * * * * * * Independent MPI-IO to file system from all compute nodes - shared file * * * * * * * * * * * * * * * * * *
		 */

		double te[2];

		//mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;

		MPI_File_open (MPI_COMM_WORLD, fileNameFS, mode, MPI_INFO_NULL, &fileHandle);
		double t2 = MPI_Wtime();
		totalBytes[1][1] += writeFile(datum, count, 1);
		totalBytes[1][1] += writeFile(datum, count, 1);
		totalBytes[1][1] += writeFile(datum, count, 1);
		totalBytes[1][1] += writeFile(datum, count, 1);
		te[1] = MPI_Wtime() - t2;
		MPI_File_close (&fileHandle);

		MPI_File_open (MPI_COMM_WORLD, fileNameFSBN, mode, MPI_INFO_NULL, &fileHandle);
		double t1 = MPI_Wtime();
		totalBytes[1][0] += writeFile(datum, count, 0);
		totalBytes[1][0] += writeFile(datum, count, 0);
		totalBytes[1][0] += writeFile(datum, count, 0);
		totalBytes[1][0] += writeFile(datum, count, 0);
		te[0] = MPI_Wtime() - t1;
		MPI_File_close (&fileHandle);

		/* free buffer */
		datum->freeElement (1);

		//* * * * * * * * * * * * * * * * * * * * * * * * * * * End of IO * * * * * * * * * * * * * * * * * * * * * * * * * * *

		double tEnd = MPI_Wtime();
		
		double max[4];

		MPI_Reduce(&tION[0], &max[0], 1, MPI_DOUBLE, MPI_MAX, rootps, MPI_COMM_WORLD);
		MPI_Reduce(&tION[1], &max[1], 1, MPI_DOUBLE, MPI_MAX, rootps, MPI_COMM_WORLD);
		MPI_Reduce(&te[0], &max[2], 1, MPI_DOUBLE, MPI_MAX, rootps, MPI_COMM_WORLD);
		MPI_Reduce(&te[1], &max[3], 1, MPI_DOUBLE, MPI_MAX, rootps, MPI_COMM_WORLD);

		MPI_Finalize ();

		if (myrank == rootps) {
			printf ("Times: %d %d %d # %d | %lf %lf | %lf %lf\n", commsize, ppn, omp_get_num_threads(), 4*8*fileSize, max[0], max[1], max[2], max[3]);
			//printf ("Times: %d %4.2f %4.2f %4.2f %d %6.3f\n", myrank, totalBytes[0][0]*1.0/oneMB, tION_elapsed_0, tION_elapsed_1, bridgeNodeInfo[0], tEnd-tStart);   // rank, MB, MB, sec..
		}

		return 0;

}

