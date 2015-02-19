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

int *bridgeNodeAll; 	//[MidplaneSize*2];				//2 integers per rank
/*
uint8_t **depthInfo; //[numBridgeNodes][MidplaneSize];
bool *visited, *processed;		//[MidplaneSize], processed[MidplaneSize];
uint8_t **revisit;	//[MidplaneSize][2];
int *newBridgeNode;	//[MidplaneSize];
*/

//int bridgeNodeAll[MidplaneSize*2];				//2 integers per rank
uint8_t depthInfo[numBridgeNodes][MidplaneSize];
bool visited[MidplaneSize], processed[MidplaneSize];
uint8_t revisit[MidplaneSize][2];
int newBridgeNode[MidplaneSize];


/*
 *  Independent MPI-IO
 */
int writeFile(dataBlock *datum, int count, int level, int all) {

		double start=0.0, end=0.0, test1=0.0, test2=0.0;
		int nbytes, totalBytes = 0, i;

		/* write to ION or storage */
		//Should I write myself?
		//if (hopsToBridgeNode > 1)
		int result;
		if (all == 0) {
			tIOStart = MPI_Wtime();

			if (newBridgeNode[myrank] != -1 && bridgeNodeInfo[1] > 1) {	
//				MPI_Request req[4];
//				MPI_Isend (datum->getAlphaBuffer(), count, MPI_DOUBLE, bridgeRanks[newBridgeNode[myrank]], bridgeRanks[newBridgeNode[myrank]], MPI_COMM_WORLD, &req[0]);	
//				MPI_Isend (datum->getAlphaBuffer(), count, MPI_DOUBLE, bridgeRanks[newBridgeNode[myrank]], bridgeRanks[newBridgeNode[myrank]], MPI_COMM_WORLD, &req[1]);	
//				MPI_Isend (datum->getAlphaBuffer(), count, MPI_DOUBLE, bridgeRanks[newBridgeNode[myrank]], bridgeRanks[newBridgeNode[myrank]], MPI_COMM_WORLD, &req[2]);	
//				MPI_Isend (datum->getAlphaBuffer(), count, MPI_DOUBLE, bridgeRanks[newBridgeNode[myrank]], bridgeRanks[newBridgeNode[myrank]], MPI_COMM_WORLD, &req[3]);	
//				MPI_Waitall (1, req, stat);
					MPI_Send (datum->getAlphaBuffer(), count, MPI_DOUBLE, bridgeRanks[newBridgeNode[myrank]], bridgeRanks[newBridgeNode[myrank]], MPI_COMM_WORLD);	
			}
			//if i am bridge node, mpi_recv from all my sender nodes
			if (bridgeNodeInfo[1] == 1) {

				int arrayLength = myWeight;//myWeight*4

				MPI_Request req[arrayLength], wrequest[myWeight+1];
				MPI_Status stat, wstatus[myWeight+1];

				double shuffledNodesData[arrayLength][count];

			//TODO file seek //file set view??	
				//result = MPI_File_write_shared (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);

				//MPI_Offset disp = myrank *count*sizeof(double); 
				//MPI_File_seek (fileHandle, disp, MPI_SEEK_SET);
				//result = MPI_File_write(fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				//if (NONBLOCKING == 1)
				//	result = MPI_File_iwrite(fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &wrequest[myWeight]);

				result = MPI_File_write_at (fileHandle, myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				//result = MPI_File_iwrite_at (fileHandle, myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &wrequest[myWeight]);

#pragma omp parallel for
				//get data from my senders
				for (i=0; i<myWeight ; i++) {
				//	if (shuffledNodes[i] == -1) break;
					MPI_Irecv (shuffledNodesData[i], count, MPI_DOUBLE, shuffledNodes[i], myrank, MPI_COMM_WORLD, &req[i]); 
				}

          test1 = MPI_Wtime();
#pragma omp parallel for
				for (i=0; i<myWeight ; i++) {
					int idx;
					MPI_Waitany (myWeight, req, &idx, &stat);
					//result = MPI_File_iwrite_at (fileHandle, shuffledNodes[idx]*count*sizeof(double), shuffledNodesData[idx], count, MPI_DOUBLE, &wrequest[idx]);
					//result = MPI_File_write_shared (fileHandle, shuffledNodesData[idx], count, MPI_DOUBLE, &status);
					
//					disp = shuffledNodes[idx]*count*sizeof(double); 
//					MPI_File_seek (fileHandle, disp, MPI_SEEK_SET);
//					result = MPI_File_write (fileHandle, shuffledNodesData[idx], count, MPI_DOUBLE, &status);
					if (NONBLOCKING == 1)
						result = MPI_File_iwrite (fileHandle, shuffledNodesData[idx], count, MPI_DOUBLE, &wrequest[i]);
					result = MPI_File_write_at (fileHandle, shuffledNodes[idx]*count*sizeof(double), shuffledNodesData[idx], count, MPI_DOUBLE, &status);
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
				printf("%d %d BN write time %lf test2=%lf\n", myrank, level, tIOEnd, test2);
			}

		}
		else {
			 
				MPI_Request wrequest;
				start = MPI_Wtime();
				//result = MPI_File_write_shared (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				result = MPI_File_write_at (fileHandle, myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
				//result = MPI_File_write_all (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);

				//result = MPI_File_iwrite_at (fileHandle, myrank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &wrequest);

	//			MPI_Offset disp = myrank *count*sizeof(double); 
	//			MPI_File_seek (fileHandle, disp, MPI_SEEK_SET); 
	//			result = MPI_File_write(fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);

//				result = MPI_File_iwrite(fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &wrequest);
//				MPI_Wait(&wrequest, &status);

				if ( result != MPI_SUCCESS) {
					prnerror (result, "MPI_File_write_at Error:");
				}
				start = MPI_Wtime()-start;
				printf("%d %d all write time %lf\n", myrank, level, start);

		}


		/* Worse performance
		MPI_Offset disp = myrank*count*sizeof(double); 
		int result = MPI_File_set_view (fileHandle, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
		result = MPI_File_write_all (fileHandle, datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
		*/

/*
		if ( result != MPI_SUCCESS) {
			prnerror (result, "MPI_File_write_at Error:");
		}*/
		//MPI_Get_elements( &status, MPI_CHAR, &nbytes );
		//totalBytes += nbytes;

		return totalBytes;

}


int build5DTorus (int root) {

	double tStart = MPI_Wtime();
	for (int rank = root ; rank < root+MidplaneSize ; rank++) {
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

		for(int diff=0; diff<hopDiff ;diff++) {

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
#endif
	if(!nodeList.empty()) printf("%d: queue size %d\n", myrank, nodeList.size());

	//if (childNum != 3) printf ("%d: ERROR: childNum not 3 but %d\n", myrank, childNum);
int BAG = 256;
#ifdef VESTA
	BAG = 64;
#endif

	if (numNodes > BAG) {
	//	numNodes = 0;
		//no more expansion for this node?
		while(!nodeList.empty())	nodeList.pop();	
		printf("%d: queue size %d\n", myrank, nodeList.size());
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
		else
			printf("%d: nodeList empty\n");
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
	int root = floor(myrank/MidplaneSize);

//	int numBridgeNodes = 8; 		// TODO Number of bridge nodes?
//	int bridgeRanks[numBridgeNodes];

	if (myrank == root) {

		MPI_Request request[numBridgeNodes], requestAll[numBridgeNodes];
		MPI_Status status[numBridgeNodes];

#ifdef DEBUG
		printf("I am the root %d\n", myrank);
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
		
			//Introduce all bridge nodes to each other
			result = MPI_Isend (bridgeRanks, numBridgeNodes, MPI_INT, bridgeRanks[i], 101, MPI_COMM_WORLD, &request[i]);
			if (result != MPI_SUCCESS) 
				prnerror (result, "MPI_Isend Error: ");

			//Send all nodes' bridgeNodeInfo to all bridge nodes//TODO hardcoded 512
			result = MPI_Isend (bridgeNodeAll, 2*MidplaneSize, MPI_INT, bridgeRanks[i], 102, MPI_COMM_WORLD, &requestAll[i]);
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
		printf("%d: send recv overhead in root process %6.3f\n", myrank, tEnd-tStart);
#endif

		MPI_Status st;	
		MPI_Recv(newBridgeNode, MidplaneSize, MPI_INT, bridgeRanks[0], bridgeRanks[0], MPI_COMM_WORLD, &st);	
	}
	
	if (bridgeNodeInfo[1] == 1) {

		MPI_Request requestSend, requestRecv, requestRecvAll;
		MPI_Status statusSend, statusRecv, statusRecvAll;

		tStart = MPI_Wtime();

		result = MPI_Isend (&myrank, 1, MPI_INT, root, 100, MPI_COMM_WORLD, &requestSend);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Isend Error: ");

		result = MPI_Irecv (bridgeRanks, numBridgeNodes, MPI_INT, root, 101, MPI_COMM_WORLD, &requestRecv);
		if (result != MPI_SUCCESS) 
			prnerror (result, "MPI_Irecv Error: ");

		result = MPI_Irecv (bridgeNodeAll, 2*MidplaneSize, MPI_INT, root, 102, MPI_COMM_WORLD, &requestRecvAll);
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

		result = build5DTorus (root);

		//	initialize MidplaneSize array
		//	all nodes but the bridge nodes are unvisited
		
		for (i=0; i<MidplaneSize ; i++) 
			  visited[i] = false, processed[i] = false;

		for (int bn=0; bn<numBridgeNodes ; bn++) { 
		  visited[bridgeRanks[bn]] = true;
		  for (i=0; i<MidplaneSize ; i++) {
				revisit[i][0] = 255, revisit[i][1] = 255;
		  	//depthInfo[bn][i] = bridgeNodeAll[i*2+1];			
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
		for (j=0; j<MidplaneSize ; j++) 
			if (newBridgeNode[j] > 0)
				printf("%d: %d (%d) is the new BN for %d\n", myrank, bridgeRanks[newBridgeNode[j]], newBridgeNode[j], j);
#endif

		for (bn=0; bn<numBridgeNodes ; bn++) 
			if (bridgeRanks[bn] == myrank) myWeight = int(avgWeight[bn]);

	  if (myrank == bridgeRanks[0])
			MPI_Send(newBridgeNode, MidplaneSize, MPI_INT, root, bridgeRanks[0], MPI_COMM_WORLD);	//TODO subcomm for the midplane
	
#ifdef DEBUG
		printf ("\n%d: myWeight = %d\n", myrank, myWeight);
#endif
		shuffledNodes = (int *) malloc (myWeight * sizeof(int));
		for (j=0; j<myWeight ; j++) 
			shuffledNodes[j] = -1;

		//MPI_Barrier (COMM_BRIDGE_NODES);	

		int k=-1;
		for (j=0; j<MidplaneSize ; j++) { 
			if (bridgeRanks[newBridgeNode[j]] == myrank && bridgeNodeAll[j*2+1]>1) {
				k++;
#ifdef DEBUG
				printf("%d: k=%d newBridgeNode[%d]=%d bn=%d %d\n", myrank, k, j, newBridgeNode[j], bridgeRanks[newBridgeNode[j]], shuffledNodes[k]); 
#endif
				shuffledNodes[k] = j;
			}
		}

#ifdef DEBUG
		if (k+1 != myWeight) printf("%d: Error in number of shuffledNodes: %d %d\n", myrank, k, myWeight);
#endif

	

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

	//bcast from root 
	MPI_Bcast(newBridgeNode, MidplaneSize, MPI_INT, root, MPI_COMM_WORLD);	
	MPI_Bcast(bridgeRanks, numBridgeNodes, MPI_INT, root, MPI_COMM_WORLD);	

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

		bridgeNodeAll = new int [MidplaneSize*2];
/*
		visited = new bool [MidplaneSize];
		processed = new bool [MidplaneSize];
		newBridgeNode = new int [MidplaneSize];

		depthInfo = new uint8_t *[numBridgeNodes];
		for (i=0; i<numBridgeNodes ; i++)
			depthInfo[i] = new uint8_t[MidplaneSize];
		revisit = new uint8_t *[MidplaneSize];
		for (i=0; i<MidplaneSize ; i++)
			revisit[i] = new uint8_t[2];
*/

		double tStart = MPI_Wtime();	//entire execution

//TODO if I am on core 0
		getPersonality(myrank);

		initTree(numBridgeNodes);

#ifdef DEBUG
		printf("%d: mybridgeNodeInfo: %d %d\n", myrank, bridgeNodeInfo[0], bridgeNodeInfo[1]);
#endif 
		//gather bridgeNodeInfo at the root
		//int bridgeNodeAll[512*2];				//2 integers per rank
		double ts = MPI_Wtime();
//TODO per partition create sub comm for core0
//TODO if I am on core 0
		MPI_Gather (bridgeNodeInfo, 2, MPI_INT, bridgeNodeAll, 2, MPI_INT, 0, MPI_COMM_WORLD);
		ts = MPI_Wtime() - ts;
#ifdef DEBUG
		if (myrank == 0)
			printf("%d: bridgeNodeInfo %d %d %d %d %lf\n", myrank, bridgeNodeAll[2], bridgeNodeAll[3], bridgeNodeAll[6], bridgeNodeAll[7], ts);
#endif

		//form intra-communicator - mainly reqd for bridge nodes
		MPI_Comm_split (MPI_COMM_WORLD, bridgeNodeInfo[1], myrank, &COMM_BRIDGE_NODES);

		double tOStart = MPI_Wtime();
		formBridgeNodesRoutes ();
		double tOEnd = MPI_Wtime();

		//printf("%d: overhead %6.3f\n", myrank, tOEnd-tOStart);

		int count = fileSize;				//weak scaling
		dataBlock *datum = new dataBlock(count);			//initializes alpha array to random double values

		//Testing BGQ compute nodes to IO nodes performance
		//Write to /dev/null
		//writeFlag = checkION();

		/*
		 * * * * * * * * * * * * * * * * * * * * * Independent MPI-IO to IO nodes from all compute nodes - shared file * * * * * * * * * * * * * * * * * *
		 */

		/* open file */
		mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
		MPI_File_open (COMM_BRIDGE_NODES, fileNameION, mode, MPI_INFO_NULL, &fileHandle);

		/* allocate buffer */
		datum->allocElement (1);

		tIOStart = MPI_Wtime();

		totalBytes[0][0] += writeFile(datum, count, 0, 0);			
		/*totalBytes[0] += writeFile(datum, count, 0);		
		totalBytes[0] += writeFile(datum, count, 0);	
		totalBytes[0] += writeFile(datum, count, 0);
		totalBytes[0] += writeFile(datum, count, 0);			//writes 
*/
		tIOEnd = MPI_Wtime();
		double tION_elapsed_0 = tIOEnd - tIOStart;

		/* close file */
		MPI_File_close (&fileHandle);

		MPI_File_open (MPI_COMM_WORLD, fileNameION, mode, MPI_INFO_NULL, &fileHandle);
		tIOStart = MPI_Wtime();
		totalBytes[0][1] += writeFile(datum, count, 0, 1);
/*
		totalBytes[0] += writeFile(datum, count, 1);			
		totalBytes[0] += writeFile(datum, count, 1);		
		totalBytes[0] += writeFile(datum, count, 1);	
		totalBytes[0] += writeFile(datum, count, 1);
		totalBytes[0] += writeFile(datum, count, 1);			//writes 
*/
		tIOEnd = MPI_Wtime();
		double tION_elapsed_1 = tIOEnd - tIOStart;

		/* close file */
		MPI_File_close (&fileHandle);

		/*
		 * * * * * * * * * * * * * * * * * * * * * Independent MPI-IO to file system from all compute nodes - shared file * * * * * * * * * * * * * * * * * *
		 */

		mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
		MPI_File_open (COMM_BRIDGE_NODES, fileNameFS, mode, MPI_INFO_NULL, &fileHandle);
		totalBytes[1][0] += writeFile(datum, count, 1, 0);
		MPI_File_close (&fileHandle);

		MPI_File_open (MPI_COMM_WORLD, fileNameFS, mode, MPI_INFO_NULL, &fileHandle);
		totalBytes[1][1] += writeFile(datum, count, 1, 1);
		MPI_File_close (&fileHandle);

		/* free buffer */
		datum->freeElement (1);

		//* * * * * * * * * * * * * * * * * * * * * * * * * * * End of IO * * * * * * * * * * * * * * * * * * * * * * * * * * *


		double tEnd = MPI_Wtime();
		MPI_Finalize ();

		printf ("Times: %d %4.2f %4.2f %4.2f %d %6.3f\n", myrank, totalBytes[0][0]*1.0/oneMB, tION_elapsed_0, tION_elapsed_1, bridgeNodeInfo[0], tEnd-tStart);   // rank, MB, MB, sec..

		return 0;

}
