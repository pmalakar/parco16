#ifndef __bgqBridgeNode__
#define __bgqBridgeNode__

#define NONBLOCKING 0

#define prnl printf("\n")

#define __STDC_FORMAT_MACROS 1 // In C++ the macros are not automatically defined just by including the file. 

using namespace std;

const int INITIALIZER = -1;

//BGQ specific

int numNodes, myWeight, myBNIdx;

//float avgWeight[numBridgeNodes];
float currentSum=0.0, currentAvg=0.0;

class Node {

	private:
		int nodeId, rootId, lastIndex;
		int depth;				//depth of the node in the tree
		int childIndex[9];
		Node *childIndexPtr[9];
		Node *parent;

	public:

		Node *leftSibling, *rightSibling;	//tree
		Node *prev, *next;					//queue
		Node() { lastIndex = -1;}
		~Node() {}

		Node (int currentNode) {
			lastIndex = -1;
			nodeId = currentNode;			
			parent = NULL;
			depth=0;
			rootId = currentNode;
			//printf("depth %d\n", depth);
		}

		Node* addChild (int id, int rootid) {
			Node *node = new Node(id);
			node->parent = this;
			node->depth = depth + 1;
			node->rootId = rootid;
			childIndexPtr[++lastIndex] = node;
			return node;
		}		

		Node* getChild (int i) {
			if (i>lastIndex) return NULL;
			return childIndexPtr[i];
		}
	
		int getChildId (int i) {
			if (i>lastIndex) return -1;
			return childIndexPtr[i]->getNodeId();
		}

		int getNodeId () 			{ return nodeId; }
		int getRootId () 			{ return rootId; }
		int getNumChildren () { return lastIndex+1; }
		int getDepth() 				{ return depth; }
		Node *getParent () 		{ return parent; }

		void printChildrenInfo (int rank) {
			int i;
			printf("%d: num children of %d = %d\n", rank, nodeId, lastIndex);
			for (i=0; i<=lastIndex; i++)
				if (rank == nodeId) printf("First level of (%d): %d[%d]\n ", nodeId, childIndexPtr[i]->getNodeId(), i);
		}

}*head, *tail, *root;

//queue <Node *> nodeList;
//queue <Node *> rootNodeList[numBridgeNodes];

//Node **bridgeNodeRootList;

int writeFlag=1; 

// Variables
int oneKB = 1024;
int oneMB = 1024*1024;
int myrank, commsize, mode, fileSize;
int totalBytes[2][3] = {0};

MPI_File fileHandle;
MPI_Status status;
MPI_Request request;
char *fileNameION = "/dev/null";
char *fileNameFS = "dummyFile";
char *fileNameFSBN = "dummyFileBN";
char *fileNameFSCO = "dummyFileCO";

double tIOStart, tIOEnd;
double tION_elapsed[2]={0.0,0.0}, tFS_elapsed[2]={0.0,0.0};

//store results
double *alphaSum;	//reduce result stored at root

#define MAXBUF (1024*32)

int hPSet;  // Punit Event Set handle
int hL2Set;  // Punit Event Set handle
int hNWSet; // Network Event Set handle
int hIOSet; // I/O Event Set handle

#define xMicroSec 1000000

// Function Declarations

// Functions

//adapted verbatim from https://wiki.alcf.anl.gov/parts/index.php/Blue_Gene/Q#Allocating_Memory
void * bgq_malloc(size_t n)
{
    void * ptr;
    size_t alignment = 32; /* 128 might be better since that ensures every heap allocation 
                            * starts on a cache-line boundary */
    posix_memalign( &ptr , alignment , n );
    return ptr;
}


double ClockSpeed()
{
    Personality_t personality;
    Kernel_GetPersonality(&personality, sizeof(Personality_t));
    return (personality.Kernel_Config.FreqMHz * 1.0e6);
}


//data used for computation, analysis, IO etc..

class dataBlock {

 private:	
	double *alpha, *beta, *gamma;
//    double U[1024][1024], V[1024][1024], W[1024][1024];
	double **U, **V, **W;
	int numElem;
	int latsize, lonsize;

 public:
	dataBlock() {}
	dataBlock(int count) {

		numElem = count;

	}

	void allocElement (int type) {
		
		if (type == 1) {
		  try {
			alpha = new double[numElem];
			for (int i=0; i<numElem ; i++) {
				alpha[i] = rand() % 100;
			}
		  }
		  catch (bad_alloc& ba) {
				cerr << "Bad allocation for alpha\n" << ba.what() << endl;
		  }
		}
	}

	void freeElement (int type) {

		if (type == 1) delete [] alpha;
	}

	double *getAlphaBuffer() {
		return &alpha[0];
	}

};


int writeFile (dataBlock *, int, int);

/* Allocate and free memory before and after use */ //Though this is bit of overhead but .. Just to ensure there is no spatial locality interference affecting the statistics
    
void alloc_free (dataBlock *datum, int type) {
    
    if (type == 1) {

		datum->allocElement(1);
		datum->allocElement(2);
		datum->allocElement(3);

    }
    else {

		datum->freeElement(1);
		datum->freeElement(2);
		datum->freeElement(3);

    }

}

void getData(int);

void prnerror (int error_code, char *string)
{
	
	char error_string[256];
	int length_of_error_string;	
	MPI_Error_string(error_code, error_string, &length_of_error_string);
	fprintf(stderr, "%3d: %s in %s\n", error_code, error_string, string);
	MPI_Finalize();
	exit(-1);
}

int min (int a, int b) {
	if (a<=b) return a;
	else return b;
}

int findNeighbours (int);

#endif
