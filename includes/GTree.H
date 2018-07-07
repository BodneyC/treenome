#ifndef _GTREE_
#define _GTREE_

#if __cplusplus < 201103L
#define nullptr 0
#endif /*c++11 check*/

#define NBASES 4

#include <iostream>
#include <string>
#include <vector>
#include "SeqRead.H"

typedef struct Node {
	long occs;
	float weight;
	short offset;
	long readNum;
	Node* subnodes[NBASES];

	Node(): occs(0), weight(0), offset(0), readNum(0)
	{
		for(int i = 0; i < NBASES; i++)
			subnodes[i] = nullptr;
	}
} Node;

// GTree Helpers
namespace GTH {
	extern std::vector<SeqRead> seqReads;

	char retLabel(int label);
	short mostOccs(Node *node);
	short countChildren(Node *node);
	void updateWeight(Node *node, char qual);
}

// GTree class used in construction and destruction of tree
class GTree {
public:
	~GTree() { deleteTree(root); }

	/* Tree creation functions */
	void createRoot(short ind);
	void addReadOne(long readNum, short offset); 
	void balanceNode(Node *node);

	/* Sequence creation functions */
	void addToSeq(long offset, std::string &sequence);
	void followPath(Node *node, short ind, std::string &sequence);

	/* Usefull non-essential */
	Node *getRoot() { if(root) return root; 
		else return nullptr; }
	void printAllPaths(short label) { 
		printAllPaths(root, 0, label);
		basePaths = occuPaths = "";
	}
private:
	Node* root;
	std::string basePaths;
	std::string occuPaths;

	/* Creation/destruction */
	void createNode(Node *node, short label, char qual);
	void deleteTreeLL(Node **node, std::string readSub);
	void deleteTree(Node *root);

	/* Public->private overloads */
	void printAllPaths(Node *node, int len, short label);
};

// Structure which holds read data and the four GTrees
typedef struct TreeTop {
	std::string filename, sequence;
	long nReads;
	int readLength;
	bool ioSuccess;
	GTree trees[NBASES];

	explicit TreeTop(const std::string &filename);
	bool readSuccess() { return ioSuccess;  }
	void createRoots();
	void processReadsOne();
	void buildSequence();
	short maxPath();
	void printSequence();
	bool rootOccsExist();
} TreeTop;

#endif /*_GTREE_*/
