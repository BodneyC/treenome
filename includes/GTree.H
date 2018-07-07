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

typedef std::vector<std::string> VecString;
typedef unsigned long ulong;

typedef struct Node {
	ulong occs;
	float weight;
	short offset;
	ulong readNum;
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
}

// GTree class used in construction and destruction of tree
class GTree {
public:
	~GTree() { deleteTree(root); }

	/* Tree creation functions */
	void createRoot(short ind,
			VecString &reads,
			VecString &quals);
	void addReadFull(ulong readNum, short offset, 
			std::string &read, std::string &qual);
	void addReadOne(ulong readNum, short offset, 
			VecString &reads,
			VecString &quals);
	void balanceNode(Node *node,
			VecString &reads,
			VecString &quals);

	/* Sequence creation functions */
	void addToSeq(ulong offset, std::string &sequence);
	void followPath(Node *node, short ind, std::string &sequence);

	/* Used in readFull...()s */
	// cleanBranches..() only needed if ..ReadFull() is used
	void cleanBranches() { cleanBranches(root); }
	Node* cleanBranchesNR(short offset, std::string &read);

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

	/* Helper function */
	short countChildren(Node *node);
	void updateWeight(Node *node, char qual);

	/* Public->private overloads */
	void cleanBranches(Node *node);
	void printAllPaths(Node *node, int len, short label);
};

// Structure which holds read data and the four GTrees
typedef struct TreeTop {
	std::string filename, sequence;
	ulong nReads;
	int readLength;
	bool ioSuccess;
	VecString reads; 
	VecString quals;
	GTree trees[NBASES];

	TreeTop(const std::string &filename);
	bool readSuccess() { return ioSuccess;  }
	void createRoots();

	void TMPTMP();

	void processReadsOne();
	void processReadsFull();
	void processReadsFullClean();
	void processReadsFullCleanNR();
	void processReadsFullCleanNRBalanced();

	void buildSequence();
	short maxPath();
} TreeTop;

#endif /*_GTREE_*/
