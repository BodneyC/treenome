#ifndef _GTREE_
#define _GTREE_

#if __cplusplus < 201103L
#define nullptr 0
#endif /*c++11 check*/

#define NBASES 4

// Little bithack converting ACTG -> 0123
#define BASE_IND(x) (((x) & 0xF) >> 1)

#include <iostream>
#include <string>
#include <vector>

typedef unsigned long ulong;
typedef unsigned int uint;

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

// GTree class used in construction and destruction of tree
class GTree {
public:
	~GTree() { deleteTree(root); }

	Node *getRoot() { if(root) return root; 
		else return nullptr; }
	void createRoot(short ind,
			std::vector<std::string> &reads,
			std::vector<std::string> &quals);

	void addReadFull(ulong readNum, short offset, 
			std::string &read, std::string &qual);
	void addReadOne(ulong readNum, short offset, 
			std::vector<std::string> &reads,
			std::vector<std::string> &quals);

	// cleanBranches..() only needed if ..ReadFull() is used
	void cleanBranches() { cleanBranches(root); }
	Node* cleanBranchesNR(short offset, std::string &read);
	void balanceNode(Node *node,
			std::vector<std::string> &reads,
			std::vector<std::string> &quals);

	void printAllPaths(short label) { 
		printAllPaths(root, 0, label);
		basePaths = occuPaths = "";
	}
private:
	Node* root;
	std::string basePaths;
	std::string occuPaths;

	void createNode(Node *node, short label, char qual);
	void updateWeight(Node *node, char qual);
	short countChildren(Node *node);
	void deleteTree(Node *root);
	void deleteTreeLL(Node **node, std::string readSub);
	void printAllPaths(Node *node, int len, short label);
	void cleanBranches(Node *node);
};

// Structure which holds read data and the four GTrees
typedef struct TreeTop {
	std::string filename;
	ulong nReads;
	int readLength;
	bool ioSuccess;
	std::vector<std::string> reads; 
	std::vector<std::string> quals;
	GTree trees[NBASES];

	TreeTop(const std::string &filename);
	bool readSuccess() { return ioSuccess;  }
	void createRoots();
	short requestBase(ulong readNum, short offset);

	void processReadsOne();
	void processReadsFull();
	void processReadsFullClean();
	void processReadsFullCleanNR();
	void processReadsFullCleanNRBalanced();
} TreeTop;

#endif /*_GTREE_*/
