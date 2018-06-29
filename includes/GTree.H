#ifndef _GTREE_
#define _GTREE_

#if __cplusplus < 201103L
#define nullptr 0
#error fish
#endif /*c++11 check*/

#define NBASES 4

#include <iostream>
#include <string>
#include <vector>

typedef unsigned long ulong;
typedef unsigned int uint;

typedef struct Node {
	ulong occurences;
	Node* subnodes[NBASES];

// The ifdeffing of ...ReadOne() is because it adds ten bytes per node
#ifdef READSONE
	short offset;
	ulong readNumber;

	Node(): occurences(0), offset(0), readNumber(0)
#else
	Node(): occurences(0) 
#endif /*READSONE*/
	{
		for(int i = 0; i < NBASES; i++)
			subnodes[i] = nullptr;
	}
} Node;

// GTree class used in construction and destruction of tree
class GTree {
public:
	GTree() { root = new Node; }
	GTree(const GTree &orig) { root = new Node; }
	~GTree() { deleteTree(root); }

	Node *getRoot() { return root; }
	void addReadFull(std::string &read);
	void addReadOne();
	// cleanBranches..() only needed if ..FullClean is used
	void cleanBranches() { cleanBranches(root); }
	void cleanBranchesNR(std::string &read);
	void deleteTree() { deleteTree(root); }
	void printAllPaths(short label) { 
		printAllPaths(root, 0, label);
		basePaths = occuPaths = "";
	}
private:
	Node* root;
	std::string basePaths;
	std::string occuPaths;

	void deleteTree(Node* root);
	void deleteTreeNR(Node **node, std::string readSub);
	void printAllPaths(Node *node, int len, short label);
	void cleanBranches(Node *node);
};

// Structure which holds read data and the four GTrees
typedef struct TreeTop {
	long nReads, readLength;
	std::vector<std::string> reads; 
	std::vector<std::string> qualities;
	GTree trees[NBASES];

	TreeTop(): nReads(0), readLength(0) {}	
	void processReadsFull();
	void processReadsFullClean();
	void processReadsFullCleanNR();
	void processReadsOne();
} TreeTop;

#endif /*_GTREE_*/
