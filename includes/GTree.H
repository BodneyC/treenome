#ifndef _GTREE_
#define _GTREE_

#if _cplusplus <= 201103L
#define nullptr 0
#endif /*c++11 check*/

#include <iostream>
#include <string>
#include <vector>

#define NBASES 4

typedef unsigned long ulong;

typedef struct Node {
	ulong occurences;
	Node* subnodes[NBASES];

	Node(): occurences(0) 
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
	// Only needed if ..FullClean is used
	void cleanBranches() { cleanBranches(root); }
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
	void processReadsOne();
} TreeTop;

#endif /*_GTREE_*/
