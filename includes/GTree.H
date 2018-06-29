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
	float weight;
	short offset;
	ulong readNum;
	Node* subnodes[NBASES];

	Node(): occurences(0), weight(0), offset(0), readNum(0)
	{
		for(int i = 0; i < NBASES; i++)
			subnodes[i] = nullptr;
	}
} Node;

// GTree class used in construction and destruction of tree
class GTree {
public:
	~GTree() { deleteTree(root); }

	void createRoot(std::vector<std::string> &reads, int label);
	void addReadFull(ulong readNum, short offset, std::string &read, std::string &qual);
	void addReadOne(std::string &read, std::string &qual);
	// cleanBranches..() only needed if ..FullClean is used
	void cleanBranches() { cleanBranches(root); }
	void cleanBranchesNR(std::string &read);
	void deleteTree() { deleteTree(root); }

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

	void deleteTree(Node *root);
	void deleteTreeNR(Node **node, std::string readSub);
	void printAllPaths(Node *node, int len, short label);
	void cleanBranches(Node *node);
};

// Structure which holds read data and the four GTrees
typedef struct TreeTop {
	ulong nReads;
	int readLength;
	std::vector<std::string> reads; 
	std::vector<std::string> qualities;
	GTree trees[NBASES];
	bool ioSuccess;

	TreeTop(std::string filename);
	bool readSuccess() { return ioSuccess;  }
	void createRoots();
	void processReadsFull();
	void processReadsFullClean();
	void processReadsFullCleanNR();
	void processReadsFullCleanNRBalanced();
	void processReadsOne();
} TreeTop;

#endif /*_GTREE_*/
