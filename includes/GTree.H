#ifndef _GTREE_
#define _GTREE_

#include <iostream>
#include <string>
#include <vector>

typedef unsigned long ulong;

typedef struct Node {
	ulong occurences;
	std::vector<Node*> subnodes;

	Node(): occurences(0), subnodes(4) 
	{
		for(int i = 0; i < 4; i++)
			subnodes[i] = NULL;
	}
} Node;

class GTree {
public:
	GTree() { root = new Node; }
	GTree(const GTree &orig) { root = new Node; }
	~GTree() { deleteTree(root); }

	void addReadFULL(std::string &read);
	void addReadONE();
	void printAllPaths(short label) { printAllPaths(root, 0, label); }
private:
	Node* root;
	std::string basePaths;
	std::string occuPaths;

	void deleteTree(Node* root);
	void printAllPaths(Node *node, int len, short label);
};

typedef struct TreeTop {
	long nReads, readLength;
	std::vector<std::string> reads; 
	std::vector<std::string> qualities;
	std::vector<GTree> trees;

	TreeTop(): nReads(0), readLength(0), trees(4) {}	
} TreeTop;

#endif /*_GTREE_*/
