/********************************************************************
 * Filename: GTree.H [C++ header code]
 *
 * Description: Declaration of TreeTop, GTree, and Node classes/
 *				structs
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#ifndef _GTREE_
#define _GTREE_

#if __cplusplus < 201103L
#define nullptr 0
#endif /*c++11 check*/

#define NBASES 4
#define RES 500

#include <iostream>
#include <string>
#include <vector>
#include "SeqRead.H"

// The node type for the GTrees 
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
	GTree();
	~GTree() {  }

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
	unsigned long head, nodesCnt;
	std::vector< std::vector<Node> > nodes;
	std::string basePaths;
	std::string occuPaths;

	/* Creation/destruction */
	void createNode(Node *node, short label, char qual);

	/* Public->private overloads */
	void printAllPaths(Node *node, int len, short label);
};

// Structure which holds read data and the four GTrees
typedef struct TreeTop {
	std::string sequence;
	long nReads;
	int readLength;
	GTree trees[NBASES];

	TreeTop();

	void createRoots();
	void processReadsOne();
	void buildSequence();
	short maxPath();
	bool rootOccsExist();
	void printSequence();
	void printTrees();
} TreeTop;

#endif /*_GTREE_*/
