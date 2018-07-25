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
#error This code requires C++11 for threading
#endif /*c++11 check*/

#define NBASES 4
#define RES 50
#define NUM_THREADS 4

#include <atomic>
#include <mutex>
#include <omp.h>
#include <iostream>
#include <string>
#include <vector>
#include "SeqRead.H"

// The node type for the GTrees 
typedef struct Node {
	long occs;
	std::atomic<double> weight;
	short offset;
	long readNum;
	Node* subnodes[NBASES];
	omp_lock_t lock;
#ifdef PARENT
	Node* parent;
#endif

	Node();
	//Copy and assignment functions needed for atomics
	Node(const Node& tmpNode);
	Node& operator=(const Node& tmpNode);
} Node;

// GTree Helpers
namespace GTH {
	extern std::vector<SeqRead> seqReads;

	float getNewWeight(float curWeight, long occs, char qual);
	void updateWeight(Node* node, char qual);
	char retLabel(int label);
	short mostOccs(Node *node);
	short countChildren(Node *node);
}

// GTree class used in construction and destruction of tree
class GTree {
public:
	GTree();
	~GTree() {  }

	/* Tree creation functions */
	void createRoot(short ind);
	void addReadOne(long readNum, short offset); 
	bool balanceNode(Node *node, bool mode);

	/* Sequence creation functions */
	void addToSeq(long offset, std::string &sequence);
	void followPath(Node *node, short ind, std::string &sequence);

	/* Usefull non-essential */
	Node* getRoot() { 
		if(root) 
			return root; 
		else 
			return nullptr; 
	}
	void printAllPaths(short label) { 
		printAllPaths(root, 0, label);
		basePaths = occuPaths = "";
	}
	std::string storeTree(short label);
private:
	Node* root;
	std::atomic< unsigned long > head, nodesCnt;
	std::vector< std::vector< Node > > nodes;
	std::string basePaths, occuPaths, treeString;
	omp_lock_t lock;

	/* Creation/destruction */
	void createNode(Node *node, short label, char qual, long rN, int offset);

	/* Public->private overloads */
	void printAllPaths(Node *node, int len, short label);
	void storeTree(Node* node);
};

// Structure which holds read data and the four GTrees
typedef struct TreeTop {
	std::string sequence, treeStrings[NBASES];
	long nReads;
	int readLength;
	GTree trees[NBASES];

	TreeTop();

	void createRoots();
	void threadFunc(unsigned long i);
	void processReadsOne();
	void buildSequence();
	short maxPath();
	bool rootOccsExist();
	void printSequence();
	void printTrees();
	void storeTrees();
} TreeTop;

#endif /*_GTREE_*/
