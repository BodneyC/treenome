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
#define RES 500
#define NUM_THREADS 2

#include <atomic>
#include <mutex>
#include <thread>
#include <iostream>
#include <string>
#include <vector>
#include "SeqRead.H"

// The node type for the GTrees 
typedef struct Node {
	std::atomic<long> occs;
	std::atomic<float> weight;
	short offset;
	long readNum;
	Node* subnodes[NBASES];
	std::mutex nodeMut;

	Node();
	//Copy and assignment functions needed for atomics
	Node(const Node& tmpNode);
	Node& operator=(const Node& tmpNode);
} Node;

// GTree Helpers
namespace GTH {
	extern std::vector<SeqRead> seqReads;

	char retLabel(int label);
	short mostOccs(Node *node);
	short countChildren(Node *node);
	float updateWeight(Node *node, char qual);
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
	std::atomic<Node*> root;
	unsigned long head, nodesCnt;
	std::vector< std::vector<Node> > nodes;
	std::string basePaths;
	std::string occuPaths;
	std::mutex gtMut[2];

	/* Creation/destruction */
	void createNode(Node *node, short label, char qual, long rN, int offset);

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
	void threadFunc(unsigned long i);
	void processReadsOne();
	void buildSequence();
	short maxPath();
	bool rootOccsExist();
	void printSequence();
	void printTrees();
} TreeTop;

#endif /*_GTREE_*/
