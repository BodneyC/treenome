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

#include <iostream>
#include <string>
#include <vector>
#include "SeqRead.H"
#include "Node.H"

// GTree Helpers
namespace GTH {
	extern std::vector<SeqRead> seqReads;

	float getNewWeight(float curWeight, uint64_t occs, char qual);
	void updateWeight(Node* node, char qual);
	char retLabel(int label);
	signed short mostOccs(Node* node);
	short countChildren(Node* node);
}

// GTree class used in construction and destruction of tree
class GTree {
public:
	GTree();
	~GTree() {  }

	/* Tree creation functions */
	void createRoot(short ind);
	void addReadOne(uint64_t readNum, short offset); 
	bool balanceNode(Node* node, bool mode);

	/* Sequence creation functions */
	void addToSeq(uint64_t offset, std::string &sequence);
	void followPath(Node* node, short ind, std::string &sequence);

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
	std::atomic<uint64_t> head, nodesCnt, nNodes;
	std::vector< std::vector< Node > > nodes;
	std::string basePaths, occuPaths, treeString;
	omp_lock_t lock;

	/* Creation/destruction */
	void createNode(Node* node, short label, char qual, uint64_t rN, int offset);

	/* Public->private overloads */
	void printAllPaths(Node* node, int len, short label);
	void storeTree(Node* node);
};

class GTreefReads : public GTree {
public:

private:

};

class GTreefFile : public GTree {

};

#endif /*_GTREE_*/

