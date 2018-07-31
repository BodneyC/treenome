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
 * Copyright ( c ) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#ifndef _GTREE_
#define _GTREE_

////TMP
#include <mutex>

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <limits>
#include "SeqRead.H"
#include "Node.H"
#include "Phred33.H"

// GTree Helpers
namespace GTH {
	extern std::vector<SeqRead> seqReads;

	char retLabel( int label );
	float getNewWeight( double curWeight, int64_t occs, char qual );
	void removeDoubleEnding( std::string& doubleString );
	template <typename T>
	std::string valToString( T& val );
}

// GTree class used in construction and destruction of tree
template <typename T>
class GTree {
public:
	GTree();
	~GTree() {  }

	/* Sequence creation functions */
	void addToSeq( uint64_t offset, std::string &sequence );
	void followPath( T* node, short ind, std::string &sequence );

	/* Usefull non-essential */
	T* getRoot() { 
		if( root ) 
			return root; 
		else 
			return nullptr; 
	}
	void printAllPaths( short label ) { 
		printAllPaths( root, 0, label );
		basePaths = occuPaths = "";
	}
	std::string storeTree( short label );
protected:
	T* root;
	std::atomic<uint64_t> head, nNodes;
	std::string basePaths, occuPaths, treeString;

	/* Helpers requiring templates */
	signed short mostOccs( T* node );
	short countChildren( T* node );

	/* Public->private overloads */
	void printAllPaths( T* node, int len, short label );
	void storeTree( T* node );
};

class GTreefReads : public GTree<Node> {
public:
	GTreefReads();

	/* Tree creation functions */
	void createRoot( short ind );
	void addReadOne( uint64_t readNum, short offset ); 
	void balanceNode( Node* node );
private:
	std::vector< std::vector< Node > > nodes;
	void createNode( Node* node, short label, char qual, uint64_t rN, int offset );
	std::atomic<uint64_t> nodesCnt;
	omp_lock_t lock;
	////TMP
	std::mutex gtMut;

	void potAddNode( Node* node );
	void updateWeight( Node* node, char qual );
	void updateWeightAndOccs( Node* node, char qual );
};

class GTreefFile : public GTree<pNode> {
public:
	GTreefFile(): tmpNode( nullptr ) {  }

	void processSString( std::stringstream& ss );

private:
	pNode* tmpNode;
	std::vector< pNode > pNodes;

	void createRoot( std::stringstream& ss );
	void createNode( struct NodeInfo& nInf );
	void getNextNode( struct NodeInfo& nInf, std::stringstream& ss );

};

#endif /*_GTREE_*/

