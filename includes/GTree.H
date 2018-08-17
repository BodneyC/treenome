/********************************************************************
 * Filename: GTree.H [C++ header code]
 *
 * Description: Declaration of GTree class and NodeInfo struct
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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <deque>
#include <limits>
#include "SeqRead.H"
#include "Node.H"
#include "Phred33.H"

#define USAGE_ERROR -1
#define IN_FILE_ERROR -2
#define OUT_FILE_ERROR -3
#define THREAD_ERROR -4

#define RES 1000

// GTree Helpers
namespace GTH {
	extern std::vector<SeqRead> seqReads;
	extern double thresh;
	extern double startWeights[NBASES];
	extern int64_t startOccs[NBASES];

	char retLabel( int label );
	float getNewWeight( double curWeight, int64_t occs, char qual );
	void removeDoubleEnding( std::string& doubleString );
	template <typename T>
	std::string valToString( T& val );
}

struct NodeInfo {
	int ind, comCnt;
	int64_t occs;
	double weight;

	NodeInfo(): ind( 0 ), comCnt( 0 ), occs( 0 ), weight( 0 ) {  }
};

// GTree class used in construction and destruction of tree
class GTree {
public:
	GTree(): root( nullptr ), head( 0 ), nNodes( 0 ), 
		basePaths( "" ), occuPaths( "" ), treeString( "" ),
		tmpNode( nullptr ) {  }
	~GTree() {  }

	/* Sequence creation functions */
	void addToSeq( uint64_t offset, std::string& sequence );
	void followPath( Node* node, short ind, std::string& sequence );
	/* Usefull non-essential */
	Node* getRoot();
	int64_t getDNodeSize();
	void printAllPaths( short label );
	std::string storeTree( short label );
	void writeTreeToFile( std::ofstream& storeFile );
	/* Tree creation functions */
	void init();
	void createRoot( short ind );
	void addReadOne( uint64_t readNum, short offset ); 
	void balanceNode( Node* node );
	/* Tree Recreation Functions */
	void resizeVector( int64_t tmp64 );
	void writeVector( std::ifstream& inFile );
	void processSString( std::stringstream& ss );
private:
	Node* root;
	std::atomic<int64_t> head, nNodes;
	std::string basePaths, occuPaths, treeString;
	omp_lock_t lock;
	std::deque< Node > dNodes;
	// TMP
	Node* tmpNode;

	/* Helper functions */
	signed short highestThresh( Node* node );
	short countChildren( Node* node );
	/* Tree creation functions */
	void potAddNode( Node* node );
	void updateWeight( Node* node, char qual );
	void updateWeightAndOccs( Node* node, char qual );
	void createNode( Node* node, short label, char qual, uint64_t rN, int offset );
	/* Tree Recreation Functions */
	void createRoot( std::stringstream& ss );
	void createNode( struct NodeInfo& nInf );
	void getNextNode( struct NodeInfo& nInf, std::stringstream& ss );
	/* Public->private overloads */
	void printAllPaths( Node* node, int len, short label );
	void storeTree( Node* node );
};

#endif /*_GTREE_*/

