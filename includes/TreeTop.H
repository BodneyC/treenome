/********************************************************************
 * Filename: TreeTop.H [C++ header code]
 *
 * Description: Header for forest class containing GTrees
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#ifndef _TREETOP_
#define _TREETOP_

#include "GTree.H"

// Structure which holds read data and the four GTrees
class TreeTop {
public:
	std::string treeStrings[NBASES];

	TreeTop(): sequence( "" ), nReads( 0 ), readLength( 0 ) {  }

	/* Common */
	void buildSequence();
	void printTrees();
	signed int storeTrees( std::string& sFilename );
	void printSequence();
	signed int storeSequence( std::string& oFilename );
	void analyseTrees();
	void printAnalysis();
	/* FReads */
	void processReadsOne();
	/* FFile */
	signed int reconstructTrees( std::string& iFilename );
	
private:
	std::string sequence;
	int64_t nReads;
	int readLength;
	GTree trees[NBASES];
	int32_t nodeCnt[NBASES], maxDepth[NBASES], maxDepthAThresh[NBASES];

	short maxPath();
	bool rootOccsAboveThresh();
	bool rootOccsExist();
	/* FReads */
	void createRoots();
	void threadFunc( uint64_t i );
};

#endif /*_TREETOP_*/
