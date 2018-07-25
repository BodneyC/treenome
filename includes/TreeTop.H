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
typedef struct TreeTop {
	std::string sequence, treeStrings[NBASES];
	int64_t nReads;
	int readLength;
	GTree trees[NBASES];

	TreeTop();

	void createRoots();
	void threadFunc(uint64_t i);
	void processReadsOne();
	void buildSequence();
	short maxPath();
	bool rootOccsExist();
	void printSequence();
	void printTrees();
	void storeTrees();
} TreeTop;

#endif /*_TREETOP_*/
