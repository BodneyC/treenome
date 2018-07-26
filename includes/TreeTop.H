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

#include <sstream>
#include "GTree.H"

// Structure which holds read data and the four GTrees
template <typename T>
class TreeTop {
public:
	std::string treeStrings[NBASES];

	TreeTop();

	/* From reads */
	void processReadsOne();

	/* From file */
	void reconstructTree(std::stringstream ss[NBASES]);

	/* Common */
	void buildSequence();
	void printTrees();
	void storeTrees();
	void printSequence();
private:
	std::string sequence;
	int64_t nReads;
	int readLength;
	T trees[NBASES];

	void createRoots();
	short maxPath();
	bool rootOccsExist();
	void threadFunc(uint64_t i);
};

#endif /*_TREETOP_*/
