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

#include <fstream>
#include <sstream>
#include "GTree.H"

// Structure which holds read data and the four GTrees
template <class T>
class TreeTop {
public:
	std::string treeStrings[NBASES];

	TreeTop();

	/* Common */
	void buildSequence();
	void printTrees();
	void storeTrees();
	void printSequence();
private:
	std::string sequence;
protected:
	int64_t nReads;
	int readLength;
	T trees[NBASES];

	short maxPath();
	bool rootOccsExist();
};

// Class for processing reads
class TreeTopfReads : public TreeTop<GTreefReads> {
public:
	TreeTopfReads();
	void processReadsOne();
private:
	void createRoots();
	void threadFunc(uint64_t i);
};

// Class for processing outfiles
class TreeTopfFile : public TreeTop<GTreefFile> {
public:
	explicit TreeTopfFile(const std::string& _iFilename): 
		iFilename(_iFilename) {  }

	void reconstructTrees();
private:
	std::stringstream ss[NBASES];
	std::string iFilename;
};

#endif /*_TREETOP_*/
