/********************************************************************
 * Filename: TreeTopfFile.C [C++ source code]
 *
 * Description: Implementation of TreeTop subclass, for building the 
 *		tree from a file outputted from the program
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/TreeTop.H"

/** ------------ Tree Reconstruction --------------- **/
void TreeTopfFile::reconstructTrees()
{
	std::ifstream inFile(iFilename);
	std::stringstream ss[NBASES];
	std::string line;

	int i = 0;
	while(std::getline(inFile, line)) {
		ss[i / 2] << line;
		if(!(i % 2))
			ss[i / 2] << '\n';
		i++;
	}

	int nNodes[NBASES];

	for(int i = 0; i < NBASES; i++) {
		std::getline(ss[i], line, '\n');
		nNodes[i] = std::stol(line);
		std::getline(ss[i], line, '\n');
	}

	for(int i = 0; i < NBASES; i++) {
		for(int j = 0; j < nNodes[i]; j++) {

		}
	}
}

