/********************************************************************
 * Filename: InputFile.H [C++ header code]
 *
 * Description: Declaration of InputFile class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#ifndef _INP_FILE_
#define _INP_FILE_

#include <fstream>
#include <vector>
#include <string>
#include "GTree.H"

#define FILE_ERROR -2

// Though sparse, this class may be expanded upon for processing/changing
// of storage mnethod
class InputFile {
public:
	long nReads, readLength;

	explicit InputFile(const std::string &_filename):
		nReads(0), readLength(0), filename(_filename) {}

	bool readFastQ();

private:
	std::string filename;
};

#endif /*_INP_FILE_*/
