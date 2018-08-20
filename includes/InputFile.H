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

#include "GTree.H"

#define FILE_ERROR -2

// Though sparse, this class may be expanded upon for processing/changing
// of storage method
class InputFile {
public:
	long nReads, readLength;
	int phredBase;

	InputFile( const std::string &_filename, int _pB ):	nReads( 0 ), 
		readLength( 0 ), phredBase( _pB ), filename( _filename ) {  }

	bool readFastQ();

private:
	std::string filename;
};

#endif /*_INP_FILE_*/
