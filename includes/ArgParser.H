/********************************************************************
 * Filename: ArgParser.H [C++ header code]
 *
 * Description: Simple class for command line argument parsing
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright ( c ) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#ifndef _ARGPARSER_
#define _ARGPARSER_

#include <string>
#include <vector>
#include <fstream>

#define USAGE_ERROR -1
#define IN_FILE_ERROR -2
#define OUT_FILE_ERROR -3

struct CMDArgs {
	std::string iFilename, oFilename, lFilename;
	bool printToScreen, storeToFile, loadFile;
	int phredBase;

	CMDArgs(): 
		iFilename( "" ), oFilename( "" ), lFilename( "" ),
		printToScreen( 0 ), storeToFile( 0 ), loadFile( 0 ),
		phredBase( 33 ) {  }
};

class ArgParser {
public:
	ArgParser( int& argc, char** argv );
	int fillArgs( struct CMDArgs &argList );
	bool argExists( const std::string &arg );
	const std::string stringOption( const std::string &arg );
	int inFileCheck( std::string filename );
	int outFileCheck( std::string filename );
private:
	int argc;
	std::vector<std::string> args;
};

#endif /*_ARGPARSER_*/
