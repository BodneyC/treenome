/********************************************************************
 * Filename: TreeNome.C [C++ source code]
 *
 * Description: Entry point for TreeNome de novo DNA assembler
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (xcx ) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
/* TODO:
 * - Double accuracy?
 * - Reconstruct tree from file ( possible second program )
 * - Quality vs. occurrences in sequence creation
 */
#include "../includes/tclap/CmdLine.h"
#include "../includes/InputFile.H"
#include "../includes/TreeTop.H"

#define USAGE_ERROR -1
#define IN_FILE_ERROR -2
#define THREAD_ERROR -3

struct CMDArgs {
	std::string fFilename, sFilename, lFilename;
	bool printToScreen, storeToFile, loadFile;
	int phredBase;

	CMDArgs(): 
		fFilename( "" ), sFilename( "" ), lFilename( "" ),
		printToScreen( 0 ), storeToFile( 0 ), loadFile( 0 ),
		phredBase( 33 ) {  }
};

signed int inFileCheck( std::string filename )
{
	std::ifstream testFile( filename );
	if( testFile.good() )
		return 0;
	else
		return IN_FILE_ERROR;
}

signed int returnArgs( int argc, char** argv, struct CMDArgs& argList ) {
	try {
		TCLAP::CmdLine cmd( "Tree based de novo DNA assembler", ' ', "1.04" );
		TCLAP::ValueArg<std::string> fFileArg( "f", "fastqfile", "Input file in fastq format", false, "", "string" );
		TCLAP::ValueArg<std::string> lFileArg( "l", "trenomefile", 
				"Input file which was previously outputted from TreeNome", false, "", "string" );
		TCLAP::ValueArg<std::string> sFileArg( "s", "storefile", 
				"File in which to store the tree data", false, "", "string" );
		TCLAP::ValueArg<int> thrArg( "t", "threads", "Number of threads to use", false, /*1*/ 8, "int" );
		std::vector<int> pAllowed = { 33, 64 };
		TCLAP::ValuesConstraint<int> pAllowedVC( pAllowed );
		TCLAP::ValueArg<int> phredArg( "", "phred", 
				"Phred base of qualities in fastq", false, 33, &pAllowedVC );
		TCLAP::SwitchArg oSwitch( "o", "stdout", "Print to stdout", cmd, 0 );

		cmd.xorAdd( fFileArg, lFileArg );
		cmd.add( sFileArg );
		cmd.add( thrArg );
		cmd.add( phredArg );
		cmd.parse( argc, argv );

		if( !fFileArg.isSet() && !lFileArg.isSet() )
			return USAGE_ERROR;
		argList.fFilename = fFileArg.getValue();
		if( fFileArg.isSet() && inFileCheck( argList.fFilename ) )
				return IN_FILE_ERROR;
		argList.lFilename = lFileArg.getValue();
		if( lFileArg.isSet() ) {
			argList.loadFile = 1;
			if( inFileCheck( argList.lFilename ) )
				return IN_FILE_ERROR;
		}
		argList.sFilename = sFileArg.getValue();
		if( sFileArg.isSet() )
			argList.storeToFile = 1;
		NUM_THREADS = thrArg.getValue();
		if( NUM_THREADS > omp_get_max_threads() )
			return THREAD_ERROR;
		argList.phredBase = phredArg.getValue();
		argList.printToScreen = oSwitch.getValue();
	} catch (TCLAP::ArgException &e) {
		std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	return 0;
}

void writeTreesToDisk( std::string sFilename, TreeTop<GTreefReads>& treeTop )
{
	std::ofstream outFile( sFilename );
	for( int i = 0; i < NBASES; i++ )
		outFile << treeTop.treeStrings[i].c_str() << std::endl;
}

signed int createTreeFromReads( struct CMDArgs& argList )
{
	InputFile inpFile( argList.fFilename, argList.phredBase );
	if( !inpFile.readFastQ() ) {
		std::cout << "[ERR]: Input file could not be opened" << std::endl;
		return FILE_ERROR;
	} else {
		std::cout << inpFile.nReads << " records found\n" << std::endl;
	}

	TreeTopfReads treeTop;
	treeTop.processReadsOne();
	std::cout << "\n----------" << std::endl;
	if( argList.printToScreen )
		treeTop.printTrees();
	if( argList.storeToFile ) {
		treeTop.storeTrees();
		writeTreesToDisk( argList.sFilename, treeTop );
	}

	//treeTop.buildSequence();
	//treeTop.printSequence();

	return 0;
}

signed int loadTreeFromFile( struct CMDArgs& argList )
{
	TreeTopfFile treeTop( argList.lFilename );
	treeTop.reconstructTrees();
	if( argList.printToScreen )
		treeTop.printTrees();

	return 0;
}

int main( int argc, char** argv )
{
	struct CMDArgs argList;
	signed int progFail = returnArgs( argc, argv, argList );

	switch ( progFail ) {
	case USAGE_ERROR:
		std::cout << "[ERR]: Incorrect command usage" << std::endl;
		break;
	case IN_FILE_ERROR:
		std::cout << "[ERR]: Input file not supplied or not present" << std::endl;
		break;
	case THREAD_ERROR:
		std::cout << "[ERR]: " << NUM_THREADS << " threads requested, " << omp_get_max_threads() << " available."<< std::endl;
		break;
	}
	if( progFail ) {
		return progFail;
	}

	omp_set_num_threads( NUM_THREADS );

	std::cout << NUM_THREADS << std::endl;

	if( argList.loadFile ) {
		progFail = loadTreeFromFile( argList );
	} else {
		progFail = createTreeFromReads( argList );
	}
	if( progFail )
		return progFail;

	return 0;
}
