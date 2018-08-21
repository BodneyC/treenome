/********************************************************************
 * Filename: TreeNome.C [C++ source code]
 *
 * Description: Entry point for TreeNome de novo DNA assembler
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
/* TODO:
 * - Analysis of trees
 */
#include "../includes/tclap/CmdLine.h"
#include "../includes/InputFile.H"
#include "../includes/TreeTop.H"
#include "sys/types.h"
#include "sys/sysinfo.h"

#define OUT

typedef struct CMDArgs {
	std::string fFilename, stFilename, lFilename, ssFilename;
	bool printToScreen, storeToFile, loadFile, analyse;
	int phredBase;

	CMDArgs(): 
		fFilename( "" ), stFilename( "" ), lFilename( "" ), ssFilename( "" ),
		printToScreen( 0 ), storeToFile( 0 ), loadFile( 0 ), analyse( 0 ),
		phredBase( 33 ) {  }
} CMDArgs;

signed int inFileCheck( std::string filename )
{
	std::ifstream testFile( filename );
	if( testFile.good() )
		return 0;
	else
		return IN_FILE_ERROR;
}

signed int returnArgs( int argc, char** argv, CMDArgs& argList ) 
{
	try {
		TCLAP::CmdLine cmd( "Tree based de novo DNA assembler", ' ', "1.04" );
		TCLAP::ValueArg<std::string> fFileArg( "f", "fastqfile", "Input file in fastq format", false, "", "string" );
		TCLAP::ValueArg<std::string> lFileArg( "l", "load-file", 
				"Input file which was previously outputted from TreeNome", false, "", "string" );
		TCLAP::ValueArg<std::string> stFileArg( "", "store-tree", 
				"File in which to store the tree data", false, "", "string" );
		TCLAP::ValueArg<std::string> ssFileArg( "s", "store-sequence", 
				"File in which to store the sequence data", true, "", "string" );
		TCLAP::ValueArg<int> threadArg( "t", "threads", "Number of threads to use", false, 1, "int" );
		TCLAP::ValueArg<double> threshArg( "", "thresh", "Threashold value", false, 0.3f, "0 to 1" );
		std::vector<std::string> pAllowed = { "Phred+33", "Phred+64", "Solexa+64" };
		TCLAP::ValuesConstraint<std::string> pAllowedVC( pAllowed );
		TCLAP::ValueArg<std::string> scoreArg( "", "score-sys", 
				"Scoring system. Phred+33 default (Sanger)", false, "Phred+33", &pAllowedVC );
		TCLAP::SwitchArg oSwitch( "o", "stdout", "Print to stdout", cmd, 0 );
		TCLAP::SwitchArg aSwitch( "a", "analyse", "Perform tree analysis", cmd, 0 );

		cmd.xorAdd( fFileArg, lFileArg );
		cmd.add( stFileArg );
		cmd.add( ssFileArg );
		cmd.add( threadArg );
		cmd.add( threshArg );
		cmd.add( scoreArg );
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

		argList.stFilename = stFileArg.getValue();
		if( stFileArg.isSet() )
			argList.storeToFile = 1;

		if( ssFileArg.isSet() ) {
			argList.ssFilename = ssFileArg.getValue();
			if( argList.ssFilename.substr( argList.ssFilename.find_last_of( "." ) + 1 ) != "gno" )
				argList.ssFilename += ".gno";
		}

		NUM_THREADS = threadArg.getValue();
		if( NUM_THREADS > omp_get_max_threads() )
			return THREAD_ERROR;
		
		GTH::scoreSys = 0;
		if( scoreArg.getValue() == "Phred+33" )
			argList.phredBase = 33;
		else if( scoreArg.getValue() == "Phred+64" )
			argList.phredBase = 64;
		else if( scoreArg.getValue() == "Solexa+64" ) {
			argList.phredBase = 59;
			GTH::scoreSys = 1;
		}

		argList.printToScreen = oSwitch.getValue();
		argList.analyse = aSwitch.getValue();
		GTH::thresh = threshArg.getValue();
	} catch ( TCLAP::ArgException &e ) {
		std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	return 0;
}

void writeTreesToDisk( std::string stFilename, TreeTop* treeTop )
{
	std::ofstream outFile( stFilename );
	for( int i = 0; i < NBASES; i++ )
		outFile << treeTop->treeStrings[i].c_str() << std::endl;
}

TreeTop* createTreeFromReads( CMDArgs& argList, OUT double& t2 )
{
	InputFile inpFile( argList.fFilename, argList.phredBase );
	if( !inpFile.readFastQ() ) {
		std::cout << "[ERR]: Input file could not be opened" << std::endl;
		return nullptr;
	} else {
		std::cout << inpFile.nReads << " records found\n" << std::endl;
	}

	TreeTop* treeTop = new TreeTop;

	double t1 = omp_get_wtime();
	treeTop->processReadsOne();
	t2 = omp_get_wtime() - t1;

	return treeTop;
}

TreeTop* loadTreeFromFile( CMDArgs& argList )
{
	TreeTop* treeTop = new TreeTop;
	treeTop->reconstructTrees( argList.lFilename );

	return treeTop;
}

int64_t parseLine( std::string& line ) {
	size_t posOfK = line.find_last_of( 'k' );
	if( posOfK == std::string::npos )
		return -1;

    int j = 0;
	while ( line[j] < '0' || line[j] > '9' ) 
		j++;

	std::string val = line.substr( j, posOfK - j );

	int64_t retVal = std::stoll( val, nullptr, 10 );

    return retVal;
}

int64_t getMemInUse() { 
	std::ifstream self( "/proc/self/status" );
	if( !self.is_open() )
		return IN_FILE_ERROR;

    int result = -1;
	std::string line;

    while( getline( self, line, '\n' ) ) {
        if( line.substr( 0, 7 ) == "VmSize:" ) {
            result = parseLine( line );
            break;
        }
    }

    self.close();
    return result;
}

// https://stackoverflow.com/questions/7276826/c-format-number-with-commas
template< typename T>
std::string valWCommas( T val )
{
	std::string st = std::to_string( val );
	int i = st.find_last_of( "." );
	if( i == -1 ) {
		i = st.length() - 3;
	} else {
		i -= 3;
		st.erase( st.find_last_not_of( '0' ) + 1, std::string::npos );
		if( st[st.length() - 1] == '.' )
			st = st.substr( 0, st.length() - 1 );
	}

	while( i > 0 ) {
		st.insert( i, "," );
		i -=3;
	}

	return st;
}

signed int analysis( TreeTop* treeTop, double timeToConstructTrees, double timeToBuildSeq )
{
	std::cout << "------------------------------------------\nTiming Information:\n" << std::endl;
	std::cout << "Number of threads in use       : " << NUM_THREADS << std::endl;
	std::cout << "Time to construct trees  ( s ) : " << timeToConstructTrees << std::endl;
	std::cout << "Time to build sequence   ( s ) : " << timeToBuildSeq << '\n' << std::endl;

	struct sysinfo memInf;
	sysinfo ( &memInf );

	int64_t totVirtMem = memInf.totalram;
	totVirtMem *= memInf.mem_unit;
	double totVirtMemKB = static_cast<double>( totVirtMem ) / 1024.0;

	int64_t vMemInUse = getMemInUse();
	if( vMemInUse == -1 )
		return vMemInUse;

	std::cout << "------------------------------------------\nMemory Information:\n" << std::endl;
	std::cout << "Virtual memory in use   ( KB ) : " << valWCommas( vMemInUse ) << std::endl;
	std::cout << "Total available         ( KB ) : " << valWCommas( totVirtMemKB ) << std::endl;
	std::cout << "Percent used             ( % ) : " << static_cast<double>( vMemInUse ) / totVirtMemKB << std::endl;

	return 0;
}

int main( int argc, char** argv )
{
	CMDArgs argList;
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

	double timeToConstructTrees;

	TreeTop* treeTop;
	if( argList.loadFile ) {
		treeTop = loadTreeFromFile( argList );
	} else {
		treeTop = createTreeFromReads( argList, timeToConstructTrees );
	}
	if( !treeTop )
		return FILE_ERROR;

	if( argList.storeToFile ) {
		progFail = treeTop->storeTrees( argList.stFilename );
		if( progFail )
			return progFail;
	}

	// Has to be before buildSequence()
	if( argList.printToScreen )
		treeTop->printTrees();

	// Has to be before buildSequence()
	treeTop->analyseTrees();

	double t1 = omp_get_wtime();
	treeTop->buildSequence();
	double timeToBuildSeq = omp_get_wtime() - t1;

	progFail = treeTop->storeSequence( argList.ssFilename );
	if( progFail )
		return progFail;

	if( argList.printToScreen ) {
		treeTop->printSequence();
	}

	if( argList.analyse ) {
		analysis( treeTop, timeToConstructTrees, timeToBuildSeq );
		treeTop->printAnalysis();
	}

	delete treeTop;

	return 0;
}
