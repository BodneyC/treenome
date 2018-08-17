/********************************************************************
 * Filename: TreeTop.C [C++ source code]
 *
 * Description: Implementation of TreeTop class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright ( c ) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/TreeTop.H"

#define MER_LEN 3

int NUM_THREADS;

/** ------------- Sequence Generation -------------- **/
// Needed?
bool TreeTop::rootOccsExist()
{
	short ret = 0;

	for( int i = 0; i < NBASES; i++ )
		if( trees[i].getRoot()->occs == 0 )
			ret++;

	return ret == NBASES ? 0 : 1;
}

bool TreeTop::rootOccsAboveThresh()
{
	short ret = 0;

	for( int i = 0; i < NBASES; i++ ) {
		if( GTH::startOccs[i] ) {
			if( GTH::startWeights[i] / static_cast<double>(GTH::startOccs[i]) < GTH::thresh) {
				ret++;
			}
		} else if( !GTH::startOccs[i] ) {
			ret++;
		}
	}

	return ret == NBASES ? 0 : 1;
}

short TreeTop::maxPath()
{
	int start = 0;
	double maxRat = GTH::thresh;

	for( short i = 0; i < NBASES; i++ ) {
		if( GTH::startOccs[i] ) {
			double curRat = GTH::startWeights[i] / static_cast<double>(GTH::startOccs[i]);
			if( curRat >= maxRat ) {
				maxRat = curRat;
				start = i;
			}
		}
	}

	GTH::startWeights[start]--;
	GTH::startOccs[start]--;
	trees[start].followPath( trees[start].getRoot(), start, sequence );

	return start;
}

void TreeTop::buildSequence()
{
	maxPath();

	uint64_t offset = 1;
//#pragma omp parallel num_threads( NUM_THREADS ) shared( offset )
//{
	while( rootOccsAboveThresh() ) {

		// Possibly make is a tighter gap as its working from single letters
		// ( this would actually be the k-mer match )
//#pragma omp single
//{
		if( offset >= sequence.length() - MER_LEN ) {
			sequence += 'N';
			maxPath();
			offset += MER_LEN + 1;
		}

		if( trees[BASE_IND( sequence[offset] )].getRoot()->getRatio() >= GTH::thresh )
			trees[BASE_IND( sequence[offset] )].addToSeq( offset, sequence );
		offset++;
//}
	}
//}
}

/** ------------ Tree Reconstruction --------------- **/
signed int TreeTop::storeTrees( std::string& sFilename )
{
	std::ofstream storeFile( sFilename, std::ios::binary );
	if( !storeFile.is_open() )
		return OUT_FILE_ERROR;

	for( int i = 0; i < NBASES; i++ ) {
		int64_t tmp64 = trees[i].getDNodeSize();
		storeFile.write( ( char* ) &tmp64, sizeof( int64_t ) );
	}

	for( int i = 0; i < NBASES; i++ )
		trees[i].writeTreeToFile( storeFile );

	//for( int i = 0; i < NBASES; i++ )
	//	treeStrings[i] = trees[i].storeTree( i );
	
	storeFile.close();

	return 0;
}

signed int TreeTop::reconstructTrees( std::string& iFilename )
{
	std::ifstream inFile( iFilename, std::ios::binary );
	if( !inFile.is_open() )
		return IN_FILE_ERROR;

	for(int i = 0; i < NBASES; i++) {
		int64_t tmp64;
		inFile.read( ( char* ) &tmp64, sizeof( int64_t ) );
		std::cout << tmp64 << std::endl;
		trees[i].resizeVector( tmp64 );
	}

	for(int i = 0; i < NBASES; i++) {
		trees[i].writeVector( inFile );
	}

	//std::ifstream inFile( iFilename );
	//std::string line;

	//int i = 0;
	//while( std::getline( inFile, line ) ) {
	//	ss[i / 2] << line;
	//	if( !( i % 2 ) )
	//		ss[i / 2] << '\n';
	//	i++;
	//}

	//for( int i = 0; i < NBASES; i++ ) {
	//	trees[i].processSString( ss[i] );
	//}
	
	return 0;
}

/** --------------- Read Processing ---------------- **/
void TreeTop::threadFunc( uint64_t i )
{
	for( short j = 0; j < GTH::seqReads[i].size(); j++ )
		//if( GTH::seqReads[i].getBaseInd( j ) == 1 )
			trees[GTH::seqReads[i].getBaseInd( j )].addReadOne( i, j );
}

void TreeTop::processReadsOne()
{
	for( int i = 0; i < NBASES; i++ ) {
		trees[i].init();
		trees[i].createRoot( i );
	}

#pragma omp parallel num_threads( NUM_THREADS )
{
	for( uint64_t i = 0; i < GTH::seqReads.size(); i += NUM_THREADS ) {
		#pragma omp for schedule( static, 1 )
		for( int j = 0; j < NUM_THREADS; j++ ) 
			if( i + j < GTH::seqReads.size() )
				threadFunc( i + j );
	}
}
}

/** --------------- Misc Functions ----------------- **/
void TreeTop::printTrees()
{
	for( int i = 0; i < NBASES; i++ )
		trees[i].printAllPaths( i );
}

void TreeTop::printSequence()
{
	// 80 for terminal width's sake
	unsigned short TWIDTH = 120;
	uint64_t i = 0;

	if( sequence.length() > TWIDTH )
		for( ; i < sequence.length() - TWIDTH; i += TWIDTH )
			std::cout << sequence.substr( i, TWIDTH ) << std::endl;

	std::cout << sequence.substr( i, sequence.length() - i ) << std::endl;
}

