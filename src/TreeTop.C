/********************************************************************
 * Filename: TreeTop.C [C++ source code]
 *
 * Description: Implementation of TreeTop class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Details: Where read processing and the building of the sequence 
 *		takes places; also the analysis printing functions
 *
 *******************************************************************/
#include "../includes/TreeTop.H"

#define MER_LEN 3

int NUM_THREADS;

/** ------------- Sequence Generation -------------- **/
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
	trees[start].followBranch( trees[start].getRoot(), start, sequence );

	return start;
}

void TreeTop::buildSequence()
{
	maxPath();

	uint64_t offset = 1;
	while( rootOccsAboveThresh() ) {

		// Possibly make is a tighter gap as its working from single letters
		// ( this would actually be the k-mer match )
		if( offset >= sequence.length() - MER_LEN ) {
			sequence += 'N';
			maxPath();
			offset += MER_LEN + 1;
		}

		if( trees[BASE_IND( sequence[offset] )].getRoot()->getRatio() >= GTH::thresh )
			trees[BASE_IND( sequence[offset] )].addToSeq( offset, sequence );
		offset++;
	}
}

signed int TreeTop::storeSequence( std::string& oFilename )
{
	std::ofstream storeFile( oFilename );
	if( !storeFile.is_open() )
		return OUT_FILE_ERROR;

	storeFile << sequence;

	storeFile.close();

	return 0;
}

/** ------------ Tree Reconstruction --------------- **/
signed int TreeTop::storeTrees( std::string& sFilename )
{
	std::ofstream storeFile( sFilename, std::ios::binary );
	if( !storeFile.is_open() )
		return OUT_FILE_ERROR;

	for( int i = 0; i < NBASES; i++ ) {
		storeFile.write( ( char* ) &GTH::startOccs[i], sizeof( int64_t ) );
		storeFile.write( ( char* ) &GTH::startWeights[i], sizeof( double ) );
	}

	for( int i = 0; i < NBASES; i++ ) {
		int64_t tmp64 = trees[i].getNNodes();
		storeFile.write( ( char* ) &tmp64, sizeof( int64_t ) );
	}

	for( int i = 0; i < NBASES; i++ )
		trees[i].writeTreeToFile( storeFile );

	storeFile.close();

	return 0;
}

signed int TreeTop::reconstructTrees( std::string& iFilename )
{
	std::ifstream inFile( iFilename, std::ios::binary );
	if( !inFile.is_open() )
		return IN_FILE_ERROR;

	for( int i = 0; i < NBASES; i++ ) {
		inFile.read( ( char* ) &GTH::startOccs[i], sizeof( int64_t ) );
		inFile.read( ( char* ) &GTH::startWeights[i], sizeof( double ) );
	}

	for(int i = 0; i < NBASES; i++) {
		int64_t tmp64;
		inFile.read( ( char* ) &tmp64, sizeof( int64_t ) );
		trees[i].resizeVector( tmp64 );
	}

	for(int i = 0; i < NBASES; i++)
		trees[i].writeVector( inFile );

	inFile.close();
	
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
	for( int i = 0; i < NBASES; i++ )
		trees[i].init( i );

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
char retBase( int index )
{
	char base;
	switch( index ) {
	case 0:
		base = 'A'; 
		break;
	case 1:
		base = 'C';
		break;
	case 2:
		base = 'G';
		break;
	case 3:
		base = 'T';
		break;
	default:
		base = 'N';
	}
	return base;
}

void TreeTop::analyseTrees()
{
	for( int i = 0; i < NBASES; i++ ) {
		nodeCnt[i] = trees[i].getNNodes();
		maxDepth[i] = trees[i].maxDepth( trees[i].getRoot() );
		maxDepthAThresh[i] = trees[i].maxDepthAThresh( trees[i].getRoot() );
	}
}

void TreeTop::printAnalysis()
{
	std::cout << "\n------------------------------------------\nTree Information:\n" << std::endl;

	for( int i = 0; i < NBASES; i++ ) {
		std::cout << "- " << retBase( i ) << "-tree:" << std::endl;
		std::cout << "  Node count                   : " << nodeCnt[i] << std::endl;
		std::cout << "  Maximum depth of tree        : " << maxDepth[i] << std::endl;
		std::cout << "  Max depth above threshold    : " << maxDepthAThresh[i] << '\n' << std::endl;
	}
}

void TreeTop::printTrees()
{
	for( int i = 0; i < NBASES; i++ )
		trees[i].printAllPaths( i );
}

void TreeTop::printSequence()
{
	// 120 for terminal width's sake
	unsigned short TWIDTH = 120;
	uint64_t i = 0;

	if( sequence.length() > TWIDTH )
		for( ; i < sequence.length() - TWIDTH; i += TWIDTH )
			std::cout << sequence.substr( i, TWIDTH ) << std::endl;

	std::cout << sequence.substr( i, sequence.length() - i ) << std::endl;
}

