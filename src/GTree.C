/********************************************************************
 * Filename: GTree.C [C++ source code]
 *
 * Description: Implementation of GTree class, templatized for 
 *		Node/pNode purposes
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/GTree.H"

#ifdef __MINGW64__
#include <sstream>
namespace mingw_fix {
	{
		std::ostringstream ss;
		ss << val ;
		return ss.str() ;
	}
}
#endif /*__MINGW32__*/

/** --------------- Helper Functions --------------- **/
namespace GTH {
	std::vector<SeqRead> seqReads;
	double thresh;
	double startWeights[NBASES];
	int64_t startOccs[NBASES];
	int scoreSys;

	char retLabel( int label )
	{
		switch( label ) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'T';
		case 3:
			return 'G';
		default: // 7 but it shouldn't
			return 'P';
		}
	}

	void removeDoubleEnding( std::string& doubleString )
	{
		doubleString.erase( doubleString.find_last_not_of('0') + 1, std::string::npos );
		doubleString.erase( doubleString.find_last_not_of('.') + 1, std::string::npos );
	}
	
	template <typename T>
	std::string valToString( T& val )
	{
#ifdef __MINGW64__
		return mingw_fix::to_string( val );
#else
		return std::to_string( val );
#endif /*__MINGW64__*/
	}
}

/** ------------- GTree Cons and Dees -------------- **/
GTree::GTree(): 
	root( nullptr ), head( 0 ), nNodes( 0 ), 
	basePaths( "" ), occuPaths( "" ), treeString( "" ),
	tmpNode( nullptr )
{ 
}

/** ------------- - Helper Functions --------------- **/
Node* GTree::getRoot()
{
	if( root ) 
		return root; 
	else 
		return nullptr; 
}

short GTree::countChildren( Node* node )
{
	short children = 0;

	for( int i = 0; i < NBASES; i++ )
		if( node->subnodes[i] )
			children++;

	return children;
}

signed short GTree::highestThresh( Node* node )
{
	short ind = -1;
	//double maxRat = std::numeric_limits<double>::lowest();
	double maxRat = GTH::thresh;

	for( short i = 0; i < NBASES; i++ ) {
		Node* tmpNode = node->subnodes[i];
		if( tmpNode && tmpNode->occs ) {
			double curRat = tmpNode->getRatio();
			if( curRat > maxRat ) {
				maxRat = curRat;
				ind = i;
			}
		}
	}

	return ind;
}

void GTree::followPath( Node* node, short ind, std::string &sequence )
{
	short children;
	do {
		children = countChildren( node );
		sequence += GTH::retLabel( ind );
		ind = highestThresh( node );
		node->occs--;
		node->weight = node->weight.load() - 1;
		if( ind == -1 )
			return;
		node = node->subnodes[ind];
	} while( children );
}

/** --------------- Genome Creation ---------------- **/
void GTree::addToSeq( uint64_t offset, std::string &sequence )
{
	Node* node = root;
	short ind = 0;
	std::vector< Node* > nPath;

	// End of current sequence	
	for( uint32_t i = offset + 1; i < sequence.length(); i++ ) {
		ind = BASE_IND( sequence[i] );
		if( node->subnodes[ind] && countChildren( node->subnodes[ind] ) ) {
			nPath.push_back( node );
			node = node->subnodes[ind];
		} else {
			return;
		}
	}

	ind = highestThresh( node );
	if( ind != -1 ) {
		followPath( node->subnodes[ind], ind, sequence );
		
		// Only if something is contributed to the sequence should the 
		// occurences be lowered
		for( uint32_t i = 0; i < nPath.size(); i++ ) {
			nPath[i]->occs--;
			nPath[i]->weight = nPath[i]->weight.load() - 1;
		}
	}
}

/** ---------------- Tree Storage ------------------ **/
std::string GTree::storeTree( short label ) {
	treeString += GTH::valToString( nNodes ) + '\n';
	std::string val[2] = {
		GTH::valToString( root->occs ),
		GTH::valToString( root->weight )
	};
	GTH::removeDoubleEnding( val[1] );
	treeString += GTH::retLabel( label );
	treeString += ':' + val[0] + ':' + val[1] + ';';
	storeTree( root );
	return treeString;
}

void GTree::storeTree( Node* node )
{
	for( int i = 0; i < NBASES; i++ ) {
		if( node->subnodes[i] ) {
			treeString += GTH::retLabel( i );
			std::string val[2] = {
				GTH::valToString( node->subnodes[i]->occs ),
				GTH::valToString( node->subnodes[i]->weight )
			};
			GTH::removeDoubleEnding( val[1] );
			treeString += ':' + val[0] + ':' + val[1] + ';';
			storeTree( node->subnodes[i] );
		}
	}
	treeString += ',';
}

/** ---------------- Path Printing ----------------- **/
void GTree::printAllPaths( short label ) { 
	printAllPaths( root, 0, label );
	basePaths = occuPaths = "";
}

void GTree::printAllPaths( Node* node, int len, short label )
{
    if( !node )
        return;

	occuPaths.erase( len, occuPaths.length() );
	std::string val = GTH::valToString( node->occs );
    occuPaths += val;
	occuPaths += "-";
	basePaths.erase( len, basePaths.length() );
	basePaths += GTH::retLabel( label );
	for( uint32_t i = 0; i < val.length(); i++ )
		basePaths += '-';
	len += val.length() + 1;
	bool check = 0;
	for( int i = 0; i < NBASES; i++ )
		if( node->subnodes[i] )
			check = true;
	std::cout << "WEIGHT: " << node->weight << '\n';
    if( !check ) {
		occuPaths.erase( occuPaths.length() - 1 );
		std::cout << occuPaths << ": EOS" << std::endl;
		basePaths.erase( basePaths.length() - 1 );
		std::cout << basePaths << ": EOS\n" << std::endl;
        return;
    }
	for( int i = 0; i < NBASES; i++ )
		printAllPaths( node->subnodes[i], len, i );
}

