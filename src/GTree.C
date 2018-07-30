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
	template< typename T > std::string to_string( const T &val )
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

	//float getNewWeight( double curWeight, int64_t occs, char qual )
	//{
	//	// Cumulative average: A_{n+1} = ((A_n + x_{n+1}) - A_n) / n + 1
	//	return curWeight += ( ( static_cast<double>( qual ) - curWeight ) /
	//			static_cast<double>( occs ) );
	//}

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
template <typename T>
GTree<T>::GTree(): 
	root( nullptr ), head( 0 ), nNodes( 0 ), 
	basePaths( "" ), occuPaths( "" ), treeString( "" )
{ 
}

/** ------------- - Helper Functions --------------- **/
template <typename T>
short GTree<T>::countChildren( T* node )
{
	short children = 0;

	for( int i = 0; i < NBASES; i++ )
		if( node->subnodes[i] )
			children++;

	return children;
}

template <typename T>
signed short GTree<T>::mostOccs( T* node )
{
	short ind = -1;
	double maxRat = std::numeric_limits<double>::lowest();

	for( short i = 0; i < NBASES; i++ ) {
		T* tmpNode = node->subnodes[i];
		if( tmpNode && tmpNode->occs ) {
			double curRat = tmpNode->weight / static_cast<double>(tmpNode->occs);
			if( curRat > maxRat ) {
				maxRat = curRat;
				ind = i;
			}
		}
	}

	return ind;
}

template <typename T>
void GTree<T>::followPath( T* node, short ind, std::string &sequence )
{
	short children;
	do {
		children = countChildren( node );
		sequence += GTH::retLabel( ind );
		ind = mostOccs( node );
		node->occs--;
		double tmpD = node->weight - 1;
		node->weight = tmpD;
		if( ind == -1 )
			return;
		node = node->subnodes[ind];
	} while( children );
}

/** --------------- Genome Creation ---------------- **/
template <typename T>
void GTree<T>::addToSeq( uint64_t offset, std::string &sequence )
{
	T* node = root;
	short ind = 0;
	uint32_t i, j, seqLength = sequence.length() - offset;
	T** nPath = new T*[seqLength];

	for( i = 0; i < seqLength; i++ )
		nPath[i] = nullptr;

	// End of current sequence	
	for( i = offset + 1, j = 0; i < sequence.length(); i++, j++ ) {
		ind = BASE_IND( sequence[i] );
		if( node->subnodes[ind] && countChildren( node->subnodes[ind] ) ) {
			nPath[j] = node;
			node = node->subnodes[ind];
		} else {
			delete[] nPath;
			return;
		}
	}

	ind = mostOccs( node );
	if( ind != -1 ) {
		followPath( node->subnodes[ind], ind, sequence );
		
		// Only if something is contributed to the sequence should the 
		// occurences be lowered
		for( i = 0; i < seqLength; i++ )
			if( nPath[i] ) {
				nPath[i]->occs--;
				double tmpD = nPath[i]->weight - 1;
				nPath[i]->weight = tmpD;
			}
	}

	delete[] nPath;
}

/** ---------------- Tree Storage ------------------ **/
template <typename T>
std::string GTree<T>::storeTree( short label ) {
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

template <typename T>
void GTree<T>::storeTree( T* node )
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
template <typename T>
void GTree<T>::printAllPaths( T* node, int len, short label )
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
	std::cout << "\nWEIGHT: " << node->weight << '\n';
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

template class GTree<Node>;
template class GTree<pNode>;
