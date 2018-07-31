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

/** ------------- Sequence Generation -------------- **/
template <class T>
bool TreeTop<T>::rootOccsExist()
{
	short ret = 0;

	for( int i = 0; i < NBASES; i++ )
		if( trees[i].getRoot()->occs == 0 )
			ret++;

	return ret == NBASES ? 0 : 1;
}

template <class T>
short TreeTop<T>::maxPath()
{
	int start = 0;
	double maxRat = std::numeric_limits<double>::lowest();

	for( short i = 0; i < NBASES; i++ ) {
		if( trees[i].getRoot()->occs ) {
			double curRat = trees[i].getRoot()->weight / static_cast<double>(trees[i].getRoot()->occs);
			if( curRat > maxRat ) {
				maxRat = curRat;
				start = i;
			}
		}
	}

	trees[start].followPath( trees[start].getRoot(), start, sequence );

	return start;
}

template <class T>
void TreeTop<T>::buildSequence()
{
	maxPath();

	uint64_t offset = 1;
	while( rootOccsExist() ) {

		// Possibly make is a tighter gap as its working from single letters
		// ( this would actually be the k-mer match )
		if( offset == sequence.length() - 1 ) {
			sequence += 'N';
			maxPath();
			offset += 2;
		}

		if( trees[BASE_IND( sequence[offset] )].getRoot()->occs > 0 )
			trees[BASE_IND( sequence[offset] )].addToSeq( offset, sequence );
		offset++;
	}
}

/** --------------- Misc Functions ----------------- **/
template <class T>
void TreeTop<T>::storeTrees()
{
	for( int i = 0; i < NBASES; i++ )
		treeStrings[i] = trees[i].storeTree( i );
}

template <class T>
void TreeTop<T>::printTrees()
{
	for( int i = 0; i < NBASES; i++ )
		trees[i].printAllPaths( i );
}

template <class T>
void TreeTop<T>::printSequence()
{
	// 80 for terminal width's sake
	unsigned short TWIDTH = 120;
	uint64_t i = 0;

	if( sequence.length() > TWIDTH )
		for( ; i < sequence.length() - TWIDTH; i += TWIDTH )
			std::cout << sequence.substr( i, TWIDTH ) << std::endl;

	std::cout << sequence.substr( i, sequence.length() - i ) << std::endl;
}

template class TreeTop<GTreefReads>;
template class TreeTop<GTreefFile>;
