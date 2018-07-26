/********************************************************************
 * Filename: TreeTopfReads.C [C++ source code]
 *
 * Description: Implementation of TreeTop subclass, for building the 
 *		tree from a set of reads
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/TreeTop.H"

TreeTopfReads::TreeTopfReads()
{
	for( int i = 0; i < NBASES; i++ )
		trees[i].createRoot( i );
}

/** --------------- Read Processing ---------------- **/
void TreeTopfReads::threadFunc( uint64_t i )
{
	for( short j = 0; j < GTH::seqReads[i].size(); j++ )
		//if( GTH::seqReads[i].getBaseInd( j ) == 3 )
			trees[GTH::seqReads[i].getBaseInd( j )].addReadOne( i, j );
}

void TreeTopfReads::processReadsOne()
{
#pragma omp parallel num_threads( NUM_THREADS )
{
	for( uint64_t i = 0; i < GTH::seqReads.size(); i += NUM_THREADS ) {
	//std::cout << i << std::endl;
#pragma omp for schedule( static, 1 )
	for( int j = 0; j < NUM_THREADS; j++ ) 
		if( i + j < GTH::seqReads.size() )
			threadFunc( i + j );
	}
}
}
