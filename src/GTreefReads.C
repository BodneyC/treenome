/********************************************************************
 * Filename: GTreefReads.C [C++ source code]
 *
 * Description: Implementation of GTree subclass, for building the 
 *		tree from a set of reads
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/GTree.H"

GTreefReads::GTreefReads(): nodesCnt( 0 )
{
	// Because reallocation of std::vectors is assured if its class is a 
	// std::vector, a vector of vectors has been used with a 1-by-1 
	// .push_back()
	nodes.resize( 1 );
	nodes[nodesCnt].resize( RES );

	omp_init_lock( &lock );
}
/** ---------------- Tree Creation ----------------- **/
void GTreefReads::updateWeight( Node* node, char qual ) {
	double newWeight, curWeight = node->weight;
	double pBQual = GTH::phredQuals[static_cast<int>( qual )];
	do {
		newWeight = curWeight + pBQual;
	} while( !( node->weight.compare_exchange_weak( curWeight, newWeight ) ) );
}

void GTreefReads::updateWeightAndOccs( Node* node, char qual )
{
	node->occs++;
	updateWeight(node, qual);
}

void GTreefReads::createRoot( short ind )
{
	uint32_t i;
	char qual;
	root = &( nodes[nodesCnt][head] );
	nNodes++;

	for( i = 0; i < GTH::seqReads.size(); i++ ) {
		int offset = -1;
		// Find algorithm
		for( int j = 0; j < GTH::seqReads[i].size(); j++ ) {
			if( GTH::seqReads[i].getBaseInd(j) == ind ) {
				offset = j;
				break;
			}
		}
		if( offset != -1 ) {
			qual = GTH::seqReads[i].getQual( offset );
			root->readNum = i;
			root->offset = offset;
			root->weight = GTH::phredQuals[static_cast<int>( qual )];
			root->occs = 1;
			break;
		}
	}
	// Create second for balancing purposes
	//offset++;
	//if( offset && offset != GTH::seqReads[i].size() ) {
	//	ind = GTH::seqReads[i].getBaseInd( offset );
	//	qual = GTH::seqReads[i].getQual( offset );
	//	createNode( root, ind, qual, i, offset );
	//}
}

void GTreefReads::createNode( Node* node, short ind, char qual, uint64_t rN, int offset )
{
	omp_set_lock( &lock );
	head++;
	if( head == RES ) {
		std::vector<Node> tmpVec( RES );
		nodes.push_back( tmpVec );
		nodesCnt++;
		head = 0;
	}
	node->subnodes[ind] = &( nodes[nodesCnt][head] );

	nNodes++;
	omp_unset_lock( &lock );

	// Doesn't set occs as addRead...() will do that
	Node* tmpNode = node->subnodes[ind];
	tmpNode->readNum = rN;
	tmpNode->offset = offset;
	tmpNode->weight = GTH::phredQuals[static_cast<int>( qual )];
	tmpNode->occs = 1;
	for( int i = 0; i < NBASES; i++ )
		tmpNode->subnodes[i] = nullptr;
}

/** --------------- Read Processing ---------------- **/
void GTreefReads::addReadOne( uint64_t readNum, short offset ) 
{
	//std::lock_guard<std::mutex> cncn(gtMut);
	std::vector<Node*> paths;

	Node* node = root;
	SeqRead* read = &GTH::seqReads[readNum];
	bool retBool = 0, clearBool = 0;

	if( root->offset == offset && ( uint64_t )root->readNum == readNum ) {
		return;
	}

	for( int i = offset + 1; i < GTH::seqReads[readNum].size(); i++ ) {
		short ind = ( *read ).getBaseInd( i );
		paths.push_back( node );

		omp_set_lock( &node->lock );
		if( !node->subnodes[ind] ) {
			//std::cout << ( *read ).getQual( i ) << std::endl;
			createNode( node, ind, ( *read ).getQual( i ), readNum, i );
			retBool = clearBool = 1;
			if( countChildren( node ) == 1 ) {
				balanceNode( node, 1 );
			}
		}
		//if( i + 1 == GTH::seqReads[readNum].size() &&
		//		node->subnodes[ind] &&
		//		!countChildren( node->subnodes[ind]) ) {
		//	balanceNode( node->subnodes[ind], 0 );
		//}
		omp_unset_lock( &node->lock );

		if( clearBool ) {
			for( unsigned int j = 0, k = offset; j < paths.size(); j++, k++ )
				updateWeightAndOccs( paths[j], ( *read ).getQual(k) );
		}
		if( retBool ){
			//std::cout << "THR: " << omp_get_thread_num() << std::endl;
			//printAllPaths(1);
			paths.clear();
			// If EOS is reached, occurrences should be increased
			break;
		}

		node = node->subnodes[ind];
	}
	//std::cout << "FISHF" <<std::endl;
}

bool GTreefReads::balanceNode( Node* node, bool mode )
{
	// Get the offset and read before overiding/updating
	int64_t lReadNum = node->readNum;
	SeqRead* lRead = &GTH::seqReads[lReadNum];
	short lOffset = node->offset + 1;

	// If there is nothing to balance with:
	// (will obviously cause imbalanced weights/occs)
	if( lOffset == ( *lRead ).size() ) 
		return 0;

	short lInd = ( *lRead ).getBaseInd( lOffset );
	char lQual = ( *lRead ).getQual( lOffset );
	
	if(node->subnodes[lInd] && 
			lOffset == node->subnodes[lInd]->offset && 
			lReadNum == node->subnodes[lInd]->readNum)
		return 0;

	// If the paths are different:
	if( !node->subnodes[lInd] ) {
		createNode( node, lInd, lQual, lReadNum, lOffset );
		return mode;
	}
	node = node->subnodes[lInd];
	updateWeightAndOccs(node, lQual);

	// If the paths follow the same route:
	int64_t rReadNum = node->readNum;
	SeqRead* rRead = &GTH::seqReads[rReadNum];
	short rOffset = node->offset + 1;
	lOffset++;

	std::vector<short> lIndPath;
	std::vector<short> rIndPath;
	short tmpLOffset = lOffset, tmpROffset = rOffset;
	short tmpLOffset2 = lOffset, tmpROffset2 = rOffset;

	while( lOffset < ( *lRead ).size() && rOffset < ( *rRead ).size() ) {
		lIndPath.push_back( ( *lRead ).getBaseInd( lOffset ) );
		rIndPath.push_back( ( *rRead ).getBaseInd( rOffset ) );

		if(( *rRead ).getBaseInd( rOffset ) != ( *lRead ).getBaseInd( lOffset )) {
			for( unsigned short i = 0; i < lIndPath.size() - 1; i++ ) {
				createNode( node, lIndPath[i], ( *lRead ).getQual( tmpLOffset ), lReadNum, tmpLOffset );
				node = node->subnodes[lIndPath[i]];
				updateWeightAndOccs( node, ( *rRead ).getQual( tmpROffset ) );
				tmpLOffset++;
				tmpROffset++;
			}
			createNode( node, lIndPath.back(), ( *lRead ).getQual( tmpLOffset ), lReadNum, tmpLOffset );
			createNode( node, rIndPath.back(), ( *rRead ).getQual( tmpROffset ), rReadNum, tmpROffset );
			return 1;
		}

		lOffset++;
		rOffset++;
	} 

	if( lOffset < ( *lRead ).size() ) {
		lInd = ( *lRead ).getBaseInd( tmpLOffset2 );
		lQual = ( *lRead ).getQual( tmpLOffset2 );
		createNode( node, lInd, lQual, lReadNum, tmpLOffset2 );
	} 
	if( rOffset < ( *rRead ).size() ) {
		short rInd = ( *rRead ).getBaseInd( tmpROffset2 );
		char rQual = ( *rRead ).getQual(tmpROffset2 );
		createNode( node, rInd, rQual, rReadNum, tmpROffset2 );
	}

	return 0;
}
