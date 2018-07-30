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
	int offset = -1;
	uint32_t i;
	char qual;
	root = &( nodes[nodesCnt][head] );
	nNodes++;

	for( i = 0; i < GTH::seqReads.size(); i++ ) {
		offset = -1;
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

	// Doesn't set occs as addRead...() will do that
	Node* tmpNode = node->subnodes[ind];
	tmpNode->readNum = rN;
	tmpNode->offset = offset;
	tmpNode->weight = GTH::phredQuals[static_cast<int>( qual )];
	tmpNode->occs = 1;
	for( int i = 0; i < NBASES; i++ )
		tmpNode->subnodes[i] = nullptr;
	omp_unset_lock( &lock );
}

/** --------------- Read Processing ---------------- **/
void GTreefReads::addReadOne( uint64_t readNum, short offset ) 
{
	std::vector<Node*> paths;

	Node* node = root;
	SeqRead* read = &GTH::seqReads[readNum];
	bool retBool = 0, clearBool = 0;

	if( root->offset == offset && root->readNum == readNum ) {
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
		if( i + 1 == GTH::seqReads[readNum].size() &&
				node->subnodes[ind] &&
				!countChildren( node->subnodes[ind]) ) {
			balanceNode( node->subnodes[ind], 0 );
		}

		if( clearBool ) {
			for( int j = 0, k = offset; j < paths.size(); j++, k++ )
				updateWeightAndOccs( paths[j], ( *read ).getQual(k) );
		}
		omp_unset_lock( &node->lock );
		if( retBool ){
			paths.clear();
			// If EOS is reached, occurrences should be increased
			break;
		}
		node = node->subnodes[ind];
	}
}

bool GTreefReads::balanceNode( Node* node, bool mode )
{
	// Get the offset and read before overiding/updating
	int64_t lReadNum = node->readNum;
	SeqRead* lRead = &GTH::seqReads[lReadNum];
	short lOffset = node->offset + 1;
	short tmpOff = lOffset;

	// If there is nothing to balance with:
	// (will obviously cause imbalanced weights/occs)
	if( lOffset == ( *lRead ).size() ) 
		return 0;

	short lInd = ( *lRead ).getBaseInd(lOffset );
	char lQual = ( *lRead ).getQual(lOffset );
	
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

	// If the paths follow the same route:
	int64_t rReadNum = node->readNum;
	SeqRead* rRead = &GTH::seqReads[rReadNum];
	short rOffset = node->offset + 1;
	short rInd = 0;
	char rQual;
	lOffset++;

	std::vector<Node*> paths;
	std::vector<short> indPath;
	std::vector<int64_t> rnPath;
	short tmpLOffset = lOffset, tmpROffset = rOffset;
	//std::string tmpLString = "", tmpRString = "";
	//for( ; tmpLOffset < ( *lRead ).size(); tmpLOffset++)
	//	tmpLString += ( *lRead ).getCharBase( tmpLOffset );

	while( lOffset < ( *lRead ).size() && rOffset < ( *rRead ).size() ) {
		lInd = ( *lRead ).getBaseInd( lOffset );
		lQual = ( *lRead ).getQual( lOffset );
		rInd = ( *rRead ).getBaseInd( rOffset );
		rQual = ( *rRead ).getQual( rOffset );

		paths.push_back( node );


		if(( *rRead ).getBaseInd( rOffset ) != ( *lRead ).getBaseInd( lOffset )) {
			createNode( node, lInd, lQual, lReadNum, lOffset );
			createNode( node, rInd, rQual, rReadNum, rOffset );

			for( unsigned short j = 0, k = tmpOff; j < paths.size(); j++, k++ )
				updateWeightAndOccs( paths[j], ( *lRead ).getQual( k ) );

			return 1;
		}

		bool side = 0;
		if( rReadNum < lReadNum || ( rReadNum == lReadNum && rOffset < lOffset ))
			side = 1;
		if( side ) {
			createNode( node, lInd, lQual, lReadNum, lOffset );
			//updateWeightAndOccs(node->subnodes[rInd], rQual);
		} else {
			createNode( node, rInd, rQual, rReadNum, rOffset );
			//updateWeightAndOccs(node->subnodes[lInd], lQual);
		}

		node = node->subnodes[lInd];

		lOffset++;
		rOffset++;
	} 


	// There will be imbalances caused by running out of read length
	if( lOffset < ( *lRead ).size() ) {
		//updateWeight(node, lQual);
		lInd = ( *lRead ).getBaseInd( lOffset );
		lQual = ( *lRead ).getQual( lOffset );
		createNode( node, lInd, lQual, lReadNum, lOffset );
	} 
	if( rOffset < ( *rRead ).size() ) {
		//updateWeight(node, rQual);
		rInd = ( *rRead ).getBaseInd( rOffset );
		rQual = ( *rRead ).getQual( rOffset );
		createNode( node, rInd, rQual, rReadNum, rOffset );
	}

	return 0;
}
