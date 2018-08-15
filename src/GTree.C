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

int64_t GTree::getDNodeSize()
{
	return dNodes.size();
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
void GTree::addToSeq( uint32_t offset, std::string &sequence )
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
void GTree::writeTreeToFile( std::ofstream& storeFile )
{
	for( int64_t i = 0; i < dNodes.size(); i++ )
		storeFile.write( ( char* )&dNodes[i], sizeof( Node ) );
}

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

/** --------------- Tree from Reads ---------------- **/
void GTree::init() {
	// Because reallocation of std::vectors is assured if its class is a 
	// std::vector, a vector of vectors has been used with a 1-by-1 
	// .push_back()
	omp_init_lock( &lock );
	dNodes.resize( RES );
}

/** ---------------- Tree Creation ----------------- **/
void GTree::updateWeight( Node* node, char qual ) {
	double newWeight, curWeight = node->weight;
	double pBQual = GTH::phredQuals[GTH::scoreSys][static_cast<int>( qual )];
	do {
		newWeight = curWeight + pBQual;
	} while( !( node->weight.compare_exchange_weak( curWeight, newWeight ) ) );
}

void GTree::updateWeightAndOccs( Node* node, char qual )
{
	node->occs++;
	updateWeight(node, qual);
}

void GTree::createRoot( short ind )
{
	uint32_t i;
	char qual;
	root = &( dNodes[head] );
	nNodes++;

	for( i = 0; i < GTH::seqReads.size(); i++ ) {
		int offset = -1;
		// Find algorithm in bool vec
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
			root->weight = GTH::phredQuals[GTH::scoreSys][static_cast<int>( qual )];
			root->occs = 1;
			break;
		}
	}
}

void GTree::createNode( Node* node, short ind, char qual, uint32_t rN, int offset )
{
	omp_set_lock( &lock );
	head++;
	if( ( unsigned ) head == dNodes.size() ) {
		dNodes.resize( dNodes.size() + RES );
	}
	node->subnodes[ind] = &( dNodes[head] );
	omp_unset_lock( &lock );
	nNodes++;

	// Doesn't set occs as addRead...() will do that
	Node* tmpNode = node->subnodes[ind];
	tmpNode->readNum = rN;
	tmpNode->offset = offset;
	tmpNode->weight = GTH::phredQuals[GTH::scoreSys][static_cast<int>( qual )];
	tmpNode->occs = 1;
	for( int i = 0; i < NBASES; i++ )
		tmpNode->subnodes[i] = nullptr;
}

/** --------------- Read Processing ---------------- **/
void GTree::addReadOne( uint32_t readNum, short offset ) 
{
	//std::lock_guard<std::mutex> cncn(gtMut);
	std::vector<Node*> paths;

	Node* node = root;
	SeqRead* read = &GTH::seqReads[readNum];
	bool returnBool = 0, updateBool = 0;

	// Edge case
	if( root->offset == offset && ( uint32_t )root->readNum == readNum )
		return;

	for( int i = offset + 1; i < GTH::seqReads[readNum].size(); i++ ) {
		short ind = ( *read ).getBaseInd( i );
		paths.push_back( node );

		omp_set_lock( &node->lock );
		if( !node->subnodes[ind] ) {
			//std::cout << ( *read ).getQual( i ) << std::endl;
			createNode( node, ind, ( *read ).getQual( i ), readNum, i );
			returnBool = updateBool = 1;
			if( countChildren( node ) == 1 )
				// Even if it doesn't balance it, a node was created and so
				// the path needs updating
				balanceNode( node );
		}
		if( i + 1 == GTH::seqReads[readNum].size() ) {
			// If the end is reached, they still count as occurrences of the path...
			updateBool = 1;
			// If a subnode exists and it wasn't created above
			if( node->subnodes[ind] && !returnBool)
				paths.push_back( node->subnodes[ind] );
			// For balancing purposes
			if( node->subnodes[ind] && !countChildren( node->subnodes[ind]) )
				potAddNode( node->subnodes[ind] );
		}
		omp_unset_lock( &node->lock );

		// It will always enter updateBool, just needs to know when
		if( updateBool )
			for( unsigned int j = 0, k = offset; j < paths.size(); j++, k++ )
				updateWeightAndOccs( paths[j], ( *read ).getQual(k) );
		if( returnBool )
			break;

		node = node->subnodes[ind];
	}
	paths.clear();
}

void GTree::potAddNode(Node* node)
{
	int32_t rN = node->readNum;
	SeqRead* tRead = &GTH::seqReads[rN];
	short offset = node->offset + 1;

	if( offset == ( *tRead ).size() )
		return;

	short ind = ( *tRead ).getBaseInd( offset );
	char qual = ( *tRead ).getQual( offset );

	createNode( node, ind, qual, rN, offset );
}

void GTree::balanceNode( Node* node )
{
	// Get the offset and read before overiding/updating
	int32_t lReadNum = node->readNum;
	SeqRead* lRead = &GTH::seqReads[lReadNum];
	short lOffset = node->offset + 1;
	bool clearBool = 0;

	// If there is nothing to balance with:
	// (will obviously cause imbalanced weights/occs)
	if( lOffset == ( *lRead ).size() )
		return;

	short lInd = ( *lRead ).getBaseInd( lOffset );
	char lQual = ( *lRead ).getQual( lOffset );
	
	if(node->subnodes[lInd] && 
			lOffset == node->subnodes[lInd]->offset && 
			lReadNum == node->subnodes[lInd]->readNum)
		return;

	// If the paths are different:
	if( !node->subnodes[lInd] ) {
		createNode( node, lInd, lQual, lReadNum, lOffset );
		return;
	}
	node = node->subnodes[lInd];
	updateWeightAndOccs(node, lQual);

	// If the paths follow the same route:
	int32_t rReadNum = node->readNum;
	SeqRead* rRead = &GTH::seqReads[rReadNum];
	short rOffset = node->offset + 1;
	lOffset++;

	std::vector<short> lIndPath;
	std::vector<short> rIndPath;
	short tmpLOffset = lOffset, tmpROffset = rOffset;

	while( lOffset < ( *lRead ).size() && rOffset < ( *rRead ).size() ) {
		lIndPath.push_back( ( *lRead ).getBaseInd( lOffset ) );
		rIndPath.push_back( ( *rRead ).getBaseInd( rOffset ) );
		if(( *rRead ).getBaseInd( rOffset ) != ( *lRead ).getBaseInd( lOffset )) {
			clearBool = 1;
			break;
		}
		lOffset++;
		rOffset++;
	} 

	// Update to the end of shared bases
	for( unsigned short i = 0; i < lIndPath.size() - clearBool; i++ ) {
		createNode( node, lIndPath[i], ( *lRead ).getQual( tmpLOffset ), lReadNum, tmpLOffset );
		node = node->subnodes[lIndPath[i]];
		updateWeightAndOccs( node, ( *rRead ).getQual( tmpROffset ) );
		tmpLOffset++;
		tmpROffset++;
	}

	// If they differ, create a node each
	if( clearBool ) {
		createNode( node, lIndPath.back(), ( *lRead ).getQual( tmpLOffset ), lReadNum, tmpLOffset );
		createNode( node, rIndPath.back(), ( *rRead ).getQual( tmpROffset ), rReadNum, tmpROffset );
		return;
	}

	// If one ends, create the relevant node
	if( lOffset < ( *lRead ).size() ) {
		lInd = ( *lRead ).getBaseInd( tmpLOffset );
		lQual = ( *lRead ).getQual( tmpLOffset );
		createNode( node, lInd, lQual, lReadNum, tmpLOffset );
	} 
	if( rOffset < ( *rRead ).size() ) {
		short rInd = ( *rRead ).getBaseInd( tmpROffset );
		char rQual = ( *rRead ).getQual(tmpROffset );
		createNode( node, rInd, rQual, rReadNum, tmpROffset );
	}

	return;
}

/** --------------- Tree From File ----------------- **/
void GTree::resizeDeque( int64_t tmp64 )
{
	nNodes = tmp64;
	dNodes.resize( nNodes );
}

void GTree::writeDeque( std::ifstream& inFile )
{
	inFile.read( ( char* ) &dNodes, nNodes * sizeof( Node ) );
}


void GTree::getNextNode( struct NodeInfo& nInf, std::stringstream& ss )
{
	std::string line;
	char tmpChar;
	while( 1 ) {
		if( !(ss >> tmpChar) )
			return;
		if( tmpChar != ',' )
			break;
		nInf.comCnt++;
	}
	nInf.ind = BASE_IND( tmpChar );
	ss >> tmpChar;
	std::getline( ss, line, ':' );
	nInf.occs = std::stol( line );
	std::getline( ss, line, ';' );
	nInf.weight = std::stod( line );
}

void GTree::createNode( struct NodeInfo& nInf )
{
	for( int i = 0; i < nInf.comCnt; i++ )
		tmpNode = tmpNode->parent;
	
	tmpNode->subnodes[nInf.ind] = &( dNodes[head] );
	head++;

	tmpNode->subnodes[nInf.ind]->parent = tmpNode;
	tmpNode = tmpNode->subnodes[nInf.ind];

	tmpNode->occs = nInf.occs;
	tmpNode->weight = nInf.weight;
}

void GTree::createRoot( std::stringstream& ss )
{
	struct NodeInfo nInf;
	getNextNode( nInf, ss );
	root = &( dNodes[0] );
	tmpNode = root;
	head++;

	tmpNode->occs = nInf.occs;
	tmpNode->weight = nInf.weight;
}

void GTree::processSString( std::stringstream& ss )
{
	std::string line;
	std::getline( ss, line, '\n' );

	nNodes = std::stol( line );
	dNodes.resize( nNodes );

	createRoot( ss );

	for( unsigned int i = 0; i < nNodes - 1; i++ ) {
		struct NodeInfo nInf;
		getNextNode( nInf, ss );
		createNode( nInf );
	}
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

