/********************************************************************
 * Filename: GTree.C [C++ source code]
 *
 * Description: Implementation of GTree class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "GTree.H"

#ifdef __MINGW64__
#include <sstream>
namespace mingw_fix {
	template< typename T > std::string to_string(const T &val)
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

	char retLabel(int label)
	{
		switch(label) {
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

	signed short mostOccs(Node *node)
	{
		short ind = -1;
		int64_t occCnt = 0;
		for(int i = 0; i < NBASES; i++) {
			Node* tmpNode = node->subnodes[i];
			if(tmpNode && tmpNode->occs > occCnt) {
				occCnt = tmpNode->occs;
				ind = i;
			}
		}
		return ind;
	}

	short countChildren(Node *node)
	{
		short children = 0;

		for(int i = 0; i < NBASES; i++)
			if(node->subnodes[i])
				children++;

		return children;
	}

	float getNewWeight(double curWeight, int64_t occs, char qual)
	{
		// Cumulative average: A_{n+1} = ((A_n + x_{n+1}) - A_n) / n + 1
		return curWeight += ((static_cast<double>(qual) - curWeight) /
				static_cast<double>(occs));
	}

	void updateWeight(Node* node, char qual)
	{
		double newWeight, curWeight = node->weight;
#pragma omp atomic
		node->occs++;
		do {
			newWeight = GTH::getNewWeight(curWeight, node->occs, qual);
		} while(!(node->weight.compare_exchange_weak(curWeight, newWeight)));
	}

	void removeDoubleEnding(std::string& doubleString)
	{
		doubleString.erase(doubleString.find_last_not_of('0') + 1, std::string::npos);
		doubleString.erase(doubleString.find_last_not_of('.') + 1, std::string::npos);
	}
}

/** ------------- GTree Cons and Dees -------------- **/
GTree::GTree(): 
	root(nullptr), head(0), nodesCnt(0), basePaths(""), occuPaths(""), treeString("")
{ 
	// Because reallocation of std::vectors is assured if its class is a 
	// std::vector, a vector of vectors has been used with a 1-by-1 
	// .push_back()
	nodes.resize(1);
	nodes[nodesCnt].resize(RES);
	omp_init_lock(&lock);
}

/** ---------------- Tree Creation ----------------- **/
void GTree::createRoot(short ind)
{
	int offset = -1;
	uint32_t i;
	char qual;
	root = &(nodes[nodesCnt][head]);

	for(i = 0; i < GTH::seqReads.size(); i++) {
		offset = -1;
		// Find algorithm
		for(int j = 0; j < GTH::seqReads[i].size(); j++) {
			if(GTH::seqReads[i].getBaseInd(j) == ind) {
				offset = j;
				break;
			}
		}
		if(offset != -1) {
			qual = GTH::seqReads[i].getQual(offset);
			root->readNum = i;
			root->offset = offset;
			root->weight = qual;
			root->occs = 1;
			break;
		}
	}
	// Create second for balancing purposes
	offset++;
	if(offset && offset != GTH::seqReads[i].size()) {
		ind = GTH::seqReads[i].getBaseInd(offset);
		qual = GTH::seqReads[i].getQual(offset);
		createNode(root, ind, qual, i, offset);
	}
}

void GTree::createNode(Node *node, short ind, char qual, uint64_t rN, int offset)
{
	omp_set_lock(&lock);
	head++;
	if(head == RES) {
		std::vector<Node> tmpVec(RES);
		nodes.push_back(tmpVec);
		nodesCnt++;
		head = 0;
	}
	node->subnodes[ind] = &(nodes[nodesCnt][head]);
	omp_unset_lock(&lock);

	// Doesn't set occs as addRead...() will do that
	Node *tmpNode = node->subnodes[ind];
	tmpNode->readNum = rN;
	tmpNode->offset = offset;
	tmpNode->weight = qual;
	tmpNode->occs = 1;
	for(int i = 0; i < NBASES; i++)
		tmpNode->subnodes[i] = nullptr;
}

/** --------------- Read Processing ---------------- **/
void GTree::addReadOne(uint64_t readNum, short offset) 
{
	std::vector<Node*> paths;

	Node *node = root;
	SeqRead *read = &GTH::seqReads[readNum];
	bool retBool = 0, clearBool = 0;

	if(root->offset == offset && root->readNum == readNum)
		return;

	for(int i = offset + 1; i < GTH::seqReads[readNum].size(); i++) {
		short ind = (*read).getBaseInd(i);
		paths.push_back(node);

		omp_set_lock(&node->lock);
		if(!node->subnodes[ind]) {
			createNode(node, ind, (*read).getQual(i), readNum, i);
			retBool = clearBool = 1;
			if(GTH::countChildren(node) == 1) {
				clearBool = balanceNode(node, 1);
			}
		}
		if(i + 1 == GTH::seqReads[readNum].size() &&
				node->subnodes[ind] &&
				!GTH::countChildren(node->subnodes[ind])) {
			balanceNode(node->subnodes[ind], 0);
		}
		omp_unset_lock(&node->lock);

		if(clearBool) {
			for(int j = 0, k = offset; j < paths.size(); j++, k++)
				GTH::updateWeight(paths[j], (*read).getQual(k));
		}
		node = node->subnodes[ind];
		if(retBool){
			paths.clear();
			// If EOS is reached, occurrences should be increased
			break;
		}
	}
}

bool GTree::balanceNode(Node *node, bool mode)
{
	std::vector<Node*> paths;
	// Get the offset and read before overiding/updating
	int64_t lReadNum = node->readNum;
	SeqRead *lRead = &GTH::seqReads[lReadNum];
	short lOffset = node->offset + 1;
	short tmpOff = lOffset;

	// If there is nothing to balance with:
	// (will obviously cause imbalanced weights/occs)
	if(lOffset == (*lRead).size()) 
		return 0;

	short lInd = (*lRead).getBaseInd(lOffset);
	char lQual = (*lRead).getQual(lOffset);
	
	if(node->subnodes[lInd] && 
			lOffset == node->subnodes[lInd]->offset && 
			lReadNum == node->subnodes[lInd]->readNum)
		return 0;

	// If the paths are different:
	if(!node->subnodes[lInd]) {
		createNode(node, lInd, lQual, lReadNum, lOffset);
		return mode;
	}
	node = node->subnodes[lInd];

	// If the paths follow the same route:
	int64_t rReadNum = node->readNum;
	SeqRead *rRead = &GTH::seqReads[rReadNum];
	short rOffset = node->offset + 1;
	short rInd = 0;
	char rQual;
	lOffset++;
	while(lOffset < (*lRead).size() && rOffset < (*rRead).size()) {
		lInd = (*lRead).getBaseInd(lOffset);
		rInd = (*rRead).getBaseInd(rOffset);
		lQual = (*lRead).getQual(lOffset);
		rQual = (*rRead).getQual(rOffset);

		paths.push_back(node);

		if((*rRead).getBaseInd(rOffset) != (*lRead).getBaseInd(lOffset)) {
			createNode(node, lInd, lQual, lReadNum, lOffset);
			createNode(node, rInd, rQual, rReadNum, rOffset);

			for(unsigned short j = 0, k = tmpOff; j < paths.size(); j++, k++)
				GTH::updateWeight(paths[j], (*lRead).getQual(k));

			return 1;
		}

		// Because we need to save the balance info of the same node regardless
		// of which read is being processed
		createNode(node, rInd, rQual, rReadNum, rOffset);

		node = node->subnodes[lInd];

		lOffset++;
		rOffset++;
	} 

	// There will be imbalances caused by running out of read length
	if(lOffset < (*lRead).size()) {
		lInd = (*lRead).getBaseInd(lOffset);
		lQual = (*lRead).getQual(lOffset);
		createNode(node, lInd, lQual, lReadNum, lOffset);
	} 
	if(rOffset < (*rRead).size()) {
		rInd = (*rRead).getBaseInd(rOffset);
		rQual = (*rRead).getQual(rOffset);
		createNode(node, rInd, rQual, rReadNum, rOffset);
	}

	return 0;
}

/** --------------- Genome Creation ---------------- **/
void GTree::followPath(Node *node, short ind, std::string &sequence)
{
	short children;
	do {
		children = GTH::countChildren(node);
		sequence += GTH::retLabel(ind);
		ind = GTH::mostOccs(node);
		node->occs--;
		if(ind == -1)
			return;
		node = node->subnodes[ind];
	} while(children);
}

void GTree::addToSeq(uint64_t offset, std::string &sequence)
{
	Node *node = root;
	short ind = 0;
	uint32_t i, j, seqLength = sequence.length() - offset;
	Node **path = new Node*[seqLength];

	for(i = 0; i < seqLength; i++)
		path[i] = nullptr;

	// End of current sequence	
	for(i = offset + 1, j = 0; i < sequence.length(); i++, j++) {
		ind = BASE_IND(sequence[i]);
		if(node->subnodes[ind] && GTH::countChildren(node->subnodes[ind])) {
			path[j] = node;
			node = node->subnodes[ind];
		} else {
			delete[] path;
			return;
		}
	}

	ind = GTH::mostOccs(node);
	if(ind != -1) {
		followPath(node->subnodes[ind], ind, sequence);
		
		// Only if something is contributed to the sequence should the 
		// occurences be lowered
		for(i = 0; i < seqLength; i++)
			if(path[i]) 
				path[i]->occs--;
	}

	delete[] path;
}

/** ---------------- Tree Storage ------------------ **/
std::string GTree::storeTree(short label) {
	treeString += GTH::retLabel(label);
#ifdef __MINGW64__
	std::string val[2] = {
		mingw_fix::to_string(root->occs),
		mingw_fix::to_string(root->weight)
	};
#else
	std::string val[2] = {
		std::to_string(root->occs),
		std::to_string(root->weight)
	};
#endif /*__MINGW32__*/
	GTH::removeDoubleEnding(val[1]);
	treeString += ':' + val[0] + ':' + val[1] + ';';
	storeTree(root);
	return treeString;
}

void GTree::storeTree(Node* node)
{
	for(int i = 0; i < NBASES; i++) {
		if(node->subnodes[i]) {
			treeString += GTH::retLabel(i);
#ifdef __MINGW64__
			std::string val[2] = {
				mingw_fix::to_string(node->subnodes[i]->occs),
				mingw_fix::to_string(node->subnodes[i]->weight)
			};
#else
			std::string val[2] = {
				std::to_string(node->subnodes[i]->occs),
				std::to_string(node->subnodes[i]->weight)
			};
#endif /*__MINGW32__*/
			GTH::removeDoubleEnding(val[1]);
			treeString += ':' + val[0] + ':' + val[1] + ';';
			storeTree(node->subnodes[i]);
		}
	}
	treeString += ',';
}

/** ---------------- Path Printing ----------------- **/
void GTree::printAllPaths(Node* node, int len, short label)
{
    if(!node)
        return;

	occuPaths.erase(len, occuPaths.length());
#ifdef __MINGW64__
	std::string val = mingw_fix::to_string(node->occs);
#else
	std::string val = std::to_string(node->occs);
#endif /*__MINGW32__*/
    occuPaths += val;
	occuPaths += "-";
	basePaths.erase(len, basePaths.length());
	basePaths += GTH::retLabel(label);
	for(uint32_t i = 0; i < val.length(); i++)
		basePaths += '-';
	len += val.length() + 1;
	bool check = 0;
	for(int i = 0; i < NBASES; i++)
		if(node->subnodes[i])
			check = true;
    if(!check) {
		occuPaths.erase(occuPaths.length() - 1);
		std::cout << occuPaths << ": EOS" << std::endl;
		basePaths.erase(basePaths.length() - 1);
		std::cout << basePaths << ": EOS\n" << std::endl;
        return;
    }
	for(int i = 0; i < NBASES; i++)
		printAllPaths(node->subnodes[i], len, i);
}
