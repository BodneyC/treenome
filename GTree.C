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
#include "includes/GTree.H"


//////TMP
#include<stdio.h>

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

	short mostOccs(Node *node)
	{
		short ind = -1;
		long occCnt = 0;
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

	float getNewWeight(float curWeight, long occs, char qual)
	{
		// Cumulative average: A_{n+1} = ((A_n + x_{n+1}) - A_n) / n + 1
		return curWeight += ((static_cast<float>(qual) - curWeight) /
				static_cast<double>(occs));
	}

	void updateWeight(Node* node, char qual)
	{
		float newWeight, curWeight = node->weight;
#pragma omp atomic
		node->occs++;
		do {
			newWeight = GTH::getNewWeight(curWeight, node->occs, qual);
		} while(!(node->weight.compare_exchange_weak(curWeight, newWeight)));
	}
}

/** ------------- GTree Cons and Dees -------------- **/
GTree::GTree(): 
	root(nullptr), head(0), nodesCnt(0), basePaths(""), occuPaths("")
{ 
	// Because reallocation of std::vectors is assured if its class is a 
	// std::vector, a vector of vectors has been used with a 1-by-1 
	// .push_back()
	nodes.resize(1);
	nodes[nodesCnt].resize(RES);
}

/** ---------------- Tree Creation ----------------- **/
void GTree::createRoot(short ind)
{
	int offset = -1;
	uint i;
	char qual;
	root = &(nodes[nodesCnt][head]);
	head++;

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

void GTree::createNode(Node *node, short ind, char qual, long rN, int offset)
{
	{
		//std::lock_guard<std::mutex> cnLock(gtMut);
		// node (LeafNode) becomes Node, 
		// node->subnodes[ind]
		if(nodes[nodesCnt].size() == head) {
			std::vector<Node> tmpVec(RES);
			nodes.push_back(tmpVec);
			nodesCnt++;
			head = 0;
		}
		node->subnodes[ind] = &(nodes[nodesCnt][head]);
		head++;
	}

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
void GTree::addReadOne(long readNum, short offset) 
{
	std::vector<Node*> paths;

	std::lock_guard<std::mutex> rpLock(gtMut);
	Node *node = root;
	SeqRead *read = &GTH::seqReads[readNum];
	bool retBool = 0, clearBool = 0;

	for(int i = offset + 1; i < GTH::seqReads[readNum].size(); i++) {
		short ind = (*read).getBaseInd(i);
		std::cout << "\nROOTOCCS: " << (*read).getCharBase(i)<<", thr: "<< omp_get_thread_num()<<", "<<root<<","<<node <<",";
		paths.push_back(node);
		//GTH::updateWeight(node, (*read).getQual(i-1));
		//{
			//omp_set_lock(&node->lock);
			if(!node->subnodes[ind]) {
				createNode(node, ind, (*read).getQual(i), readNum, i);
				retBool = clearBool = 1;
				if(GTH::countChildren(node) == 1) {
					//std::cout << "NEGROS" << std::endl;
					clearBool = balanceNode(node);
				}
			}
			//// ??????????????
			//if(i + 1 == GTH::seqReads[readNum].size() &&
			//		node->subnodes[ind] ) {
			//	balanceNode(node->subnodes[ind]);
			//	retBool = 1;
			//}
			//omp_unset_lock(&node->lock);
		//}
		if(clearBool) {
			for(int j = 0, k = offset; j < paths.size(); j++, k++)
				GTH::updateWeight(paths[j], (*read).getQual(k));
			paths.clear();
		}
		if(retBool){
			// If EOS is reached, occurrences should be increased

			printf("\nTHR: %d, RN: %ld\n", omp_get_thread_num(), readNum);
				printAllPaths(2);
			break;
		}
		node = node->subnodes[ind];
	}

	//std::cout << "Final occs: " << root->occs << std::endl;
}

bool GTree::balanceNode(Node *node)
{
	std::vector<Node*> paths;
	// Get the offset and read before overiding/updating
	long lReadNum = node->readNum;
	SeqRead *lRead = &GTH::seqReads[lReadNum];
	short lOffset = node->offset + 1;
	short tmpOff = lOffset;

	// If there is nothing to balance with:
	// (will obviously cause imbalanced weights/occs)
	if(lOffset == (*lRead).size()) 
		return 0;
		//std::cout << "\nRootOccs:"<<root->occs<<"\nlOffset "<<lOffset << "lRead size " << (*lRead).size() <<", lReadNum "<<lReadNum<<"\n";

	short lInd = (*lRead).getBaseInd(lOffset);
	char lQual = (*lRead).getQual(lOffset);
	
	// If the paths are different:
	if(!node->subnodes[lInd]) {
		createNode(node, lInd, lQual, lReadNum, lOffset);
		return 1;
	}
	node = node->subnodes[lInd];

	// If the paths are literally the same:
	if(lOffset == node->offset && lReadNum == node->readNum)
		return 0;

	// If the paths follow the same route:
	long rReadNum = node->readNum;
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

		//std::cout << "Node: Root: "<<node<<","<<root<<"lOffset "<<lOffset<<", lReadNum "<<lReadNum<<"\nrOffset "<<rOffset<<", rReadNum "<<rReadNum<<std::endl;

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
			//GTH::updateWeight(node, (*lRead).getQual(lOffset));

		node = node->subnodes[lInd];

		lOffset++;
		rOffset++;
	} 

	// There will be imbalances caused by running out of read length
	if(lOffset < (*lRead).size()) {
		//GTH::updateWeight(node, (*lRead).getQual(rOffset));
		createNode(node, lInd, lQual, lReadNum, lOffset);
	} 
	if(rOffset < (*rRead).size()) {
		//GTH::updateWeight(node, (*rRead).getQual(rOffset));
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

void GTree::addToSeq(long offset, std::string &sequence)
{
	Node *node = root;
	short ind = 0;
	uint i, j, seqLength = sequence.length() - offset;
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

/** ---------------- Path Printing ----------------- **/
void GTree::printAllPaths(Node *node, int len, short label)
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
	for(uint i = 0; i < val.length(); i++)
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
