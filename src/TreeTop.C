/********************************************************************
 * Filename: TreeTop.C [C++ source code]
 *
 * Description: Implementation of TreeTop class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "includes/GTree.H"
#include "includes/InputFile.H"

TreeTop::TreeTop(): 
	sequence(""), nReads(0), readLength(0)
{
	for(int i = 0; i < NBASES; i++)
		trees[i].createRoot(i);
}

/** --------------- Read Processing ---------------- **/
void TreeTop::threadFunc(unsigned long i)
{
	for(short j = 0; j < GTH::seqReads[i].size(); j++)
		//if(GTH::seqReads[i].getBaseInd(j) == 3)
			trees[GTH::seqReads[i].getBaseInd(j)].addReadOne(i, j);
}

void TreeTop::processReadsOne()
{
#pragma omp parallel num_threads(NUM_THREADS)
{
	for(unsigned long i = 0; i < GTH::seqReads.size(); i += NUM_THREADS) {
	//std::cout << i << std::endl;
#pragma omp for schedule(static, 1)
	for(int j = 0; j < NUM_THREADS; j++) 
		if(i + j < GTH::seqReads.size())
			threadFunc(i + j);
	}
}
}

/** ------------- Sequence Generation -------------- **/
bool TreeTop::rootOccsExist()
{
	short ret = 0;

	for(int i = 0; i < NBASES; i++)
		if(trees[i].getRoot()->occs == 0)
			ret++;

	return ret == NBASES ? 0 : 1;
}

short TreeTop::maxPath()
{
	int start = 0;
	long maxOccs = 0;
	Node* treeRoots[NBASES];

	for(short i = 0; i < NBASES; i++) {
		treeRoots[i] = trees[i].getRoot();
		if(treeRoots[i]->occs > maxOccs) {
			maxOccs = treeRoots[i]->occs;
			start = i;
		}
	}

	trees[start].followPath(trees[start].getRoot(), start, sequence);

	return start;
}

void TreeTop::buildSequence()
{
	maxPath();

	long offset = 1;
	while(rootOccsExist()) {

		// Possibly make is a tighter gap as its working from single letters
		// (this would actually be the k-mer match)
		if((unsigned) offset == sequence.length() - 1) {
			sequence += 'N';
			maxPath();
			offset += 2;
		}

		if(trees[BASE_IND(sequence[offset])].getRoot()->occs > 0)
			trees[BASE_IND(sequence[offset])].addToSeq(offset, sequence);
		offset++;
	}
}

/** --------------- Misc Functions ----------------- **/
void TreeTop::storeTrees()
{
	for(int i = 0; i < NBASES; i++)
		treeStrings[i] = trees[i].storeTree(i);
}

void TreeTop::printTrees()
{
	for(int i = 0; i < NBASES; i++)
		trees[i].printAllPaths(i);
}

void TreeTop::printSequence()
{
	// 80 for terminal width's sake
	unsigned short TWIDTH = 80;
	unsigned long i = 0;

	if(sequence.length() > TWIDTH)
		for(; i < sequence.length() - TWIDTH; i += TWIDTH)
			std::cout << sequence.substr(i, TWIDTH) << std::endl;

	std::cout << sequence.substr(i, sequence.length() - i) << std::endl;
}

