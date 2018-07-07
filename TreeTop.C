#include "includes/GTree.H"
#include "includes/InputFile.H"

TreeTop::TreeTop(const std::string &_filename): 
	filename(_filename), sequence(""), nReads(0), readLength(0), ioSuccess(0)
{
	InputFile iFile(filename);
	ioSuccess = iFile.readFastQ();
	nReads = iFile.nReads;
	readLength = iFile.readLength;

	for(int i = 0; i < NBASES; i++)
		trees[i].createRoot(i);
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
		if((unsigned) offset == sequence.length() - 1) {
			sequence += 'N';
			maxPath();

			offset++;
		}

		if(trees[BASE_IND(sequence[offset])].getRoot()->occs > 0)
			trees[BASE_IND(sequence[offset])].addToSeq(offset, sequence);
		offset++;
	}
	std::cout << sequence << std::endl;
}

/** --------------- Read Processing ---------------- **/
void TreeTop::processReadsOne()
{
	for(unsigned long i = 0; i < GTH::seqReads.size(); i++)
		for(short j = 0; j < GTH::seqReads[i].size(); j++)
			trees[GTH::seqReads[i].getBaseInd(j)].addReadOne(i, j);
}

/** --------------- Misc Functions ----------------- **/
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

