#include "includes/GTree.H"
#include "includes/InputFile.H"

TreeTop::TreeTop(const std::string &_filename): 
	filename(_filename), sequence(""), nReads(0), readLength(0), ioSuccess(0)
{
	InputFile iFile(filename);
	ioSuccess = iFile.readFastQ(reads, quals);
	nReads = iFile.nReads;
	readLength = iFile.readLength;

	for(int i = 0; i < NBASES; i++)
		trees[i].createRoot(i, reads, quals);
}

short TreeTop::maxPath()
{
	int start = 0;
	ulong maxOccs = 0;
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

	ulong offset = 1;
	while(offset < 20 /*trees[start].getRoot()->occs*/) {

		// Possibly make is a tighter gap as its working from single letters
		if(offset == sequence.length() - 1) {
			sequence += 'N';
			maxPath();
			offset++;
		}

		std::cout << "T/O/L: " << sequence[offset] << "/" << 
			offset << "/" << sequence.length() << std::endl;
		std::cout << "On Sequence: " << sequence << std::endl;

		trees[BASE_IND(sequence[offset])].addToSeq(offset, sequence);
		offset++;
	}
}

void TreeTop::processReadsOne()
{
	for(ulong i = 0; i < nReads; i++)
		for(short j = 0; j < (short)reads[i].length(); j++)
			trees[BASE_IND(reads[i][j])].addReadOne(i, j, reads, quals);
}

// The following is by no means DRY but each exists for experimental purpose,
//	realistically only one would exist and be called processReads()
void TreeTop::processReadsFullCleanNRBalanced()
{
	for(ulong i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		std::string curQual = quals[i];

		// Add the full string to the tree tip to toe
		short curLength = curRead.length();
		for(short j = 0; j < curLength; j++)
			trees[BASE_IND(curRead[j])].addReadFull(i, j, curRead, curQual);

		for(short j = 0; j < curLength; j++) {
			// Delete linked list of ones for memories sake
			Node *tmpNode = trees[BASE_IND(curRead[j])].cleanBranchesNR(j, curRead);
			// Because the string causing the previous node (tmpNode) to be
			//	created could have come from anywhere both reads and 
			//	quals are needed
			trees[BASE_IND(curRead[j])].balanceNode(tmpNode, reads, quals);
		}
	}
}

void TreeTop::processReadsFullCleanNR()
{
	for(ulong i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		std::string curQual = quals[i];

		short curLength = curRead.length();
		for(short j = 0; j < curLength; j++)
			trees[BASE_IND(curRead[j])].addReadFull(i, j, curRead, curQual);

		for(short j = 0; j < curLength; j++)
			trees[BASE_IND(curRead[j])].cleanBranchesNR(j, curRead);
	}
}

void TreeTop::processReadsFullClean()
{
	for(ulong i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		std::string curQual = quals[i];

		short curLength = curRead.length();
		for(short j = 0; j < curLength; j++)
			trees[BASE_IND(curRead[j])].addReadFull(i, j, curRead, curQual);

		for(short j = 0; j < NBASES; j++)
			trees[j].cleanBranches();
	}
}

void TreeTop::processReadsFull()
{
	for(ulong i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		std::string curQual = quals[i];
		//curRead = curRead.substr(0, curRead.find('N'));
		//curQual = curQual.substr(0, curRead.find('N'));
		short curLength = curRead.length();
		for(short j = 0; j < curLength; j++)
			trees[BASE_IND(curRead[j])].addReadFull(i, j, curRead, curQual);
	}
}
