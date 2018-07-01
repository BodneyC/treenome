#include "includes/GTree.H"
#include "includes/InputFile.H"

TreeTop::TreeTop(const std::string &_filename): 
	filename(_filename), nReads(0), readLength(0), ioSuccess(0)
{
	InputFile iFile(filename);
	ioSuccess = iFile.readFastQ(reads, quals);
	nReads = iFile.nReads;
	readLength = iFile.readLength;

	for(int i = 0; i < NBASES; i++)
		trees[i].createRoot(i, reads, quals);
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
			if(i > 0 && tmpNode) {
				// Because the string causing the previous node (tmpNode) to be
				//	created could have come from anywhere both reads and 
				//	quals are needed
				//trees[BASE_IND(curRead[j])].balanceNode(tmpNode, reads, quals);
			}
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
