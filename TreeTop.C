#include "includes/GTree.H"
#include "includes/InputFile.H"

TreeTop::TreeTop(std::string filename): nReads(0), readLength(0), ioSuccess(0)
{
	InputFile iFile(filename);
	ioSuccess = iFile.readFastQ(reads, qualities);
	nReads = iFile.nReads;
	readLength = iFile.readLength;

	for(int i = 0; i < NBASES; i++)
		trees[i].createRoot(reads, i);
}

void TreeTop::processReadsOne()
{

}

// The following is by no means DRY but each exists for experimental purpose,
//		realistically only one would exist and be called processReads()
void TreeTop::processReadsFullCleanNRBalanced()
{
	for(ulong i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		std::string curQual = qualities[i];
		curRead = curRead.substr(0, curRead.find('N'));
		curQual = curQual.substr(0, curRead.find('N'));

		int curLength = curRead.length();
		for(int j = 0; j < curLength; j++)
			// Little bithack converting ACTG -> 0123
			trees[(curRead[j] & 0xF) >> 1].addReadFull(i, j, curRead, curQual);

		for(int j = 0; j < curLength; j++) {
			trees[(curRead[0] & 0xF) >> 1].cleanBranchesNR(curRead);
			curRead.erase(0, 1);
		}

		if(i > 0) {
			// Attempt to balance tree
		}
	}
}

void TreeTop::processReadsFullCleanNR()
{
	for(ulong i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		std::string curQual = qualities[i];
		curRead = curRead.substr(0, curRead.find('N'));
		curQual = curQual.substr(0, curRead.find('N'));

		int curLength = curRead.length();
		for(int j = 0; j < curLength; j++)
			// Little bithack converting ACTG -> 0123
			trees[(curRead[j] & 0xF) >> 1].addReadFull(i, j, curRead, curQual);

		for(int j = 0; j < curLength; j++) {
			trees[(curRead[0] & 0xF) >> 1].cleanBranchesNR(curRead);
			curRead.erase(0, 1);
		}
	}
}

void TreeTop::processReadsFullClean()
{
	for(ulong i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		std::string curQual = qualities[i];
		curRead = curRead.substr(0, curRead.find('N'));
		curQual = curQual.substr(0, curRead.find('N'));

		int curLength = curRead.length();
		for(int j = 0; j < curLength; j++)
			// Little bithack converting ACTG -> 0123
			trees[(curRead[j] & 0xF) >> 1].addReadFull(i, j, curRead, curQual);

		for(int j = 0; j < NBASES; j++)
			trees[j].cleanBranches();
	}
}

void TreeTop::processReadsFull()
{
	for(ulong i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		std::string curQual = qualities[i];
		curRead = curRead.substr(0, curRead.find('N'));
		curQual = curQual.substr(0, curRead.find('N'));
		int curLength = curRead.length();
		for(int j = 0; j < curLength; j++)
			// Little bithack converting ACTG -> 0123
			trees[(curRead[j] & 0xF) >> 1].addReadFull(i, j, curRead, curQual);
	}
}
