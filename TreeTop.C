#include "includes/GTree.H"

void TreeTop::processReadsFullClean()
{
	for(long i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		curRead = curRead.substr(0, curRead.find('N'));
		int curLength = curRead.length();
		for(int j = 0; j < curLength; j++) {
			// Little bithack converting ACTG -> 0123
			trees[(curRead[0] & 0xF) >> 1].addReadFull(curRead);
			curRead.erase(0, 1);
		}
		for(int j = 0; j < NBASES; j++)
			trees[j].cleanBranches();
	}
}

void TreeTop::processReadsFull()
{
	for(long i = 0; i < nReads; i++) {
		std::string curRead = reads[i];
		curRead = curRead.substr(0, curRead.find('N'));
		int curLength = curRead.length();
		for(int j = 0; j < curLength; j++) {
			// Little bithack converting ACTG -> 0123
			trees[(curRead[0] & 0xF) >> 1].addReadFull(curRead);
			curRead.erase(0, 1);
		}
	}
}

void processReadsOne()
{

}
