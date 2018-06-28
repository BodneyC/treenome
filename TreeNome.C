#include <fstream>
#include "includes/ArgParser.H"
#include "includes/InputFile.H"
#include "includes/GTree.H"

int processReadsFullClean(TreeTop &treeTop)
{
	for(long i = 0; i < treeTop.nReads; i++) {
		std::string curRead = treeTop.reads[i];
		curRead = curRead.substr(0, curRead.find('N'));
		int curLength = curRead.length();
		for(int j = 0; j < curLength; j++) {
			// Little bithack converting ACTG -> 0123
			treeTop.trees[(curRead[0] & 0xF) >> 1].addReadFULL(curRead);
			curRead.erase(0, 1);
		}
	}

	return 0;
}

int processReadsFull(TreeTop &treeTop)
{
	for(long i = 0; i < treeTop.nReads; i++) {
		std::string curRead = treeTop.reads[i];
		curRead = curRead.substr(0, curRead.find('N'));
		int curLength = curRead.length();
		for(int j = 0; j < curLength; j++) {
			// Little bithack converting ACTG -> 0123
			treeTop.trees[(curRead[0] & 0xF) >> 1].addReadFULL(curRead);
			curRead.erase(0, 1);
		}
	}

	return 0;
}

int main(int argc, char** argv)
{
	std::string iFilename = "";
	ArgParser argParser(argc, argv);

	if(argParser.argExists("-f")) {
		iFilename = argParser.stringOption("-f");
	} else {
		std::cout << "Input file not supplied" << std::endl;
		return FILE_ERROR;
	}

	InputFile iFile(iFilename);
	TreeTop treeTop;

	iFile.readFastQ(treeTop);
	if(iFile.nReads == FILE_ERROR ) {
		std::cout << "Input file could not be opened" << std::endl;
		return FILE_ERROR;
	} else {
		std::cout << iFile.nReads << " records found\n" << std::endl;
	}


	processReadsFullClean(treeTop);
	treeTop.trees[0].printAllPaths(0);


	return 0;
}
