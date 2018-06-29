#include <fstream>
#include "includes/ArgParser.H"
#include "includes/InputFile.H"
#include "includes/GTree.H"

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

	treeTop.processReadsFullCleanNR();
	treeTop.trees[0].printAllPaths(0);

	return 0;
}
