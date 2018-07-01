#include <fstream>
#include "includes/ArgParser.H"
#include "includes/InputFile.H"
#include "includes/GTree.H"

#define USAGE_ERROR -1

void argHelp()
{
	std::cout << "Command usage: \n\n"
		"\t""./TreeNome [-f ./path/to/file]"
		"Where:\n\n"
		"\t""-f <string>\n"
		"\t Path to input file (fastq format)" << std::endl;
}

int main(int argc, char** argv)
{
	ArgParser argParser(argc, argv);

	if(argParser.argExists("-h")) {
		argHelp();
		return USAGE_ERROR;
	}
	std::string iFilename = "";
	if(argParser.argExists("-f")) {
		iFilename = argParser.stringOption("-f");
	} else {
		std::cout << "Input file not supplied" << std::endl;
		argHelp();
		return FILE_ERROR;
	}

	TreeTop treeTop(iFilename);
	bool readSuccess = treeTop.readSuccess();
	if(!readSuccess) {
		std::cout << "Input file could not be opened" << std::endl;
		return FILE_ERROR;
	} else {
		std::cout << treeTop.nReads << " records found\n" << std::endl;
	}

	treeTop.processReadsOne();
	treeTop.trees[0].printAllPaths(0);

	return 0;
}
