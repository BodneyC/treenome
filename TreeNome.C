#include <fstream>
#include "includes/ArgParser.H"
#include "includes/InputFile.H"
#include "includes/GTree.H"

#define USAGE_ERROR -1

struct CMDArgs {
	std::string iFilename;
};

void argHelp()
{
	std::cout << "\nCommand usage: \n\n"
		"\t""./TreeNome [-f ./path/to/file]"
		"Where:\n\n"
		"\t""-f <string>\n"
		"\t Path to input file (fastq format)\n" << std::endl;
}

int main(int argc, char** argv)
{
	ArgParser argParser(argc, argv);
	struct CMDArgs argList;

	if(argParser.argExists("-h")) {
		argHelp();
		return USAGE_ERROR;
	}
	if(argParser.argExists("-f")) {
		argList.iFilename = argParser.stringOption("-f");
	} else {
		std::cout << "[ERR]: Input file not supplied" << std::endl;
		argHelp();
		return FILE_ERROR;
	}

	TreeTop treeTop(argList.iFilename);
	bool readSuccess = treeTop.readSuccess();
	if(!readSuccess) {
		std::cout << "[ERR]: Input file could not be opened" << std::endl;
		return FILE_ERROR;
	} else {
		std::cout << treeTop.nReads << " records found\n" << std::endl;
	}

	for(int i = 0; i < NBASES; i++)
		treeTop.trees[i].printAllPaths(i);
	treeTop.processReadsOne();
	for(int i = 0; i < NBASES; i++)
		treeTop.trees[i].printAllPaths(i);
	treeTop.buildSequence();

	return 0;
}
