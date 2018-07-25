/********************************************************************
 * Filename: TreeNome.C [C++ source code]
 *
 * Description: Entry point for TreeNome de novo DNA assembler
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
/* TODO:
 * - Double accuracy?
 * - Reconstruct tree from file (possible second program)
 * - Quality vs. occurrences in sequence creation
 */
#include "../includes/ArgParser.H"
#include "../includes/InputFile.H"
#include "../includes/TreeTop.H"

void argHelp()
{
	std::cout << "\nCommand usage: \n\n"
		"\t  ./TreeNome [-f ./path/to/inFile] [-o] [-s ./path/to/outFile]"
		"\nWhere:\n\n"
		"\t  -f <string>\n"
		"\t   Path to fastq input file\n\n"
		"\t  -h\n"
		"\t   Print command line help\n\n"
		"\t  -l <string>\n"
		"\t   Load tree from file (as outputted by TreeNome)\n\n"
		"\t  -o\n"
		"\t   Output to screen\n\n"
		"\t  -s <string>\n"
		"\t   Store tree to file\n"
		<< std::endl;
}

void writeTreesToDisk(std::string oFilename, TreeTop& treeTop)
{
	std::ofstream outFile(oFilename);
	for(int i = 0; i < NBASES; i++)
		outFile << treeTop.treeStrings[i].c_str() << std::endl;
}

signed int createTreeBuildSequence(struct CMDArgs& argList)
{
	InputFile inpFile(argList.iFilename);
	if(!inpFile.readFastQ()) {
		std::cout << "[ERR]: Input file could not be opened" << std::endl;
		return FILE_ERROR;
	} else {
		std::cout << inpFile.nReads << " records found\n" << std::endl;
	}

	TreeTop treeTop;
	treeTop.processReadsOne();
	std::cout << "\n----------" << std::endl;
	if(argList.printToScreen)
		treeTop.printTrees();
	if(argList.storeToFile) {
		treeTop.storeTrees();
		writeTreesToDisk(argList.oFilename, treeTop);
	}
	treeTop.buildSequence();
	treeTop.printSequence();

	return 0;
}

signed int loadTreeFromFile(struct CMDArgs& argList)
{

	return USAGE_ERROR;
}

int main(int argc, char** argv)
{
	omp_set_num_threads(NUM_THREADS);
	ArgParser argParser(argc, argv);
	struct CMDArgs argList;
	int argSucc = argParser.fillArgs(argList), progSucc;

	if(argSucc)
		argHelp();
	switch (argSucc) {
	case USAGE_ERROR:
		std::cout << "[ERR]: Incorrect command usage" << std::endl;
		return USAGE_ERROR;
	case IN_FILE_ERROR:
		std::cout << "[ERR]: Input file not supplied or not present" << std::endl;
		return IN_FILE_ERROR;
	case OUT_FILE_ERROR:
		std::cout << "[ERR]: Output file not supplied or not present" << std::endl;
		return OUT_FILE_ERROR;
	}

	if(argList.loadFile)
		progSucc = loadTreeFromFile(argList);
	else
		progSucc = createTreeBuildSequence(argList);

	if(progSucc)
		return progSucc;

	return 0;
}
