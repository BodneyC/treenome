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
 * - Pre-processing of reads, consider qualities in this
 * - Either move to TCLAP or update ArgParser a little to fill CMDArgs
 * - Quality vs. occurrences in sequence creation
 * - put() function in SeqRead
 * - Output file for printTrees() and printSequence()
 */
#include <fstream>
#include "includes/ArgParser.H"
#include "includes/InputFile.H"
#include "includes/GTree.H"

#define USAGE_ERROR -1

struct CMDArgs {
	std::string iFilename;
	bool printToScreen;
	bool preProcess;
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
	if(argParser.argExists("-p"))
		argList.printToScreen = 1;
	if(argParser.argExists("-c"))\
		argList.preProcess = 1;

	InputFile inpFile(argList.iFilename);
	if(!inpFile.readFastQ()) {
		std::cout << "[ERR]: Input file could not be opened" << std::endl;
		return FILE_ERROR;
	} else {
		std::cout << inpFile.nReads << " records found\n" << std::endl;
	}

	TreeTop treeTop;
	//if(argList.preProcess)
	//	treeTop.preProcess();
	treeTop.processReadsOne();
	if(argList.printToScreen)
		treeTop.printTrees();
	treeTop.buildSequence();
	treeTop.printSequence();
	//if(argList.printToScreen)
	//treeTop.printTrees();

	return 0;
}
