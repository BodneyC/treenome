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
 * - Quality vs. occurrences in sequence creation
 * - put() function in SeqRead
 * - Output file for printTrees() and printSequence()
 */
#include "includes/ArgParser.H"
#include "includes/InputFile.H"
#include "includes/GTree.H"

void argHelp()
{
	std::cout << "\nCommand usage: \n\n"
		"\t""./TreeNome [-f ./path/to/file] [-o] [-p]"
		"\nWhere:\n\n"
		"\t""-f <string>\n"
		"\t   Path to input file (fastq format)\n\n"
		"\t""-o\n"
		"\t   Output to screen\n\n"
		"\t""-p\n"
		"\t   Perform pre-processing\n" << std::endl;
}

int main(int argc, char** argv)
{
	ArgParser argParser(argc, argv);
	struct CMDArgs argList;
	int argSucc = argParser.fillArgs(argList);

	if(argSucc == USAGE_ERROR) {
		argHelp();
		return USAGE_ERROR;
	}
	if(argSucc == FILE_ERROR) {
		std::cout << "[ERR]: Input file not supplied" << std::endl;
		argHelp();
		return USAGE_ERROR;
	}

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
