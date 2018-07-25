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
 * - Data race or algorithm?
 * - Derive leafNode from Node (somehow)
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
		"\t  ./TreeNome [-f ./path/to/inFile] [-o] [-s ./path/to/outFile]"
		"\nWhere:\n\n"
		"\t  -f <string>\n"
		"\t   Path to input file (fastq format)\n\n"
		"\t  -o\n"
		"\t   Output to screen\n\n"
		"\t  -s <string>\n"
		"\t   Store tree to file\n" << std::endl;
}

void writeTreesToDisk(std::string oFilename, TreeTop& treeTop)
{
	std::ofstream outFile(oFilename);
	for(int i = 0; i < NBASES; i++)
		outFile << treeTop.treeStrings[i].c_str() << std::endl;
}

int main(int argc, char** argv)
{
	omp_set_num_threads(NUM_THREADS);
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
