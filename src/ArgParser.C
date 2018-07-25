/********************************************************************
 * Filename: ArgParser.C [C++ header code]
 *
 * Description: Implementation for command line argument parsing class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/ArgParser.H"

ArgParser::ArgParser(int& argc, char** argv) 
{
	this->argc = argc;
	for(int i = 1; i < argc; i++)
		args.push_back(std::string(argv[i]));
	this->argc--;
}
int ArgParser::fillArgs(struct CMDArgs &argList)
{
	for(int i = 0; i < argc; i++) {
		if(args[i] == "-f") {
			i++;
			if(inFileCheck(args[i]))
				return IN_FILE_ERROR;
			argList.iFilename = args[i];
		} else if(args[i] == "-s") {
			i++;
			if(outFileCheck(args[i]))
				return OUT_FILE_ERROR;
			argList.oFilename = args[i];
			argList.storeToFile = 1;
		} else if (args[i] == "-l") {
			i++;
			if(inFileCheck(args[i]))
				return IN_FILE_ERROR;
			argList.lFilename = args[i];
			argList.loadFile = 1;
		} else if(args[i] == "-o"){
			argList.printToScreen = 1;
		} else {
			return USAGE_ERROR;
		}
	}
	return 0;
}
bool ArgParser::argExists(const std::string &arg)
{
	for(int i = 0; i < argc; i++)
		if(args[i] == arg)
			return true;
	return false;
}
const std::string ArgParser::stringOption(const std::string &arg)
{
	for(int i = 0; i < argc; i++) 
		if(args[i] == arg && i + 1 < argc)
			return args[i + 1];
	return "NONE";
}
int ArgParser::inFileCheck(std::string filename)
{
	std::ifstream testFile(filename);
	if(testFile.good())
		return 0;
	else
		return IN_FILE_ERROR;
}
int ArgParser::outFileCheck(std::string filename)
{
	std::ofstream testFile(filename);
	if(testFile.good())
		return 0;
	else
		return OUT_FILE_ERROR;
}
