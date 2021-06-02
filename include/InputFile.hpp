/********************************************************************
 * Filename: InputFile.hpp [C++ header code]
 *
 * Description: Declaration of InputFile class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 *******************************************************************/
#ifndef _INP_FILE_
#define _INP_FILE_

#include "GTree.hpp"

#define FILE_ERROR -2

// Though sparse, this class may be expanded upon for processing/changing
// of storage method
class InputFile {
public:
  long n_reads, read_length;
  int phred_base;

  InputFile(const std::string &_filename, int _pB)
      : n_reads(0), read_length(0), phred_base(_pB), filename(_filename) {}

  bool read_fastq();

private:
  std::string filename;
};

#endif /*_INP_FILE_*/
