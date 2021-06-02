/********************************************************************
 * Filename: InputFile.cpp [C++ source code]
 *
 * Description: Implementation of InputFile class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Details: Serial FASTQ file reader, not at all fast, not the
 *		point of the project
 *
 *******************************************************************/
#include "InputFile.hpp"

bool InputFile::read_fastq() {
  std::string seq_line = "";
  std::string qual_line = "";
  int64_t i = 3, j = 1;
  short len = 0;

  std::ifstream in_file(filename.c_str());
  if (!in_file) {
    return 0;
  }

  // Bit of a hacky way of doing this...
  for (std::string line; std::getline(in_file, line); i++, j++) {
    if (i == 3)
      GTH::read_len = read_length = line.length();
    if (!(i % 4)) {
      len = line.find('N');
      line = line.substr(0, len);

      seq_line = line;
    }
    if (!(j % 4)) {
      line = line.substr(0, len);

      qual_line = line;
      GTH::seq_reads.push_back(SeqRead(seq_line, qual_line, phred_base));

      // NEW
      GTH::start_occs[BASE_IND(seq_line[0])]++;
      GTH::start_weights[BASE_IND(seq_line[0])] +=
          GTH::phred_qualities[GTH::score_sys][qual_line[0] - phred_base];
    }
  }

  n_reads = (i - 3) / 4;

  return 1;
}
