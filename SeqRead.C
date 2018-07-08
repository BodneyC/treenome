/********************************************************************
 * Filename: SeqRead.C [C++ source code]
 *
 * Description: Implementation of SeqRead class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "includes/SeqRead.H"

// Order important here, matches BASE_IND()
char SeqRead::ind2base[4] = { 'A', 'C', 'T', 'G' };

SeqRead::SeqRead(const std::string &read, std::string qual)
{
	qualities = qual;

	for(uint i = 0; i < read.length(); i++)
		pushBase(BASE_IND(read[i]));
}

void SeqRead::pushBase(char base)
{
	sequence.push_back(base & 0x02);
	sequence.push_back(base & 0x01);
}

short SeqRead::getBaseInd(short offset)
{
	offset *= 2;
	short retVal = sequence[offset];
	retVal <<= 1;
	retVal += sequence[offset + 1];

	return retVal;
}

char SeqRead::getCharBase(short ind)
{
	ind = getBaseInd(ind);
	return ind2base[ind];
}

char SeqRead::getQual(short offset)
{
	return qualities[offset];
}

short SeqRead::getInd(short ind)
{
	return ind2base[ind];
}
