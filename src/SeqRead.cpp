/********************************************************************
 * Filename: SeqRead.cpp [C++ source code]
 *
 * Description: Implementation of SeqRead class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Details: Tightly-packed 2-bit-array representing a read where:
 *		| B. | Represents |
 *		| 00 |      A     |
 *		| 01 |      C     |
 *		| 02 |      T     |
 *		| 03 |      G     |
 *
 *******************************************************************/
#include "SeqRead.hpp"

// Order important here, matches BASE_IND()
char SeqRead::idx_to_base[4] = { 'A', 'C', 'T', 'G' };

SeqRead::SeqRead( const std::string &read, std::string qual, int pB )
{
	for( uint32_t i = 0; i < read.length(); i++ ) {
		push_base( BASE_IND( read[i] ) );
		qual[i] -= pB;
	}

	qualities = qual;
}

void SeqRead::push_base( char base )
{
	sequence.push_back( base & 0x02 );
	sequence.push_back( base & 0x01 );
}

short SeqRead::get_base_idx( short offset )
{
	offset *= 2;
	short retVal = sequence[offset];
	retVal <<= 1;
	retVal += sequence[offset + 1];

	return retVal;
}

char SeqRead::get_char_base( short ind )
{
	ind = get_base_idx( ind );
	return idx_to_base[ind];
}

char SeqRead::get_quality( short offset )
{
	return qualities[offset];
}

