/********************************************************************
 * Filename: SeqRead.hpp [C++ header code]
 *
 * Description: Declaration of SeqRead class for 2-bit storage of 
 *				bases
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 *******************************************************************/
#ifndef _SEQ_READ_
#define _SEQ_READ_

// Little bithack converting ACTG -> 0123
#define BASE_IND( x ) ((( x ) & 0xF ) >> 1 )

#include <vector>
#include <string>

class SeqRead {
public:
	SeqRead(): qualities( "" ) {  }
	explicit SeqRead( const std::string &read, std::string qual, int pB );

	short size() { return sequence.size() / 2; }

	void push_base( char base );
	short get_base_idx( short offset );
	char get_char_base( short ind );

	char get_quality( short offset );

private:
	std::vector<bool> sequence;
	std::string qualities;
	static char idx_to_base[4];
};

#endif /*_SEQ_READ_*/
