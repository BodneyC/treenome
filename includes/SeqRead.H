/********************************************************************
 * Filename: SeqRead.H [C++ header code]
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

	void pushBase( char base );
	short getBaseInd( short offset );
	char getCharBase( short ind );

	char getQual( short offset );

private:
	std::vector<bool> sequence;
	std::string qualities;
	static char ind2base[4];
};

#endif /*_SEQ_READ_*/
