/********************************************************************
 * Filename: InputFile.C [C++ source code]
 *
 * Description: Implementation of InputFile class
 *
 * Author: Primary - Benjamin Carrington
 *		   Secondary - Dr. Ben Mora
 *
 * Organisation: Swansea University
 * Copyright (c) 2018, Benjamin Carrington, all rights reserved
 *
 *******************************************************************/
#include "../includes/InputFile.H"

bool InputFile::readFastQ()
{
	std::string seqLine = "";
	std::string qualLine = "";
	int64_t i = 3, j = 1;
	short len = 0;

	std::ifstream inpFile( filename.c_str() );
	if( !inpFile ) {
		return 0;
	}

	// Bit of a hacky way of doing this...
	for( std::string line; std::getline( inpFile, line ); i++, j++ ) {
		if( i == 3 )
			readLength = line.length();
		if( !( i % 4 ) ) {
			len = line.find( 'N' );
			line = line.substr( 0, len );

			seqLine = line;
		}
		if( !(j % 4) ) {
			line = line.substr( 0, len );

			qualLine = line;
			GTH::seqReads.push_back( SeqRead( seqLine, qualLine, phredBase ) );

			// NEW
			GTH::startOccs[BASE_IND( seqLine[0] )]++;
			GTH::startWeights[BASE_IND( seqLine[0] )] += GTH::phredQuals[GTH::scoreSys][qualLine[0] - phredBase];
		}
	}

	nReads = ( i - 3 ) / 4;

	return 1;
}
