#ifndef GisBinFile_H
#define GisBinFile_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
//#include <cstring>
using namespace std;

bool GisUncompress (int inSize, unsigned char* inChars, int outSize, unsigned char* outChars );

bool Gis_is_little_endian();

class GisBinFile
{
public:
	GisBinFile( const string& name, const char* mode = "r");

	virtual ~GisBinFile(){};

	void setEndian (const char* mode = "little");

	void setWordSize (int wordSize)
	{ wordSize_ = wordSize; }

	void setIsInteger (bool trueFalse)
	{ isInteger_ = trueFalse; }

	bool readRow(int row, float* floatValues);
	bool readCompressdRow(int row, float* floatValues);

	bool read (float* floatValue);	//IEEE 32 bits

	bool read (int* intValue);		// 32 bits

	bool read (short* shortValue); //16 bits

	bool read4Bytes(char* char4Values);

	bool readNChar (char* charValues, int nValues);

	bool good()
	{ return file_.good(); }

	int nRows ()
	{ return nRows_; }

	int nCols ()
	{ return nCols_; }

	void nRows (int nrows)
	{ nRows_ = nrows; }

	void nCols (int ncols)
	{ nCols_ = ncols; }

	bool isCompressed()
	{ return ( compressed_ == 1); }

	void isCompressed( bool compFlag )
	{ compressed_ = compFlag; }

	void gotoPos(long newPos);

protected:
// -- File pointer
	fstream file_;
	string fileName_;
	string mode_;
	bool swappMode_;
	int wordSize_;
	bool isInteger_;
	char tempStore_[4];
	int nRows_;
	int nCols_;
	bool compressed_;

private:
// No copy allowed
	GisBinFile(const GisBinFile&);
	GisBinFile& operator=(const GisBinFile&);
};
#endif
