#ifdef WIN32
#pragma warning ( disable: 4786 )
#endif
#include "../header/GisBinFile.h"
#include "zlib.h"
//#include <cstring>

// -- Constructor

#include <cstring>
GisBinFile::GisBinFile ( const string& name, const char* mode ):
	fileName_ ( name ), mode_ ( mode ), swappMode_ (false),
	wordSize_ (2), isInteger_ (true), compressed_ (false)
{
	const char* filename = fileName_.c_str();
	if ( mode_ == string("r") )
		file_.open ( filename, ios::in | ios::binary );
	else
		file_.open ( filename, ios::out | ios::binary );
//	if ( ! file_.good() )
}

void
GisBinFile::setEndian (const char* mode)
{
	swappMode_ = true;
	if ( mode == string("little") )
	{
		if ( Gis_is_little_endian() ) // Is hardware little endian?
			swappMode_ = false;
	}
	else if ( ! Gis_is_little_endian() )// mode not Big and system
		swappMode_ = false;				// Big Endian
}

bool
GisBinFile::read(short* shortValue)
{
	if ( file_.good() )
	{
		if ( ( wordSize_ == 2 ) && isInteger_ )
		{
			if ( file_.read((char*)shortValue,wordSize_) )
			{
				if (swappMode_)
					*shortValue = ( ((*shortValue >> 8) & 0x00FF) +
								   ((*shortValue << 8) & 0xFF00) );
			}
			return true;
		}
	}
	return false;
}

bool
GisBinFile::read4Bytes(char* char4Values)
{
	if (swappMode_)
	{
		if ( file_.get(tempStore_[3]) )
			if ( file_.get(tempStore_[2]) )
				if ( file_.get(tempStore_[1]) )
					if ( file_.get(tempStore_[0]) )
					{
						memcpy (char4Values, tempStore_, 4);
						return true;
					}
	}
	else if ( file_.read(char4Values,4) )
		return true;
	return false;
}

bool
GisBinFile::read(int* intValue)
{
	if ( file_.good() )
		if ( wordSize_ == 4 )
			return this->read4Bytes((char*) intValue);
	return false;
}

bool
GisBinFile::read(float* floatValue)
{
	if ( file_.good() )
		if ( wordSize_ == 4 )
			return this->read4Bytes((char*) floatValue);
	return false;
}

bool
GisBinFile::readNChar (char* charValues, int nValues)
{
	if ( file_.good() )
		if ( file_.read ( charValues, nValues ) )
			return true;
	return false;
}

void
GisBinFile::gotoPos(long newPos)
{
	file_.seekg( newPos );
}

bool
GisBinFile::readRow(int row, float* floatValues)
{
	if ( file_.good() )
	{
		this->gotoPos( 0L );
		if (this->isCompressed())
			return this->readCompressdRow(row, floatValues);
		else
		{
			char nBytes = 4;
			int totBytes = this->nCols()*nBytes;
			this->gotoPos( row*totBytes );
			unsigned char* readValues = new unsigned char[totBytes];
			this->readNChar ((char *)readValues, totBytes);

			int i = 0;
			int j = 0;
			while ( i < totBytes )
			{
				char fvalchar[4];
				for ( int k = 0; k < 4 ; k++)
				{
					if (swappMode_)
						fvalchar[3-k] = readValues[i++];
					else
						fvalchar[k] = readValues[i++];
				}
				floatValues[j++] = *((float*)&fvalchar[0]);
			}
			delete[] readValues;
			return true;
		}
	}
	return false;
}

bool
GisBinFile::readCompressdRow(int row, float* floatValues)
{
	char nBytes;
	if ( file_.get(nBytes) )
	{
		if ( (int)nBytes == 4 )
		{
			int totBytes = this->nCols()*nBytes;
			this->gotoPos( row*nBytes + 1L );
			int rowPtr, nextRowPtr;
			this->read (&rowPtr);
			this->read (&nextRowPtr);
			this->gotoPos( rowPtr );
			unsigned char* expandedValues = new unsigned char[totBytes];
			char compressFlag;
			file_.get(compressFlag);
			if ( compressFlag == 0x31)
			{
				int rowSize = nextRowPtr - rowPtr -1;
				unsigned char* charValues = new unsigned char[rowSize];
				this->readNChar ((char *)charValues, rowSize);
				if ( ! GisUncompress ( rowSize, charValues, totBytes, expandedValues ) )
				{
					delete[] charValues;
					return false;
				}
				delete[] charValues;
			}
			else
				this->readNChar ((char *)expandedValues, totBytes);
			int i = 0;
			int j = 0;
			while ( i < totBytes )
			{
				char fvalchar[4];
				for ( int k = 0; k < 4; k++)
				{
					if (swappMode_)
						fvalchar[3-k] = expandedValues[i++];
					else
						fvalchar[k] = expandedValues[i++];
				}
				floatValues[j++] = *((float*)&fvalchar[0]);
			}
			delete[] expandedValues;
			return true;
		}
	}
	return false;
}

bool
GisUncompress ( int inSize, unsigned char* inChars, int outSize, unsigned char* outChars )
{
	z_stream streamctrl;
	streamctrl.zalloc = (alloc_func)0;
	streamctrl.zfree  = (free_func)0;
	streamctrl.opaque = (voidpf)0;
	streamctrl.avail_in  = inSize;
	streamctrl.next_in   = inChars;
	int decompSize = outSize;
	streamctrl.avail_out = decompSize;
	streamctrl.next_out  = outChars;

	int err = inflateInit (&streamctrl); //zlib init
	if (err == Z_OK)
	{
		err = inflate (&streamctrl, Z_FINISH);
		int availBytes = decompSize - streamctrl.avail_out;
		if (!(err == Z_STREAM_END || err == Z_OK))
		{	//Some error
			if (!(err == Z_BUF_ERROR && availBytes == decompSize))
				inflateEnd (&streamctrl);
		}
		else
		{ //No error
            inflateEnd (&streamctrl);
			return true;
		}
	}
	return false;
}

bool Gis_is_little_endian()
{
	union
	{
		int anyInt;
		char anyChar[sizeof(int)];
	} testUnion;
    testUnion.anyInt = 1;
	if (testUnion.anyChar[0] == 1)
		return true;
	else
		return false;
}

