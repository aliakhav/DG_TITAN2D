#ifdef WIN32
#pragma warning ( disable: 4786 )
#endif
#include "../header/GisAscFile.h"

// -- Constructor

GisAscFile::GisAscFile ( const string& name, const char* mode ):
	mode_ ( mode )
{
	const char* filename = name.c_str();
	long rw_flag = 0;
	if ( mode_ == string("r") )
		//file_.open ( filename, ios_base::in );
		file_.open ( filename, ios::in );
	else
		//file_.open ( filename, ios_base::out );
		file_.open ( filename, ios::out );
//	if ( ! file_.good() )
}

void
GisAscFile::rewind()
{
	if ( file_.good() )
	{
		if ( mode_ == string("r") ) //valid only for reading ??
			file_.seekg(0, ios::beg );
	}
}

