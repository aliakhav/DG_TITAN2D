#ifndef GisAscFile_H
#define GisAscFile_H
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

class GisAscFile {
public:

	GisAscFile(const string& name, const char* mode = "r");
	
	virtual ~GisAscFile(){} 

	bool good()
	{ return file_.good(); }

	void getline(char *outChar, int nChar, char termChar = '\n') //'\n' = 0x0A
	{ file_.getline(outChar, nChar, termChar); }

	void getAscInt(int &intValue)
	{ file_ >> intValue; }

	void getAscDouble(double &doubleValue)
	{ file_ >> doubleValue; }

	void rewind();

	string mode() { return mode_; }

protected:
	
// -- File pointer

	fstream file_;
	string mode_;

private:
	
// No copy allowed
	GisAscFile(const GisAscFile&);
	GisAscFile& operator=(const GisAscFile&);
};

#endif
