#ifndef GisRasterHdr_H
#define GisRasterHdr_H

#include <string>
#include <iostream>
#include <fstream>
using namespace std;

class GisRasterHdr
{
public:

	GisRasterHdr(const string& name);
	
	virtual ~GisRasterHdr(){} 

	bool isCompressed()
	{ return (_compressed == 1); }

	int Rows()
	{ return _rows; }

	int Cols()
	{ return _cols; }

	double XRes()
	{ return _ewresol;}

	double YRes()
	{ return _nsresol;}

	double North()
	{ return _north;}

	double South()
	{ return _south;}

	double East()
	{ return _east;}

	double West()
	{ return _west;}

	bool good()
	{ return _status; }
	void print();

protected:

	int _projId; //0=XY,1=UTM,2=SP(spherical?),3=LL(LatLong),99=other
	int _zoneId; //UTM zone: from (-180 + (zone-1)*6) to (-180 + zone*6)
	double _north;
	double _south;
	double _east;
	double _west;
	int _cols;
	int _rows;
	double _ewresol; //E-W direction resolution - Cell size
	double _nsresol; //N-S direction resolution - Cell size
	int _formatId;   // Don't care
	int _compressed; // 0 = uncompressed, 1 = compressed
	bool _status; // true = OK, 1false = some error

private:
	
// No copy allowed
	GisRasterHdr(const GisRasterHdr&);
	GisRasterHdr& operator=(const GisRasterHdr&);
};

#endif
