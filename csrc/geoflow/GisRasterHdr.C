#include "../header/GisRasterHdr.h"
#include "../header/GisAscFile.h"
#include <cstring>

GisRasterHdr::GisRasterHdr ( const string& name )
{
	_status = false;
	GisAscFile headerFile( name );
	if (headerFile.good())
	{
		char charText[20];
		headerFile.getline(charText, 20, ':');
		if ( strcmp(charText, "proj") == 0)
		{ //At least test first line!!!
			headerFile.getAscInt ( _projId );
			headerFile.getline(charText, 20, ':');
			headerFile.getAscInt ( _zoneId );
			headerFile.getline(charText, 20, ':');
			headerFile.getAscDouble (_north);
			headerFile.getline(charText, 20, ':');
			headerFile.getAscDouble (_south);
			headerFile.getline(charText, 20, ':');
			headerFile.getAscDouble (_east);
			headerFile.getline(charText, 20, ':');
			headerFile.getAscDouble (_west);
			headerFile.getline(charText, 20, ':');
			headerFile.getAscInt (_cols);
			headerFile.getline(charText, 20, ':');
			headerFile.getAscInt (_rows);
			headerFile.getline(charText, 20, ':');
			headerFile.getAscDouble (_ewresol);
			headerFile.getline(charText, 20, ':');
			headerFile.getAscDouble (_nsresol);
			headerFile.getline(charText, 20, ':');
			headerFile.getAscInt (_formatId);
			headerFile.getline(charText, 20, ':');
			headerFile.getAscInt (_compressed);
			_status = true;
		}
	}
}

void
GisRasterHdr::print()
{
	cout << "projId " << _projId << "\n";
	cout << "zoneId " << _zoneId << "\n";
	cout << "north " << _north << "\n";
	cout << "south " << _south << "\n";
	cout << "east " << _east << "\n";
	cout << "west " << _west << "\n";
	cout << "cols " << _cols << "\n";
	cout << "rows " << _rows << "\n";
	cout << "ewres " << _ewresol << "\n";
	cout << "nsres " << _nsresol << "\n";
	cout << "formatId " << _formatId << "\n";
	cout << "compressed " << _compressed << "\n";
}

