#include <iostream>
#include <cstdlib>

#include "FastFDK.hpp"

/*
Distance det: 150 mm
Distance src: 137 mm
img width   : 100
img height  : 100
img depth   : 360
pixel size  : 0.33 mm
vol width   : 100
vol height  : 100
vol depth   : 100
voxel size  : 0.15 mm
*/

using namespace std;



//! display message with command help and quit program 
void printHelpAndQuit( char const * msg )
{
	if( msg != NULL )
		cerr << msg << endl << endl;
	cerr << "Command manual:" << endl << endl;
	cerr << "  FastFDK <images> <pxlSize> <detDist> <srcDist> <volW> <volH> <volume>" << endl << endl;
	cerr << "  images : input ims image stack : -ln(I/I°)" << endl;
	cerr << "  pxlSize: pixel size (mm)" << endl;
	cerr << "  detDist: distance from origin to dector (mm)" << endl;
	cerr << "  srcDist: distance from origin to source (mm)" << endl;
	cerr << "  volW   : volume slice width (vxl)" << endl;
	cerr << "  volH   : number of volume slices" << endl;
	exit(1);
}


int main( int argc, char* argv[] )
{
	try
	{
		// check if help requested
		if( argc == 1 || strcmp( argv[1], "-h" ) == 0 || 
			strcmp( argv[1], "--help" ) == 0 )
			printHelpAndQuit( NULL );
	
		// check number of arguments
		if( argc != 8 )
			printHelpAndQuit( "Invalid number of arguments" );
		
		// load projection images
		ims_float in( argv[1] );

		// get parameter values
		Real pxlSize = atof(argv[2]);
		Real detDist = atof(argv[3]);
		Real srcDist = atof(argv[4]);
		int vw = atoi(argv[5]);
		int vh = atoi(argv[6]);
		
		// compute optimal voxel size
		Real vs = (Real)(2*srcDist*sin(atan((in.width()/2)*pxlSize/
			(detDist+srcDist)))/vw);
		
		cout << "Voxel size: " << vs << " mm" << endl;
		
		// instantiate fdk engine
		FastFDK fdk( in.width(), in.height(), in.depth(), pxlSize, 
			detDist, srcDist, vw, vh, vs, 360 ); 

		// reconstruct and save volume
		fdk.reconstruct( in ).save( argv[7] );
		
	}
	catch( std::exception& e )
	{
		cout << "Exception: " << e.what() << endl;
	}
	catch( ... )
	{
		cout << "Exception" << endl;
	}


	return 0;
}
