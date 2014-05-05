#ifndef FGENERATEDISTRIBUTION_HPP
#define FGENERATEDISTRIBUTION_HPP



#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
//
#include "Utils/FMath.hpp"
#include "Utils/FParameters.hpp"

/**  return a random number between 0 and 1 */
//double getRandom() {
//	return drand48();
//} ;
void initRandom() {
	srand48( static_cast<long int>(time(0))) ;
} ;
FReal getRandom() {
	return static_cast<FReal>(drand48());
	//return static_cast<FReal>(rand()/FReal(RAND_MAX));
} ;
//!  \fn   unifRandonPointsOnUnitCube(const int N , FReal * points)

//! \brief Generate N points uniformly distributed on the unit cube

//!
//! \param N the number of points uniformly randomly sample on the unit cube
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//! \example  generateDistributions.cpp
void unifRandonPointsOnUnitCube(const int N , FReal * points) {
	//
	initRandom() ;
	int j = 0;
	for (int i = 0 ; i< N ; ++i, j+=4)  {
		//
		points[j]	  =	getRandom()  ;
		points[j+1] =	getRandom()  ;
		points[j+2] =	getRandom()  ;
		//
	}
};
//!  \fn   unifRandonPointsOnCube(const int N , FReal * points)

//! \brief Generate N points uniformly distributed on the cube of length R

//!
//! \param N the number of points uniformly randomly sample on the unit cube
//! \param Lx the the X-length of the  cube
//! \param Ly the the Y-length of the  cube
//! \param Lz the the Z-length of the  cube
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//! \example  generateDistributions.cpp
void unifRandonPointsOnCube(const int N , const FReal Lx,  const FReal Ly,  const FReal Lz, FReal * points) {
	//
	unifRandonPointsOnUnitCube(N , points) ;
	int j =0 ;
	for (int i = 0 ; i< N ; ++i, j+=4)  {
		points[j]	   *= Lx ;
		points[j+1]  *= Ly ;
		points[j+2]  *= Lz ;
	}
};
//!  \fn   unifRandonPointsOnUnitSphere(const int N , FReal * points)

//! \brief Generate N points uniformly distributed on the unit sphere

//!
//! \param N the number of points uniformly randomly sample on the unit sphere
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//! \example  generateDistributions.cpp
void unifRandonPointsOnUnitSphere(const int N , FReal * points) {
	FReal u, v, theta, phi, sinPhi ;
	//
	initRandom() ;
	int j = 0 ;
	for (int i = 0 ; i< N ; ++i, j+=4)  {
		//
		u = getRandom() ;  v = getRandom() ;
		theta  = FMath::FTwoPi*u ;
		phi     = FMath::ACos(2*v-1);
		sinPhi = FMath::Sin(phi);
		//
		points[j]	  =	FMath::Cos(theta)*sinPhi ;
		points[j+1] =	FMath::Sin(theta)*sinPhi ;
		points[j+2] =	2*v-1 ;
		//
	}
};
//!  \fn  nonunifRandonPointsOnElipsoid(const int N , const FReal &a, const FReal &b, const FReal &c, FReal * points)

//! \brief  Generate N points non uniformly distributed on the ellipsoid of  aspect ratio a:b:c

//!
//! \param N the number of points
//! \param a  the x  semi-axe length
//! \param b  the y  semi-axe length
//! \param c  the z  semi-axe length
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//!
void nonunifRandonPointsOnElipsoid(const int N , const FReal &a, const FReal &b, const FReal &c, FReal * points) {
	//
	FReal u, v , cosu ;
	int j =0 ;
	for (int i = 0 ; i< N ; ++i, j+=4)  {
		u = getRandom() ;  v = getRandom() ;
		u  = FMath::FPi*u - FMath::FPiDiv2;   v   = FMath::FTwoPi*v - FMath::FPi;
		cosu = FMath::Cos(u) ;
		points[j]	   = a*cosu*FMath::Cos(v)  ;
		points[j+1]  = b*cosu*FMath::Sin(v)  ;
		points[j+2]  = c*FMath::Sin(u)  ;
	}
};
//!  \fn  nonunifRandonPointsOnElipsoid(const int N , const FReal &a, const FReal &c, FReal * points)

//! \brief  Generate N points uniformly distributed on the ellipsoid of  aspect ratio a:a:c

//!
//! \param N the number of points
//! \param a  the x  semi-axe length
//! \param c  the z  semi-axe length
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//!
void unifRandonPointsOnProlate(const int N , const FReal &a, const FReal &c, FReal * points){
	//
	FReal u, w,v ,ksi ;
	FReal e = (a*a*a*a)/(c*c*c*c) ;
	bool isgood = false;
	int j =0 , cpt =0 ;
	//
	for (int i = 0 ; i< N ; ++i, j+=4)  {
		// Select a random point on the prolate
		do {
			cpt++	;
			u = getRandom() ;  v = getRandom() ;
			u  = 2.0*u - 1.0;   v   = FMath::FTwoPi*v;
			w =FMath::Sqrt(1-u*u) ;
			points[j]	   = a*w*FMath::Cos(v)  ;
			points[j+1]  = a*w*FMath::Sin(v)  ;
			points[j+2]  = c*u ;
			// Accept the position ?
			ksi = a*getRandom()  ;
			//			std::cout << "Gradf  "<<  points[j]*points[j] + points[j+1] *points[j+1]  +e*points[j+2] *points[j+2]  << std::endl;
			isgood = (points[j]*points[j] + points[j+1] *points[j+1]  +e*points[j+2] *points[j+2]  < ksi*ksi );
		} while (isgood);
	}
	std::cout.precision(4);
	std::cout << "Total tested points: "<< cpt << " % of rejected points: "<<100*static_cast<FReal>(cpt-N)/cpt << " %" <<std::endl;

} ;

//!  \fn  unifRandonPointsOnSphere(const int N , const FReal R, FReal * points)

//! \brief Generate N points uniformly distributed on the sphere of radius R

//!
//! \param N the number of points uniformly randomly sample on the sphere
//! \param R the radius of the sphere
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//!
void unifRandonPointsOnSphere(const int N , const FReal R, FReal * points) {
	//
	unifRandonPointsOnUnitSphere(N , points) ;
	int j =0 ;
	for (int i = 0 ; i< N ; ++i, j+=4)  {
		points[j]	   *= R ;
		points[j+1]  *= R ;
		points[j+2]  *= R ;
	}
};
//!  \fn FReal plummerDist(int & cpt, const FReal &R)

//! \brief   Radial Plummer distribution

//!
//! \param cpt : counter to know how many random selections we need to obtain a radius less than R
//! \param R    : Radius of the sphere that contains the particles
//! @return Return the radius according to the Plummer distribution either double type or float type
//!
FReal plummerDist(int & cpt, const FReal &R) {
	//
	FReal radius ,u ;
	do  {
		//radius =  getRandom() ;
		u        = FMath::pow (getRandom() , 2.0/3.0) ;
		radius = FMath::Sqrt (u/(1.0-u));
		cpt++;
		if(radius  <=R){
			//			std::cout << radius << "    "  <<std::endl;
			return static_cast<FReal>(radius);
		}
	} while (true);
}
//! \fn void unifRandonPlummer(const int N , const FReal R, const FReal M, FReal * points)

//! \brief  Build N points following the Plummer distribution

//! First we construct N points uniformly distributed on the unit sphere. Then the radius in construct according to the Plummr distribution.
//!
//! \param N the number of points following the Plummer distribution
//! \param R the radius of the sphere that contains all the points
//! \param M the total mass of all the particles inside the Sphere or radius R
//! \param points array of size 4*N and stores data as follow x,y,z,0,x,y,z,0....
//!
void unifRandonPlummer(const int N , const FReal R, const FReal M, FReal * points) {
	//
	unifRandonPointsOnUnitSphere(N , points) ;
	//
	FReal r , rm= 0.0;
	//	FReal Coeff =  3.0*M/(4.0*FMath::FPi*R*R*R) ;
	//am1 = 0 ;//1/FMath::pow(1+R*R,2.5);
	int cpt = 0 ;
	for (int i = 0,j=0 ; i< N ; ++i, j+=4)  {
		// u \in []
		r = plummerDist(cpt,R) ;
		rm = std::max(rm, r);
		points[j]	   *= r ;
		points[j+1]  *= r ;
		points[j+2]  *= r ;
	}

	std::cout << "Total tested points: "<< cpt << " % of rejected points: "
			<<100*static_cast<FReal>(cpt-N)/cpt << " %" <<std::endl;

} ;
//! \fn void exportCVS(std::ofstream& file, const int N, const FReal * particles )

//! \brief  Export particles in CVS Format
//!
//! Export particles in CVS Format as follow
//!      x ,  y  , z , physicalValue
//! It is useful to plot the distribution with paraView
//!
void exportCVS(std::ofstream& file, const int N, const FReal * particles ){
	int j = 0;
	file << " x ,  y , z, q " <<std::endl;
	for(int i = 0 ; i< N; ++i, j+=4){
		file <<    particles[j]    << " , "    <<   particles[j+1]    << " , "   <<   particles[j+2]    << " , "   <<   particles[j+3]   <<std::endl;
	}
}
//
void exportCOSMOS(std::ofstream& file, const int N, const FReal * particles ){
	int j = 0;
	file << " x ,  y , z, q " <<std::endl;
	for(int i = 0 ; i< N; ++i, j+=4){
		file <<    particles[j]    << "  "    <<   particles[j+1]    << "  "   <<   particles[j+2]    << "  0.0 0.0 0.0  "   <<   particles[j+3]   <<"  " << i << std::endl;
	}
}
//
void exportVTK(std::ofstream& VTKfile, const int N, const FReal * particles ){
	int j = 0;
	//---------------------------
	// print generic information
	//---------------------------
	VTKfile << "# vtk DataFile Version 3.0" << "\n";
	VTKfile << "#  Generated bt exportVTK" << "\n";

	VTKfile << "ASCII" << "\n";
	VTKfile << "DATASET POLYDATA" << "\n";
	//
	//---------------------------------
	// print nodes ordered by their TAG
	//---------------------------------
	VTKfile << "POINTS " << N << "  float" << "\n";
	//
	for(int i = 0 ; i< N; ++i, j+=4){
		VTKfile <<    particles[j]    << "  "    <<   particles[j+1]    << "   "   <<   particles[j+2]      <<std::endl;
	}
	// ------------------------------------------
	VTKfile << "\n";
	VTKfile << "VERTICES  " <<  N << " " << 2*N << "\n";
	for(int i = 0 ; i< N; ++i, j+=4){
		VTKfile <<    "  1 "    << " "    <<i<<std::endl;
	}
	VTKfile << "POINT_DATA  " <<  N << "\n";
	VTKfile << "SCALARS PhysicalValue  float 1" << "\n"
			<< "LOOKUP_TABLE default" << "\n" ;
	j = 0 ;
	for(int i = 0 ; i< N; ++i, j+=4){
		VTKfile <<    particles[j+3]    << " "    <<std::endl;
	}
	VTKfile << "\n";
};

void exportVTKxml(std::ofstream& VTKfile, const int N, const FReal * particles ){
	int j = 0;

	VTKfile << "<?xml version=\"1.0\"?>" <<std::endl
			<< "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"> "<<std::endl
			<< "<PolyData>"<<std::endl
			<< "<Piece NumberOfPoints=\" " << N << " \"  NumberOfVerts=\" "<<N <<" \" NumberOfLines=\" 0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">"<<std::endl
			<< "<Points>"<<std::endl
			<< "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"> "<<std::endl ;
	j = 0 ;
	for(int i = 0 ; i< N; ++i, j+=4){
		VTKfile <<    particles[j]    << "  "    <<   particles[j+1]    << "   "   <<   particles[j+2]      << "   "   ;
	}
	VTKfile <<std::endl<< "</DataArray> "<<std::endl
			<< "</Points> "<<std::endl
			<< "<PointData Scalars=\"PhysicalValue\" > "<<std::endl
			<< "<DataArray type=\"Float64\" Name=\"PhysicalValue\"  format=\"ascii\">"<<std::endl ;
	j = 0 ;
	for(int i = 0 ; i< N; ++i, j+=4){
		VTKfile <<    particles[j+3]    << " "   ;
	}
	VTKfile <<std::endl << "</DataArray>"<<std::endl
			<< "	</PointData>"<<std::endl
			<< "	<CellData>"<<" </CellData>"<<std::endl
			<< "	<Verts>"<<std::endl
			<< "	<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<std::endl ;
	for(int i = 0 ; i< N; ++i){
		VTKfile <<   i   << " "   ;
	}
	VTKfile<<std::endl << "</DataArray>" <<std::endl
			<< "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<std::endl ;
	for(int i = 1 ; i< N+1; ++i){
		VTKfile <<   i   << " "   ;
	}
	VTKfile<<std::endl  << "</DataArray>"<<std::endl
			<< "	</Verts>"<<std::endl
			<< "<Lines></Lines>"<<std::endl
			<< "<Strips></Strips>"<<std::endl
			<< "<Polys></Polys>"<<std::endl
			<< "</Piece>"<<std::endl
			<< "</PolyData>"<<std::endl
			<< "</VTKFile>"<<std::endl;
} ;
//

#endif
