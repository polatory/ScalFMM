// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
// author Berenger Bramas and Olivier Coulaud
//
#ifndef FFmaGenericLoader_HPP
#define FFmaGenericLoader_HPP

#include <ios>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
//
#include "Utils/FGlobal.hpp"
#include "Utils/FAssert.hpp"
#include "FAbstractLoader.hpp"
#include "Utils/FPoint.hpp"


//! \class  FmaRWParticle
//!
//! \brief The Particle class used in FMA loader and writer
//!
//! In this class, we use template parameters in order to let the user choose how many datas he wish to store.
//! The datas are all FReal.
//! The first ones are : PosX,PosY,PosZ,physicalValue ==> Required
//! The following ones are : Potential,forceX,forceY,forceZ == optionnal
//!
template<int READ, int WRITE>
class FmaRWParticle {
    //Data stored
    FReal data[WRITE];
public:
    FmaRWParticle(){
	if(WRITE<4){
	    std::cout << "Cannot create FmaRWParticle with less than 4 as value for WRITE" << std::endl;
	    exit(0);
	}
    }
    
    //Get a FPoint from the position
    FPoint getPosition() const{
	return FPoint(data[0],data[1],data[2]);
    }
    
    //Set the position from a FPoint
    void setPosition(FPoint & inPoint){
	data[0] = inPoint.getX();
	data[1] = inPoint.getY();
	data[2] = inPoint.getZ();
    }

    //Get a FReal from the physicalValue
    FReal getPhysicalValue() const{
	return data[3];
    }
    
    //Get a ptr to be able to set the physicalValue
    FReal* setPhysicalValue() {
	return &data[3];
    }

    //Get a FReal from the potential
    FReal getPotential() const{
	FAssertLF(READ>4,"Cannot access to Potential with READ<=4");
	return data[4];
    }
    
    //Get a ptr to be able to set the potential
    FReal* setPotential() {
	FAssertLF(WRITE>4,"Cannot set Potential with WRITE<=4");
	return &data[4];
    }
    
    //Get a ptr to read the forces
    FReal* getForces() {
	FAssertLF(READ>7,"Cannot access to forces[] with READ<=8");
	return &data[5];
    }
    
    //Get a ptr to write the forces
    FReal* setForces() {
	FAssertLF(WRITE>7,"Cannot set Forces[] with WRITE<=7");
	return &data[5];
    }

    //Get directly a ptr to the data
    FReal  * getPtrFirstData(){
	return data;
    }
    //Same as above with const qualifier
    const  FReal * getPtrFirstData() const{
	return data;
    }
    
    //Get READ
    unsigned int getReadDataNumber() const{
	return (unsigned int) (READ);
    }
    
    //Get WRITE
    unsigned int getWriteDataNumber() const{
	return WRITE;
    }

    //Get size of Class Particle
    unsigned int getWriteDataSize() const { 
	return sizeof(FmaRWParticle<READ,WRITE>);
    }
    //Get Size of array (should be same as above...)
    long unsigned int getClassSize() const {
	return WRITE*sizeof(FReal);
    }
};


//! \class  FmaR4W4Particle
//!
//! \brief Basic  Particle class used in FMA loader and writer
//!
//!  Here we consider only the position and the physical value in the structure
//!  We read (R4) and write  (W4)  the four values: the position, the physical value.
//!  This class is used in the generateDistributions example.
//!
// class FmaR4W4Particle {
// public:
// 	FPoint position;            ///< position of the particle
// 	FReal  physicalValue;    ///<  its physical value
// 	/**
// 	 *  return a pointer on the first value of the structure
// 	 */
// 	FReal  * getPtrFirstData()
// 	{return position.getDataValue() ;}
// 	const  FReal * getPtrFirstData() const
// 	{return position.getDataValue() ;}
// 	/**
// 	 *  return The number of data to read. it is used in FMAGenericLoader  (here 4)
// 	 */
// 	int getReadDataNumber()
// 	{ return 4;}
// 	/**
// 	 *  return The number of data to write. it is used in FMAGenericWriter (here 4)
// 	 */
// 	int getWriteDataNumber() const
// 	{ return 4;}
// 	/**
// 	 *  return size in byte of the structure. it is used in FMAGenericWriter
// 	 */
// 	unsigned int getWriteDataSize() const
// 	{ return sizeof(FmaR4W4Particle);}
// };

//! \class  FmaR4W8Particle
//!
//! \brief The Particle class used in FMA loader and writer
//!
//! In this class we consider the position, the physical value, the potential and the force in the structure
//! but we read (R4) only the four first values: the position and the physical value and we write (W8) all the data (8 values) in a file.
//! This class can be used if you read a file  generated by generateDistributions and you want to store the result of a direct computations or a FFM one..
//!
// class FmaR4W8Particle {
// public:
// 	FPoint position;            ///< position of the particle
// 	FReal  physicalValue;    ///< its physical value (mass or charge
// 	FReal  potential;           ///< the potential
// 	FReal  forces[3];            ///<the force

// 	/**
// 	 *  return a pointer on the first value of the structure
// 	 */
// 	FReal * getPtrFirstData()
// 	{return position.getDataValue() ;}
// 	const  FReal * getPtrFirstData() const
// 	{return position.getDataValue() ;}
// 	/**
// 	 *  return The number of data to read. it is used in FMAGenericLoader  (here 4)
// 	 */
// 	int getReadDataNumber()
// 	{ return 4;}
// 	/**
// 	 *  return The number of data to write. it is used in FMAGenericWriter (here 8)
// 	 */
// 	int getWriteDataNumber() const
// 	{ return 8;}
// 	/**
// 	 *  return size in byte of the structure. it is used in FMAGenericWriter
// 	 */
// 	unsigned int getWriteDataSize() const
// 	{ return sizeof(FmaR4W8Particle);}
// };
//
//! \class  FmaR8W8Particle
//!
//! \brief The Particle class used in FMA loader and writer
//!
//! Same as FmaR4W8Particle class but here we want to read  (R8) all the data (8 values) and want to store also all the data.
//! This class is used if you create an array of particles  from a file generate by a previous computation (directComputation or a FMM computation)
//! and you want to compare your result with the read values .

//!
// class FmaR8W8Particle : public FmaR4W8Particle {
// public:
// 	/**
// 	 *  return The number of data to read. it is used in FMAGenericLoader  (here 8)
// 	 *  Useful to read the result of the output of DirectComputation and to compare the potential and the force
// 	 *  with the FMM computation
// 	 */
// 	int getReadDataNumber()
// 	{ return 8;}
// };
//!\typedef FmaBasicParticle an alias of FmaR4W4Particle
//!  Particle contains 4 values of type FReal and we read and write the four values
//!
typedef FmaRWParticle<4,4> FmaBasicParticle  ;
typedef FmaRWParticle<4,8> FmaRParticle        ;
typedef FmaRWParticle<8,8> FmaParticle          ;
//
//! \class  FFmaGenericLoader
//!
//! \brief Read a set of particles in FMA format
//!
//! The FMA format is a simplest format to store the particles in a file.
//!
//!  It is organized as follow<br>
//!   DatatypeSise Number_of_record_per_line  <br>
//!    NB_particles half_Box_width Center_X Center_Y Center_Z <br>
//!    Particle Values
//!
//!  if  Number_of_record_per_line is  <br>
//!      4    the Particle Values  are X Y Z Q <br>
//!      8   the Particle Values  are X X Y Z Q  P FX FY FZ<br>
//!
//! There is 3 methods to read the data from the file <br>
//!    -# 	fillParticle(FPoint*, FReal*);<br>
//!    -#    fillParticle(FReal*, int);<br>
//!    -#    fillParticle(PartClass);<br>
//!
//!  \code
//!     FFmaGenericLoader  loader("../Data/unitCubeXYZQ20k.fma");    // extension fma --> ascii format
//!
//!	 FmaRParticle * const particles = new FmaRParticle[nbParticles];
//!	  memset(particles, 0, sizeof(FmaRParticle) * nbParticles) ;
//!
//!    FPoint position ;
//!    Freal   physicalValue ;
//!    for(int idx = 0 ; idx < loader.getNumberOfParticles() ; ++idx){
//!              loader.fillParticle(particles[idx]);
//!    }
//! \endcode
//!
//! \warning It works only in shared memory (doesn't work with MPI)
//!
//
class FFmaGenericLoader : public FAbstractLoader  {
protected:
	std::fstream *file;                   ///< the stream used to read the file
	bool binaryFile  ;                    ///< if true the file to read is in binary mode
	FPoint     centerOfBox;            ///< The center of box (read from file)
	FReal boxWidth;                     ///< the box width (read from file)
	FSize    nbParticles;                ///< the number of particles (read from file)
	unsigned int typeData[2];      ///< {Size of the data to read, number of data on 1 line}
private:
	FReal   *tmpVal ;                         ///  Temporary array to read data
	 int  otherDataToRead  ;       ///< <<number of other data (>4)to read in a particle record
public:
	/**
	 * The constructor need the file name
	 * @param filename the name of the file to open
	 * @param binary   true if the file to open is in binary mode
	 *
	 *  This function also read the header of the fma file (Two first lines)
	 *  This means that we can obtain after the call the number of particles, ...
	 * you can test if file is successfully open by calling hasNotFinished()
	 */
	FFmaGenericLoader(const std::string & filename,const bool binary ) :file(nullptr),binaryFile(binary),
	centerOfBox(0.0,0.0,0.0),boxWidth(0.0),nbParticles(0),tmpVal(nullptr),otherDataToRead(0){
		if(binary) {
			this->file = new std::fstream (filename.c_str(),std::ifstream::in| std::ios::binary);
		}
		else {
			this->file = new std::fstream(filename.c_str(),std::ifstream::in) ;
		}
		// test if open
		if(! this->file->is_open()){
			std::cerr << "File "<< filename<<" not opened! " <<std::endl;
			std::exit( EXIT_FAILURE);
		}
		this->readHeader();
	}
	/**
	 * The constructor need the file name
	 * @param filename the name with the extension .fma, .bfma of the file to open. We check the extension to know if the file is in ASCII mode or binary mode
	 *
	 *  This function also read the header of the  file (Two first lines)
	 *  This means that we can obtain after the call the number of particles, ...
	 * you can test if file is successfully open by calling hasNotFinished()
	 */
	FFmaGenericLoader(const std::string & filename) : file(nullptr),binaryFile(false),
			centerOfBox(0.0,0.0,0.0),boxWidth(0.0),nbParticles(0),tmpVal(nullptr),otherDataToRead(0) {
		std::string ext(".bfma");
		// open particle file
		if(filename.find(ext) != std::string::npos) {
			binaryFile = true;
			this->file = new std::fstream (filename.c_str(),std::ifstream::in| std::ios::binary);
		}
		else if(filename.find(".fma")!=std::string::npos ) {
			this->file = new std::fstream(filename.c_str(),std::ifstream::in) ;
		}
		else  {
			std::cout << "Input file not allowed only .fma or .bfma extensions" <<std::endl;
			std::exit ( EXIT_FAILURE) ;
		}
		// test if open
		if(! this->file->is_open()){
			std::cerr << "File "<< filename<<" not opened! " <<std::endl;
			std::exit( EXIT_FAILURE);
		}
		this->readHeader();
	}
	/**
	 * Default destructor, simply close the file
	 */
	virtual ~FFmaGenericLoader(){
		file->close();
		delete file ;
		delete tmpVal;
	}

	/**
	 * To know if file is open and ready to read
	 * @return true if loader can work
	 */
	bool isOpen() const{
		return this->file->is_open() && !this->file->eof();
	}

	/**
	 * To get the number of particles from this loader
	 * @param the number of particles the loader can fill
	 */
	FSize getNumberOfParticles() const{
		return this->nbParticles;
	}

	/**
	 * The center of the box from the simulation file opened by the loader
	 * @return box center
	 */
	FPoint getCenterOfBox() const{
		return this->centerOfBox;
	}
	/**
	 * The box width from the simulation file opened by the loader
	 * @return box width
	 */
	FReal getBoxWidth() const{
		return this->boxWidth;
	}
	/**
	 * The box width from the simulation file opened by the loader
	 * @return the number of data per record (Particle)
	 */
	unsigned int getNbRecordPerline(){
		return typeData[1]; }
	/**
	 * To know if the data are in float or in double type
	 * @return the type of the values float (4) or double (8)
	 */
	unsigned int getDataType(){
		return typeData[0]; }

	/**
	 * Fill a particle form the current position in the file
	 *
	 * @param outParticlePositions the position of particle to fill (FPoint class)
	 * @param outPhysicalValue the physical value of particle to fill (FReal)
	 *
	 */
	void fillParticle(FPoint*const outParticlePositions, FReal*const outPhysicalValue){
	    if(binaryFile){
			file->read((char*)(outParticlePositions), sizeof(FReal)*3);
			file->read((char*)(outPhysicalValue), sizeof(FReal));
			if(otherDataToRead> 0){
				file->read((char*)(this->tmpVal), sizeof(FReal)*otherDataToRead);
			}
		}
		else{
			FReal x,y,z;
			//			(*this->file)  >> &outParticlePositions>> &outPhysicalValue;
			(*this->file)  >> x >> y >> z >> (*outPhysicalValue);
			outParticlePositions->setPosition(x,y,z);
			
			if(otherDataToRead> 0){
				for (int 	i = 0 ; i <otherDataToRead; ++i){
					(*this->file) >> x ;
				}
			}
		}
		//		std::cout <<  "X Y Z Q " << *outParticlePositions<<" "<<*outPhysicalValue <<std::endl;

	}
	/**
	 * Fill a particle form the current position in the file
	 * @param dataToRead is an array of type FReal. It contains all the values of a particles (for instance X,Y,Z,Q, ..
	 * @param nbDataToRead number of value to read (I.e. size of the array)
	 */
    void fillParticle(FReal* dataToRead, const unsigned int nbDataToRead){
		if(binaryFile){
			file->read((char*)(dataToRead), sizeof(FReal)*nbDataToRead);
			if(nbDataToRead< typeData[1]){
				file->read((char*)(this->tmpVal), sizeof(FReal)*(typeData[1]-nbDataToRead));
			}
		}
		else{

			for (unsigned int i = 0 ; i <nbDataToRead; ++i){
				(*this->file)  >>dataToRead[i];
			}
			if(nbDataToRead< typeData[1]){
				FReal x;
				for (unsigned int 	i = 0 ; i <typeData[1]-nbDataToRead; ++i){
					(*this->file) >> x ;
				}
			}
		}
	}

	/**
	 * Fill a particle form the current position in the file
	 * @param dataToRead is the particle  class. If the class is different from FmaBasicParticle, FmaRParticle or FmaParticle the we have to implement the following method
	 *   getPtrFirstData(), getReadDataNumber() see FmaBasicParticle class for more details.
	 */
	template <class dataPart>
	void fillParticle(dataPart &dataToRead){
	    int otherDataRead = typeData[1] - dataToRead.getReadDataNumber() ;
		if(binaryFile){
			file->read((char*)(dataToRead.getPtrFirstData()), sizeof(FReal)*(dataToRead.getReadDataNumber()));
			if( otherDataRead > 0){
				file->read((char*)(this->tmpVal), sizeof(FReal)*(otherDataRead));
			}
		}
		else{
			FReal * val = dataToRead.getPtrFirstData();
			for (unsigned int i = 0 ; i <dataToRead.getReadDataNumber(); ++i){
				(*this->file)  >>*val;
				++val;
			}
			if( otherDataRead > 0){
				FReal x;
				for (int i = 0 ; i <otherDataRead ;++i){
					(*this->file)  >>x;
				}
			}
		}
	}
	/**
	 * Fill a set of particles form the current position in the file.
	 *  If the file is a binary file and we read all record per particle then we read and fill the array in one instruction
	 * @param dataToRead is an array of the particle  class. If the class is different from FmaBasicParticle, FmaRParticle or FmaParticle the we have to implement the following method
	 *   getPtrFirstData(), getReadDataNumber() see FmaBasicParticle class for more details.
	 */
	template <class dataPart>
	void fillParticle(dataPart *dataToRead, const int N){
	    int otherDataRead = typeData[1] - (*dataToRead).getReadDataNumber() ;
		if(binaryFile && otherDataRead == 0 ){
		    file->read((char*)((*dataToRead).getPtrFirstData()), sizeof(FReal)*(N*(*dataToRead).getReadDataNumber()));
		}
		else {
			for (int i = 0 ; i <N; ++i){
				this->fillParticle(dataToRead[i]) ;
			}
		}
	}
private:
	void readHeader() {
		if(this->binaryFile){
			this->readBinaryHeader();
		}
		else {
			this->readAscciHeader();
		}

		std::cout << "   nbParticles: " <<this->nbParticles << std::endl
				<< "   Box width:   " <<this->boxWidth << std::endl
				<< "   Center:        " << this->centerOfBox << std::endl;
	}
	void readAscciHeader() {
		std::cout << " File open in ASCII mode "<< std::endl ;
		FReal x,y,z;
		(*this->file) >> typeData[0]>> typeData[1];
		std::cout << "   Datatype "<< typeData[0] << " "<< typeData[1] << std::endl;
		(*this->file) >> this->nbParticles >> this->boxWidth >> x >> y >> z;
		this->centerOfBox.setPosition(x,y,z);
		this->boxWidth *= 2;
		otherDataToRead = typeData[1] -4;
	};
	void readBinaryHeader(){
		std::cout << " File open in binary mode "<< std::endl;
		file->seekg (std::ios::beg);
		file->read((char*)&typeData,2*sizeof(unsigned int));
		std::cout << "   Datatype "<< typeData[0] << " "<< typeData[1] << std::endl;
		if(typeData[0] != sizeof(FReal)){
			std::cerr << "Size of elements in part file " << typeData[0] << " is different from size of FReal " << sizeof(FReal)<<std::endl;
			std::exit( EXIT_FAILURE);
		}
		else{
			file->read( (char*)&(this->nbParticles), sizeof(FSize) );
			file->read( (char*)&(this->boxWidth) ,sizeof(this->boxWidth) );
			this->boxWidth *= 2;

			FReal x[3];
			file->read( (char*)x,sizeof(FReal)*3);
			this->centerOfBox.setPosition(x[0],x[1],x[2]);
		}
		otherDataToRead = typeData[1] -4;
		if(otherDataToRead>0){
			tmpVal = new FReal[otherDataToRead];
		}
	}

};
//
//! \class  FFmaGenericWriter
//!
//! \brief Write a set of particles in FMA format (ASCII or binary mode)
//!
//! The FMA format is a simplest format to store the particles in a file.
//!
//!  It is organized as follow<br>
//!    DatatypeSise Number_of_record_per_line
//!    NB_particles half_Box_width Center_X Center_Y Center_Z <br>
//!    X Y Z PhysicalValue // one particle per line in ascii mode
//!
//!  DatatypeSise = 4 float ; 8 double
//!
//!  Number_of_record_per_line 4 (X,Y,Z,Q)  ; 8 (X,Y,Z,Q,P,FX,FY,FZ  )
//!  \code
//!     FFmaGenericWriter writer ("data.bfma");    // open binary fma file (extension .bfma)
//!
//!      // write the header of the file
//!     writer.writeHeader( loader.getCenterOfBox(), loader.getBoxWidth() , NbPoints, sizeof(FReal), nbData) ;
//!
//!      // write the data here particles is an array of float or double and a particle has nbData value
//!     writer.writeArrayOfReal(particles, nbData, NbPoints);
//!
//! \endcode
//! \warning It works only in shared memory (doesn't work with MPI)
//!

class FFmaGenericWriter {

protected:
	std::fstream *file;                   ///< the stream used to read the file
	bool binaryFile  ;                    ///< if true the file to read is in binary mode

public:
	/**
	 * The constructor need the file name
	 * @param filename the name of the file to open
	 * @param binary   true if the file to open is in binary mode
	 *
	 *  This function also read the header of the fma file (First line)
	 *  This means that we can obtain after the call the number of particles, ...
	 * you can test if file is successfully open by calling hasNotFinished()
	 */
	FFmaGenericWriter(const std::string & filename): binaryFile(false){
		std::string ext(".bfma");
		// open particle file
		if(filename.find(ext) !=std::string::npos) {
			binaryFile = true;
			this->file = new std::fstream (filename.c_str(),std::ifstream::out| std::ios::binary);
		}
		else if(filename.find(".fma")!=std::string::npos ) {
			this->file = new std::fstream(filename.c_str(),std::ifstream::out) ;
			this->file->precision(10);
		}
		else  {
			std::cout << "Input file not allowed only .fma or .bfma extensions" <<std::endl;
			std::exit ( EXIT_FAILURE) ;
		}
		// test if open
		if(! this->file->is_open()){
			std::cerr << "File "<< filename<<" not opened! " <<std::endl;
			std::exit( EXIT_FAILURE);
		}
	}
	/**
	 * The constructor need the file name
	 * @param filename the name of the file to open
	 * @param binary   true if the file to open is in binary mode
	 *
	 *  This function also read the header of the fma file (First line)
	 *  This means that we can obtain after the call the number of particles, ...
	 * you can test if file is successfully open by calling hasNotFinished()
	 */
	FFmaGenericWriter(const std::string & filename, const bool binary ) : file(nullptr), binaryFile(binary)
	{
		if(binary) {
			this->file = new std::fstream (filename.c_str(),std::ifstream::out| std::ios::binary);
		}
		else {
			this->file = new std::fstream(filename.c_str(),std::ifstream::out) ;
			this->file->precision(10);
		}
		// test if open
		if(! this->file->is_open()){
			std::cerr << "File "<< filename<<" not opened! " <<std::endl;
			std::exit( EXIT_FAILURE);
		}
	}
	/**
	 * Default destructor, simply close the file
	 */
	virtual ~FFmaGenericWriter(){
		file->close();
		delete file ;
	}

	/**
	 * To know if file is open and ready to read
	 * @return true if loader can work
	 */
	bool isOpen() const{
		return this->file->is_open() && !this->file->eof();
	}
	//!
	//! Write the header of FMA file
	//! \warning All values inside typePart should be of the same type {float or double}
	//!
	//!@param centerOfBox  The centre of the Box (FPoint class)
	//!@param boxWidth       The width of the box
	//!@param nbParticles     Number of particles in the box (or to save)
	//!@param data               Data type of the particle class (FmaBasicParticle, FmaRParticle or FmaParticle)
	//!
	template <class typePart>
	void writeHeader(const FPoint &centerOfBox,const FReal &boxWidth, const FSize &nbParticles, const typePart  data) {
		unsigned int typeFReal[2]  = {sizeof(FReal) , sizeof(typePart) / sizeof(FReal) };
		const unsigned int ndata  = data.getWriteDataNumber();
		std::cout <<"    WriteHeader: typeFReal: " << typeFReal[0]  << "  nb Elts: " << typeFReal[1]  <<"   NData to write "<< ndata<< "\n";
		if (ndata != typeFReal[1]){
			typeFReal[1] = ndata;
		}
		FReal x = boxWidth *0.5;
		if(this->binaryFile) {
			this->writerBinaryHeader(centerOfBox,x,nbParticles,typeFReal);
		}
		else {
			this->writerAscciHeader(centerOfBox,x,nbParticles,typeFReal);
		}
	}
	//!
	//! Write the header of FMA file. Should be used if we write the particles with writeArrayOfReal method
	//!
	//!
	//!@param centerOfBox      The center of the Box (FPoint class)
	//!@param boxWidth            The width of the box
	//!@param nbParticles          Number of particles in the box (or to save)
	//!@param dataType             Data type of the values in particle
	//!@param nbDataPerRecord  Number of record/value per particle
	//!
	void writeHeader(const FPoint &centerOfBox,const FReal &boxWidth, const FSize &nbParticles,
			const unsigned int  dataType, const unsigned int  nbDataPerRecord) {
		unsigned int typeFReal[2]  = {dataType , nbDataPerRecord };
		FReal x = boxWidth *0.5;
		if(this->binaryFile) {
			this->writerBinaryHeader(centerOfBox,x,nbParticles,typeFReal);
		}
		else {
			this->writerAscciHeader(centerOfBox,x,nbParticles,typeFReal);
		}
	}

	//!
	//!@warning the type dataPart should be FmaBasicParticle, FmaRParticle or FmaParticle  class or should implement all methods inside FmaBasicParticle.
	//!
	//!@param dataToWrite array of particles of type dataPart
	//!@param N number of element in the array
	//!
	//!  example 1
	//!\code
	//!	FmaRParticle *  particles = new FmaRParticle[nbParticles];
	//!     memset(particles, 0, sizeof(FmaRParticle) * nbParticles) ;
	//!    ...
	//! 	FFmaGenericWriter writer(filenameOut) ;
	//! 	Fwriter.writeHeader(Centre,BoxWith, nbParticles,*particles) ;
	//! 	Fwriter.writeArrayOfParticles(particles, nbParticles);
	//! \endcode
	//!
	//!  example2
	//!\code
	//!	FReal *  particles = new FReal[4*NbPoints] ; // store 4 data per particle
	//!     memset(particles, 0, sizeof(FmaRParticle) * nbParticles) ;
	//!    ...
	//!    FmaBasicParticle *ppart = (FmaBasicParticle*)(&particles[0]);
	//! 	FFmaGenericWriter writer(filenameOut) ;
	//! 	Fwriter.writeHeader(Centre,BoxWith, nbParticles,*particles) ;
	//! 	Fwriter.writeArrayOfParticles(particles, nbParticles);
	//! \endcode

	template <class dataPart>
	void writeArrayOfParticles(const dataPart *dataToWrite, const FSize N){
		//		std::cout << "NB points to write: "<< N <<std::endl;
		if(binaryFile){
			unsigned int recordSize=  dataToWrite[0].getWriteDataSize() ;
			unsigned int typeFReal[2]      = {sizeof(FReal) , sizeof(dataPart) / sizeof(FReal) };
			//			std::cout << "typeData "<< typeFReal[0] << " "<< typeFReal[1] <<"  "<< std::endl;
			//
			if (sizeof(dataPart) == recordSize){
				//				std::cout << "Size to write:  "<<N*dataToWrite[0].getWriteDataSize() <<std::endl;
				file->write((char*)(dataToWrite[0].getPtrFirstData()), N*recordSize);
			}
			else {
				file->write((char* )&typeFReal[0],2*sizeof(unsigned int));
				//				std::cout << "Size to write:   N* "<<typeFReal[0] *typeFReal[1]  <<std::endl;
				for (int i = 0 ; i <N ; ++i){
					file->write((char*)(dataToWrite[i].getPtrFirstData()), recordSize );
					//					const FReal * val = dataToWrite[i].getPtrFirstData() ;
					//					std::cout << i <<"   ";
					//					for( int j=0; j<typeFReal[1] ; ++j){
					//						std::cout << *val << "   ";++val;
					//					}
					//					std::cout <<std::endl;
				}
			}
		}
		else{ // ASCII part
			const int ndata = dataToWrite[0].getWriteDataNumber();
//			std::cout << "typeData "<< sizeof(FReal) << " "<<ndata << std::endl;
			this->file->precision(10);

			for (int i = 0 ; i <N ; ++i){
				const FReal * val = dataToWrite[i].getPtrFirstData() ;
				for (int j= 0 ; j <ndata ; ++j){
					(*this->file)  << *val << "    "; ++val;
				}
				(*this->file)  <<std::endl;
			}
		}
	}
	//!
	//! Write an array of data in a file Fill
	//!
	//!@param dataToWrite array of particles of type FReal
	//!@param nbData number of data per particle
	//!@param N number of particles
	//!
	//!  The size of the array is N*nbData
	//!
	/*	//!  example
	//!\code
	//!	FmaRParticle * const particles = new FmaRParticle[nbParticles];
	//!     memset(particles, 0, sizeof(FmaRParticle) * nbParticles) ;
	//!    ...
	//! 	FFmaGenericWriter writer(filenameOut) ;
	//! 	Fwriter.writeHeader(Centre,BoxWith, nbParticles,*particles) ;
	//! 	Fwriter.writeArrayOfReal(particles, nbParticles);
	//! \endcode
	 */
	void writeArrayOfReal(const FReal *dataToWrite, const unsigned int nbData, const FSize N){
		if(binaryFile){
			file->write((char*)(dataToWrite), N*nbData*sizeof(FReal));
		}
		else{
			this->file->precision(10);
			//			std::cout << "N "<< N << " nbData "<< nbData<<std::endl;
			//			exit(-1);
			int k = 0;
			for (int i = 0 ; i <N ; ++i){
				//				std::cout << "i "<< i << "  ";
				for (unsigned int jj= 0 ; jj<nbData ; ++jj, ++k){
					(*this->file)  << dataToWrite[k] << "    ";
					//					std::cout      << dataToWrite[k]<< "  ";
				}
				(*this->file)  <<std::endl;
				//				std::cout <<std::endl;
			}
			//			std::cout << "END"<<std::endl;
		}
	}
private:
	void writerAscciHeader( const FPoint &centerOfBox,const FReal &boxWidth,
			const FSize &nbParticles, const unsigned int *typeFReal) {
		this->file->precision(10);
		(*this->file) << typeFReal[0] <<"   "<<typeFReal[1]<<std::endl;
		(*this->file) << nbParticles << "   "<<  boxWidth << "   "
				<<  centerOfBox.getX()  << "  " << centerOfBox.getY() << " "<<centerOfBox.getZ()
				<< std::endl;
	}
	void writerBinaryHeader(const FPoint &centerOfBox,const FReal &boxWidth,
			const FSize &nbParticles, const unsigned int *typeFReal) {
		file->seekg (std::ios::beg);
		file->write((char*)typeFReal,2*sizeof(unsigned int));
		if(typeFReal[0]  != sizeof(FReal)){
			std::cout << "Size of elements in part file " << typeFReal[0] << " is different from size of FReal " << sizeof(FReal)<<std::endl;
			std::exit( EXIT_FAILURE);
		}
		else{
			file->write( (char*)&(nbParticles), sizeof(FSize) );
			//			std::cout << "nbParticles "<< nbParticles<<std::endl;
			file->write( (char*)&(boxWidth) ,sizeof(boxWidth) );
			file->write( (char*)(centerOfBox.getDataValue()),sizeof(FReal)*3);
		}
	}

};


#endif //FFmaGenericLoader_HPP


