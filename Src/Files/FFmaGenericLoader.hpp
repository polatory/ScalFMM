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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//
#include "Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "Utils/FPoint.hpp"

//
//! \class  FFmaGenericLoader
//!
//! \brief Read a set of particles in FMA format
//!
//! The FMA format is a simplest format to store the particles in a file.
//!
//!  It is organized as follow<br>
//!    NB_particles half_Box_width Center_X Center_Y Center_Z <br>
//!    X Y Z PhysicalValue // one particle by line
//!
//!  \code
//!     FFmaGenericLoader  loader("../ADir/Tests/particles.basic.txt");    // default ascii format
//!     if(!loader.isOpen()){
//!         std::cout << "Loader Error\n";
//!        return 1;
//!    }
//!
//!   FOctree<FBasicParticle, TestCell, FSimpleLeaf> tree(loader.getBoxWidth(),loader.getCenterOfBox());
//!
//!    FPoint position ;
//!    Freal   physicalValue ;
//!    for(int idx = 0 ; idx < loader.getNumberOfParticles() ; ++idx){
//!        loader.fillParticle(&position, &physicalValue);
//!        tree.insert(position, physicalValue);
//!    }
//! \endcode
//!
//!
//!
 //
class FFmaGenericLoader : public FAbstractLoader {
protected:
	std::fstream *file;                   ///< the stream used to read the file
	bool binaryFile  ;                    ///< if true the file to read is in binary mode
	FPoint     centerOfBox;            ///< The center of box (read from file)
	FReal boxWidth;                     ///< the box width (read from file)
	int nbParticles;                       ///< the number of particles (read from file)
  
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
	FFmaGenericLoader(const std::string & filename,const bool binary = false) :file(nullptr),binaryFile(binary),
	centerOfBox(0.0,0.0,0.0),boxWidth(0.0),nbParticles(0){
		if(binary) {
			this->file = new std::fstream (filename,std::ifstream::in| std::ios::binary);
		}
		else {
			this->file = new std::fstream(filename,std::ifstream::in) ;
		}
		// test if open
		if(this->file->is_open()){
			if(binaryFile){
			    this->readBinaryHeader();
			}
			else {
				this->readAscciHeader();
			}
		}
		else {
			this->boxWidth = 0;
			this->nbParticles = -10;
		}
		std::cout << " nbParticles: " <<this->nbParticles << std::endl
				<< " Box width:   " <<this->boxWidth << std::endl
				<< " Center:        " << this->centerOfBox << std::endl;
	}

	/**
	 * Default destructor, simply close the file
	 */
	virtual ~FFmaGenericLoader(){
		file->close();
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
		return FSize(this->nbParticles);
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
	 * Fill a particle form th curent position in the file
	 * @warning to work with the loader, particles has to expose a setPosition method
	 * @param the particle to fill
	 */
	void fillParticle(FPoint*const outParticlePositions, FReal*const outPhysicalValue){
		if(binaryFile){
			file->read((char*)(outParticlePositions), sizeof(FReal)*3);
			file->read((char*)(outPhysicalValue), sizeof(FReal));
		}
		else{
			FReal x,y,z;
//			(*this->file)  >> &outParticlePositions>> &outPhysicalValue;
			(*this->file)  >> x >> y >> z >> (*outPhysicalValue);
			outParticlePositions->setPosition(x,y,z);
		}
//		std::cout <<  "X Y Z Q " << *outParticlePositions<<" "<<*outPhysicalValue <<std::endl;

	}
private:
	void readAscciHeader() {
		FReal x,y,z;
		(*this->file) >> this->nbParticles >> this->boxWidth >> x >> y >> z;
		this->centerOfBox.setPosition(x,y,z);
		this->boxWidth *= 2;
	}
	void readBinaryHeader() {
	  int sizeOfElement;
	  file->seekg (std::ios::beg);
	  file->read((char*)&sizeOfElement,sizeof(int));
	  if(sizeOfElement != sizeof(FReal)){
	    std::cout << "Size of elements in part file " << sizeOfElement << " is different from size of FReal" << std::endl;
	    exit(0);
	  }
	  else{
	    file->read( (char*)&(this->nbParticles), sizeof(FSize) );
	    printf("NbPart found %d \n",this->nbParticles);
	    file->read( (char*)&(this->boxWidth) ,sizeof(this->boxWidth) );
	    this->boxWidth *= 2;
	    FReal x[3];
	    file->read( (char*)x,sizeof(FReal)*3);
	    this->centerOfBox.setPosition(x[0],x[1],x[2]);
	  }
	}

};


#endif //FFmaGenericLoader_HPP


