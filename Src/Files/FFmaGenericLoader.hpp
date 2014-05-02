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
#ifndef FFmaGenericLoader_HPP
#define FFmaGenericLoader_HPP


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/FPoint.hpp"

/**
 * @author Berenger Bramas and Olivier Coulaud
 * @class FFmaLoader
 * Please read the license
 *
 * Load a file with a format like :
 * NB_particles Box_width Center_X Center_Y Center_Z // init
 * X Y Z PhysicalValue // one particle by line
 * ....
 * @code
 *    FFmaLoader<FBasicParticle> loader("../ADir/Tests/particles.basic.txt"); <br>
 *    if(!loader.isOpen()){ <br>
 *        std::cout << "Loader Error\n"; <br>
 *        return 1; <br>
 *    } <br>
 * <br>
 *    FOctree<FBasicParticle, TestCell, FSimpleLeaf> tree(loader.getBoxWidth(),loader.getCenterOfBox()); <br>
 * <br>
 *    for(int idx = 0 ; idx < loader.getNumberOfParticles() ; ++idx){ <br>
 *        FBasicParticle* const part = new FBasicParticle(); <br>
 *        loader.fillParticle(part); <br>
 *        tree.insert(part); <br>
 *    } <br>
 * @endcode
 *
 * Particle has to extend {FExtendPhysicalValue,FExtendPosition}
 */
class FFmaGenericLoader : public FAbstractLoader {
protected:
	std::fstream *file;                   //< The file to read
	bool binaryFile  ;                     //< if true the file to read is in binary mode
	FPoint     centerOfBox;     //< The center of box read from file
	FReal boxWidth;                     //< the box width read from file
	int nbParticles;                       //< the number of particles read from file

public:
	/**
	 * The constructor need the file name
	 * @param filename the name of the file to open
	 * you can test if file is successfully open by calling hasNotFinished()
	 */
//	FFmaGenericLoader(const char* const filename,const bool binary = false) :file(nullptr),binaryFile(binary),
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
	void readAscciHeader() {
		FReal x,y,z;
		(*this->file) >> this->nbParticles >> this->boxWidth >> x >> y >> z;
		this->centerOfBox.setPosition(x,y,z);
		this->boxWidth *= 2;
	}
	void readBinaryHeader() {
		file->seekg (std::ios::beg);
		file->read( (char*)&(this->nbParticles), sizeof(int) );
		file->read( (char*)&(this->boxWidth) ,sizeof(this->boxWidth) );
		this->boxWidth *= 2;
		FReal x[3];
		file->read( (char*)x,sizeof(FReal)*3);
		this->centerOfBox.setPosition(x[0],x[1],x[2]);
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
	 * Fill a particle
	 * @warning to work with the loader, particles has to expose a setPosition method
	 * @param the particle to fill
	 */
	void fillParticle(FPoint*const outParticlePositions, FReal*const outPhysicalValue){
		FReal x,y,z,data;
		if(binaryFile){
			file->read((char*)(outParticlePositions), sizeof(FReal)*3);
			file->read((char*)(outPhysicalValue), sizeof(FReal));
		}
		else{
			(*this->file)  >> x >> y >> z >> data;
			outParticlePositions->setPosition(x,y,z);
			(*outPhysicalValue) = data;
		}
//		std::cout <<  "X Y Z Q " << *outParticlePositions<<" "<<*outPhysicalValue <<std::endl;

	}
};


#endif //FFmaGenericLoader_HPP

