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
#ifndef FFMALOADER_HPP
#define FFMALOADER_HPP
#error("FFMALOADER_HPP ")

#include <iostream>
#include <fstream>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/FPoint.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
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
class FFmaLoader : public FAbstractLoader {
protected:
    std::ifstream file;                   //< The file to read
    FPoint           centerOfBox;     //< The center of box read from file
    FReal boxWidth;                     //< the box width read from file
    int nbParticles;                       //< the number of particles read from file

public:
    /**
    * The constructor need the file name
    * @param filename the name of the file to open
    * you can test if file is successfuly open by calling hasNotFinished()
    */
    FFmaLoader(const char* const filename): file(filename,std::ifstream::in){
        // test if open
        if(this->file.is_open()){
            FReal x,y,z;
            this->file >> this->nbParticles >> this->boxWidth >> x >> y >> z;
            this->centerOfBox.setPosition(x,y,z);
            this->boxWidth *= 2;
        }
        else {
             this->boxWidth = 0;
             this->nbParticles = 0;
        }
    }

    /**
    * Default destructor, simply close the file
    */
    virtual ~FFmaLoader(){
        file.close();
    }

    /**
      * To know if file is open and ready to read
      * @return true if loader can work
      */
    bool isOpen() const{
        return this->file.is_open() && !this->file.eof();
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
    void fillParticle(FPoint*const inParticlePositions, FReal*const inPhysicalValue){
        FReal x,y,z,data;
        this->file >> x >> y >> z >> data;
        inParticlePositions->setPosition(x,y,z);
        (*inPhysicalValue) = data;
    }

};


#endif //FFmaLoader_HPP

