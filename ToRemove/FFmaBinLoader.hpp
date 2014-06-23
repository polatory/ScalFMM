// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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
#ifndef FFMABINLOADER_HPP
#define FFMABINLOADER_HPP


#include <cstdio>

#include "../Utils/FGlobal.hpp"
#error(" FFmaBinLoader.hpp ")
#include "FAbstractLoader.hpp"
#include "../Utils/FPoint.hpp"
#include "../Utils/FLog.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmaBinLoader
* Please read the license
*
* Load a file with a format like :
* NB_particles Box_width Center_X Center_Y Center_Z // init
* X Y Z PhysicalValue // one particle by line
* ....
* <code>
*    FFmaBinLoader<FBasicParticle> loader("../Adir/Tests/particles.basic.txt"); <br>
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
* </code>
*
* Particle has to extend {FExtendPhysicalValue,FExtendPosition}
*/
class FFmaBinLoader : public FAbstractLoader {
protected:
    FILE* const file;            //< The file to read
    FPoint centerOfBox;     //< The center of box read from file
    FReal boxWidth;              //< the box width read from file
    FSize nbParticles;             //< the number of particles read from file

    size_t removeWarning;

public:
    /**
    * The constructor need the file name
    * @param filename the name of the file to open
    * you can test if file is successfuly open by calling hasNotFinished()
    */
    FFmaBinLoader(const char* const filename): file(fopen(filename, "rb")), removeWarning(0) {
        // test if open
        if(this->file != NULL) {
            int sizeOfElement(0);
            removeWarning += fread(&sizeOfElement, sizeof(int), 1, file);
            FLOG(if(sizeOfElement != int(sizeof(FReal)) ){)
                FLOG( FLog::Controller.writeFromLine("Warning type size between file and FReal are differents\n", __LINE__, __FILE__); )
                    printf("%d sizeofelement\n",sizeOfElement);
            FLOG(})
            removeWarning += fread(&this->nbParticles, sizeof(FSize), 1, file);

            removeWarning += fread(&this->boxWidth, sizeof(FReal), 1, file);
            this->boxWidth *= 2;

            FReal x,y,z;
            removeWarning += fread(&x, sizeof(FReal), 1, file);
            removeWarning += fread(&y, sizeof(FReal), 1, file);
            removeWarning += fread(&z, sizeof(FReal), 1, file);
            this->centerOfBox.setPosition(x,y,z);
        }
        else {
             this->boxWidth = 0;
             this->nbParticles = 0;
        }
    }

    /**
    * Default destructor, simply close the file
    */
    virtual ~FFmaBinLoader(){
        if(file) fclose(file);
    }

    /**
      * To know if file is open and ready to read
      * @return true if loader can work
      */
    bool isOpen() const{
        return this->file != NULL;
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
      * Fill a particle
      * @warning to work with the loader, particles has to expose a setPosition method
      * @param the particle to fill
      */
    void fillParticle(FPoint*const inParticlePosition, FReal*const physicalValue){
        FReal x,y,z,data;

        removeWarning += fread(&x, sizeof(FReal), 1, file);
        removeWarning += fread(&y, sizeof(FReal), 1, file);
        removeWarning += fread(&z, sizeof(FReal), 1, file);
        removeWarning += fread(&data, sizeof(FReal), 1, file);

        inParticlePosition->setPosition(x,y,z);
        (*physicalValue) = data;
    }

};


#endif //FFmaBinLoader_HPP


