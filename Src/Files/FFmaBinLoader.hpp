// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FFMABINLOADER_HPP
#define FFMABINLOADER_HPP


#include <cstdio>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/FPoint.hpp"
#include "../Utils/FDebug.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmaBinLoader
* Please read the license
*
* Load a file with a format like :
* NB_particles Box_width Box_X Box_Y Box_Z // init
* X Y Z // one particle by line
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
template <class ParticleClass>
class FFmaBinLoader : public FAbstractLoader<ParticleClass> {
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
            FDEBUG(if(sizeOfElement != int(sizeof(FReal)) ){)
                FDEBUG( FDebug::Controller.writeFromLine("Warning type size between file and FReal are differents\n", __LINE__, __FILE__); )
                    printf("%d sizeofelement\n",sizeOfElement);
            FDEBUG(})
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
    void fillParticle(ParticleClass& inParticle){
        FReal x,y,z,data;

        removeWarning += fread(&x, sizeof(FReal), 1, file);
        removeWarning += fread(&y, sizeof(FReal), 1, file);
        removeWarning += fread(&z, sizeof(FReal), 1, file);
        removeWarning += fread(&data, sizeof(FReal), 1, file);

        inParticle.setPosition(x,y,z);
        inParticle.setPhysicalValue(data);
    }

};


#endif //FFmaBinLoader_HPP


