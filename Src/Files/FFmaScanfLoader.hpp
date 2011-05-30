#ifndef FFMASCANFLOADER_HPP
#define FFMASCANFLOADER_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmaScanfLoader
* Please read the license
*
* Load a file with a format like :
* NB_particles Box_width Box_X Box_Y Box_Z // init
* X Y Z // one particle by line
* ....
* <code>
*    FFmaScanfLoader<FBasicParticle> loader("../FMB++/Tests/particles.basic.txt"); <br>
*    if(!loader.isValide()){ <br>
*        std::cout << "Loader Error\n"; <br>
*        return 1; <br>
*    } <br>
* <br>
*    FOctree<FBasicParticle, TestCell, FSimpleLeaf, 10, 3> tree(loader.getBoxWidth(),loader.getCenterOfBox()); <br>
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
class FFmaScanfLoader : public FAbstractLoader<ParticleClass> {
protected:
    FILE* file;                 //< The file to read
    F3DPosition centerOfBox;    //< The center of box read from file
    FReal boxWidth;             //< the box width read from file
    int nbParticles;            //< the number of particles read from file

public:
    /**
    * The constructor need the file name
    * @param filename the name of the file to open
    * you can test if file is successfuly open by calling isValide()
    */
    FFmaScanfLoader(const char* const filename): file(0){
        file = fopen(filename,"r");
        // test if open
        if(this->file){
            float x,y,z, fBoxWidth;
            const int nbReadElements = fscanf(file,"%d %f %f %f %f",&this->nbParticles,&fBoxWidth,&x,&y,&z);
            if(nbReadElements == 5){
                this->boxWidth = fBoxWidth;
                this->centerOfBox.setPosition(x,y,z);
                this->boxWidth *= 2;
            }
            else{
                fclose(file);
                file = NULL;
            }
        }
        else {
             this->boxWidth = 0;
             this->nbParticles = 0;
        }
    }

    /**
    * Default destructor, simply close the file
    */
    virtual ~FFmaScanfLoader(){
        fclose(file);
    }

    /**
      * To know if file is open and ready to read
      * @return true if loader can work
      */
    bool isValide() const{
        return this->file != NULL;
    }

    /**
      * To get the number of particles from this loader
      * @param the number of particles the loader can fill
      */
    long getNumberOfParticles() const{
        return this->nbParticles;
    }

    /**
      * The center of the box from the simulation file opened by the loader
      * @return box center
      */
    F3DPosition getCenterOfBox() const{
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
        if(this->file){
            float x,y,z,data;
            const int nbReadElements = fscanf(this->file,"%f %f %f %f",&x,&y,&z,&data);
            if(nbReadElements == 4){
                inParticle.setPosition(x,y,z);
                inParticle.setPhysicalValue(data);
            }
            else{
                fclose(this->file);
                this->file = NULL;
            }
        }
    }

};


#endif //FFMASCANFLOADER_HPP

// [--LICENSE--]
