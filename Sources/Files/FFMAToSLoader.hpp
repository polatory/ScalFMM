#ifndef FFMATOSLOADER_HPP
#define FFMATOSLOADER_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFMAToSLoader
* Please read the license
*
* Load a file with a format like :
* NB_particules Box_width Box_X Box_Y Box_Z // init
* X Y Z // one particule by line
* ....
* <code>
*    FFMAToSLoader<FBasicParticule> loader("../FMB++/Tests/particules.basic.txt"); <br>
*    if(!loader.isValide()){ <br>
*        std::cout << "Loader Error\n"; <br>
*        return 1; <br>
*    } <br>
* <br>
*    FOctree<FBasicParticule, TestCell, FSimpleLeaf, 10, 3> tree(loader.getBoxWidth(),loader.getCenterOfBox()); <br>
* <br>
*    for(int idx = 0 ; idx < loader.getNumberOfParticules() ; ++idx){ <br>
*        FBasicParticule* const part = new FBasicParticule(); <br>
*        loader.fillParticule(part); <br>
*        tree.insert(part); <br>
*    } <br>
* </code>
*
* Particule has to extend {FExtendPhysicalValue,FExtendPosition}
*/
template <class ParticuleClass>
class FFMAToSLoader : public FAbstractLoader<ParticuleClass> {
protected:
    std::ifstream file;         //< The file to read
    F3DPosition centerOfBox;    //< The center of box read from file
    FReal boxWidth;             //< the box width read from file
    int nbParticules;           //< the number of particules read from file

public:
    /**
    * The constructor need the file name
    * @param filename the name of the file to open
    * you can test if file is successfuly open by calling isValide()
    */
    FFMAToSLoader(const char* const filename): file(filename,std::ifstream::in){
        // test if open
        if(this->file.is_open()){
            FReal x,y,z;
            this->file >> this->nbParticules >> this->boxWidth >> x >> y >> z;
            this->centerOfBox.setPosition(x,y,z);
            this->boxWidth *= 2;
        }
        else {
             this->boxWidth = 0;
             this->nbParticules = 0;
        }
    }

    /**
    * Default destructor, simply close the file
    */
    virtual ~FFMAToSLoader(){
        file.close();
    }

    /**
      * To know if file is open and ready to read
      * @return true if loader can work
      */
    bool isValide() const{
        return this->file.is_open() && !this->file.eof();
    }

    /**
      * To get the number of particules from this loader
      * @param the number of particules the loader can fill
      */
    long getNumberOfParticules() const{
        return this->nbParticules;
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
      * Fill a particule
      * @warning to work with the loader, particules has to expose a setPosition method
      * @param the particule to fill
      */
    void fillParticule(ParticuleClass* const inParticule){
        FReal x,y,z,data;
        int isTarget;
        this->file >> x >> y >> z >> data >> isTarget;
        inParticule->setPosition(x,y,z);
        inParticule->setPhysicalValue(data);
        if(isTarget) inParticule->setAsTarget();
        else inParticule->setAsSource();
    }

};


#endif //FFMATOSLOADER_HPP

// [--LICENSE--]
