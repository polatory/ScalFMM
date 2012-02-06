// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================
#ifndef FHLOADER_HPP
#define FHLOADER_HPP


#include <iostream>
#include <fstream>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FHLoader
* Please read the license
*
* Load a file with a format like :
* NB_particles Box_width Box_X Box_Y Box_Z // init
* X Y Z // one particle by line
* ....
* <code>
*    FHLoader<FBasicParticle> loader("../FMB++/Tests/particles.basic.txt"); <br>
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
class FHLoader : public FAbstractLoader<ParticleClass> {
protected:
    FILE* file;                 //< The file to read
    F3DPosition centerOfBox;    //< The center of box read from file
    FReal boxWidth;             //< the box width read from file
    int nbParticles;            //< the number of particles read from file

public:
    /**
    * The constructor need the file name
    * @param filename the name of the file to open
    * you can test if file is successfuly open by calling hasNotFinished()
    */
    FHLoader(const char* const filename): file(0){
        file = fopen(filename,"r");
        // test if open
        if(this->file){
            float x,y,z, fBoxWidth;
            const int nbReadElements = fscanf(file,"%d %f %f %f %f%*c",&this->nbParticles,&fBoxWidth,&x,&y,&z);
            if(nbReadElements == 5){
                this->boxWidth = fBoxWidth;
                this->centerOfBox.setPosition(x,y,z);

                char buff[512];
                char * ret = fgets(buff, 512, file);
                if(!ret){}
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
    virtual ~FHLoader(){
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
        return FSize(this->nbParticles);
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
            char buff[128];
            float x,y,z;
            const int nbReadElements = fscanf(this->file,"%s %f %f %f",buff,&x,&y,&z);
            if(nbReadElements == 4){
                inParticle.setPosition(x,y,z);
                inParticle.setData(buff[0]);
            }
            else{
                fclose(this->file);
                this->file = NULL;
            }
        }
    }

};


#endif //FHLOADER_HPP


