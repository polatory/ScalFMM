#ifndef FFMALOADER_HPP
#define FFMALOADER_HPP
// [--License--]

#include <iostream>
#include <fstream>
#include <limits>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FEwalLoader
* Please read the license
* Particle has to extend {FExtendPhysicalValue,FExtendPosition}
*/
template <class ParticleClass>
class FEwalLoader : public FAbstractLoader<ParticleClass> {
protected:
    std::ifstream file;         //< The file to read
    F3DPosition centerOfBox;    //< The center of box read from file
    FReal boxWidth;             //< the box width read from file
    int nbParticles;            //< the number of particles read from file

public:
    /**
    * The constructor need the file name
    * @param filename the name of the file to open
    * you can test if file is successfuly open by calling hasNotFinished()
        Box SPC water from DLPOLY TEST17
         2         1       417  -591626.141968
     17.200000000000      0.000000000000      0.000000000000
      0.000000000000     17.200000000000      0.000000000000
      0.000000000000      0.000000000000     17.200000000000
    */
    FEwalLoader(const char* const filename): file(filename,std::ifstream::in){
        // test if open
        if(this->file.is_open()){
            const int bufferSize = 512;
            char buffer[bufferSize];
            file.getline(buffer, bufferSize);

            int tempi(0);
            FReal tempf(0);
            file >> tempi >> tempi >> this->nbParticles >> tempf;

            FReal widthx, widthy, widthz;
            file >> widthx >> tempf >> tempf;
            file >> tempf >> widthy >> tempf;
            file >> tempf >> tempf >> widthz;

            this->centerOfBox.setPosition(0.0,0.0,0.0);
            this->boxWidth = widthx;
        }
        else {
             this->boxWidth = 0;
             this->nbParticles = 0;
        }
    }

    /**
    * Default destructor, simply close the file
    */
    virtual ~FEwalLoader(){
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
        OW               1
             5.447823189       -0.4124521286        -3.845403447
           7.64746800518      -1.34490700206      -2.81036521708
          -4406.48579000       6815.52906417       10340.2577024
      */
    void fillParticle(ParticleClass& inParticle){
        FReal x, y, z, fx, fy, fz, vx, vy, vz;
        int index;
        char type[2];

        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        file.read(type, 2);
        file >> index;
        file >> x >> y >> z;
        file >> vx >> vy >> vz;
        file >> fx >> fy >> fz;

        inParticle.setPosition(x,y,z);
        inParticle.setForces(fx,fy,fz);
        inParticle.setIndex(index);

        if( strncmp(type, "OW", 2) == 0){
            inParticle.setPhysicalValue(FReal(-0.82));
            inParticle.setType(ParticleClass::OW);
        }
        else{
            inParticle.setPhysicalValue(FReal(0.41));
            inParticle.setType(ParticleClass::HW);
        }
    }

};


#endif //FEwalLoader_HPP

// [--END--]
