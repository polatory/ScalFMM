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
#ifndef FEWALLOADER_HPP
#define FEWALLOADER_HPP


#include <iostream>
#include <fstream>
#include <limits>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/FPoint.hpp"


template <class BaseClass>
class FEwalParticle : public BaseClass {
public:
    // Type of particle
    enum Type{
        OW,
        HW,
        Undefined,
    };

private:
    Type type; //< current type
    int index; //< current index in array
    int indexInFile; //< current index in array

public:
    // Basic constructor
    FEwalParticle() : type(Undefined), index(-1), indexInFile(-1) {
    }

    Type getType() const{
        return type;
    }

    void setType(const Type inType) {
        type = inType;
    }

    int getIndex() const{
        return index;
    }

    void setIndex( const int inIndex ){
        index = inIndex;
    }

    int getIndexInFile() const{
        return indexInFile;
    }

    void setIndexInFile( const int inIndex ){
        indexInFile = inIndex;
    }
};



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
    FPoint centerOfBox;    //< The center of box read from file
    FReal boxWidth;             //< the box width read from file
    int nbParticles;            //< the number of particles read from file
    int levcfg  ;         //< DL_POLY CONFIG file key. 0,1 or 2
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

            int imcon ;
            //int tempi(0);
            FReal tempf(0);
            file >> levcfg >> imcon >> this->nbParticles;
            // Periodic case
            if( imcon > 0 ) {
                FReal widthx, widthy, widthz;
                file >> widthx >> tempf >> tempf;
                file >> tempf >> widthy >> tempf;
                file >> tempf >> tempf >> widthz;

                this->boxWidth = widthx;
            }
            // Non periodic case
            else{
                file >> this->boxWidth;
            }
            this->centerOfBox.setPosition(0.0,0.0,0.0);
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
        OW               1
             5.447823189       -0.4124521286        -3.845403447
           7.64746800518      -1.34490700206      -2.81036521708
          -4406.48579000       6815.52906417       10340.2577024
      */
    void fillParticle(ParticleClass& inParticle){
        FReal x, y, z, fx, fy, fz, vx, vy, vz;
        int index;
        char type[2];
        std::string line;
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        file.read(type, 2);
        file >> index;
        std::getline(file, line); // needed to skip the end of the line in non periodic case
        if ( levcfg == 0) {
           file >> x >> y >> z;
        }else if ( levcfg == 1) {
            file >> x >> y >> z;
            file >> vx >> vy >> vz;
        }else {
            file >> x >> y >> z;
            file >> vx >> vy >> vz;
            file >> fx >> fy >> fz;
        }


        //	std::cout << " x >> y >> z: " << x<< " " <<y<< " " <<z <<std::endl;
        inParticle.setPosition(x,y,z);
        inParticle.setForces(fx,fy,fz);
        //inParticle.setForces(vx,vy,vz);
        inParticle.setIndexInFile(index);

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


template <class ParticleClass>
class FEwalBinLoader : public FAbstractLoader<ParticleClass> {
protected:
    FILE* const file;         //< The file to read
    FPoint centerOfBox;    //< The center of box read from file
    double boxWidth;             //< the box width read from file
    int nbParticles;            //< the number of particles read from file
    double energy;
    int removeWarning;

    template<class Type>
    Type readValue(){
        int sizeBefore, sizeAfter;
        Type value;
        removeWarning = fread(&sizeBefore, sizeof(int), 1, file);
        removeWarning = fread(&value, sizeof(Type), 1, file);
        removeWarning = fread(&sizeAfter, sizeof(int), 1, file);
        if( sizeBefore != sizeof(Type) ) printf("Error in loader ewal Size before %d should be %d\n", sizeBefore, sizeof(Type));
        if( sizeAfter != sizeof(Type) ) printf("Error in loader ewal Size after %d should be %d\n", sizeAfter, sizeof(Type));
        return value;
    }

    template<class Type>
    Type* readArray(Type array[], const int size){
        int sizeBefore, sizeAfter;
        removeWarning = fread(&sizeBefore, sizeof(int), 1, file);
        removeWarning = fread(array, sizeof(Type), size, file);
        removeWarning = fread(&sizeAfter, sizeof(int), 1, file);
        if( sizeBefore != int(sizeof(Type) * size) ) printf("Error in loader ewal Size before %d should be %d\n", sizeBefore, size*sizeof(Type));
        if( sizeAfter != int(sizeof(Type) * size) ) printf("Error in loader ewal Size after %d should be %d\n", sizeAfter, size*sizeof(Type));
        return array;
    }

public:
    /**
    * The constructor need the file name
    * @param filename the name of the file to open
    * you can test if file is successfuly open by calling hasNotFinished()
        energy box size nb particles
        [index charge x y z fx fy fz]
        int double double ...
    */
    FEwalBinLoader(const char* const filename): file(fopen(filename, "rb")) {
        // test if open
        if(this->file != NULL){
            energy = readValue<double>();
            double boxDim[3];
            boxWidth = readArray<double>(boxDim,3)[0];
            nbParticles = readValue<int>();

            centerOfBox.setPosition(0.0,0.0,0.0);
        }
        else {
            this->boxWidth = 0;
            this->nbParticles = 0;
        }
    }
    /**
    * Default destructor, simply close the file
    */
    virtual ~FEwalBinLoader(){
        fclose(file);
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

    FReal getEnergy() const{
        return this->energy;
    }

    /**
      * Fill a particle
      * @warning to work with the loader, particles has to expose a setPosition method
      * @param the particle to fill
        [index charge x y z fx fy fz]
      */
    void fillParticle(ParticleClass& inParticle){
        double x, y, z, fx, fy, fz, charge;
        int index;

        int size;
        removeWarning = fread(&size, sizeof(int), 1, file);
        if(size != 60) printf("Error in loader ewal Size %d should be %d\n", size, 60);

        removeWarning = fread(&index, sizeof(int), 1, file);
        removeWarning = fread(&charge, sizeof(double), 1, file);

        removeWarning = fread(&x, sizeof(double), 1, file);
        removeWarning = fread(&y, sizeof(double), 1, file);
        removeWarning = fread(&z, sizeof(double), 1, file);

        removeWarning = fread(&fx, sizeof(double), 1, file);
        removeWarning = fread(&fy, sizeof(double), 1, file);
        removeWarning = fread(&fz, sizeof(double), 1, file);

        removeWarning = fread(&size, sizeof(int), 1, file);
        if(size != 60) printf("Error in loader ewal Size %d should be %d\n", size, 60);

        inParticle.setPosition(x,y,z);
        inParticle.setForces(fx,fy,fz);
        inParticle.setIndexInFile(index);
        inParticle.setPhysicalValue(charge);
    }

};


#endif //FEwalLoader_HPP


