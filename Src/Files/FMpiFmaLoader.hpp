#ifndef FMPIFMALOADER_HPP
#define FMPIFMALOADER_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/F3DPosition.hpp"
#include "../Utils/FMpi.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FMpiFmaLoader
* Please read the license
*
* Load a file with a format like :
* NB_particles Box_width Box_X Box_Y Box_Z // init
* X Y Z // one particle by line
* ....
* <code>
*    FMpiFmaLoader<FBasicParticle> loader("../FMB++/Tests/particles.basic.txt"); <br>
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
class FMpiFmaLoader : public FAbstractLoader<ParticleClass> {
protected:
    F3DPosition centerOfBox;    //< The center of box read from file
    FReal boxWidth;             //< the box width read from file
    int totalNbParticles;       //< the number of particles read from file
    int nbParticles;            //< the number of particles read from file
    bool isOpenFlag;            //< to knwo if the file is open now
    FReal* particles;           //< the particles loaded from the binary file
    MPI_Offset idxParticles;    //< to iterate on the particles array

public:
    /**
    * The constructor need the file name
    * @param filename the name of the file to open
    * you can test if file is successfuly open by calling hasNotFinished()
    */
    FMpiFmaLoader(const char* const filename, FMpi& app)
            : boxWidth(0), totalNbParticles(0), nbParticles(0), isOpenFlag(false), particles(0), idxParticles(0) {
        char nonConstFilename[512];
        strcpy(nonConstFilename,filename);
        MPI_File file;
        if(MPI_File_open(MPI::COMM_WORLD, nonConstFilename, MPI::MODE_RDONLY, MPI::INFO_NULL, &file) == MPI_SUCCESS){
            int sizeOfElement(0);
            FReal xyzBoxWidth[4];

            MPI_Status status;
            if( MPI_File_read(file, &sizeOfElement, 1, MPI_INT, &status) == MPI_SUCCESS
                && MPI_File_read(file, &this->totalNbParticles, 1, MPI_INT, &status) == MPI_SUCCESS
                && MPI_File_read(file, xyzBoxWidth, 4, MPI_FLOAT, &status) == MPI_SUCCESS ){

                FDEBUG(if(sizeOfElement != sizeof(FReal)){)
                    FDEBUG( FDebug::Controller.writeFromLine("Warning type size between file and FReal are differents\n", __LINE__, __FILE__); )
                FDEBUG(})

                this->boxWidth = xyzBoxWidth[3];
                this->centerOfBox.setPosition(xyzBoxWidth[0],xyzBoxWidth[1],xyzBoxWidth[2]);
                this->boxWidth *= 2;
                this->isOpenFlag = true;

                // load my particles
                MPI_Offset headDataOffSet;
                MPI_File_get_position(file, &headDataOffSet);

                MPI_Offset filesize(0);
                MPI_File_get_size(file, &filesize); /* in bytes */
                filesize = (filesize - headDataOffSet) / sizeof(FReal);
                if(filesize/4 != this->totalNbParticles){
                    printf("Error fileSize %lld, nbPart %d\n",filesize/4, this->totalNbParticles);
                }
                // in number of floats
                const long startPart = app.getLeft(this->totalNbParticles);
                const long endPart = app.getRight(this->totalNbParticles);
                nbParticles = (endPart - startPart);
                const int bufsize = nbParticles * 4;
                // local number to read
                particles = new FReal[bufsize];

                MPI_File_read_at(file, headDataOffSet + startPart * 4 * sizeof(FReal), particles, bufsize, MPI_FLOAT, &status);

                // check if needed
                int count(0);
                MPI_Get_count(&status, MPI_INT, &count);
                FDEBUG(if(count  / 4 != this->nbParticles){)
                    FDEBUG( FDebug::Controller<< "Error read " << count << " data, nbPart is " << this->nbParticles << __LINE__ << " " << __FILE__ << "\n"; )
                FDEBUG(})
            }
            else{
                this->totalNbParticles = 0;
            }
            MPI_File_close(&file);
        }
        else {
             this->boxWidth = 0;
             this->totalNbParticles = 0;
        }
    }

    /**
    * Default destructor, simply close the file
    */
    virtual ~FMpiFmaLoader(){
        if(isOpen()){
            delete [] particles;
        }
    }

    /**
      * To know if file is open and ready to read
      * @return true if loader can work
      */
    bool isOpen() const{
        return this->isOpenFlag;
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
        inParticle.setPosition(particles[idxParticles],particles[idxParticles+1],particles[idxParticles+2]);
        inParticle.setPhysicalValue(particles[idxParticles+3]);
        idxParticles += 4;
    }

};


#endif //FMPIFMALOADER_HPP

// [--LICENSE--]
