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
#ifndef FMPIFMALOADER_HPP
#define FMPIFMALOADER_HPP


#include <iostream>
#include <fstream>

#include "../Utils/FGlobal.hpp"
#include "FAbstractLoader.hpp"
#include "../Utils/FPoint.hpp"
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
* @code
*    FMpiFmaLoader<FBasicParticle> loader("../ADir/Tests/particles.basic.txt"); <br>
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
class FMpiFmaLoader : public FAbstractLoader {
protected:
    FPoint centerOfBox;    //< The center of box read from file
    FReal boxWidth;             //< the box width read from file
    FSize totalNbParticles;     //< the number of particles read from file
    FSize nbParticles;          //< the number of particles read from file
    bool isOpenFlag;            //< to knwo if the file is open now
    FReal* particles;           //< the particles loaded from the binary file
    MPI_Offset idxParticles;    //< to iterate on the particles array

public:
    /**
    * The constructor need the file name
    * @param filename the name of the file to open
    * you can test if file is successfuly open by calling hasNotFinished()
    */
    FMpiFmaLoader(const char* const filename, const FMpi::FComm& comm, const bool useMpiIO = false)
            : boxWidth(0), totalNbParticles(0), nbParticles(0), isOpenFlag(false), particles(0), idxParticles(0) {
        if( useMpiIO ){
            char nonConstFilename[512];
            strcpy(nonConstFilename,filename);
            MPI_File file;
            if(MPI_File_open(comm.getComm(), nonConstFilename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file) == MPI_SUCCESS){
                int sizeOfElement(0);
                FReal xyzBoxWidth[4];

                MPI_Status status;
                if( MPI_File_read(file, &sizeOfElement, 1, MPI_INT, &status) == MPI_SUCCESS
                    && MPI_File_read(file, &this->totalNbParticles, 1, MPI_LONG_LONG, &status) == MPI_SUCCESS
                    && MPI_File_read(file, xyzBoxWidth, 4, MPI_FLOAT, &status) == MPI_SUCCESS ){

                    FLOG(if(sizeOfElement != sizeof(FReal)){)
                        FLOG( FLog::Controller.writeFromLine("Warning type size between file and FReal are differents\n", __LINE__, __FILE__); )
                    FLOG(})

                    this->boxWidth = xyzBoxWidth[3];
                    this->centerOfBox.setPosition(xyzBoxWidth[0],xyzBoxWidth[1],xyzBoxWidth[2]);
                    this->boxWidth *= 2;
                    this->isOpenFlag = true;

                    // load my particles
                    MPI_Offset headDataOffSet(0);
                    MPI_File_get_position(file, &headDataOffSet);

                    MPI_Offset filesize(0);
                    MPI_File_get_size(file, &filesize); // in bytes
                    filesize = (filesize - headDataOffSet) / sizeof(FReal);
                    if(filesize/4 != this->totalNbParticles){
                        printf("Error fileSize %lld, nbPart %lld\n",filesize/4, this->totalNbParticles);
                    }
                    // in number of floats
                    const FSize startPart = comm.getLeft(this->totalNbParticles);
                    const FSize endPart   = comm.getRight(this->totalNbParticles);
                    nbParticles = (endPart - startPart);
                    const FSize bufsize = nbParticles * 4;
                    // local number to read
                    particles = new FReal[bufsize];

                    if( sizeof(FReal) == sizeof(float) ){
                        MPI_File_read_at(file, headDataOffSet + startPart * 4 * sizeof(FReal), particles, int(bufsize), MPI_FLOAT, &status);
                    }
                    else{
                        MPI_File_read_at(file, headDataOffSet + startPart * 4 * sizeof(FReal), particles, int(bufsize), MPI_DOUBLE, &status);
                    }


                    // check if needed
                    int count(0);
                    MPI_Get_count(&status, MPI_INT, &count);
                    FLOG(if(count  / 4 != this->nbParticles){)
                        FLOG( FLog::Controller<< "Error read " << count << " data, nbPart is " << this->nbParticles << __LINE__ << " " << __FILE__ << "\n"; )
                    FLOG(})
                }
                else{
                    this->totalNbParticles = 0;
                }
                MPI_File_close(&file);
            }
        }
        else {
            FILE* file(fopen(filename, "rb"));
            size_t removeWarning(0);
            // test if open
            if(file != NULL) {
                int sizeOfElement(0);
                removeWarning += fread(&sizeOfElement, sizeof(int), 1, file);
                FLOG(if(sizeOfElement != int(sizeof(FReal)) ){)
                    FLOG( FLog::Controller.writeFromLine("Warning type size between file and FReal are differents\n", __LINE__, __FILE__); )
                FLOG(})
                removeWarning += fread(&this->totalNbParticles, sizeof(FSize), 1, file);

                removeWarning += fread(&this->boxWidth, sizeof(FReal), 1, file);
                this->boxWidth *= 2;

                FReal x,y,z;
                removeWarning += fread(&x, sizeof(FReal), 1, file);
                removeWarning += fread(&y, sizeof(FReal), 1, file);
                removeWarning += fread(&z, sizeof(FReal), 1, file);
                this->centerOfBox.setPosition(x,y,z);

                this->isOpenFlag = true;

                const long int headDataOffSet = ftell(file);
                fseek(file, 0L, SEEK_END);
                const long int filesize = (ftell(file) - headDataOffSet) / sizeof(FReal);

                if(filesize/4 != this->totalNbParticles){
                    printf("Error fileSize %ld, nbPart %lld\n", filesize/4, this->totalNbParticles);
                }

                // in number of floats
                const FSize startPart = comm.getLeft(this->totalNbParticles);
                const FSize endPart   = comm.getRight(this->totalNbParticles);
                nbParticles = (endPart - startPart);
                const FSize bufsize = nbParticles * 4;
                // local number to read
                particles = new FReal[bufsize];

                fseek(file, long(headDataOffSet + startPart * 4 * sizeof(FReal)), SEEK_SET);

                removeWarning += fread(particles, sizeof(FReal), int(bufsize), file);

                fclose(file);
            }
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
    void fillParticle(FPoint*const inParticlePositions, FReal*const inPhysicalValue){
        inParticlePositions->setPosition(particles[idxParticles],particles[idxParticles+1],particles[idxParticles+2]);
        (*inPhysicalValue) = (particles[idxParticles+3]);
        idxParticles += 4;
    }

};


#endif //FMPIFMALOADER_HPP


