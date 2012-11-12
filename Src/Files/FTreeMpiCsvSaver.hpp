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
#ifndef FTREEMPICSVSAVER_HPP
#define FTREEMPICSVSAVER_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/FMpi.hpp"

#include <cstring>
#include <iostream>
#include <fstream>

/** This class is to export a tree in csv file
  *
  */
template <class OctreeClass, class ContainerClass , class ParticleClass>
class FTreeMpiCsvSaver {
    FMpi::FComm comm;           //< Communicator
    const bool includeHeader;   //< To include a line of header
    int nbFrames;               //< The current frame

    char basefile[512];         //< The base file name like "~/OUT/simulation%d.csv"
    char lineFormat[512];       //< The line format
    char header[512];           //< The header line
    char defaultLine[512];      //< The default line

    static const int SizeOfHeader = 8;
    static const int SizeOfLine = 56;

public:
	/** Constructor
	 * @param inBasefile is the output file name, you must put %d in it
	 * @param communicator Mpi communicator
	 * @param inIncludeHeader tells if header must be included
	 */
    FTreeMpiCsvSaver(const char inBasefile[], const FMpi::FComm& communicator, const bool inIncludeHeader = false)
        : comm(communicator.getComm()), includeHeader(inIncludeHeader), nbFrames(0) {
        strcpy(basefile, inBasefile);
        strcpy(lineFormat, "%+12.6e,%+12.6e,%+12.6e,%+12.6e\n"); // 56 = 4 x 14 chars
        strcpy(header, "x,y,z,v\n"); // 8 chars
    }

    /** Virtual destructor
      */
    virtual ~FTreeMpiCsvSaver(){
    }

    /** to know how many frame has been saved
      */
    int getNbFrames() const {
        return nbFrames;
    }

    /** export a tree
      */
    void exportTree(OctreeClass*const tree){
        char currentFilename[512];
        sprintf(currentFilename, basefile, nbFrames++);

        // Erase existing file
        if(comm.processId() == 0){
            FILE* fd = fopen(currentFilename, "w");
            if(fd) fclose(fd);
        }

        MPI_Barrier(comm.getComm());

        // All procs open the file
        MPI_File file;
        const int fileIsOpen = MPI_File_open( comm.getComm(), currentFilename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file );

        // Is it open?
        if(fileIsOpen != MPI_SUCCESS){
            std::cout << "Cannot create parallel file, exportTree abort." << std::endl;
            return;
        }

        MPI_Offset disp = 0;
        {
            // Each procs count its particles
            int mypart = 0;
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            do{
                mypart += octreeIterator.getCurrentListTargets()->getSize();
                const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());
                if( isUsingTsm ) mypart += octreeIterator.getCurrentListSrc()->getSize();
            } while(octreeIterator.moveRight());
            // Gather particles number
            int*const particlesPerProcs = new int[comm.processCount()];
            MPI_Allgather(&mypart, 1, MPI_INT, particlesPerProcs, 1, MPI_INT, comm.getComm());
            int previousPart = 0;
            for(int idxProc = 0 ; idxProc < comm.processId() ; ++idxProc){
                previousPart += particlesPerProcs[idxProc];
            }
            delete[] particlesPerProcs;
            // How many particle before me * size of particle in file
            disp = previousPart * SizeOfLine;
            if(includeHeader && comm.processId() != 0){
                disp += SizeOfHeader;
            }
        }

        // Set view
        const char* const native = "native";
        MPI_File_set_view( file, disp, MPI_BYTE, MPI_BYTE, const_cast<char*>(native), MPI_INFO_NULL);

        // Write header if needed
        MPI_Offset offset = 0;
        if(includeHeader && comm.processId() == 0){
            MPI_File_write_at( file, 0, header, SizeOfHeader, MPI_BYTE, MPI_STATUS_IGNORE);
            offset += SizeOfHeader;
        }

        // Write particle
        char line[SizeOfLine+1];
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            while( iter.hasNotFinished() ){
                sprintf(line, lineFormat, iter.data().getPosition().getX(),iter.data().getPosition().getY(),
                        iter.data().getPosition().getZ(),getValue(&iter.data()));

                MPI_File_write_at( file, offset, line, SizeOfLine, MPI_BYTE, MPI_STATUS_IGNORE);
                offset += SizeOfLine;

                iter.gotoNext();
            }

            const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());
            if( isUsingTsm ){
                typename ContainerClass::BasicIterator iterSource(*octreeIterator.getCurrentListSrc());

                while( iterSource.hasNotFinished() ){
                    sprintf(line, lineFormat, iterSource.data().getPosition().getX(),iterSource.data().getPosition().getY(),
                            iterSource.data().getPosition().getZ(),getValue(&iterSource.data()));

                    MPI_File_write_at( file, offset, line, SizeOfLine, MPI_BYTE, MPI_STATUS_IGNORE);
                    offset += SizeOfLine;

                    iterSource.gotoNext();
                }
            }
        } while(octreeIterator.moveRight());

        MPI_File_close( &file );
    }

    /** Inherit from this class and customize this function if you need it
      */
    virtual FReal getValue(ParticleClass*const part){
        return FReal(0);
    }
};


#endif // FTREEMPICSVSAVER_HPP
