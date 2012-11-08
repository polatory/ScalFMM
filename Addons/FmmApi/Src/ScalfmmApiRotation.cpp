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
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Kernels/Rotation/FRotationCell.hpp"
#include "../../Src/Kernels/Rotation/FRotationParticle.hpp"
#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"

#include "FmmApi.h"

////////////////////// Opérateurs FMM Kernel : //////////////////////////

class KernelCell : public FBasicCell {
    FComplexe* multipole;
    FComplexe* local;
public:
    KernelCell() : multipole(0), local(0){
    }
    void attachArrays(FComplexe inMultipole[], FComplexe inLocal[]){
        multipole = inMultipole;
        local = inLocal;
    }

    const FComplexe* getMultipole() const{
        return multipole;
    }

    const FComplexe* getLocal() const{
        return local;
    }

    FComplexe* getMultipole(){
        return multipole;
    }

    FComplexe* getLocal(){
        return local;
    }
};


static const int P = 5;

typedef FRotationParticle      KernelParticleClass;
typedef FVector<KernelParticleClass>        KernelContainerClass;
typedef KernelCell               KernelCellClass;
typedef FRotationKernel<KernelParticleClass, KernelCellClass, KernelContainerClass, P >         KernelClass;

struct ScalFmmKernelHandle {
    KernelClass** kernel;
    int potentialDataSize;
    int fieldDataSize;
    int nbthread;
};

/******* Allocation : ******/
int FmmKernel_init(void *fmmCore, void **fmmKernel){
    ScalFmmKernelHandle* kernelhandle = new ScalFmmKernelHandle;
    memset(kernelhandle, 0, sizeof(ScalFmmKernelHandle));

    int NbLevels;
    FmmCore_getParameter(fmmCore, FMMCORE_TREE_HEIGHT, &NbLevels);
    FReal boxWidth;
    FmmCore_getParameter(fmmCore, FMMCORE_ROOT_BOX_WIDTH, &boxWidth);
    FReal centerOfBox[3];
    FmmCore_getParameter(fmmCore, FMMCORE_ROOT_BOX_CENTER, centerOfBox);

    FmmCore_getParameter(fmmCore, FMMCORE_THREADS_NUMBER, &kernelhandle->nbthread);

    KernelClass original( NbLevels, boxWidth, FPoint(centerOfBox) );
    kernelhandle->kernel = new KernelClass*[kernelhandle->nbthread];
    for(int idxThread = 0 ; idxThread < kernelhandle->nbthread ; ++idxThread){
        kernelhandle->kernel[idxThread] = new KernelClass(original);
    }

    kernelhandle->potentialDataSize = 1;
    kernelhandle->fieldDataSize = 4;

    *fmmKernel = kernelhandle;

    return FMMAPI_NO_ERROR;
}/* : alloue et initialise le FmmKernel */
int FmmKernel_free(void *fmmKernel){
    ScalFmmKernelHandle* kernelhandle = (ScalFmmKernelHandle*) fmmKernel;

    for(int idxThread = 0 ; idxThread < kernelhandle->nbthread ; ++idxThread){
        delete kernelhandle->kernel[idxThread];
    }

    delete[] kernelhandle->kernel;
    delete kernelhandle;

    return FMMAPI_NO_ERROR;
} /* libére le FmmKernel */




/******* Configuration : ***/

int FmmKernel_isParameterUsed(void * /*fmm*/, int *name, int *flag){
    switch( *name ){
    case FMMKERNEL_POTENTIAL_DATA_SIZE:
    case FMMKERNEL_FIELD_DATA_SIZE:
    case  FMMKERNEL_HANDLES_P2P:
        *flag = FMMAPI_SUPPORTED_PARAMETER;
        break;
    case FMMKERNEL_ACCURACY :
        *flag = FMMAPI_UNSUPPORTED_PARAMETER;
        break;
    default:
        *flag = FMMAPI_UNKNOWN_PARAMETER;
    }

    return FMMAPI_NO_ERROR;
}

int FmmKernel_setParameter(void *fmmKernel, int *name, void*value){
    /*ScalFmmKernelHandle* kernelhandle = (ScalFmmKernelHandle*) fmmKernel;*/
    int flag;

    FmmKernel_isParameterUsed(fmmKernel, name, &flag);
    if( flag != FMMAPI_SUPPORTED_PARAMETER){
        return flag;
    }

    switch( *name ){
    case FMMKERNEL_POTENTIAL_DATA_SIZE :
    case FMMKERNEL_FIELD_DATA_SIZE :
        return FMMAPI_SUPPORTED_PARAMETER;
    default:
        return FMMAPI_UNKNOWN_PARAMETER;
    }

    return FMMAPI_NO_ERROR;
}

int FmmKernel_setParameter(void *fmmKernel, int name, void*value){
    return FmmKernel_setParameter( fmmKernel, &name, value);
}

int FmmKernel_getParameter(void *fmmKernel, int *name, void*value){
    ScalFmmKernelHandle* kernelhandle = (ScalFmmKernelHandle*) fmmKernel;
    int flag;

    FmmKernel_isParameterUsed(fmmKernel, name, &flag);
    if( flag != FMMAPI_SUPPORTED_PARAMETER){
        return flag;
    }

    switch( *name ){
    case FMMKERNEL_POTENTIAL_DATA_SIZE :
        *(int*)value = kernelhandle->potentialDataSize*sizeof(FReal);
        break;
    case FMMKERNEL_FIELD_DATA_SIZE :
        *(int*)value = kernelhandle->fieldDataSize*sizeof(FReal);
        break;
    default:
        return FMMAPI_UNKNOWN_PARAMETER;
    }

    return FMMAPI_NO_ERROR;
}

int FmmKernel_getParameter(void *fmmKernel, int name, void*value){
    return FmmKernel_getParameter(fmmKernel, &name, value);
}


/****** Données FMM : *****/
int FmmKernel_getMultipoleArraySize(void */*fmmCore*/, int *size) {
    *size = ((P+2)*(P+1))/2 * sizeof(FComplexe);
    return FMMAPI_NO_ERROR;
} /* Renvoie dans size la taille (en octets) de l'expansion multipôle associée à la boîte boxId */

int FmmKernel_getLocalArraySize(void */*fmmCore*/, int *size){
    *size = ((P+2)*(P+1))/2 * sizeof(FComplexe);
    return FMMAPI_NO_ERROR;
} /* Renvoie dans size la taille (en octets) de l'expansion locale associée à la boîte boxId*/


/******* Opérateurs FMM : **/
int FmmKernel_P2M(void *fmmCore, void* boxId){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(corehandle, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplexe* multipole;
    FmmCore_getMultipoleArray(fmmCore, boxId, (void**)&multipole);

    KernelCellClass cell;
    cell.attachArrays(multipole, 0);
    int coord[3];
    FmmCore_getCoord(fmmCore, boxId, coord);
    cell.setCoordinate(coord[0], coord[1], coord[2]);

    FReal* positions;
    FReal* potentials;
    int number;
    FmmCore_getSource(fmmCore, boxId, &positions, (void**)&potentials, &number);

    KernelContainerClass sources;
    KernelParticleClass part;
    for(int idxPart = 0 ; idxPart < number ; ++idxPart){
        part.setPosition(FPoint(&positions[idxPart*3]));
        part.setPhysicalValue(potentials[idxPart]);
        sources.push(part);
    }

    kernelhandle->kernel[threadId]->P2M(&cell, &sources);

    FmmCore_releaseSource(fmmCore, boxId, potentials, positions);

    return FMMAPI_NO_ERROR;
}

int FmmKernel_L2P(void *fmmCore, void* boxId){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(corehandle, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplexe* local;
    FmmCore_getLocalArray(fmmCore, boxId, (void**)&local);

    KernelCellClass cell;
    cell.attachArrays(0,local);
    int coord[3];
    FmmCore_getCoord(fmmCore, boxId, coord);
    cell.setCoordinate(coord[0], coord[1], coord[2]);

    FReal* potentials;
    FReal* positions;
    int number;
    //FmmCore_getTargetPoints(fmmCore, boxId, &positions, &number);
    FmmCore_getSource(fmmCore, boxId, &positions, (void**)&potentials, &number);

    FReal* fields = new FReal[number*kernelhandle->fieldDataSize];
    FmmCore_getTargetField(fmmCore, boxId, fields);

    KernelContainerClass targets;
    KernelParticleClass part;
    for(int idxPart = 0 ; idxPart < number ; ++idxPart){
        part.setPosition(FPoint(&positions[idxPart*3]));
        part.setForces(fields[idxPart*kernelhandle->fieldDataSize],fields[idxPart*kernelhandle->fieldDataSize+1],
                       fields[idxPart*kernelhandle->fieldDataSize+2]);
        part.setPotential(fields[idxPart*kernelhandle->fieldDataSize+3]);
        part.setPhysicalValue(potentials[idxPart]);
        targets.push(part);
    }

    kernelhandle->kernel[threadId]->L2P(&cell, &targets);

    for(int idxPart = 0 ; idxPart < number ; ++idxPart){
        fields[idxPart*kernelhandle->fieldDataSize] = targets[idxPart].getForces().getX();
        fields[idxPart*kernelhandle->fieldDataSize+1] = targets[idxPart].getForces().getY();
        fields[idxPart*kernelhandle->fieldDataSize+2] = targets[idxPart].getForces().getZ();
        fields[idxPart*kernelhandle->fieldDataSize+3] = targets[idxPart].getPotential();
    }

    FmmCore_releaseTargetPoints(fmmCore, boxId, positions);
    FmmCore_setTargetField(fmmCore, boxId, fields);
    delete[] fields;

    return FMMAPI_NO_ERROR;
}

int FmmKernel_M2M(void *fmmCore, void *boxIdFather, void *boxIdSon){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(corehandle, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplexe* multipole;
    FmmCore_getMultipoleArray(fmmCore, boxIdFather, (void**)&multipole);

    KernelCellClass cellFather;
    cellFather.attachArrays(multipole, 0);
    int coordFather[3];
    FmmCore_getCoord(fmmCore, boxIdFather, coordFather);
    cellFather.setCoordinate(coordFather[0], coordFather[1], coordFather[2]);

    FmmCore_getMultipoleArray(fmmCore, boxIdSon, (void**)&multipole);

    KernelCellClass cellSon;
    cellSon.attachArrays(multipole, 0);
    int coordChild[3];
    FmmCore_getCoord(fmmCore, boxIdSon, coordChild);
    cellSon.setCoordinate(coordChild[0], coordChild[1], coordChild[2]);

    int level;
    FmmCore_getLevel(fmmCore,boxIdFather, &level);

    const KernelCellClass* children[8];
    memset(children, 0, sizeof(KernelCellClass*)*8);
    const int mindex = ((coordChild[0]&1) * 2 + (coordChild[1]&1)) * 2 + (coordChild[2]&1);
    children[mindex] = &cellSon;

    kernelhandle->kernel[threadId]->M2M(&cellFather, children, level);

    return FMMAPI_NO_ERROR;
}

int FmmKernel_L2L(void *fmmCore, void *boxIdFather, void *boxIdSon){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(corehandle, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplexe* local;
    FmmCore_getLocalArray(fmmCore, boxIdFather, (void**)&local);

    KernelCellClass cellFather;
    cellFather.attachArrays(0, local);
    int coordFather[3];
    FmmCore_getCoord(fmmCore, boxIdFather, coordFather);
    cellFather.setCoordinate(coordFather[0], coordFather[1], coordFather[2]);

    FmmCore_getLocalArray(fmmCore, boxIdSon, (void**)&local);

    KernelCellClass cellSon;
    cellSon.attachArrays(0, local);
    int coordChild[3];
    FmmCore_getCoord(fmmCore, boxIdSon, coordChild);
    cellSon.setCoordinate(coordChild[0], coordChild[1], coordChild[2]);

    int level;
    FmmCore_getLevel(fmmCore,boxIdFather, &level);

    KernelCellClass* children[8];
    memset(children, 0, sizeof(KernelCellClass*)*8);
    const int mindex = ((coordChild[0]&1) * 2 + (coordChild[1]&1)) * 2 + (coordChild[2]&1);
    children[mindex] = &cellSon;

    kernelhandle->kernel[threadId]->L2L(&cellFather, children, level);

    return FMMAPI_NO_ERROR;
}

int FmmKernel_M2L(void *fmmCore, void *boxIdSrc, void *boxIdDest){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(corehandle, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FComplexe* multipole;
    FmmCore_getMultipoleArray(fmmCore, boxIdSrc, (void**)&multipole);
    KernelCellClass cellSrc;
    cellSrc.attachArrays(multipole,0);
    int coord[3];
    FmmCore_getCoord(fmmCore, boxIdSrc, coord);
    cellSrc.setCoordinate(coord[0], coord[1], coord[2]);

    FComplexe* local;
    FmmCore_getLocalArray(fmmCore, boxIdDest, (void**)&local);
    KernelCellClass cellDst;
    cellDst.attachArrays(0, local);
    FmmCore_getCoord(fmmCore, boxIdDest, coord);
    cellDst.setCoordinate(coord[0], coord[1], coord[2]);

    int level;
    FmmCore_getLevel(fmmCore, boxIdDest, &level);

    const int xdiff = cellSrc.getCoordinate().getX() - cellDst.getCoordinate().getX();
    const int ydiff = cellSrc.getCoordinate().getY() - cellDst.getCoordinate().getY();
    const int zdiff = cellSrc.getCoordinate().getZ() - cellDst.getCoordinate().getZ();
    const int index = (((xdiff+3) * 7) + (ydiff+3)) * 7 + zdiff + 3;

    const KernelCellClass* inter[343];
    memset(inter, 0, sizeof(KernelCellClass*)*343);
    inter[index] = &cellSrc;

    kernelhandle->kernel[threadId]->M2L(&cellDst, inter, 1, level);

    return FMMAPI_NO_ERROR;
}

int FmmKernel_P2P(void *fmmCore, void *boxIdSrc, void *boxIdDest){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)fmmCore;
    ScalFmmKernelHandle* kernelhandle;
    FmmCore_getKernelData(corehandle, (void**)&kernelhandle);
    int threadId;
    FmmCore_getParameter(fmmCore, FMMCORE_THREAD_ID, &threadId);

    FReal* positionsTargets;
    FReal* potentialsTargets;
    int numberTargets;
    //FmmCore_getTargetPoints(fmmCore, boxIdDest, &positionsTargets, &numberTargets);
    FmmCore_getSource(fmmCore, boxIdDest, &positionsTargets, (void**)&potentialsTargets, &numberTargets);

    FReal* fieldsTargets = new FReal[numberTargets*kernelhandle->fieldDataSize];
    FmmCore_getTargetField(fmmCore, boxIdDest, fieldsTargets);

    int coordTargets[3];
    FmmCore_getCoord(fmmCore, boxIdDest, coordTargets);

    KernelContainerClass targets;
    KernelParticleClass part;
    for(int idxPart = 0 ; idxPart < numberTargets ; ++idxPart){
        part.setPosition(FPoint(&positionsTargets[idxPart*3]));
        part.setForces(fieldsTargets[idxPart*kernelhandle->fieldDataSize],
                       fieldsTargets[idxPart*kernelhandle->fieldDataSize+1],
                       fieldsTargets[idxPart*kernelhandle->fieldDataSize+2]);
        part.setPotential(fieldsTargets[idxPart*kernelhandle->fieldDataSize+3]);
        part.setPhysicalValue(potentialsTargets[idxPart]);
        targets.push(part);
    }

    FReal* positionsSrc;
    FReal* potentialsSrc;
    int numberSources;
    FmmCore_getSource(fmmCore, boxIdSrc, &positionsSrc, (void**)&potentialsSrc, &numberSources);

    KernelContainerClass sources;
    for(int idxPart = 0 ; idxPart < numberSources ; ++idxPart){
        part.setPosition(FPoint(&positionsSrc[idxPart*3]));
        part.setPhysicalValue(potentialsSrc[idxPart]);
        sources.push(part);
    }

    //KernelContainerClass* neigh[27];
    //memset(neigh, 0, sizeof(KernelContainerClass*)*27);
    //neigh[0] = &sources;
    //kernelhandle->kernel->P2PRemote(FTreeCoordinate(coordTargets[0],coordTargets[1],coordTargets[2]),
    //                          &targets, &sources, neigh, 1);
    for(int idxTarget = 0 ; idxTarget < numberTargets ; ++idxTarget){
        for(int idxSource = 0 ; idxSource < numberSources ; ++idxSource){
            FReal dx = sources[idxSource].getPosition().getX() - targets[idxTarget].getPosition().getX();
            FReal dy = sources[idxSource].getPosition().getY() - targets[idxTarget].getPosition().getY();
            FReal dz = sources[idxSource].getPosition().getZ() - targets[idxTarget].getPosition().getZ();
            const FReal distSquare = (dx*dx + dy*dy + dz*dz);

            if(distSquare > 10E-6*10E-6){
                kernelhandle->kernel[threadId]->particlesInteraction(&targets[idxTarget],sources[idxSource]);
            }
        }
    }

    for(int idxPart = 0 ; idxPart < numberTargets ; ++idxPart){
        fieldsTargets[idxPart*kernelhandle->fieldDataSize] = targets[idxPart].getForces().getX();
        fieldsTargets[idxPart*kernelhandle->fieldDataSize+1] = targets[idxPart].getForces().getY();
        fieldsTargets[idxPart*kernelhandle->fieldDataSize+2] = targets[idxPart].getForces().getZ();
        fieldsTargets[idxPart*kernelhandle->fieldDataSize+3] = targets[idxPart].getPotential();
    }

    FmmCore_releaseSource(fmmCore, boxIdDest, potentialsSrc, positionsSrc);
    FmmCore_releaseTargetPoints(fmmCore, boxIdDest, positionsTargets);
    FmmCore_setTargetField(fmmCore, boxIdDest, fieldsTargets);
    delete[] fieldsTargets;

    return FMMAPI_NO_ERROR;
} /* pas mutuel, i.e. on fait seulement dans 1 sens. */


