// ===================================================================================
// Copyright ScalFmm 2014 I
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


/**
 * @file This file contains a class that inherits from FScalFMMEngine,
 * and will implement the API functions for a user defined kernel.
 */
#ifndef FUSERKERNELENGINE_HPP
#define FUSERKERNELENGINE_HPP

#include "FScalFMMEngine.hpp"


/**
 * @brief CoreCell : Cell used to store User datas
 */
class CoreCell : public FBasicCell {
    // Mutable in order to work with the API
    mutable void* userData;

    //Static members to be initialised before octree creation
    static Scalfmm_Cell_Descriptor user_cell_descriptor;

public:
    static void Init(Scalfmm_Cell_Descriptor cell_descriptor){
        user_cell_descriptor=cell_descriptor;
    }

    static Callback_init_cell GetInit(){
        return user_cell_descriptor.user_init_cell;
    }

    static Callback_free_cell GetFree(){
        return user_cell_descriptor.user_free_cell;
    }

    CoreCell() : userData(nullptr) {
    }

    //We free the cells here
    ~CoreCell(){
        if(userData){
            this->user_cell_descriptor.user_free_cell(userData);
        }
    }

    /**
     * @brief setContainer store the ptr to the user data inside our
     * struct
     */
    void setContainer(void* inContainer) const {
        userData = inContainer;
    }

    /**
     * @brief getContainer : return the user datas (in order to give
     * it back to the user defined kernel function)
     */
    void* getContainer() const {
        return userData;
    }

};

/**
 * Define here static member
 */
Scalfmm_Cell_Descriptor CoreCell::user_cell_descriptor;

/**
 * This class simply call the function pointers from Scalfmm_Kernel_Descriptor.
 * If not pointer is set the calls are skipped.
 * The userData is given at any calls.
 */
template< class CellClass, class ContainerClass>
class CoreKernel : public FAbstractKernels<CellClass,ContainerClass> {
    Scalfmm_Kernel_Descriptor kernel;
    void* userData;

public:
    CoreKernel(Scalfmm_Kernel_Descriptor inKernel, void* inUserData) : kernel(inKernel) , userData(inUserData){
    }

    /** Default destructor */
    virtual ~CoreKernel(){
    }

    /** Do nothing */
    virtual void P2M(CellClass* const cell, const ContainerClass* const container) {
        if(kernel.p2m) kernel.p2m(cell->getContainer(), container->getNbParticles(), container->getIndexes().data(), userData);
    }

    /** Do nothing */
    virtual void M2M(CellClass* const FRestrict cell, const CellClass*const FRestrict *const FRestrict children, const int level) {
        if(kernel.m2m){
            for(int idx = 0 ; idx < 8 ; ++idx){
                if( children[idx] ){
                    kernel.m2m(level, cell->getContainer(), idx, children[idx]->getContainer(), userData);
                }
            }
        }
    }

    /** Do nothing */
    virtual void M2L(CellClass* const FRestrict cell, const CellClass* interactions[], const int , const int level) {
        if(kernel.m2l_full){//all 343 interactions will be computed directly
            //First, copy the fmm cell inside an array of user cells
            void * userCellArray[343];
            for(int i=0 ; i<343 ; ++i){
                if(interactions[i] != nullptr){
                    userCellArray[i] = interactions[i]->getContainer();
                }
                else{
                    userCellArray[i] = nullptr;
                }
            }
            kernel.m2l_full(level,cell->getContainer(),userCellArray,userData);
        }
        else{
            if(kernel.m2l){
                for(int idx = 0 ; idx < 343 ; ++idx){
                    if( interactions[idx] ){
                        kernel.m2l(level, cell->getContainer(), idx, interactions[idx]->getContainer(), userData);
                    }
                }
            }
        }
    }

    /** Do nothing */
    virtual void L2L(const CellClass* const FRestrict cell, CellClass* FRestrict *const FRestrict children, const int level) {
        if(kernel.l2l){
            for(int idx = 0 ; idx < 8 ; ++idx){
                if( children[idx] ){
                    kernel.l2l(level, cell->getContainer(), idx, children[idx]->getContainer(), userData);
                }
            }
        }
    }

    /** Do nothing */
    virtual void L2P(const CellClass* const cell, ContainerClass* const container){
        if(kernel.l2p) kernel.l2p(cell->getContainer(), container->getNbParticles(), container->getIndexes().data(), userData);
    }


    /** Do nothing */
    virtual void P2P(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict /*sources*/,
                     ContainerClass* const neighbors[27], const int ){
        if(kernel.p2pinner) kernel.p2pinner(targets->getNbParticles(), targets->getIndexes().data(), userData);

        if(kernel.p2p_full){
            //Create the arrays of size and indexes
            int nbPartPerNeighbors[27];
            const int * indicesPerNeighbors[27];
            for(int idx=0 ; idx<27 ; ++idx){
                if(neighbors[idx]){
                    nbPartPerNeighbors[idx] = neighbors[idx]->getNbParticles();
                    indicesPerNeighbors[idx] = neighbors[idx]->getIndexes().data();
                }
                else{
                    nbPartPerNeighbors[idx] = 0;
                    indicesPerNeighbors[idx] = nullptr;
                }
            }
            kernel.p2p_full(targets->getNbParticles(),targets->getIndexes().data(),indicesPerNeighbors,nbPartPerNeighbors,userData);
        }
        if(kernel.p2p){
            for(int idx = 0 ; idx < 27 ; ++idx){
                if( neighbors[idx] ){
                    kernel.p2p(targets->getNbParticles(), targets->getIndexes().data(),
                               neighbors[idx]->getNbParticles(), neighbors[idx]->getIndexes().data(), userData);
                }
            }
        }
    }

    /** Do nothing */
    virtual void P2PRemote(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     ContainerClass* const [27], const int ){
    }

};


class FUserKernelEngine : public FScalFMMEngine{

private:

    //Typedefs :
    typedef FP2PParticleContainerIndexed<>           ContainerClass;
    typedef FSimpleLeaf<ContainerClass>                   LeafClass;
    typedef FOctree<CoreCell,ContainerClass,LeafClass>  OctreeClass;
    typedef CoreKernel<CoreCell,ContainerClass>     CoreKernelClass;

    //For arranger classes
    typedef FBasicParticleContainerIndexedMover<OctreeClass, ContainerClass> MoverClass;
    typedef FOctreeArranger<OctreeClass, ContainerClass, MoverClass> ArrangerClass;
    typedef FArrangerPeriodic<OctreeClass, ContainerClass, MoverClass> ArrangerClassPeriodic;


    //Attributes
    OctreeClass * octree;
    CoreKernelClass * kernel;
    ArrangerClass * arranger;


public:
    FUserKernelEngine(/*int TreeHeight, double BoxWidth , double * BoxCenter, */scalfmm_kernel_type KernelType) :
        octree(nullptr), kernel(nullptr), arranger(nullptr){
        //        octree = new OctreeClass(TreeHeight,FMath::Min(3,TreeHeight-1),BoxWidth,FPoint(BoxCenter));
        kernelType = KernelType;
        //Kernel is not set now because the user must provide a
        //Scalfmm_Kernel_descriptor
    }


    ~FUserKernelEngine(){
        delete octree;
        octree=nullptr;
        if(arranger){
            delete arranger;
            arranger=nullptr;
        }
        if(kernel){
            delete kernel;
            kernel=nullptr;
        }
    }

    void user_kernel_config( Scalfmm_Kernel_Descriptor userKernel, void * userDatas){
        if(!kernel){
            kernel = new CoreKernelClass(userKernel,userDatas);
        }
    }

    void build_tree(int TreeHeight,double BoxWidth,double* BoxCenter,Scalfmm_Cell_Descriptor user_cell_descriptor){
        CoreCell::Init(user_cell_descriptor);
        this->octree = new OctreeClass(TreeHeight,FMath::Min(3,TreeHeight-1),BoxWidth,FPoint(BoxCenter));
    }

    void tree_insert_particles( int NbPositions, double * arrayX, double * arrayY, double * arrayZ){
        for(int idPart = 0; idPart<NbPositions ; ++idPart){
            octree->insert(FPoint(arrayX[idPart],arrayY[idPart],arrayZ[idPart]),idPart);
        }
        nbPart += NbPositions;
        this->init_cell();
    }

    void tree_insert_particles_xyz( int NbPositions, double * XYZ){
        for(int idPart = 0; idPart<NbPositions ; ++idPart){
            octree->insert(FPoint(&XYZ[3*idPart]),idPart);
        }
        nbPart += NbPositions;
        this->init_cell();
    }

    /**
     * To retrieve the positions, in order to move the parts
     */
    void get_positions_xyz(int NbPositions, double * positionsToFill){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    positionsToFill[indexes[idxPart]*3+0] = sources->getPositions()[0][idxPart];
                    positionsToFill[indexes[idxPart]*3+1] = sources->getPositions()[1][idxPart];
                    positionsToFill[indexes[idxPart]*3+2] = sources->getPositions()[2][idxPart];
                }
            });
    }

    void get_positions_xyz_npart(int NbPositions, int * idxOfParticles, double * positionsToFill){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < NbPositions && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            positionsToFill[indexes[idxPart]*3+0] =  sources->getPositions()[0][idxPart];
                            positionsToFill[indexes[idxPart]*3+1] =  sources->getPositions()[1][idxPart];
                            positionsToFill[indexes[idxPart]*3+2] =  sources->getPositions()[2][idxPart];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }

    void get_positions( int NbPositions, double * X, double * Y , double * Z){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    X[indexes[idxPart]] = sources->getPositions()[0][idxPart];
                    Y[indexes[idxPart]] = sources->getPositions()[1][idxPart];
                    Z[indexes[idxPart]] = sources->getPositions()[2][idxPart];
                }
            });
    }

    void get_positions_npart(int NbPositions, int * idxOfParticles,double * X, double * Y , double * Z){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < NbPositions && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            X[indexes[idxPart]] =  sources->getPositions()[0][idxPart];
                            Y[indexes[idxPart]] =  sources->getPositions()[1][idxPart];
                            Z[indexes[idxPart]] =  sources->getPositions()[2][idxPart];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }



    //Arranger parts : following function provide a way to move parts
    //inside the tree
    void add_to_positions_xyz(int NbPositions,double * updatedXYZ){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getWPositions()[0][idxPart] += updatedXYZ[indexes[idxPart]*3+0];
                    sources->getWPositions()[1][idxPart] += updatedXYZ[indexes[idxPart]*3+1];
                    sources->getWPositions()[2][idxPart] += updatedXYZ[indexes[idxPart]*3+2];
                }
            });
    }

    void add_to_positions(int NbPositions,double * X, double * Y , double * Z){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getWPositions()[0][idxPart] += X[indexes[idxPart]];
                    sources->getWPositions()[1][idxPart] += Y[indexes[idxPart]];
                    sources->getWPositions()[2][idxPart] += Z[indexes[idxPart]];
                }
            });
    }


    void set_positions_xyz(int NbPositions, double * updatedXYZ){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getWPositions()[0][idxPart] = updatedXYZ[indexes[idxPart]*3+0];
                    sources->getWPositions()[1][idxPart] = updatedXYZ[indexes[idxPart]*3+1];
                    sources->getWPositions()[2][idxPart] = updatedXYZ[indexes[idxPart]*3+2];
                }
            });
    }

    void set_positions(int NbPositions, double * X, double * Y, double * Z){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getWPositions()[0][idxPart] = X[indexes[idxPart]];
                    sources->getWPositions()[1][idxPart] = Y[indexes[idxPart]];
                    sources->getWPositions()[2][idxPart] = Z[indexes[idxPart]];
                }
            });
    }


    void update_tree(){
        if(arranger){
            arranger->rearrange();
            //then, we need to re-allocate cells user data for the
            //cells created during the process and free user datas for
            //the cells removed during the process
            init_cell();
        }
        else{
            if(Algorithm == 2){ //case in wich the periodic algorithm is used
                arranger = new ArrangerClassPeriodic(octree);
                arranger->rearrange();
            }
            else{
                arranger = new ArrangerClass(octree);
                arranger->rearrange();
                init_cell();
            }
        }
    }

    /*
     * Call the user allocator on userDatas member field of each cell
     */
    void init_cell(){

        double boxwidth = octree->getBoxWidth();
        FPoint BoxCenter = octree->getBoxCenter();
        double boxCorner[3];
        boxCorner[0] = BoxCenter.getX() - boxwidth/2.0;
        boxCorner[1] = BoxCenter.getY() - boxwidth/2.0;
        boxCorner[2] = BoxCenter.getZ() - boxwidth/2.0;
        //apply user function on each cell
        octree->forEachCellWithLevel([&](CoreCell * currCell,const int currLevel){
                if(!(currCell->getContainer())){
                    FTreeCoordinate currCoord = currCell->getCoordinate();
                    int arrayCoord[3] = {currCoord.getX(),currCoord.getY(),currCoord.getZ()};
                    MortonIndex    currMorton = currCoord.getMortonIndex(currLevel);
                    double position[3];
                    position[0] = boxCorner[0] + currCoord.getX()*boxwidth/double(1<<currLevel);
                    position[1] = boxCorner[1] + currCoord.getY()*boxwidth/double(1<<currLevel);
                    position[2] = boxCorner[2] + currCoord.getZ()*boxwidth/double(1<<currLevel);
                    currCell->setContainer(CoreCell::GetInit()(currLevel,currMorton,arrayCoord,position));
                }
            });
    }


    void free_cell(Callback_free_cell user_cell_deallocator){
        octree->forEachCell([&](CoreCell * currCell){
                if(currCell->getContainer()){
                    user_cell_deallocator(currCell->getContainer());
                    currCell->setContainer(nullptr);
                }
            });
    }

    void execute_fmm(){
        FAssertLF(kernel,"No kernel set, please use scalfmm_user_kernel_config before calling the execute routine ... Exiting \n");
        switch(Algorithm){
        case 0:
            {
                typedef FFmmAlgorithm<OctreeClass,CoreCell,ContainerClass,CoreKernelClass,LeafClass> AlgoClassSeq;
                AlgoClassSeq algoSeq(octree,kernel);
                algoSeq.execute();
                break;
            }
        case 1:
            {
                typedef FFmmAlgorithmThread<OctreeClass,CoreCell,ContainerClass,CoreKernelClass,LeafClass> AlgoClassThread;
                AlgoClassThread algoThread(octree,kernel);
                algoThread.execute(FFmmNearField);
                break;
            }
        case 2:
            {
                typedef FFmmAlgorithmPeriodic<OctreeClass,CoreCell,ContainerClass,CoreKernelClass,LeafClass> AlgoClassPeriodic;
                AlgoClassPeriodic algoPeriod(octree,2);
                algoPeriod.setKernel(kernel);
                algoPeriod.execute();
                break;
            }
        default :
            std::cout<< "No algorithm found (probably for strange reasons) : "<< Algorithm <<" exiting" << std::endl;
        }

    }

    void intern_dealloc_handle(Callback_free_cell userDeallocator){
        free_cell(userDeallocator);
    }
};


#endif //FUSERKERNELENGINE_HPP
