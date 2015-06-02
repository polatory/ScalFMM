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
class CoreCell : public FBasicCell, public FExtendCellType {
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
            FSize nbPartPerNeighbors[27];
            const FSize * indicesPerNeighbors[27];
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

    //Getter
    void * getUserKernelDatas(){
        return userData;
    }
    //Getter
    Scalfmm_Kernel_Descriptor getKernelFct() const {
        return kernel;
    }

    void M2L_Extended(CellClass * src, CellClass * tgt, const FTreeCoordinate transfer, const int level){
        if(kernel.m2l_ext){
            int array[3] = {transfer.getX(),transfer.getY(),transfer.getZ()};
            kernel.m2l_ext(level,tgt->getContainer(),src->getContainer(),array,userData);
        }
    }

};

template<class FReal,class LeafClass>
class FUserKernelEngine : public FScalFMMEngine<FReal>{

private:
    //Typedefs
    typedef FP2PParticleContainerIndexed<FReal>           ContainerClass;

    //Typedefs :
    typedef FOctree<FReal,CoreCell,ContainerClass,LeafClass>            OctreeClass;
    typedef CoreKernel<CoreCell,ContainerClass>           CoreKernelClass;

    //For arranger classes

    //Attributes
    OctreeClass * octree;
    CoreKernelClass * kernel;
    int upperLimit;
    int treeHeight;

    // ArrangerClass * arranger;
    // ArrangerClassTyped * arrangerTyped;
    // ArrangerClassPeriodic * arrangerPeriodic;


public:
    FUserKernelEngine(/*int TreeHeight, double BoxWidth , double * BoxCenter, */scalfmm_kernel_type KernelType) :
        octree(nullptr), kernel(nullptr), upperLimit(2), treeHeight(0) /*,arranger(nullptr)*/ {
        FScalFMMEngine<FReal>::kernelType = KernelType;
    }


    ~FUserKernelEngine(){
        delete octree;
        octree=nullptr;
        // if(arranger){
        //     delete arranger;
        //     arranger=nullptr;
        // }
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
        this->treeHeight = TreeHeight;
        printf("Tree Height : %d \n",TreeHeight);
        this->octree = new OctreeClass(TreeHeight,FMath::Min(3,TreeHeight-1),BoxWidth,FPoint<FReal>(BoxCenter));
    }

    void reset_tree(Callback_reset_cell cellReset){
        double boxwidth = octree->getBoxWidth();
        FPoint<FReal> BoxCenter = octree->getBoxCenter();
        double boxCorner[3];
        boxCorner[0] = BoxCenter.getX() - boxwidth/2.0;
        boxCorner[1] = BoxCenter.getY() - boxwidth/2.0;
        boxCorner[2] = BoxCenter.getZ() - boxwidth/2.0;
        //apply user function reset on each user's cell
        octree->forEachCellWithLevel([&](CoreCell * currCell,const int currLevel){
                if(currCell->getContainer()){
                    FTreeCoordinate currCoord = currCell->getCoordinate();
                    int arrayCoord[3] = {currCoord.getX(),currCoord.getY(),currCoord.getZ()};
                    MortonIndex    currMorton = currCoord.getMortonIndex(currLevel);
                    double position[3];
                    position[0] = boxCorner[0] + currCoord.getX()*boxwidth/double(1<<currLevel);
                    position[1] = boxCorner[1] + currCoord.getY()*boxwidth/double(1<<currLevel);
                    position[2] = boxCorner[2] + currCoord.getZ()*boxwidth/double(1<<currLevel);
                    cellReset(currLevel,currMorton,arrayCoord,position,currCell->getContainer());
                }
            });
    }


    void tree_insert_particles( int NbPositions, double * X, double * Y, double * Z, PartType type){
        if(type == BOTH){
            for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                octree->insert(FPoint<FReal>(X[idPart],Y[idPart],Z[idPart]),idPart);
            }
            FScalFMMEngine<FReal>::nbPart += NbPositions;
        }else{
            if(type==SOURCE){
                for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                    octree->insert(FPoint<FReal>(X[idPart],Y[idPart],Z[idPart]),FParticleTypeSource,idPart);
                }
                FScalFMMEngine<FReal>::nbPart += NbPositions;
            }else{
                for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                    octree->insert(FPoint<FReal>(X[idPart],Y[idPart],Z[idPart]),FParticleTypeTarget,idPart);
                }
                FScalFMMEngine<FReal>::nbPart += NbPositions;
            }
        }

        this->init_cell();
    }

    void tree_insert_particles_xyz( int NbPositions, double * XYZ, PartType type){
        if(type == BOTH){
            for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                octree->insert(FPoint<FReal>(&XYZ[3*idPart]),idPart);
            }
            FScalFMMEngine<FReal>::nbPart += NbPositions;
        }else{
            if(type==SOURCE){
                for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                    octree->insert(FPoint<FReal>(&XYZ[3*idPart]),FParticleTypeSource,idPart);
                }
                FScalFMMEngine<FReal>::nbPart += NbPositions;
            }else{
                for(FSize idPart = 0; idPart<NbPositions ; ++idPart){
                    octree->insert(FPoint<FReal>(&XYZ[3*idPart]),FParticleTypeTarget,idPart);
                }
                FScalFMMEngine<FReal>::nbPart += NbPositions;
            }
        }
        this->init_cell();
    }

    /**
     * To retrieve the positions, in order to move the parts
     */
    void get_positions_xyz(int NbPositions, double * positionsToFill, PartType type){
        FScalFMMEngine<FReal>::template generic_get_positions_xyz<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,positionsToFill,type);
    }

    void get_positions_xyz_npart(int NbPositions, int * idxOfParticles, double * positionsToFill,PartType type){
        FScalFMMEngine<FReal>::template generic_get_positions_xyz_npart<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,idxOfParticles,positionsToFill,type);
    }

    void get_positions( int NbPositions, double *X, double *Y , double *Z, PartType type){
        FScalFMMEngine<FReal>::template generic_get_positions<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,X,Y,Z,type);
    }

    void get_positions_npart(int NbPositions, int * idxOfParticles,double * X, double * Y , double * Z,PartType type){
        FScalFMMEngine<FReal>::template generic_get_positions_npart<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,idxOfParticles,X,Y,Z,type);
    }



    //Arranger parts : following function provide a way to move parts
    //inside the tree
    void add_to_positions_xyz(int NbPositions,double * updatedXYZ,PartType type){
        FScalFMMEngine<FReal>::template generic_add_to_positions_xyz<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,updatedXYZ,type);
    }

    void add_to_positions(int NbPositions,double * X, double * Y , double * Z,PartType type){
        FScalFMMEngine<FReal>::template generic_add_to_positions<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,X,Y,Z,type);
    }

    void set_positions_xyz(int NbPositions, FReal * updatedXYZ, PartType type){
        FScalFMMEngine<FReal>::template generic_set_positions_xyz<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,updatedXYZ,type);
    }

    void set_positions(int NbPositions, FReal * X,FReal * Y,FReal * Z, PartType type){
        FScalFMMEngine<FReal>::template generic_set_positions<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,X,Y,Z,type);
    }

    void set_positions_xyz_npart(int NbPositions, int* idxOfParticles, FReal * updatedXYZ, PartType type){
        FScalFMMEngine<FReal>::template generic_set_positions_xyz_npart<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,idxOfParticles,updatedXYZ,type);
    }
    void set_positions_npart(int NbPositions, int* idxOfParticles, FReal * X, FReal * Y , FReal * Z, PartType type){
        FScalFMMEngine<FReal>::template generic_set_positions_npart<ContainerClass,LeafClass,CoreCell>(octree,NbPositions,idxOfParticles,X,Y,Z,type);
    }


    // void update_tree(){
    //     if(arranger){
    //         arranger->rearrange();
    //         //then, we need to re-allocate cells user data for the
    //         //cells created during the process and free user datas for
    //         //the cells removed during the process
    //         init_cell();
    //     }
    //     else{
    //         if(FScalFMMEngine<FReal>::Algorithm == 2){ //case in wich the periodic algorithm is used
    //             arranger = new ArrangerClassPeriodic(octree);
    //             arranger->rearrange();
    //             init_cell();
    //         }
    //         else{
    //             arranger = new ArrangerClass(octree);
    //             arranger->rearrange();
    //             init_cell();
    //         }
    //     }
    // }

    /*
     * Call the user allocator on userDatas member field of each cell
     */
    void init_cell(){
        void * generic_ptr = nullptr;
        if(kernel){
            generic_ptr = kernel->getUserKernelDatas();
        }
        else{
            std::cout <<"Warning, no user kernel data set, need to call kernel config first"<< std::endl;
        }
        double boxwidth = octree->getBoxWidth();
        FPoint<FReal> BoxCenter = octree->getBoxCenter();
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
                    currCell->setContainer(CoreCell::GetInit()(currLevel,currMorton,arrayCoord,position,generic_ptr));
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

    void set_upper_limit(int inUpperLimit){
        upperLimit = inUpperLimit;
    }

    /**
     * @brief This function is called if the FMM is not computed on
     * all the standards levels
     *
     */
    void internal_M2L(){
        if(this->kernel->getKernelFct().m2l_ext){
            if(upperLimit > 1){ // if upperLimit == 1, then, M2L has been
                // done at level 2, and hence all the far
                // field has been calculated.
                //Starting at the lower level where the M2L has not been done.
                typename OctreeClass::Iterator octreeIterator(octree); //lvl : 1

                while(octreeIterator.level() != upperLimit){
                    octreeIterator.moveDown();
                }

                //I'm at the upperLimit, so the lowest level where M2L has been done.
                do{
                    CoreCell * currentTgt = octreeIterator.getCurrentCell(); // This one is targeted

                    //Then, we get the interaction list at this lvl. This will provide us with lots of source cells.
                    const CoreCell * currentInteractionList[343];
                    //Get an iterator for the sources
                    typename OctreeClass::Iterator upAndDownIterator = octreeIterator;

                    {//This is supposed to be done for multiple level. You
                        //need to go up until level 2. And then, to go down
                        //until level upperLimit. I think it's possible ...
                        while(upAndDownIterator.level() >= 2){
                            upAndDownIterator.moveUp();

                            //There, we get the interaction list of all parents of tgt cell
                            const int nbInteract = octree->getInteractionNeighbors(currentInteractionList,
                                                                                   upAndDownIterator.getCurrentGlobalCoordinate(),
                                                                                   upAndDownIterator.level());
                            int currentLevel = upAndDownIterator.level();
                            if(nbInteract){
                                //Then, we do M2L for each child at level upperLimit of each 343 Interaction cells.
                                for(int idxSrc = 0; idxSrc < 343 ; ++idxSrc){
                                    if(currentInteractionList[idxSrc]){//Check if it exist
                                        const CoreCell * currentSource = currentInteractionList[idxSrc]; //For clarity, will be otpimised out, anyway
                                        MortonIndex idx = currentSource->getMortonIndex();

                                        //At this point, we instanciate
                                        //the number of child needed.
                                        //This only depends on diffenrence
                                        //between current level and
                                        //upperLimit level
                                        int totalNumberOfChild = FMath::pow(8,upperLimit-currentLevel);

                                        for(int idxChildSrc = 0; idxChildSrc < totalNumberOfChild ; ++idxChildSrc){//For all 8^{number of levels to down} children
                                            MortonIndex indexOfChild = ((idx << 3*(upperLimit-currentLevel))+idxChildSrc);
                                            CoreCell * src = octree->getCell(indexOfChild,upperLimit); //Get the cell
                                            if(src){//check if it exists
                                                FTreeCoordinate srcCoord = src->getCoordinate();
                                                FTreeCoordinate tgtCoord = currentTgt->getCoordinate();
                                                //Build tree coord translation vector
                                                FTreeCoordinate transfer;
                                                transfer.setPosition(tgtCoord.getX()-srcCoord.getX(),
                                                                     tgtCoord.getY()-srcCoord.getY(),
                                                                     tgtCoord.getZ()-srcCoord.getZ());
                                                kernel->M2L_Extended(src,currentTgt,transfer,octreeIterator.level());
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }while(octreeIterator.moveRight());
            }
            else{
                FAssertLF("No reasons to be there, seriously ...\nExiting anyway...");
            }
        }

    }

    void execute_fmm(){
        FAssertLF(kernel,"No kernel set, please use scalfmm_user_kernel_config before calling the execute routine ... Exiting \n");
        FAbstractAlgorithm * abstrct = nullptr;
        switch(FScalFMMEngine<FReal>::Algorithm){
        case 0:
            {
                typedef FFmmAlgorithm<OctreeClass,CoreCell,ContainerClass,CoreKernelClass,LeafClass> AlgoClassSeq;
                AlgoClassSeq * algoSeq = new AlgoClassSeq(octree,kernel);
                FScalFMMEngine<FReal>::algoTimer = algoSeq;
                abstrct = algoSeq;
                //algoSeq->execute(); will be done later
                break;
            }
        case 1:
            {
                typedef FFmmAlgorithmThread<OctreeClass,CoreCell,ContainerClass,CoreKernelClass,LeafClass> AlgoClassThread;
                AlgoClassThread*  algoThread = new AlgoClassThread(octree,kernel);
                FScalFMMEngine<FReal>::algoTimer = algoThread;
                abstrct = algoThread;
                //algoThread->execute(); will be done later
                break;
            }
        case 2:
            {
                typedef FFmmAlgorithmPeriodic<FReal,OctreeClass,CoreCell,ContainerClass,CoreKernelClass,LeafClass> AlgoClassPeriodic;
                AlgoClassPeriodic algoPeriod(octree,2);
                algoPeriod.setKernel(kernel);
                algoPeriod.execute();
                break;
            }
        case 3:
            {
                typedef FFmmAlgorithmThreadTsm<OctreeClass,CoreCell,ContainerClass,CoreKernelClass,LeafClass> AlgoClassTargetSource;
                AlgoClassTargetSource* algoTS = new AlgoClassTargetSource(octree,kernel);
                FScalFMMEngine<FReal>::algoTimer = algoTS;
                abstrct = algoTS;
                //algoTS->execute(); will be done later
                break;
            }
        default :
            std::cout<< "No algorithm found (probably for strange reasons) : "<< FScalFMMEngine<FReal>::Algorithm <<" exiting" << std::endl;
        }
        if (FScalFMMEngine<FReal>::Algorithm != 2){
            if(upperLimit != 2){
                abstrct->execute(FFmmP2M | FFmmM2M | FFmmM2L, upperLimit, treeHeight);
                printf("\tUpPass finished\n");
                internal_M2L();
                printf("\tStrange M2L finished\n");
                abstrct->execute(FFmmL2L | FFmmL2P | FFmmP2P, upperLimit, treeHeight);
                printf("\tDownPass finished\n");
            }
            else{
                abstrct->execute();
            }
        }
    }

    void intern_dealloc_handle(Callback_free_cell userDeallocator){
        free_cell(userDeallocator);
    }
};


#endif //FUSERKERNELENGINE_HPP
