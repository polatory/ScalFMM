#ifndef FFMMALGORITHMOMP4_HPP
#define FFMMALGORITHMOMP4_HPP

// @SCALFMM_PRIVATE

#include "../Utils/FGlobal.hpp"
#include "../Utils/FAssert.hpp"
#include "../Utils/FLog.hpp"

#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"
#include "../Utils/FAlgorithmTimers.hpp"

#include "FCoreCommon.hpp"
#include "FP2PExclusion.hpp"

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FFmmAlgorithmOmp4
 * @brief
 * Please read the license
 *
 * This class is a basic FMM algorithm
 * It just iterates on a tree and call the kernels with good arguments.
 *
 * Of course this class does not deallocate pointer given in arguements.
 */
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass, class P2PExclusionClass = FP2PMiddleExclusion>
class FFmmAlgorithmOmp4 : public FAbstractAlgorithm, public FAlgorithmTimers {

    OctreeClass* const tree;       //< The octree to work on
    KernelClass** kernels;    //< The kernels

    const int MaxThreads;

    const int OctreeHeight;

    const int leafLevelSeparationCriteria;
public:
    /** The constructor need the octree and the kernels used for computation
     * @param inTree the octree to work on
     * @param inKernels the kernels to call
     * An assert is launched if one of the arguments is null
     */
    FFmmAlgorithmOmp4(OctreeClass* const inTree, KernelClass* const inKernels, const int inLeafLevelSeperationCriteria = 1)
: tree(inTree) , kernels(nullptr),
  MaxThreads(omp_get_max_threads()), OctreeHeight(tree->getHeight()), leafLevelSeparationCriteria(inLeafLevelSeperationCriteria)
{

        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");
        FAssertLF(leafLevelSeparationCriteria < 3, "Separation criteria should be < 3");

        this->kernels = new KernelClass*[MaxThreads];
        #pragma omp parallel num_threads(MaxThreads)
        {
            #pragma omp critical (InitFFmmAlgorithmOmp4)
            {
                this->kernels[omp_get_thread_num()] = new KernelClass(*inKernels);
            }
        }

        FAbstractAlgorithm::setNbLevelsInTree(tree->getHeight());

        FLOG(FLog::Controller << "FFmmAlgorithmOmp4 (Max Thread " << omp_get_max_threads() << ")\n");
}

    /** Default destructor */
    virtual ~FFmmAlgorithmOmp4(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;
    }

protected:
    /**
     * To execute the fmm algorithm
     * Call this function to run the complete algorithm
     */
    void executeCore(const unsigned operationsToProceed) override {

        #pragma omp parallel
        {
            #pragma omp master
            {
                Timers[P2MTimer].tic();
                if(operationsToProceed & FFmmP2M)
                    bottomPass();
                Timers[P2MTimer].tac();

                Timers[M2MTimer].tic();
                if(operationsToProceed & FFmmM2M)
                    upwardPass();
                Timers[M2MTimer].tac();

                Timers[M2LTimer].tic();
                if(operationsToProceed & FFmmM2L)
                    transferPass();
                Timers[M2LTimer].tac();

                Timers[L2LTimer].tic();
                if(operationsToProceed & FFmmL2L)
                    downardPass();
                Timers[L2LTimer].tac();

                Timers[NearTimer].tic();
                if( (operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P) )
                    directPass((operationsToProceed & FFmmP2P),(operationsToProceed & FFmmL2P));
                Timers[NearTimer].tac();

                #pragma omp taskwait
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void bottomPass(){
        FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
        FLOG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            // We need the current cell that represent the leaf
            // and the list of particles
            CellClass* taskCell = octreeIterator.getCurrentCell();
            ContainerClass* taskParticles = octreeIterator.getCurrentListSrc();
            #pragma omp task firstprivate(taskCell, taskParticles) depend(inout:taskCell[0]) depend(in:taskParticles[0])
            {
                kernels[omp_get_thread_num()]->P2M( taskCell , taskParticles);
            }
        } while(octreeIterator.moveRight());


        FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        for(int idxLevel = OctreeHeight - 2 ; idxLevel > FAbstractAlgorithm::lowerWorkingLevel-1 ; --idxLevel){
            octreeIterator.moveUp();
        }

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = FMath::Min(OctreeHeight - 2, FAbstractAlgorithm::lowerWorkingLevel - 1) ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel ){
            FLOG(FTic counterTimeLevel);
            // for each cells
            do{
                // We need the current cell and the child
                // child is an array (of 8 child) that may be null
                CellClass* taskCell = octreeIterator.getCurrentCell();
                CellClass* taskChild[8];
                memcpy(taskChild, octreeIterator.getCurrentChild(), 8*sizeof(CellClass*));
#pragma omp task firstprivate(taskCell, taskChild, idxLevel) depend(inout:taskCell[0]) depend(in:taskChild[0],taskChild[1],taskChild[2],taskChild[3],taskChild[4],taskChild[5],taskChild[6],taskChild[7])
                {
                    kernels[omp_get_thread_num()]->M2M( taskCell , taskChild, idxLevel);
                }
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }


        FLOG( FLog::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Transfer
    /////////////////////////////////////////////////////////////////////////////

    /** M2L  */
    void transferPass(){
      #ifdef SCALFMM_USE_EZTRACE

      eztrace_start();
#endif
        this->transferPassWithOutFinalize() ;
#ifdef SCALFMM_USE_EZTRACE
      eztrace_stop();
#endif
        }

    void transferPassWithOutFinalize(){
        FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);
        // Goto the right level
        octreeIterator.moveDown();
        for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; ++idxLevel){
            octreeIterator.moveDown();
        }
        ////////////////////////////////////////////////////////////////
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
        //
        // for each levels
        for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < FAbstractAlgorithm::lowerWorkingLevel ; ++idxLevel ){
            FLOG(FTic counterTimeLevel);
            const int separationCriteria = (idxLevel != FAbstractAlgorithm::lowerWorkingLevel-1 ? 1 : leafLevelSeparationCriteria);
            // for each cell we apply the M2L with all cells in the implicit interaction list
            do{
                const CellClass* taskNeigh[343];
                const int counter = tree->getInteractionNeighbors(taskNeigh, octreeIterator.getCurrentGlobalCoordinate(), idxLevel, separationCriteria);
                if(counter){
                    CellClass* taskCell = octreeIterator.getCurrentCell();

    #pragma omp task firstprivate(taskCell,taskNeigh, idxLevel, counter) depend(inout:taskCell[0]) depend(in:taskNeigh[0],taskNeigh[1],taskNeigh[2],taskNeigh[3],taskNeigh[4],taskNeigh[5],taskNeigh[6],taskNeigh[7],taskNeigh[8],taskNeigh[9],taskNeigh[10],taskNeigh[11],taskNeigh[12],taskNeigh[13],taskNeigh[14],taskNeigh[15],taskNeigh[16],taskNeigh[17],taskNeigh[18],taskNeigh[19],taskNeigh[20],taskNeigh[21],taskNeigh[22],taskNeigh[23],taskNeigh[24],taskNeigh[25],taskNeigh[26],taskNeigh[27],taskNeigh[28],taskNeigh[29],taskNeigh[30],taskNeigh[31],taskNeigh[32],taskNeigh[33],taskNeigh[34],taskNeigh[35],taskNeigh[36],taskNeigh[37],taskNeigh[38],taskNeigh[39],taskNeigh[40],taskNeigh[41],taskNeigh[42],taskNeigh[43],taskNeigh[44],taskNeigh[45],taskNeigh[46],taskNeigh[47],taskNeigh[48],taskNeigh[49],taskNeigh[50],taskNeigh[51],taskNeigh[52],taskNeigh[53],taskNeigh[54],taskNeigh[55],taskNeigh[56],taskNeigh[57],taskNeigh[58],taskNeigh[59],taskNeigh[60],taskNeigh[61],taskNeigh[62],taskNeigh[63],taskNeigh[64],taskNeigh[65],taskNeigh[66],taskNeigh[67],taskNeigh[68],taskNeigh[69],taskNeigh[70],taskNeigh[71],taskNeigh[72],taskNeigh[73],taskNeigh[74],taskNeigh[75],taskNeigh[76],taskNeigh[77],taskNeigh[78],taskNeigh[79],taskNeigh[80],taskNeigh[81],taskNeigh[82],taskNeigh[83],taskNeigh[84],taskNeigh[85],taskNeigh[86],taskNeigh[87],taskNeigh[88],taskNeigh[89],taskNeigh[90],taskNeigh[91],taskNeigh[92],taskNeigh[93],taskNeigh[94],taskNeigh[95],taskNeigh[96],taskNeigh[97],taskNeigh[98],taskNeigh[99],taskNeigh[100],taskNeigh[101],taskNeigh[102],taskNeigh[103],taskNeigh[104],taskNeigh[105],taskNeigh[106],taskNeigh[107],taskNeigh[108],taskNeigh[109],taskNeigh[110],taskNeigh[111],taskNeigh[112],taskNeigh[113],taskNeigh[114],taskNeigh[115],taskNeigh[116],taskNeigh[117],taskNeigh[118],taskNeigh[119],taskNeigh[120],taskNeigh[121],taskNeigh[122],taskNeigh[123],taskNeigh[124],taskNeigh[125],taskNeigh[126],taskNeigh[127],taskNeigh[128],taskNeigh[129],taskNeigh[130],taskNeigh[131],taskNeigh[132],taskNeigh[133],taskNeigh[134],taskNeigh[135],taskNeigh[136],taskNeigh[137],taskNeigh[138],taskNeigh[139],taskNeigh[140],taskNeigh[141],taskNeigh[142],taskNeigh[143],taskNeigh[144],taskNeigh[145],taskNeigh[146],taskNeigh[147],taskNeigh[148],taskNeigh[149],taskNeigh[150],taskNeigh[151],taskNeigh[152],taskNeigh[153],taskNeigh[154],taskNeigh[155],taskNeigh[156],taskNeigh[157],taskNeigh[158],taskNeigh[159],taskNeigh[160],taskNeigh[161],taskNeigh[162],taskNeigh[163],taskNeigh[164],taskNeigh[165],taskNeigh[166],taskNeigh[167],taskNeigh[168],taskNeigh[169],taskNeigh[170],taskNeigh[171],taskNeigh[172],taskNeigh[173],taskNeigh[174],taskNeigh[175],taskNeigh[176],taskNeigh[177],taskNeigh[178],taskNeigh[179],taskNeigh[180],taskNeigh[181],taskNeigh[182],taskNeigh[183],taskNeigh[184],taskNeigh[185],taskNeigh[186],taskNeigh[187],taskNeigh[188],taskNeigh[189],taskNeigh[190],taskNeigh[191],taskNeigh[192],taskNeigh[193],taskNeigh[194],taskNeigh[195],taskNeigh[196],taskNeigh[197],taskNeigh[198],taskNeigh[199],taskNeigh[200],taskNeigh[201],taskNeigh[202],taskNeigh[203],taskNeigh[204],taskNeigh[205],taskNeigh[206],taskNeigh[207],taskNeigh[208],taskNeigh[209],taskNeigh[210],taskNeigh[211],taskNeigh[212],taskNeigh[213],taskNeigh[214],taskNeigh[215],taskNeigh[216],taskNeigh[217],taskNeigh[218],taskNeigh[219],taskNeigh[220],taskNeigh[221],taskNeigh[222],taskNeigh[223],taskNeigh[224],taskNeigh[225],taskNeigh[226],taskNeigh[227],taskNeigh[228],taskNeigh[229],taskNeigh[230],taskNeigh[231],taskNeigh[232],taskNeigh[233],taskNeigh[234],taskNeigh[235],taskNeigh[236],taskNeigh[237],taskNeigh[238],taskNeigh[239],taskNeigh[240],taskNeigh[241],taskNeigh[242],taskNeigh[243],taskNeigh[244],taskNeigh[245],taskNeigh[246],taskNeigh[247],taskNeigh[248],taskNeigh[249],taskNeigh[250],taskNeigh[251],taskNeigh[252],taskNeigh[253],taskNeigh[254],taskNeigh[255],taskNeigh[256],taskNeigh[257],taskNeigh[258],taskNeigh[259],taskNeigh[260],taskNeigh[261],taskNeigh[262],taskNeigh[263],taskNeigh[264],taskNeigh[265],taskNeigh[266],taskNeigh[267],taskNeigh[268],taskNeigh[269],taskNeigh[270],taskNeigh[271],taskNeigh[272],taskNeigh[273],taskNeigh[274],taskNeigh[275],taskNeigh[276],taskNeigh[277],taskNeigh[278],taskNeigh[279],taskNeigh[280],taskNeigh[281],taskNeigh[282],taskNeigh[283],taskNeigh[284],taskNeigh[285],taskNeigh[286],taskNeigh[287],taskNeigh[288],taskNeigh[289],taskNeigh[290],taskNeigh[291],taskNeigh[292],taskNeigh[293],taskNeigh[294],taskNeigh[295],taskNeigh[296],taskNeigh[297],taskNeigh[298],taskNeigh[299],taskNeigh[300],taskNeigh[301],taskNeigh[302],taskNeigh[303],taskNeigh[304],taskNeigh[305],taskNeigh[306],taskNeigh[307],taskNeigh[308],taskNeigh[309],taskNeigh[310],taskNeigh[311],taskNeigh[312],taskNeigh[313],taskNeigh[314],taskNeigh[315],taskNeigh[316],taskNeigh[317],taskNeigh[318],taskNeigh[319],taskNeigh[320],taskNeigh[321],taskNeigh[322],taskNeigh[323],taskNeigh[324],taskNeigh[325],taskNeigh[326],taskNeigh[327],taskNeigh[328],taskNeigh[329],taskNeigh[330],taskNeigh[331],taskNeigh[332],taskNeigh[333],taskNeigh[334],taskNeigh[335],taskNeigh[336],taskNeigh[337],taskNeigh[338],taskNeigh[339],taskNeigh[340],taskNeigh[341],taskNeigh[342])
                    {
                        if(counter){
                            kernels[omp_get_thread_num()]->M2L(  taskCell, taskNeigh, counter, idxLevel);
                        }
                    }
                }
            } while(octreeIterator.moveRight());
            ////////////////////////////////////////////////////////////////
            // move up  and goto left
            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;

            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }

        FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );

    }

    /////////////////////////////////////////////////////////////////////////////
    // Downward
    /////////////////////////////////////////////////////////////////////////////

    void downardPass(){ // second L2L
        FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();

        for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; ++idxLevel){
            octreeIterator.moveDown();
        }

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = FAbstractAlgorithm::lowerWorkingLevel - 1;
        // for each levels exepted leaf level
        for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < heightMinusOne ; ++idxLevel ){
            FLOG(FTic counterTimeLevel);
            // for each cells
            do{
                CellClass* taskCell = octreeIterator.getCurrentCell();
                CellClass* taskChild[8];
                memcpy(taskChild, octreeIterator.getCurrentChild(), 8*sizeof(CellClass*));

#pragma omp task firstprivate(taskCell, taskChild, idxLevel) depend(in:taskCell[0]) depend(inout:taskChild[0],taskChild[1],taskChild[2],taskChild[3],taskChild[4],taskChild[5],taskChild[6],taskChild[7])
                {
                    kernels[omp_get_thread_num()]->L2L( taskCell , taskChild, idxLevel);
                }

            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;

            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }

        FLOG( FLog::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }


    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////

    /** P2P */
    void directPass(const bool p2pEnabled, const bool l2pEnabled){
        FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        const int heightMinusOne = OctreeHeight - 1;

        FLOG( computationCounter.tic() );

        const int SizeShape = P2PExclusionClass::SizeShape;
        FVector<typename OctreeClass::Iterator> shapes[SizeShape];

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();

        // for each leafs
        do{
            if(l2pEnabled){
                CellClass* taskCell = octreeIterator.getCurrentCell();
                ContainerClass* taskParticles = octreeIterator.getCurrentListTargets();
                #pragma omp task firstprivate(taskCell, taskParticles) depend(in:taskCell[0]) depend(inout:taskParticles[0])
                {
                    kernels[omp_get_thread_num()]->L2P(taskCell, taskParticles);
                }
            }
            if(p2pEnabled){
                // There is a maximum of 26 neighbors
                ContainerClass* neighbors[27];
                const int counter = tree->getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),heightMinusOne);

                ContainerClass* taskParticlesTgt = octreeIterator.getCurrentListTargets();
                ContainerClass* taskParticlesSrc = octreeIterator.getCurrentListSrc();
                const FTreeCoordinate coord = octreeIterator.getCurrentGlobalCoordinate();

                if(taskParticlesTgt == taskParticlesSrc){
                    #pragma omp task firstprivate(neighbors, counter, taskParticlesTgt, coord) depend(inout:taskParticlesTgt[0],neighbors[0],neighbors[1],neighbors[2],neighbors[3],neighbors[4],neighbors[5],neighbors[6],neighbors[7],neighbors[8],neighbors[9],neighbors[10],neighbors[11],neighbors[12],neighbors[13],neighbors[14],neighbors[15],neighbors[16],neighbors[17],neighbors[18],neighbors[19],neighbors[20],neighbors[21],neighbors[22],neighbors[23],neighbors[24],neighbors[25],neighbors[26])
                    {
                        kernels[omp_get_thread_num()]->P2P(coord, taskParticlesTgt,
                                taskParticlesTgt, neighbors, counter);
                    }
                }
                else{
                    #pragma omp task firstprivate(neighbors, counter, taskParticlesTgt, taskParticlesSrc, coord) depend(inout:taskParticlesTgt[0]) depend(in:taskParticlesSrc[0],neighbors[0],neighbors[1],neighbors[2],neighbors[3],neighbors[4],neighbors[5],neighbors[6],neighbors[7],neighbors[8],neighbors[9],neighbors[10],neighbors[11],neighbors[12],neighbors[13],neighbors[14],neighbors[15],neighbors[16],neighbors[17],neighbors[18],neighbors[19],neighbors[20],neighbors[21],neighbors[22],neighbors[23],neighbors[24],neighbors[25],neighbors[26])
                    {
                        kernels[omp_get_thread_num()]->P2P(coord, taskParticlesTgt,
                                taskParticlesSrc, neighbors, counter);
                    }
                }

            }

        } while(octreeIterator.moveRight());

        FLOG( computationCounter.tac() );



        FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation L2P + P2P : " << computationCounter.cumulated() << " s\n" );
    }

};


#endif // FFMMALGORITHMOMP4_HPP

