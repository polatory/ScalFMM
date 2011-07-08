// /!\ Please, you must read the license at the bottom of this page

#include "../Src/Utils/FMpi.hpp"
#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/F3DPosition.hpp"
#include "../Src/Utils/FAbstractSendable.hpp"

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Components/FTestParticle.hpp"
#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"
#include "../Src/Extensions/FExtendPhysicalValue.hpp"

#include "../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FProcFmaLoader.hpp"

#include "../Src/Components/FBasicKernels.hpp"

#include <iostream>

#include <stdio.h>
#include <stdlib.h>


// Compile by : g++ testFmmAlgorithmProc.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmmAlgorithmProc.exe

/** This program show an example of use of
  * the fmm threaded + mpi algo
  * it also check that each particles is impacted each other particles
  */


/** Fmb class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FExtendPhysicalValue}
  */
class TestParticle : public FTestParticle, public FExtendPhysicalValue {
public:
};

class FTestCellPar : public FTestCell, public FAbstractSendable{
public :
    int bytesToSendUp() const{
        return sizeof(long);
    }
    int writeUp(void* const buffer, const int) const {
        const long tmpUp = getDataUp();
        memcpy(buffer,&tmpUp,bytesToSendUp());
        return bytesToSendUp();
    }
    int bytesToReceiveUp() const{
        return sizeof(long);
    }
    int readUp(void* const buffer, const int) {
        long tmpUp;
        memcpy(&tmpUp,buffer,bytesToSendUp());
        setDataUp(tmpUp);
        return bytesToReceiveUp();
    }

    int bytesToSendDown() const{
        return sizeof(long);
    }
    int writeDown(void* const buffer, const int) const {
        const long tmpDown = getDataDown();
        memcpy(buffer,&tmpDown,bytesToSendDown());
        return bytesToSendDown();
    }
    int bytesToReceiveDown() const{
        return sizeof(long);
    }
    int readDown(void* const buffer, const int) {
        long tmpDown;
        memcpy(&tmpDown,buffer,bytesToSendDown());
        setDataDown(tmpDown + getDataDown());
        return bytesToReceiveDown();
    }
};


/////////////////////////////////////////////////////////////////////////////
// Test function
/////////////////////////////////////////////////////////////////////////////

/** This function test the octree to be sure that the fmm algorithm
  * has worked completly.
  */
template<class OctreeClass, class ContainerClass, class FmmClassProc>
void ValidateFMMAlgoProc(OctreeClass* const badTree,
                         OctreeClass* const valideTree,
                         FmmClassProc*const fmm){
    const int OctreeHeight = badTree->getHeight();
    std::cout << "Check Result\n";
    {
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level > 0 ; --level){
            int NbLeafs = 0;
            do{
                ++NbLeafs;
            } while(octreeIterator.moveRight());
            octreeIterator.gotoLeft();

            const int startIdx = fmm->getLeft(NbLeafs);
            const int endIdx = fmm->getRight(NbLeafs);
            // Check that each particle has been summed with all other

            for(int idx = 0 ; idx < startIdx ; ++idx){
                octreeIterator.moveRight();
                octreeIteratorValide.moveRight();
            }

            for(int idx = startIdx ; idx < endIdx ; ++idx){
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Error index are not equal!" << std::endl;
                }
                else{
                    if(octreeIterator.getCurrentCell()->getDataUp() != octreeIteratorValide.getCurrentCell()->getDataUp()){
                        std::cout << "M2M error at level " << level << " up bad " << octreeIterator.getCurrentCell()->getDataUp()
                                << " good " << octreeIteratorValide.getCurrentCell()->getDataUp() << " idx " << idx << std::endl;
                    }
                    if(octreeIterator.getCurrentCell()->getDataDown() != octreeIteratorValide.getCurrentCell()->getDataDown()){
                        std::cout << "L2L error at level " << level << " down bad " << octreeIterator.getCurrentCell()->getDataDown()
                                << " good " << octreeIteratorValide.getCurrentCell()->getDataDown() << " idx " << idx << std::endl;
                    }
                }

                octreeIterator.moveRight();
                octreeIteratorValide.moveRight();
            }

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorValide.moveUp();
            octreeIteratorValide.gotoLeft();
        }
    }
    {
        int NbPart = 0;
        int NbLeafs = 0;
        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();
            do{
                NbPart += octreeIterator.getCurrentListSrc()->getSize();
                ++NbLeafs;
            } while(octreeIterator.moveRight());
            std::cout << "There is " << NbPart << " particles on " << NbLeafs << " Leafs" << std::endl;
        }
        {
            const int startIdx = fmm->getLeft(NbLeafs);
            const int endIdx = fmm->getRight(NbLeafs);
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            for(int idx = 0 ; idx < startIdx ; ++idx){
                octreeIterator.moveRight();
            }

            for(int idx = startIdx ; idx < endIdx ; ++idx){
                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

                const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());

                while( iter.hasNotFinished() ){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( (!isUsingTsm && iter.data().getDataDown() != NbPart - 1) ||
                        (isUsingTsm && iter.data().getDataDown() != NbPart) ){
                        std::cout << "Problem L2P + P2P, value on particle is : " << iter.data().getDataDown() << "\n";
                    }
                    iter.gotoNext();
                }
                octreeIterator.moveRight();
            }
        }
    }

    std::cout << "Done\n";
}

/** To print an octree
  * used to debug and understand how the values were passed
  */
template<class OctreeClass>
void print(OctreeClass* const valideTree){
    typename OctreeClass::Iterator octreeIterator(valideTree);
    for(int idxLevel = valideTree->getHeight() - 1 ; idxLevel > 1 ; --idxLevel ){
        do{
            std::cout << "[" << octreeIterator.getCurrentGlobalIndex() << "] up:" << octreeIterator.getCurrentCell()->getDataUp() << " down:" << octreeIterator.getCurrentCell()->getDataDown() << "\t";
        } while(octreeIterator.moveRight());
        std::cout << "\n";
        octreeIterator.gotoLeft();
        octreeIterator.moveDown();
    }
}

class MyFackParticle : public FBasicParticle {
    int myIndex;
public:
    MyFackParticle() : myIndex(0) {
    }
    void setIndex(const int inIndex){
        this->myIndex = inIndex;
    }
    int getIndex() const{
        return this->myIndex;
    }
};

// My leaf store the indexes of the particles it receives
// in a vector
class MyLeaf : public FAbstractLeaf<MyFackParticle, FVector<int> > {
    FVector<int> indexes;
public:
    void push(const MyFackParticle& particle){
        indexes.push( particle.getIndex() );
    }
    FVector<int>* getSrc(){
        return &indexes;
    }
    FVector<int>* getTargets(){
        return &indexes;
    }
};

struct ParticlesGroup {
    int number;
    int positionInArray;
    MortonIndex index;
    ParticlesGroup(const int inNumber = 0 , const int inPositionInArray = 0, const MortonIndex inIndex = 0)
        : number(inNumber), positionInArray(inPositionInArray), index(inIndex) {
    }
};


typedef TestParticle               ParticleClass;
typedef FTestCellPar               CellClass;
typedef FVector<ParticleClass>     ContainerClass;

typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
typedef FTestKernels<ParticleClass, CellClass, ContainerClass >         KernelClass;

typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
typedef FFmmAlgorithmThreadProc<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClassProc;


struct IndexedParticle{
    MortonIndex index;
    ParticleClass particle;
};

long getTreeCoordinate(const FReal inRelativePosition, const FReal boxWidthAtLeafLevel) {
        const FReal indexFReal = inRelativePosition / boxWidthAtLeafLevel;
        const long index = FMath::dfloor(indexFReal);
        if( index && FMath::LookEqual(inRelativePosition, boxWidthAtLeafLevel * index ) ){
                return index - 1;
        }
        return index;
}

long partition(IndexedParticle arr[], const long left, const long right) {
    long i = left, j = right;
    IndexedParticle tmp;
    IndexedParticle pivot = arr[(left + right) / 2];
    /* partition */
    while (i <= j) {
        while (arr[i].index < pivot.index)
            i++;
        while (arr[j].index > pivot.index)
            j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    }
    return j;
}
void quick_sort(IndexedParticle arr[], const long left, const long right){
    /* recursion */
    long part_index = partition(arr, left, right);
    if (left < part_index)
        quick_sort(arr, left, part_index);
    if (part_index + 1 < right)
        quick_sort(arr, part_index + 1, right);
}


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app( argc, argv);

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 9);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    char defaultFilename[] = "testLoaderFMA.fma"; //../../Data/ "testLoaderFMA.fma" "testFMAlgorithm.fma" Sphere.fma
    char* filename;
    FTic counter;

    if(argc == 1){
        std::cout << "You have to give a .fma file in argument.\n";
        std::cout << "The program will try a default file : " << defaultFilename << "\n";
        filename = defaultFilename;
    }
    else{
        filename = argv[1];
        std::cout << "Opening : " << filename << "\n";
    }

    FProcFmaLoader<ParticleClass> loader(filename,app);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }


    OctreeClass realTree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    {
        OctreeClass treeInterval(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
        int myNbParticlesCounter = 0;

        //////////////////////////////////////////////////////////////////////////////////
        // We sort our particles
        //////////////////////////////////////////////////////////////////////////////////
        {
            FVector<ParticlesGroup> groups;
            ParticleClass*const realParticles = reinterpret_cast<ParticleClass*>(new char[loader.getNumberOfParticles() * sizeof(ParticleClass)]);


            //////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////


            std::cout << "Inserting & Sorting my particles ..." << std::endl;
            counter.tic();

            {
                IndexedParticle*const realParticlesIndexed = reinterpret_cast<IndexedParticle*>(new char[loader.getNumberOfParticles() * sizeof(ParticleClass)]);
                F3DPosition boxCorner(loader.getCenterOfBox() - (loader.getBoxWidth()/2));
                FTreeCoordinate host;
                const FReal boxWidthAtLeafLevel = loader.getBoxWidth() / (2 << (NbLevels - 1) );
                for(long idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                    loader.fillParticle(realParticlesIndexed[idxPart].particle);
                    host.setX( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getX() - boxCorner.getX(), boxWidthAtLeafLevel ));
                    host.setY( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getY() - boxCorner.getY(), boxWidthAtLeafLevel ));
                    host.setZ( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getZ() - boxCorner.getZ(), boxWidthAtLeafLevel ));
                    realParticlesIndexed[idxPart].index = host.getMortonIndex(NbLevels - 1);
                }

                quick_sort(realParticlesIndexed,0,loader.getNumberOfParticles()-1);

                delete [] reinterpret_cast<char*>(realParticlesIndexed);
            }

            /*{
                OctreeClass sortingTree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

                ParticleClass particle;
                for(long idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                    loader.fillParticle(particle);
                    //printf("%f %f %f\n",particle.getPosition().getX(),particle.getPosition().getY(),particle.getPosition().getZ());
                    sortingTree.insert(particle);
                }

                //////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////
                int indexPart = 0;

                OctreeClass::Iterator octreeIterator(&sortingTree);
                octreeIterator.gotoBottomLeft();
                do{
                    ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
                    const MortonIndex indexAtThisLeaf = octreeIterator.getCurrentGlobalIndex();

                    groups.push(ParticlesGroup(octreeIterator.getCurrentListTargets()->getSize(),indexPart, indexAtThisLeaf));

                    while( iter.hasNotFinished() ){
                        realParticles[indexPart] = iter.data();

                        //std::cout << "Particles with index " << indexPart << " has a morton index of " << indexAtThisLeaf << std::endl;

                        //const F3DPosition& particlePosition = realParticles[indexPart].getPosition();
                        //std::cout << "\t The real position of this particle is (" << particlePosition.getX() << ";" << particlePosition.getY() << ";" << particlePosition.getZ() << ")" << std::endl;

                        ++indexPart;
                        iter.gotoNext();
                    }
                } while(octreeIterator.moveRight());

            }*/

            counter.tac();
            std::cout << "Done  " << "(@Creating and Inserting Temporary Particles = " << counter.elapsed() << "s)." << std::endl;


            //////////////////////////////////////////////////////////////////////////////////
            // We send the particle that do not belong to us
            //////////////////////////////////////////////////////////////////////////////////

            const MortonIndex min = app.broadcast( app.reduceMin(groups[0].index) );
            const MortonIndex max = app.broadcast( app.reduceMax(groups[groups.getSize() - 1].index) );

            const MortonIndex startIndex = app.getLeft(max - min + 1) + min;
            const MortonIndex endIndex = app.getRight(max - min + 1) + min;

            printf("%d we are going from %lld to %lld\n",app.processId(), min, max);
            printf("%d I will go from %lld to %lld\n",app.processId(), startIndex, endIndex);
            printf("There is actually %d leafs\n", groups.getSize());

            int*const needToReceive = new int[app.processCount() * app.processCount()];
            memset(needToReceive,0,app.processCount() * app.processCount() * sizeof(int));

            FMpi::Request requests[app.processCount()];
            {
                int needToSend[app.processCount()];
                memset(needToSend, 0, sizeof(int) * app.processCount());

                MortonIndex rightMortonIndex = min;
                int groudIndex = 0;
                for(int idxProc = 0 ; idxProc < app.processCount() && groudIndex < groups.getSize() ; ++idxProc){
                    rightMortonIndex = app.getOtherRight(max - min + 1, idxProc) + min;
                    printf("Working with %d, he goes to %lld \n",idxProc,rightMortonIndex);

                    if(idxProc != app.processId()){
                        int size = 0;
                        int currentGroupIndex = groudIndex;
                        while(groudIndex < groups.getSize() && groups[groudIndex].index < rightMortonIndex){
                            size += groups[groudIndex].number;
                            ++groudIndex;
                        }
                        needToSend[idxProc] = size;

                        printf("%d Send %d to %d\n",app.processId(), size, idxProc);
                        app.isendData( idxProc, sizeof(ParticleClass) * size, &realParticles[groups[currentGroupIndex].positionInArray], 1, &requests[idxProc]);
                    }
                    else{
                        needToSend[idxProc] = 0;
                        while(groudIndex < groups.getSize() && groups[groudIndex].index < rightMortonIndex){
                            const int end = groups[groudIndex].positionInArray + groups[groudIndex].number;
                            for(int idxPart = groups[groudIndex].positionInArray ; idxPart < end ; ++idxPart){
                                //std::cout << "\t I keep (" << realParticles[idxPart].getPosition().getX() << ";" << realParticles[idxPart].getPosition().getY() << ";" << realParticles[idxPart].getPosition().getZ() << ")" << std::endl;
                                treeInterval.insert(realParticles[idxPart]);
                                ++myNbParticlesCounter;
                            }
                            ++groudIndex;
                        }
                    }
                }

                app.allgather(needToSend, app.processCount(), needToReceive, app.processCount());
                for(int idxSrc = 0 ; idxSrc < app.processCount() ; ++idxSrc){
                    for(int idxTest = 0 ; idxTest < app.processCount() ; ++idxTest){
                        printf("[%d][%d] = %d\n", idxSrc, idxTest, needToReceive[idxSrc * app.processCount() + idxTest]);
                    }
                }
            }


            //////////////////////////////////////////////////////////////////////////////////
            // We receive others particles and insert them in the tree
            //////////////////////////////////////////////////////////////////////////////////
            int CounterProcToReceive(0);
            int maxPartToReceive(0);
            for(int idxProc = 0 ; idxProc < app.processCount() ; ++idxProc){
                if(idxProc != app.processId() && needToReceive[app.processCount() * idxProc + app.processId()]){
                    ++CounterProcToReceive;
                    if(maxPartToReceive < needToReceive[app.processCount() * idxProc + app.processId()]){
                        maxPartToReceive = needToReceive[app.processCount() * idxProc + app.processId()];
                    }
                    printf("needToReceive[idxProc][app.processId()] %d",needToReceive[app.processCount() * idxProc + app.processId()]);
                }
            }

            printf("maxPartToReceive %d \n",maxPartToReceive);

            ParticleClass*const iterParticles = reinterpret_cast<ParticleClass*>(new char[maxPartToReceive * sizeof(ParticleClass)]);
            // we receive message from nb proc - 1 (from every other procs
            for(int idxProc = 0 ; idxProc < CounterProcToReceive ; ++idxProc){
                int source(0);
                printf("Wait data to receive\n");
                app.waitMessage(&source);

                const int nbPartFromProc = needToReceive[app.processCount() * source + app.processId()];
                int received(0);

                printf("%d receive %d\n",source,nbPartFromProc);
                app.receiveDataFromTagAndSource(sizeof(ParticleClass) * nbPartFromProc, 1, source, iterParticles,&received);

                printf("Received %d part*partcileSize %d \n",received,sizeof(ParticleClass) * nbPartFromProc);

                printf("Insert into tree\n");
                for(int idxPart = 0 ; idxPart < nbPartFromProc ; ++idxPart){
                    //std::cout << "\t We receive a new particle (" << (*iterParticles).getPosition().getX() << ";" << (*iterParticles).getPosition().getY() << ";" << (*iterParticles).getPosition().getZ() << ")" << std::endl;
                    treeInterval.insert(iterParticles[idxPart]);
                    ++myNbParticlesCounter;
                }
            }

            printf("Wait all send\n");
            for(int idxProc = 0 ; idxProc < app.processCount() ; ++idxProc){
                if(idxProc != app.processId() && needToReceive[app.processCount() * app.processId() + idxProc ]){
                    app.iWait(&requests[idxProc]);
                }
            }

            printf("Delete particle array\n");
            delete [] reinterpret_cast<char*>(realParticles);
            delete [] needToReceive;
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Now we can build the real tree
        //////////////////////////////////////////////////////////////////////////////////


        {
            //////////////////////////////////////////////////////////////////////////////////
            // We inform the master proc about the data we have
            //////////////////////////////////////////////////////////////////////////////////
            printf("Inform other about leaves we have\n");

            FVector<ParticlesGroup> groups;
            ParticleClass*const realParticles = myNbParticlesCounter?new ParticleClass[myNbParticlesCounter]:0;

            int nbLeafs = 0;

            // we might now have any particles in our interval
            if(myNbParticlesCounter){
                OctreeClass::Iterator octreeIterator(&treeInterval);
                octreeIterator.gotoBottomLeft();
                int indexPart = 0;
                do{
                    ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
                    const MortonIndex indexAtThisLeaf = octreeIterator.getCurrentGlobalIndex();

                    groups.push(ParticlesGroup(octreeIterator.getCurrentListTargets()->getSize(),indexPart, indexAtThisLeaf));

                    while( iter.hasNotFinished() ){
                        realParticles[indexPart] = iter.data();
                        ++indexPart;
                        iter.gotoNext();
                    }

                    ++nbLeafs;
                } while(octreeIterator.moveRight());
            }

            printf("%d I have %d leafs\n",app.processId(), nbLeafs);

            // receive from left and right
            if(app.isMaster()){
                app.sendData( 1, sizeof(int), &nbLeafs, 0);
            }
            else if(app.processId() == app.processCount() - 1){
                app.sendData( app.processId() - 1, sizeof(int), &nbLeafs, 0);
            }
            // receive
            int leftLeafs = 0;
            int rightLeafs = 0;
            if(!app.isMaster() && app.processId() != app.processCount() - 1){
                for(int idxToReceive = 0 ; idxToReceive < 2 ; ++idxToReceive){
                    int source(0);
                    int temp = 0;
                    app.receiveDataFromTag(sizeof(int), 0, &temp, &source);
                    if(source < app.processId()){ // come from left
                        leftLeafs = temp;
                        temp += nbLeafs;
                        app.sendData( app.processId() + 1, sizeof(int), &temp, 0);
                    }
                    else { // come from right
                        rightLeafs = temp;
                        temp += nbLeafs;
                        app.sendData( app.processId() - 1, sizeof(int), &temp, 0);
                    }
                }
            }
            else {
                if(app.isMaster()){ // come from right
                    app.receiveDataFromTag(sizeof(int), 0, &rightLeafs);
                }
                else { // come from left
                    app.receiveDataFromTag(sizeof(int), 0, &leftLeafs);
                }
            }

            printf("%d I have %d leafs on the right and %d on the left\n",app.processId(), rightLeafs, leftLeafs);

            //////////////////////////////////////////////////////////////////////////////////
            // We balance the data
            //////////////////////////////////////////////////////////////////////////////////

            const int totalNbLeafs = (leftLeafs + nbLeafs + rightLeafs);
            const int myLeftLeaf = app.getLeft(totalNbLeafs);
            const int myRightLeaf = app.getRight(totalNbLeafs);

            const bool iNeedToSendToLeft = leftLeafs < myLeftLeaf;
            const bool iNeedToSendToRight = myRightLeaf < leftLeafs + nbLeafs;

            const bool iWillReceiveFromRight = leftLeafs + nbLeafs < myRightLeaf;
            const bool iWillReceiveFromLeft = leftLeafs > myLeftLeaf;

            const bool iDoNotHaveEnoughtToSendRight = myRightLeaf < leftLeafs;
            const bool iDoNotHaveEnoughtToSendLeft = leftLeafs + nbLeafs < myLeftLeaf;

            ParticleClass* rpart(0);
            int rpartSize(0);

            // Do I need to send to right?
            if(iNeedToSendToRight){
                int iNeedToSend = leftLeafs + nbLeafs - myRightLeaf;
                int iCanSend = nbLeafs;
                int idxSend (0);

                printf("%d I need to send %d to right\n",app.processId(), iNeedToSend);
                app.sendData( app.processId() + 1, sizeof(int), &iNeedToSend, 0);

                while(idxSend < iNeedToSend && idxSend < iCanSend){
                    app.sendData(app.processId() + 1, sizeof(int),&groups[nbLeafs - idxSend - 1].number,0);
                    app.sendData(app.processId() + 1, sizeof(ParticleClass) * groups[nbLeafs - idxSend - 1].number,
                                 &realParticles[groups[nbLeafs - idxSend - 1].positionInArray],0);
                    ++idxSend;
                }
                // I need to wait (idxSend == iCanSend && idxSend < iNeedToSend)
                if( iDoNotHaveEnoughtToSendRight ){ // I need to wait
                    int nbLeafsToRead(0);
                    app.receiveDataFromTag(sizeof(int), 1, &nbLeafsToRead);
                    for(int idxToRead = 0 ; idxToRead < nbLeafsToRead ; ++idxToRead){
                        int nbPartToRead(0);
                        app.receiveDataFromTag(sizeof(int), 1, &nbPartToRead);

                        if(rpartSize < nbPartToRead){
                            rpartSize = nbPartToRead;
                            delete [] (reinterpret_cast<char*>(rpart));
                            rpart = reinterpret_cast<ParticleClass*>(new char[nbPartToRead*sizeof(ParticleClass)]);
                        }

                        app.receiveDataFromTag(nbPartToRead*sizeof(ParticleClass), 0, rpart);
                        if(idxSend < iNeedToSend){
                            app.sendData(app.processId() + 1, sizeof(int),&nbPartToRead,0);
                            app.sendData(app.processId() + 1, sizeof(ParticleClass) * nbPartToRead, rpart,0);
                            ++idxSend;
                        }
                        else{
                            //insert into tree
                            for(int idxPart = 0 ; idxPart < nbPartToRead ; ++idxPart){
                                realTree.insert(rpart[idxPart]);
                                //const F3DPosition& particlePosition = rpart[idxPart].getPosition();
                                //std::cout << "\t I received (" << particlePosition.getX() << ";" << particlePosition.getY() << ";" << particlePosition.getZ() << ")" << std::endl;
                            }
                        }
                    }
                }
            }
            // will I receive from left
            if(iNeedToSendToLeft){
                int iNeedToSend = myLeftLeaf - leftLeafs;
                int iCanSend = nbLeafs;
                int idxSend (0);

                printf("%d I need to send %d to left",app.processId(), iNeedToSend);
                app.sendData( app.processId() - 1, sizeof(int), &iNeedToSend, 1);

                while(idxSend < iNeedToSend && idxSend < iCanSend){
                    app.sendData(app.processId() - 1, sizeof(int),&groups[idxSend].number,1);
                    app.sendData(app.processId() - 1, sizeof(ParticleClass) * groups[idxSend].number,
                                 &realParticles[groups[idxSend].positionInArray],1);
                    ++idxSend;
                }
                // Can I do it now?
                if( iDoNotHaveEnoughtToSendLeft ){
                    int nbLeafsToRead(0);
                    app.receiveDataFromTag(sizeof(int), 1, &nbLeafsToRead);
                    for(int idxToRead = 0 ; idxToRead < nbLeafsToRead ; ++idxToRead){
                        int nbPartToRead(0);
                        app.receiveDataFromTag(sizeof(int), 1, &nbPartToRead);

                        if(rpartSize < nbPartToRead){
                            rpartSize = nbPartToRead;
                            delete [] (reinterpret_cast<char*>(rpart));
                            rpart = reinterpret_cast<ParticleClass*>(new char[nbPartToRead*sizeof(ParticleClass)]);
                        }

                        app.receiveDataFromTag(nbPartToRead*sizeof(ParticleClass), 1, rpart);
                        if(idxSend < iNeedToSend){
                            app.sendData(app.processId() - 1, sizeof(int),&nbPartToRead,1);
                            app.sendData(app.processId() - 1, sizeof(ParticleClass) * nbPartToRead, rpart,1);
                            ++idxSend;
                        }
                        else{
                            for(int idxPart = 0 ; idxPart < nbPartToRead ; ++idxPart){
                                realTree.insert(rpart[idxPart]);
                                //const F3DPosition& particlePosition = rpart[idxPart].getPosition();
                                //std::cout << "\t I receive (" << particlePosition.getX() << ";" << particlePosition.getY() << ";" << particlePosition.getZ() << ")" << std::endl;
                            }
                        }
                    }
                }
            }

            // If i will receive from left and I did no already have
            if(!(iNeedToSendToRight && iDoNotHaveEnoughtToSendRight) && iWillReceiveFromLeft){
                printf("%d I need to receive from left\n",app.processId());
                int nbLeafsToRead(0);
                app.receiveDataFromTag(sizeof(int), 0, &nbLeafsToRead);
                printf("%d I will receive from left %d\n",app.processId(), nbLeafsToRead);
                for(int idxToRead = 0 ; idxToRead < nbLeafsToRead ; ++idxToRead){
                    int nbPartToRead(0);
                    app.receiveDataFromTag(sizeof(int), 0, &nbPartToRead);
                    //printf("%d I will receive %d particles\n",app.processId(), nbPartToRead);
                    if(rpartSize < nbPartToRead){
                        rpartSize = nbPartToRead;
                        delete [] (reinterpret_cast<char*>(rpart));
                        rpart = reinterpret_cast<ParticleClass*>(new char[nbPartToRead*sizeof(ParticleClass)]);
                    }

                    app.receiveDataFromTag(nbPartToRead*sizeof(ParticleClass), 0, rpart);
                    for(int idxPart = 0 ; idxPart < nbPartToRead ; ++idxPart){
                        realTree.insert(rpart[idxPart]);
                        //const F3DPosition& particlePosition = rpart[idxPart].getPosition();
                        //std::cout << "\t I received (" << particlePosition.getX() << ";" << particlePosition.getY() << ";" << particlePosition.getZ() << ")" << std::endl;
                    }
                }
            }
            // If i will receive from right and I did no already have
            if(!(iNeedToSendToLeft && iDoNotHaveEnoughtToSendLeft) && iWillReceiveFromRight){
                printf("%d I need to receive from right\n",app.processId());
                int nbLeafsToRead(0);
                app.receiveDataFromTag(sizeof(int), 1, &nbLeafsToRead);
                //printf("%d I will receive from right %d\n",app.processId(), nbLeafsToRead);
                for(int idxToRead = 0 ; idxToRead < nbLeafsToRead ; ++idxToRead){
                    int nbPartToRead(0);
                    app.receiveDataFromTag(sizeof(int), 1, &nbPartToRead);
                    printf("%d I will receive %d particles\n",app.processId(), nbPartToRead);
                    if(rpartSize < nbPartToRead){
                        rpartSize = nbPartToRead;
                        delete [] (reinterpret_cast<char*>(rpart));
                        rpart = reinterpret_cast<ParticleClass*>(new char[nbPartToRead*sizeof(ParticleClass)]);
                    }

                    app.receiveDataFromTag(nbPartToRead*sizeof(ParticleClass), 1, rpart);
                    for(int idxPart = 0 ; idxPart < nbPartToRead ; ++idxPart){
                        realTree.insert(rpart[idxPart]);
                        //const F3DPosition& particlePosition = rpart[idxPart].getPosition();
                        //std::cout << "\t I received (" << particlePosition.getX() << ";" << particlePosition.getY() << ";" << particlePosition.getZ() << ")" << std::endl;
                    }
                }
            }

            printf("Will now take my own particles from %d to %d\n",FMath::Max(myLeftLeaf-leftLeafs,0) , FMath::Min(myRightLeaf,totalNbLeafs- rightLeafs) - myLeftLeaf);
            // insert the particles we already have
            for(int idxLeafInsert = FMath::Max(myLeftLeaf-leftLeafs,0) ; idxLeafInsert < FMath::Min(myRightLeaf,totalNbLeafs- rightLeafs) - myLeftLeaf ; ++idxLeafInsert){
                for(int idxPart = 0 ; idxPart < groups[idxLeafInsert].number ; ++idxPart){
                    realTree.insert(realParticles[groups[idxLeafInsert].positionInArray + idxPart]);
                    //const F3DPosition& particlePosition = realParticles[groups[idxLeafInsert].positionInArray + idxPart].getPosition();
                    //std::cout << "\t Position is (" << particlePosition.getX() << ";" << particlePosition.getY() << ";" << particlePosition.getZ() << ")" << std::endl;
                }
            }

            delete [] reinterpret_cast<char*>(rpart);
        }

        app.processBarrier();

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Create real tree
    //////////////////////////////////////////////////////////////////////////////////

    OctreeClass treeValide(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
    {
        FFmaBinLoader<ParticleClass> loaderSeq(filename);
        ParticleClass partToInsert;
        for(int idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
            loaderSeq.fillParticle(partToInsert);
            treeValide.insert(partToInsert);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Check particles in tree
    //////////////////////////////////////////////////////////////////////////////////

    {
        int totalNbLeafs = 0;
        {

            OctreeClass::Iterator octreeIterator(&treeValide);
            octreeIterator.gotoBottomLeft();
            do{
                ++totalNbLeafs;
            }while(octreeIterator.moveRight());
        }

        const int myLeftLeaf = app.getLeft(totalNbLeafs);
        const int myRightLeaf = app.getRight(totalNbLeafs);

        printf("%d should go from %d to %d leaf (on %d total leafs)\n", app.processId(), myLeftLeaf, myRightLeaf, totalNbLeafs);

        OctreeClass::Iterator octreeIteratorValide(&treeValide);
        octreeIteratorValide.gotoBottomLeft();
        for(int idxLeaf = 0 ; idxLeaf < myLeftLeaf ; ++idxLeaf){
            if(!octreeIteratorValide.moveRight()){
                printf("Error cannot access to the left leaf %d in the valide tree\n", myLeftLeaf);
            }
        }

        OctreeClass::Iterator octreeIterator(&realTree);
        octreeIterator.gotoBottomLeft();

        for(int idxLeaf = myLeftLeaf ; idxLeaf < myRightLeaf ; ++idxLeaf){
            if(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                printf("Error index are different, valide %lld invalid %lld\n",octreeIteratorValide.getCurrentGlobalIndex(),
                       octreeIterator.getCurrentGlobalIndex());
                break;
            }
            if(octreeIteratorValide.getCurrentListSrc()->getSize() != octreeIterator.getCurrentListSrc()->getSize()){
                printf("Error leafs do not have the same number of particles, valide %d, invalide %d\n",
                       octreeIteratorValide.getCurrentListSrc()->getSize(), octreeIterator.getCurrentListSrc()->getSize() );
            }

            //printf("index %lld with %d particles\n", octreeIteratorValide.getCurrentGlobalIndex(), octreeIteratorValide.getCurrentListSrc()->getSize());

            if(!octreeIteratorValide.moveRight() && idxLeaf != myRightLeaf - 1){
                printf("Error cannot valide tree end to early, idxLeaf %d myRightLeaf %d\n", idxLeaf, myRightLeaf);
                break;
            }

            if(!octreeIterator.moveRight() && idxLeaf != myRightLeaf - 1){
                printf("Error cannot test tree end to early, idxLeaf %d myRightLeaf %d\n", idxLeaf, myRightLeaf);
                break;
            }
        }



    }
    return 0;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;

    FmmClassProc algo(app,&realTree,&kernels);
    algo.execute();

    FmmClass algoValide(&treeValide,&kernels);
    algoValide.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    ValidateFMMAlgoProc<OctreeClass,ContainerClass,FmmClassProc>(&realTree,&treeValide,&algo);

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}


// [--LICENSE--]
