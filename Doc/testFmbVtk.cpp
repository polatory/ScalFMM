/*
  In the cmake

OPTION( SCALFMM_USE_VTK "Set to On to use VTK" ON )

# VTK option has been set
if( SCALFMM_USE_VTK )
        MESSAGE( STATUS "SCALFMM_USE_VTK = ${SCALFMM_USE_VTK}" )
        FIND_PACKAGE(VTK REQUIRED)
        INCLUDE(${VTK_USE_FILE})
endif()

    In scalfmm config
///////////////////////////////////////////////////////
// VTK
///////////////////////////////////////////////////////

#cmakedefine SCALFMM_USE_VTK

   in Tests/CMakeList.txt

 target_link_libraries(${execname} vtkRendering vtkWidgets)
  */

#include "../Src/Utils/FGlobal.hpp"
#ifndef SCALFMM_USE_VTK

#include <cstdio>

int main(){
    printf("VTK has been set to unused in the compile config file.\n");
    return 0;
}

#else


#ifndef FTHREAD_HPP
#define FTHREAD_HPP

// Use Win or Posix
#ifdef WINDOWS
    #include <windows.h>
#else
    #ifndef POSIX
        #warning POSIX will be used (but you did not define it)
    #endif
    #include <pthread.h>
    #include <signal.h>
#endif

class FThread {
private:
    bool isActive;      /**< Fast bool lookup */

#ifdef WINDOWS
    HANDLE handles;     /**< Windows thread */
#else
    pthread_t handles;  /**< Posix Thread*/
#endif


    /**
    * @brief Denied equality operator
    * @param none
    */
    void operator=(const FThread &){}
    /**
    * @brief Denied copy constructor
    * @param none
    */
    FThread(const FThread &){}

    /**
    * @brief Static starter function to execute posix thread
    * @brief This function set thread->isActive to false
    */
#ifdef WINDOWS
    static DWORD WINAPI Starter(LPVOID in_thread){
#else
    static void* Starter(void* in_thread){
#endif
        FThread * thread = static_cast< FThread * >(in_thread);
        thread->run();
        thread->isActive = false;
        return 0x00;
    }

public:
    /**
    * @brief Constructor
    */
    FThread(): isActive(false) {
    #ifdef WINDOWS
        handles = 0x00;
    #else
    #endif
    }
    /**
    * @brief Destructor, Warning, it waits the end of the current thread
    */
    virtual ~FThread(){
        if(!isActive) return;
        // if we destroy the thread until it has finished
        // there is a problem in your implementation algorithm
        // So we wait before destroying the thread!
        wait();
    #ifdef WINDOWS
        CloseHandle (handles);
    #else
    #endif
    }

    /**
    * @brief start the thread
    * @return true if success else false
    */
    bool start(){
        if(isActive) return false;
        isActive = true;
#ifdef WINDOWS
        handles = CreateThread( 0x00, 0x00,FThread::Starter, static_cast< void* >(this), 0x00, 0x00);
        return handles != NULL;
#else
        return pthread_create(&handles, NULL, FThread::Starter, static_cast< void* >(this)) == 0;
#endif
    }

    /**
    * @brief Fast look up to know if a thread is running
    * @return true if running else false
    */
    bool isRunning() const{
        return isActive;
    }

    /**
    * @brief Wait the end of a thread
    * @return false in case of error, true if all right
    */
    bool wait() const{
        if(!isActive) return false;
#ifdef WINDOWS
        return WaitForSingleObject(handles,INFINITE) == 0x00000000L;
#else
        return pthread_join(handles, NULL) == 0;
#endif

    }

    /**
    * @brief the function is called when thread is starting
    * @must You must implement this methode!
    */
    virtual void run() = 0;

    /**
    * Kill the thread
    */
    bool kill(){
        if(!isActive) return false;

        isActive = false;
#ifdef WINDOWS
        bool success = TerminateThread(handles,1) && CloseHandle(handles);
        handles = 0x00;
        return success;
#else
        return pthread_kill( handles, SIGKILL) == 0;
#endif
    }

};

#endif


#ifndef FMUTEX_HPP
#define FMUTEX_HPP

// Use Window or Posix
#ifdef WINDOWS
     #include <windows.h>
#else
     #ifndef POSIX
          #warning POSIX will be used (but you did not define it)
     #endif
     #include <pthread.h>
#endif

/**
* This class represent a simple way to use mutex
*
* @example FMutex mut;
* @example mut.lock();      // lock
* @example ...
* @example mut.islocked();  // fast look up
* @example ...
* @example mut.unlock();    // unlock
* @example mut.tryLock();   // try lock
*
* Ressources : http://www.codeproject.com/KB/threads/thread_class.aspx
*
* @must You may have to change this class if you are not on Windows
* @must or Posix OS
*/

class FMutex {
private:

#ifdef WINDOWS
     CRITICAL_SECTION _mutex; /**< Window mutex */
#else
     pthread_mutex_t _mutex; /**< posix mutex */
#endif

     bool locked;           /**< Fast locked look up used for copying */

     void init(){
     #ifdef WINDOWS
          InitializeCriticalSection(&_mutex);
     #else
          pthread_mutexattr_t attr;
          pthread_mutexattr_init(&attr);
          pthread_mutexattr_settype(&attr,PTHREAD_MUTEX_RECURSIVE);
          pthread_mutex_init(&_mutex,&attr);
          pthread_mutexattr_destroy(&attr);
     #endif
          locked = false;
     }

public:

     /**
     * @brief Construct a FMutex
     * @brief Posix and Win mutex
     */
     FMutex(){
          init();
     }
     /**
     * @brief Copy Constructor a mutex (copy the locked state only)
     * @param Based mutex
     *
     */
     FMutex( const FMutex &in_mutex ) {
          init();

          if(in_mutex.locked && !locked) lock();
          else if(!in_mutex.locked && locked) unlock();
     }

     /**
     * @brief Copy a mutex (copy the locked state only)
     * @param Based mutex
     * @return Current mutex
     */
     FMutex& operator=(const FMutex &in_mutex) {
          if(in_mutex.locked && !locked) lock();
          else if(!in_mutex.locked && locked) unlock();
          return *this;
     }

     /**
     * @brief Destructor
     */
     virtual ~FMutex(){
     #ifdef WINDOWS
          DeleteCriticalSection(&_mutex);
     #else
          pthread_mutex_unlock(&_mutex);
          pthread_mutex_destroy(&_mutex);
     #endif
     }

     /**
     * @brief lock a mutex
     * @return WIN true
     * @return POSIX true if success
     */
     bool lock(){
          locked = true;
     #ifdef WINDOWS
          EnterCriticalSection(&_mutex);
          return true;
     #else
          return pthread_mutex_lock(&_mutex) == 0;
     #endif
     }

     /**
     * @brief lock a mutex
     * @return true if success else false (if busy or error)
     */
     bool tryLock(){
          locked = true;
     #ifdef WINDOWS
          return TryEnterCriticalSection(&_mutex);
     #else
          return pthread_mutex_trylock(&_mutex) == 0;
     #endif
     }

     /**
     * @brief unlock a mutex
     * @return WIN true in every cases
     * @return POSIX true if success
     */
     bool unlock(){
          locked = false;
     #ifdef WINDOWS
          LeaveCriticalSection(&_mutex);
          return true;
     #else
          return pthread_mutex_unlock(&_mutex) == 0;
     #endif
     }

     /**
     * @brief Fast locked look up
     * @return true if locked else false
     * This methode use the fast look up variable but still CONST
     * if you want to test perfectly do :
     * if(myMutex.tryLock()){
     *      myMutex.unlock();
     *      // I am sure that the mutex is not locked
     * }
     * else{
     *      // The mutex is locked
     * }
     */
     bool isLocked() const{
          return locked;
     }

};

#endif



// [--License--]



#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkActor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCommand.h"

class FVTKScene {
    vtkRenderer *render;
    vtkRenderWindow *window;
    vtkRenderWindowInteractor *interactor;

public:
    FVTKScene() : render(0), window(0), interactor(0) {
        render = vtkRenderer::New();
        window = vtkRenderWindow::New();
        window->SetSize(800,600);
        window->AddRenderer(render);

        // an interactor
        interactor = vtkRenderWindowInteractor::New();
        interactor->SetRenderWindow(window);
        interactor->Initialize();

        render->SetBackground(0,0,0); // Background color white
    }

    ~FVTKScene(){
        render->Delete();
        window->Delete();
        interactor->Delete();
    }

    void addActor(vtkActor* const actor){
        // add the actor to the scene
        render->AddActor(actor);
    }

    template<class ObjectClass>
    void addActors(ObjectClass actors[], const FSize size){
        for( int idxActor = 0 ; idxActor < size ; ++idxActor){
            addActor(actors[idxActor].getActor());
        }
    }

    void addCommand(vtkCommand* command){
          interactor->AddObserver(vtkCommand::TimerEvent, command);
          interactor->CreateRepeatingTimer(250);
    }

    void show(){
        // render an image (lights and cameras are created automatically)
        window->Render();
        // begin mouse interaction
        interactor->Start();
    }
};


class FVTKSphere {
    static const double Radius = 1.0;
    static const int ThetaResolution = 4;
    static const int PhiResolution = 4;

    vtkSphereSource *sphere;
    vtkPolyDataMapper *mapper;
    vtkActor *actor;
public:
    FVTKSphere(const F3DPosition& center = F3DPosition(0,0,0))
        : sphere(0), mapper(0), actor(0) {
       sphere = vtkSphereSource::New();
       sphere->SetRadius(Radius);
       sphere->SetThetaResolution(ThetaResolution);
       sphere->SetPhiResolution(PhiResolution);
       sphere->SetCenter(center.getX(),center.getY(),center.getZ());

       mapper = vtkPolyDataMapper::New();
       mapper->SetInput(sphere->GetOutput());

       actor = vtkActor::New();
       actor->SetMapper(mapper);
       actor->GetProperty()->SetColor(0,0,1); // sphere color blue
    }

    ~FVTKSphere(){
        sphere->Delete();
        mapper->Delete();
        actor->Delete();
    }

    void setPosition(const F3DPosition& center){
        //sphere->SetCenter(center.getX(),center.getY(),center.getZ());
        actor->SetPosition(center.getX(),center.getY(),center.getZ());
    }

    vtkActor* getActor(){
        return actor;
    }
};


class FVtkTimer : public vtkCommand {
    int size;
    F3DPosition* sharedPositions;
    FVTKSphere* actors;
    FMutex* mutex;

public:
    static FVtkTimer* New(const int inSize, F3DPosition*const inSharedPositions, FVTKSphere*const inActors, FMutex*const inMutex) {
      FVtkTimer* command = new FVtkTimer;
      command->sharedPositions = inSharedPositions;
      command->actors = inActors;
      command->size = inSize;
      command->mutex = inMutex;
      return command;
    }

    virtual void Execute(vtkObject *caller, unsigned long eventId,
                         void * vtkNotUsed(callData)) {
        mutex->lock();
        for(int idx = 0 ; idx < size ; ++idx){
            actors[idx].setPosition(sharedPositions[idx]);
        }
        mutex->unlock();

        vtkRenderWindowInteractor* const iren = vtkRenderWindowInteractor::SafeDownCast(caller);
        iren->GetRenderWindow()->Render();
    }
};

class FSceneThreaded : public FThread {
    const int nbSphere;
    F3DPosition* sharedPositions;
    FMutex mutex;

    void run(){
        FVTKScene scene;

        FVTKSphere* actors = new FVTKSphere[nbSphere];
        for(int idx = 0 ; idx < nbSphere ; ++idx){
            actors[idx].setPosition(sharedPositions[idx]);
        }

        scene.addActors((FVTKSphere*)actors, nbSphere);

        FVtkTimer* command = FVtkTimer::New(nbSphere, sharedPositions, actors, &mutex);
        scene.addCommand(command);

        scene.show();

        delete[] actors;
        command->Delete();
    }

public:
    FSceneThreaded(const int inNbSphere) : nbSphere(inNbSphere), sharedPositions(0) {
        sharedPositions = new F3DPosition[inNbSphere];
    }
    ~FSceneThreaded(){
        delete[] sharedPositions;
    }
    void setPosition(const int index, const F3DPosition& center){
        mutex.lock();
        sharedPositions[index] = center;
        mutex.unlock();
    }

    template <class OctreeClass, class ContainerClass>
    void mapToOctree(OctreeClass*const tree){
        mutex.lock();
        int idxPart = 0;
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                setPosition( idxPart++, iter.data().getPosition() );
                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());
        mutex.unlock();
    }

};


#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Kernels/FElecForcesKernels.hpp"
#include "../Src/Fmb/FFmbComponents.hpp"

#include "../Src/Extensions/FExtendVelocity.hpp"

#include "../Src/Files/FRandomLoader.hpp"
#include "../Src/Files/FFmaLoader.hpp"
#include "../Src/Arranger/FOctreeArranger.hpp"


class FmbVeloParticle : public FmbParticle, public FExtendVelocity {
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FmbVeloParticle         ParticleClass;
    typedef FmbCell                 CellClass;
    typedef FVector<ParticleClass>  ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FElecForcesKernels<ParticleClass, CellClass, ContainerClass >   KernelClass;

    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test fmb algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const long NbPart       = FParameters::getValue(argc,argv,"-nb", 30);
    const FReal DT          = FParameters::getValue(argc,argv,"-dt", FReal(0.1));
    const int DevP          = FParameters::getValue(argc,argv,"-p", 5);

    //FRandomLoader<ParticleClass> loader(NbPart, 30, F3DPosition(0.5,0.5,0.5), 1);
    FFmaLoader<ParticleClass> loader("../Data/galaxy.fma.tmp");

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    {
        ParticleClass particleToFill;
        particleToFill.setPhysicalValue(FReal(0.10));

        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particleToFill);
            tree.insert(particleToFill);
        }
    }


    // -----------------------------------------------------

    KernelClass kernels( DevP, NbLevels, loader.getBoxWidth());
    FmmClass algo( &tree, &kernels);
    FOctreeArranger<OctreeClass, ContainerClass, ParticleClass> arranger(&tree);
    FSceneThreaded scene(int(loader.getNumberOfParticles()));

    scene.mapToOctree<OctreeClass, ContainerClass>(&tree);
    scene.start();

    while(scene.isRunning()){
        algo.execute();
        { // update velocity and position
            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            do{
                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
                while( iter.hasNotFinished() ){
                    kernels.computeVelocity(&iter.data(), DT);
                    kernels.updatePosition(&iter.data(), DT);
                    iter.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }
        // update tree and vtk
        arranger.rearrange();
        scene.mapToOctree<OctreeClass, ContainerClass>(&tree);

        usleep(10000000);// TODO delete
    }

    // -----------------------------------------------------

    return 0;
}

#endif
// [--END--]
