#ifndef FOPENMPTHREAD_HPP
#define FOPENMPTHREAD_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <omp.h>

#include "FAbstractThread.hpp"
#include "FDebug.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FOpenMPThread
* Please read the license
*
* This class enable openmp threading.
*
* @warning You have to put your class name as the template when inheriting.
*
* <code>
* // Example with FOpenMPThread <br>
* class TOpen : public FOpenMPThread<TOpen>{ <br>
* public: <br>
*	void threadCallback(const int inThreadId, const int){ <br>
*		printf("I am %d\n",inThreadId); <br>
*	} <br>
* }; <br>
* // ... <br>
* TOpen open; <br>
* open.executeThreads(10); <br>
* </code>
*/
template <class Derived, int DefaultThreadsNumber = 5>
class FOpenMPThread : public FAbstractThread<Derived> {
private:
	mutable omp_lock_t mutex;	//< openmp mutex

public:
	/**
	* Constructor
	* just init mutex
	*/
        FOpenMPThread(){
		omp_init_lock(&this->mutex);
	}

	/**
	* Copy constructor
	* just init mutex
	*/
        FOpenMPThread(const FOpenMPThread&){
		omp_init_lock(&this->mutex);
	}

	/** 
	* Destructor
	* just destroy mutex
	*/
        virtual ~FOpenMPThread(){
		omp_destroy_lock(&this->mutex);
	}

	/**
	* This function is used to create inThreadsNumber threads with inCallback as the callback.
	* @param inCallback the callback (must be a object method)
	* @param inThreadsNumber the number of threads to create (default is DefaultThreadsNumber)
        * <code> open.executeThreads(&TOpen::threadCallback,10); </code>
	*/
	void executeThreads(void (Derived::*inCallback)(const int, const int), const int inThreadsNumber = DefaultThreadsNumber){
		Derived* const thisDerived = dynamic_cast<Derived*>(this);
		#pragma omp parallel num_threads(inThreadsNumber)
		{
			(thisDerived->*inCallback)(omp_get_thread_num(), omp_get_num_threads());
		}
	}

	/**
	* This function is used to create inThreadsNumber threads with Object::threadCallback as the callback.
	* @param inThreadsNumber the number of threads to create (default is DefaultThreadsNumber)
        * @warning You have to implement threadCallback
	*/
	void executeThreads(const int inThreadsNumber = DefaultThreadsNumber){
		#pragma omp parallel num_threads(inThreadsNumber)
		{
			threadCallback(omp_get_thread_num(), omp_get_num_threads());
		}
	}

	/**
	* This function executed by each thread after executeThreads has been called
	* @param inThreadId the current thread index
	* @param inThreadNumbers the number of threads started
        * @warning Must be impleted in the derived class
	*/
	virtual void threadCallback(const int inThreadId, const int inThreadNumbers){
		FDEBUG( FDebug::Controller.writeFromLine("[W] You called executeThreads() but did not implement threadCallback", __LINE__, __FILE__); )
	};

protected:
	/**
	* This function lock an openmp mutex
	*/
	void lock() const {
		omp_set_lock(&this->mutex);	
	}

	/**
	* This function unlock an openmp mutex
	*/
	void unlock() const {
		omp_unset_lock(&this->mutex);
	}

	/**
	* Barrier to sync all thread
	*/
	void barrier() const {
		#pragma omp barrier
	}
};

#endif //FOPENMPTHREAD_HPP

// [--LICENSE--]
