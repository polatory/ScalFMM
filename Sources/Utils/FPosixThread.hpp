#ifndef FPOSIXTHREAD_HPP
#define FPOSIXTHREAD_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <pthread.h>

#include "FAbstractThread.hpp"
#include "FDebug.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FPosixThread
* Please read the license
*
* This class enable posix threading.
*
* @warning You have to put your class name as the template when inheriting.
*
* <code>
* // Example with FPosixThread <br>
* class TPosix : public FPosixThread<TOpen>{ <br>
* public: <br>
*	void threadCallback(const int inThreadId, const int){ <br>
*		printf("I am %d\n",inThreadId); <br>
*	} <br>
* }; <br>
* // ... <br>
* TPosix posix; <br>
* posix.executeThreads(10); <br>
* </code>
*/
template <class Derived, int DefaultThreadsNumber = 5>
class FPosixThread : public FAbstractThread<Derived> {
private:
	mutable pthread_mutex_t mutex;		//< Posix mutex
	mutable pthread_barrier_t pbarrier;	//< To use a barrier

public:
	/**
	* Constructor
	* just init mutex
	*/
        FPosixThread(){
		initMutex();
	}

	/**
	* Copy constructor
	* just init mutex
	*/
        FPosixThread(const FPosixThread&){
		initMutex();
	}

	/** 
	* Destructor
	* just destroy mutex
	*/
        virtual ~FPosixThread(){
		pthread_mutex_unlock(&this->mutex);
        	pthread_mutex_destroy(&this->mutex);
	}

	/**
	* This function is used to create inThreadsNumber threads with inCallback as the callback.
	* @param inCallback the callback (must be a object method)
	* @param inThreadsNumber the number of threads to create (default is DefaultThreadsNumber)
        * <code> posix.executeThreads(&TPosix::threadCallback,10); </code>
	*/
	void executeThreads(void (Derived::*inCallback)(const int, const int), const int inThreadsNumber = DefaultThreadsNumber){
		// init barrier
		pthread_barrier_init( &this->pbarrier, 0, inThreadsNumber);

		// Get current object address as a Derived class
		Derived* const thisDerived = dynamic_cast<Derived*>(this);

		// One descriptor per thread
		PosixDescriptor descriptors[inThreadsNumber];

		// Create thread
		for(long indexThread = 1 ; indexThread < inThreadsNumber ; ++indexThread){
			descriptors[indexThread].target = 	thisDerived;
			descriptors[indexThread].callback = 	inCallback;
			descriptors[indexThread].threadId = 	indexThread;
			descriptors[indexThread].numThreads = 	inThreadsNumber;

			pthread_create(&descriptors[indexThread].thread, 0, PthreadCallback, (void*)&descriptors[indexThread]);
		}

		// Master thread
		(thisDerived->*inCallback)(0 , inThreadsNumber);

		// Wait each thread to finish
		for(long indexThread = 1 ; indexThread < inThreadsNumber ; ++indexThread){
			pthread_join(descriptors[indexThread].thread, 0);
		}

		// destroy barrier
		pthread_barrier_destroy( &this->pbarrier );
	}

	/**
	* This function is used to create inThreadsNumber threads with Object::threadCallback as the callback.
	* @param inThreadsNumber the number of threads to create (default is DefaultThreadsNumber)
        * @warning You have to implement threadCallback
	*/
	void executeThreads(const int inThreadsNumber = DefaultThreadsNumber){
                executeThreads( &FPosixThread::threadCallback , inThreadsNumber);
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
	* Init mutext
	* Equal this->mutex = PTHREAD_MUTEX_INITIALIZER;
	*/
	void initMutex(){
		pthread_mutexattr_t attr;
		pthread_mutexattr_init(&attr);
		pthread_mutexattr_settype(&attr,PTHREAD_MUTEX_RECURSIVE);
		pthread_mutex_init(&this->mutex,&attr);
		pthread_mutexattr_destroy(&attr);
	}

	/**
	* This function lock a posix mutex
	*/
	void lock() const {
		pthread_mutex_lock( &this->mutex );
	}

	/**
	* This function unlock a posix mutex
	*/
	void unlock() const {
		pthread_mutex_unlock( &this->mutex );
	}
	
	/**
	* Barrier to sync all thread
	*/
	void barrier() const {
		pthread_barrier_wait( &this->pbarrier );
	}

private:

	/**
	* This struct is useless for users, it is use to describe a posix thread
	*/
	typedef void (Derived::*Callback)(const int, const int);
	struct PosixDescriptor{
		Derived * target;	//< object to call
		Callback callback;	//< method to call
		int threadId;		//< current thread position
		int numThreads;		//< number of threads
		pthread_t thread;	//< posix descriptor to enable wait/join functions
	};

	/**
	* This function is the normal posix thread callback
	* @param inThreadDescriptor a PosixDescriptor of the current thread
	* @return 0
	*/
	static void* PthreadCallback(void* const inThreadDescriptor)
	{
		// Simply cast the parameter and call the function
		PosixDescriptor* const threadDescriptor = (PosixDescriptor*) inThreadDescriptor;
		(threadDescriptor->target->*threadDescriptor->callback)(threadDescriptor->threadId,threadDescriptor->numThreads);
		return 0;
	}
};

#endif //FPOSIXTHREAD_HPP

// [--LICENSE--]
