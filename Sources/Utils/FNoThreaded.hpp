#ifndef FNOTHREADED_HPP
#define FNOTHREADED_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <omp.h>

#include "FAbstractThreaded.hpp"
#include "FDebug.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FNoThreaded
* Please read the license
*
* This class do not use thread. It is used on system that do not allow openmp or posix threads,
* or to make some test without these libs.
*
* @warning You have to put your class name as the template when inheriting.
*
* <code>
* // Example with FNoThreaded <br>
* class TNo : public FNoThreaded<TNo>{ <br>
* public: <br>
*	void threadCallback(const int inThreadId, const int){ <br>
*		printf("I am %d\n",inThreadId); <br>
*	} <br>
* }; <br>
* // ... <br>
* TNo notd; <br>
* notd.executeThreads(10); <br>
* </code>
*/
template <class Derived, int DefaultThreadsNumber = 1>
class FNoThreaded : public FAbstractThreaded<Derived> {
public:
	/**
	* This function is used to create inThreadsNumber threads with inCallback as the callback.
	* @param inCallback the callback (must be a object method)
	* @param inThreadsNumber the number of threads to create (default is DefaultThreadsNumber)
        * <code> open.executeThreads(&TOpen::threadCallback,10); </code>
	*/
	void executeThreads(void (Derived::*inCallback)(const int, const int), const int inThreadsNumber = DefaultThreadsNumber){
		Derived* const thisDerived = dynamic_cast<Derived*>(this);
		(thisDerived->*inCallback)(0, 1);
	}

	/**
	* This function is used to create inThreadsNumber threads with Object::threadCallback as the callback.
	* @param inThreadsNumber the number of threads to create (default is DefaultThreadsNumber)
        * @warning You have to implement threadCallback
	*/
	void executeThreads(const int inThreadsNumber = DefaultThreadsNumber){
		threadCallback(0, 1);
	}

	/**
	* This function executed by each thread after executeThreads has been called
	* @param inThreadId the current thread index (equal 0 for no thread)
	* @param inThreadNumbers the number of threads started (equal 1 for no thread)
        * @warning Must be impleted in the derived class
	*/
	virtual void threadCallback(const int inThreadId, const int inThreadNumbers){
		FDEBUG( FDebug::Controller.writeFromLine("[W] You called executeThreads() but did not implement threadCallback", __LINE__, __FILE__); )
	};

	/** Useless Virtual Destructor */
	virtual ~FNoThreaded(){}

protected:
	/**
	* This function lock a thread-type spefic mutex - here it does nothing
	*/
	void lock() const {}

	/**
	* This function unlock a thread-type spefic mutex - here it does nothing
	*/
	void unlock() const {}

	/**
	* Barrier to sync all thread - here it does nothing
	*/
	void barrier() const {}
};

#endif //FNOTHREADED_HPP

// [--LICENSE--]
