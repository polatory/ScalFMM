#ifndef FABSTRACTTHREADED_HPP
#define FABSTRACTTHREADED_HPP
// /!\ Please, you must read the license at the bottom of this page

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractThreaded
* Please read the license
*
* This class is an interface for threaded class.
* Each class that wants to use thread must inherited from FOpenMPThreaded or FPosixThreaded.
*
* Please refere to testThread.cpp to see an example.
* <code>
* // Example with FOpenMPThreaded <br>
* class TOpen : public FOpenMPThreaded<TOpen>{ <br>
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
template <class Derived>
class FAbstractThreaded{
public:
	/**
	* This function is used to create inThreadsNumber threads with inCallback as the callback.
	* @param inCallback the callback (must be a object method)
	* @param inThreadsNumber the number of threads to create
	*/
	virtual void executeThreads(void (Derived::*inCallback)(const int,const int), const int inThreadsNumber) = 0;

	/**
	* This function is used to create inThreadsNumber threads with Object::threadCallback as the callback.
	* @param inThreadsNumber the number of threads to create (default is DefaultThreadsNumber)
        * @warning You have to implement threadCallback
	*/
	virtual void executeThreads(const int inThreadsNumber) = 0;

	/**
	* This function executed by each thread after executeThreads has been called
	* @param inThreadId the current thread index
	* @param inThreadNumbers the number of threads started
        * @warning Must be impleted in the derived class
	*/
	virtual void threadCallback(const int inThreadId,const int inThreadNumbers) = 0;

	/** Useless Virtual Destructor */
	virtual ~FAbstractThreaded(){}

protected:
	/**
	* This function lock a thread-type spefic mutex
	*/
	virtual void lock() const = 0;

	/**
	* This function unlock a thread-type spefic mutex
	*/
	virtual void unlock() const = 0;

	/**
	* Barrier to sync all thread
	*/
	virtual void barrier() const = 0;
};

#endif //FABSTRACTTHREADED_HPP

// [--LICENSE--]
