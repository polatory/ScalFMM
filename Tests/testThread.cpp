// /!\ Please, you must read the license at the bottom of this page

#include "../Sources/Utils/FAbstractThread.hpp"
#include "../Sources/Utils/FOpenMPThread.hpp"
#include "../Sources/Utils/FPosixThread.hpp"
#include "../Sources/Utils/FNoThread.hpp"

#include <stdio.h>

// Compile by g++ testThread.cpp ../Sources/Utils/FDebug.cpp -lgomp -lpthread -fopenmp -o testThread.exe

/**
* In this file we show how to use the thread module
*/


/**
* TOpen is an example of the FOpenMPThreaded implementation class
*/
class TOpen : public FOpenMPThread<TOpen>{
public:
	void threadCallback(const int inThreadId, const int inThreadsNum){
		printf("I am %d on %d, \n",inThreadId, inThreadsNum);
	}

	void threadCallbackMutex(const int inThreadId, const int inThreadsNum){
		lock();
		for(long idx = 0 ; idx < 50000000 ; ++idx) {++idx;--idx;}
		printf("I am %d on %d, \n",inThreadId, inThreadsNum);
		unlock();

		barrier();
		printf("I am %d ok\n",inThreadId);
	}
};

/**
* TPosix is an example of the FPosixThreaded implementation class
*/
class TPosix : public FPosixThread<TPosix>{
public:
	void threadCallback(const int inThreadId, const int inThreadsNum){
		printf("I am %d on %d, \n",inThreadId, inThreadsNum);
	}

	void threadCallbackMutex(const int inThreadId, const int inThreadsNum){
		lock();
		for(long idx = 0 ; idx < 50000000 ; ++idx) {++idx;--idx;}
		printf("I am %d on %d, \n",inThreadId, inThreadsNum);
		unlock();

		barrier();
		printf("I am %d ok\n",inThreadId);
	}
};


int main(void){
	// create openmp thread derived class
	TOpen open;
	open.executeThreads(10);
	open.executeThreads(&TOpen::threadCallbackMutex,5);

	// create posix thread derived class
	TPosix posix;
	posix.executeThreads(10);
	posix.executeThreads(&TPosix::threadCallbackMutex,5);

	return 0;
}


// [--LICENSE--]
