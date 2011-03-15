#ifndef FASSERTABLE_HPP
#define FASSERTABLE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <sstream>
#include <iostream>

class FAbstractApplication;

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAssertable
* Please read the license
*
* This class is an interface for managing error.
*
* Please refere to testAssert.cpp to see an example
* <code>
* </code>
*/
class FAssertable {
private:	
	static FAbstractApplication* CurrentApp;	//< You must have only one app

	/** 
	* Called by Fapplication instance to set the current app
	* @param inCurrentApp current app
	*/
        static void SetCurrentApp(FAbstractApplication* const inCurrentApp){
            CurrentApp = inCurrentApp;
	}

	/** To set CurrentApp */
	friend class FAbstractApplication;

	/** To quit current application */
	void exitApplication(const int inExitCode) const;

protected:
	/** Empty Destructor */
	virtual ~FAssertable(){}

	/**
	* to write debug data with line & file
	* @param inTest if false, application will stop
	* @param inMessage a message - from any type - to print
	* @param inLinePosition line number
	* @param inFilePosition file name
	*
        * <code> assert(toto == titi, "toto is not equal titi!", __LINE__, __FILE__); </code>
	*
	* To prevent use from multiple thread we use a ostringstream before printing
	*/
	template <class Tmess, class Tline, class Tfile>
	void assert(const bool inTest, const Tmess& inMessage, const Tline& inLinePosition, const Tfile& inFilePosition, const int inExitCode = 1) const {
		if(!inTest){
			calledBeforeExit();

			std::ostringstream oss;
			oss << "Error in " << inFilePosition << " at line " << inLinePosition <<" :\n";
			oss << inMessage << "\n";
		
			std::cerr << oss.str();
			exitApplication(inExitCode);
		}
	}

	/**
	* May be implemented in the derived class to know when app will quit
	*/
	virtual void calledBeforeExit() const {}

};

#endif //FASSERTABLE_HPP

// [--LICENSE--]
