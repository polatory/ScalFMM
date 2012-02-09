// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FASSERTABLE_HPP
#define FASSERTABLE_HPP


#include <cstdlib>
#include <sstream>
#include <iostream>

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
        * <code> fassert(toto == titi, "problem : toto is not equal titi!", __LINE__, __FILE__); </code>
	*
	* To prevent use from multiple thread we use a ostringstream before printing
	*/
	template <class Tmess, class Tline, class Tfile>
        void fassert(const bool inTest, const Tmess& inMessage, const Tline& inLinePosition, const Tfile& inFilePosition, const int inExitCode = 1) const {
		if(!inTest){
                        std::ostringstream oss;
			oss << "Error in " << inFilePosition << " at line " << inLinePosition <<" :\n";
			oss << inMessage << "\n";
		
			std::cerr << oss.str();
                        exit(inExitCode);
		}
	}

};

#endif //FASSERTABLE_HPP


