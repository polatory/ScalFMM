// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
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
	* @param inExitCode an exit code
	*
	* @code
	* fassert(toto == titi, "problem : toto is not equal titi!", __LINE__, __FILE__);
	* @endcode
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


