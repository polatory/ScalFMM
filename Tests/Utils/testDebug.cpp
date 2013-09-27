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
#include "../../Src/Utils/FDebug.hpp"

/**
* In this file we show how to use the debug module.
* Warning, in FGlobal.hpp (included in FDebug.hpp) SCALFMM_USE_DEBUG might be undefined.
*/

int main(void){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> to understand the FDebug system.\n";
    //////////////////////////////////////////////////////////////
	// Print data simply
	FLOG( FDebug::Controller << "Hello Wordl\n");

	// Print a variable (formated print)
	int i = 50;
	FLOG( FDebug::Controller.writeVariableFromLine( "i", i, __LINE__, __FILE__););

	// Write a developer information
	FLOG( FDebug::Controller.writeFromLine("Strange things are there!", __LINE__, __FILE__); )

	// Flush
	FLOG( FDebug::Controller << FDebug::Flush );

	// Change stream type
	FLOG( FDebug::Controller.writeToFile("testDebug.out.temp"); )
	FLOG( FDebug::Controller << "Hello Wordl 2 the return\n");

	return 0;
}



