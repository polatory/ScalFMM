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
	FDEBUG( FDebug::Controller << "Hello Wordl\n");

	// Print a variable (formated print)
	int i = 50;
	FDEBUG( FDebug::Controller.writeVariableFromLine( "i", i, __LINE__, __FILE__););

	// Write a developer information
	FDEBUG( FDebug::Controller.writeFromLine("Strange things are there!", __LINE__, __FILE__); )

	// Flush
	FDEBUG( FDebug::Controller << FDebug::Flush );

	// Change stream type
	FDEBUG( FDebug::Controller.writeToFile("testDebug.out.temp"); )
	FDEBUG( FDebug::Controller << "Hello Wordl 2 the return\n");

	return 0;
}



