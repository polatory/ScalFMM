// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================
#include "../Src/Utils/FDebug.hpp"

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



