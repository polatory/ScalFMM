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


