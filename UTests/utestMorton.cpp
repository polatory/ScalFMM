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
#include "FUTester.hpp"

#include "../Src/Containers/FTreeCoordinate.hpp"

// compile by g++ utestMorton.cpp -o utestMorton.exe

/**
* This file is a unit test for the FTreeCoordinate class
*/


/** this class test the list container */
class TestMorton : public FUTester<TestMorton> {
        void Morton(){
            {
                FTreeCoordinate pos(5,1,7);
                FTreeCoordinate cp;
                cp.setPositionFromMorton(pos.getMortonIndex(10),10);
                assert(pos == cp);
                assert(cp.getMortonIndex(10) == pos.getMortonIndex(10));
            }
            {
                FTreeCoordinate pos(2,8,3);
                FTreeCoordinate cp;
                cp.setPositionFromMorton(pos.getMortonIndex(10),10);
                assert(pos == cp);
                assert(cp.getMortonIndex(10) == pos.getMortonIndex(10));
            }
            {
                FTreeCoordinate pos(51,11,47);
                FTreeCoordinate cp;
                cp.setPositionFromMorton(pos.getMortonIndex(10),10);
                assert(pos == cp);
                assert(cp.getMortonIndex(10) == pos.getMortonIndex(10));
            }
	}

        void Position(){
            {
                FTreeCoordinate pos(0,0,0);
                assert(pos.getMortonIndex(1) == 0);
            }
            {
                FTreeCoordinate pos(1,1,1);
                assert(pos.getMortonIndex(1) == 7);
            }
            {
                FTreeCoordinate pos(0,1,1);
                assert(pos.getMortonIndex(1) == 3);
            }
            {
                FTreeCoordinate pos(2,2,2);
                assert(pos.getMortonIndex(2) == (7 << 3) );
            }
            {
                FTreeCoordinate pos(1,2,4);
                assert(pos.getMortonIndex(3) == 84 );// 001 010 100 =>> 001010100 => 84d
            }
	}
		
	// set test
	void SetTests(){
            AddTest(&TestMorton::Morton,"Test Morton");
            AddTest(&TestMorton::Position,"Test Position");
	}
};

// You must do this
TestClass(TestMorton)


