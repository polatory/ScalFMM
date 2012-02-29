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
                uassert(pos == cp);
                uassert(cp.getMortonIndex(10) == pos.getMortonIndex(10));
            }
            {
                FTreeCoordinate pos(2,8,3);
                FTreeCoordinate cp;
                cp.setPositionFromMorton(pos.getMortonIndex(10),10);
                uassert(pos == cp);
                uassert(cp.getMortonIndex(10) == pos.getMortonIndex(10));
            }
            {
                FTreeCoordinate pos(51,11,47);
                FTreeCoordinate cp;
                cp.setPositionFromMorton(pos.getMortonIndex(10),10);
                uassert(pos == cp);
                uassert(cp.getMortonIndex(10) == pos.getMortonIndex(10));
            }
	}

        void Position(){
            {
                FTreeCoordinate pos(0,0,0);
                uassert(pos.getMortonIndex(1) == 0);
            }
            {
                FTreeCoordinate pos(1,1,1);
                uassert(pos.getMortonIndex(1) == 7);
            }
            {
                FTreeCoordinate pos(0,1,1);
                uassert(pos.getMortonIndex(1) == 3);
            }
            {
                FTreeCoordinate pos(2,2,2);
                uassert(pos.getMortonIndex(2) == (7 << 3) );
            }
            {
                FTreeCoordinate pos(1,2,4);
                uassert(pos.getMortonIndex(3) == 84 );// 001 010 100 =>> 001010100 => 84d
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


