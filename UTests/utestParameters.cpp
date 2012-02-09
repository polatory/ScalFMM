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

#include "../Src/Utils/FParameters.hpp"

// compile by g++ utestMorton.cpp -o utestMorton.exe

/**
* This file is a unit test for the FTreeCoordinate class
*/


/** this class test the list container */
class TestParameters : public FUTester<TestParameters> {
        void Lower(){
            assert(FParameters::toLower('A') == 'a');
            assert(FParameters::toLower('Z') == 'z');
            assert(FParameters::toLower('a') == 'a');
            assert(FParameters::toLower('z') == 'z');
            assert(FParameters::toLower('m') == 'm');
            assert(FParameters::toLower(';') == ';');
        }

        void CharsEquals(){
            assert(FParameters::areCharsEquals('a','A',false));
            assert(FParameters::areCharsEquals('A','A',false));
            assert(FParameters::areCharsEquals('A','a',false));
            assert(FParameters::areCharsEquals('a','a',false));
            assert(!FParameters::areCharsEquals('a','l',false));
            assert(!FParameters::areCharsEquals('a','M',false));
            assert(!FParameters::areCharsEquals('1',';',false));


            assert(!FParameters::areCharsEquals('a','A',true));
            assert(FParameters::areCharsEquals('A','A',true));
            assert(!FParameters::areCharsEquals('A','a',true));
            assert(FParameters::areCharsEquals('a','a',true));
            assert(FParameters::areCharsEquals(';',';',true));
        }

        void StrsEquals(){
            assert(FParameters::areStrEquals("","",false));
            assert(FParameters::areStrEquals("aa","Aa",false));
            assert(!FParameters::areStrEquals("az","azz",false));
            assert(FParameters::areStrEquals("aza","azA",false));

            assert(FParameters::areStrEquals("a","a",true));
            assert(FParameters::areStrEquals("az","az",true));
            assert(!FParameters::areStrEquals("aza","azA",true));
            assert(!FParameters::areStrEquals("az","azA",true));
        }

        void FindParameter(){
            const int argc = 4;
            const char* argv[4] = {
            "toto",
            "titi",
            "TaTa",
            "tu.tU"
            };

            assert(FParameters::findParameter(argc,argv,"tOto",false) == 0);
            assert(FParameters::findParameter(argc,argv,"titi",false) == 1);
            assert(FParameters::findParameter(argc,argv,"tAtA",false) == 2);
            assert(FParameters::findParameter(argc,argv,"tu.tu",false) == 3);

            for(int idx = 0 ; idx < argc ; ++idx){
                assert(FParameters::findParameter(argc,argv,argv[idx],true) == idx);
            }

            assert(FParameters::findParameter(argc,argv,"tu.tu",true) == -1);
            assert(FParameters::findParameter(argc,argv,"TOTO",true) == -1);
        }

        void StrToOther(){
            assert(FParameters::StrToOther<int>("0") == 0);
            assert(FParameters::StrToOther<int>("a",1) == 1);
        }

        // set test
        void SetTests(){
            AddTest(&TestParameters::Lower,"Test Lower");
            AddTest(&TestParameters::CharsEquals,"Test Char equality");
            AddTest(&TestParameters::StrsEquals,"Test Str equality");
            AddTest(&TestParameters::FindParameter,"Test Str equality");
            AddTest(&TestParameters::StrToOther,"Test Convert");
        }
};

// You must do this
TestClass(TestParameters)



