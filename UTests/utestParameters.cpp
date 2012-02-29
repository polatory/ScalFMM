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
            uassert(FParameters::toLower('A') == 'a');
            uassert(FParameters::toLower('Z') == 'z');
            uassert(FParameters::toLower('a') == 'a');
            uassert(FParameters::toLower('z') == 'z');
            uassert(FParameters::toLower('m') == 'm');
            uassert(FParameters::toLower(';') == ';');
        }

        void CharsEquals(){
            uassert(FParameters::areCharsEquals('a','A',false));
            uassert(FParameters::areCharsEquals('A','A',false));
            uassert(FParameters::areCharsEquals('A','a',false));
            uassert(FParameters::areCharsEquals('a','a',false));
            uassert(!FParameters::areCharsEquals('a','l',false));
            uassert(!FParameters::areCharsEquals('a','M',false));
            uassert(!FParameters::areCharsEquals('1',';',false));


            uassert(!FParameters::areCharsEquals('a','A',true));
            uassert(FParameters::areCharsEquals('A','A',true));
            uassert(!FParameters::areCharsEquals('A','a',true));
            uassert(FParameters::areCharsEquals('a','a',true));
            uassert(FParameters::areCharsEquals(';',';',true));
        }

        void StrsEquals(){
            uassert(FParameters::areStrEquals("","",false));
            uassert(FParameters::areStrEquals("aa","Aa",false));
            uassert(!FParameters::areStrEquals("az","azz",false));
            uassert(FParameters::areStrEquals("aza","azA",false));

            uassert(FParameters::areStrEquals("a","a",true));
            uassert(FParameters::areStrEquals("az","az",true));
            uassert(!FParameters::areStrEquals("aza","azA",true));
            uassert(!FParameters::areStrEquals("az","azA",true));
        }

        void FindParameter(){
            const int argc = 4;
            const char* argv[4] = {
            "toto",
            "titi",
            "TaTa",
            "tu.tU"
            };

            uassert(FParameters::findParameter(argc,argv,"tOto",false) == 0);
            uassert(FParameters::findParameter(argc,argv,"titi",false) == 1);
            uassert(FParameters::findParameter(argc,argv,"tAtA",false) == 2);
            uassert(FParameters::findParameter(argc,argv,"tu.tu",false) == 3);

            for(int idx = 0 ; idx < argc ; ++idx){
                uassert(FParameters::findParameter(argc,argv,argv[idx],true) == idx);
            }

            uassert(FParameters::findParameter(argc,argv,"tu.tu",true) == -1);
            uassert(FParameters::findParameter(argc,argv,"TOTO",true) == -1);
        }

        void StrToOther(){
            uassert(FParameters::StrToOther<int>("0") == 0);
            uassert(FParameters::StrToOther<int>("a",1) == 1);
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



