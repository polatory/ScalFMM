#include "FUTester.hpp"

#include "../Sources/Utils/FGlobal.hpp"
#include "../Sources/Utils/F3DPosition.hpp"
#include "../Sources/Containers/FTreeCoordinate.hpp"
#include "../Sources/Utils/FConvert.hpp"
#include "../Sources/Utils/FMath.hpp"

// compile by g++ utestConvert.cpp -o utestConvert.exe

/**
* This file is a unit test for the FConvert class
*/


/** this class test the list container */
class TestConvert : public FUTester<TestConvert> {
    void ValideFMathEqual(){
        assert(FMath::LookEqual(0.0001,0.0001));
        assert(!FMath::LookEqual(0.0001,0.00011));
        assert(FMath::LookEqual(0.0001,0.0000999999999));
    }

    void ConvertMorton(){
        {
            const FTreeCoordinate pos(0,1,10);
            const F3DPosition posReal = FConvert::MortonToPosition(pos.getMortonIndex(5),5,10.0);
            assert( FMath::LookEqual(posReal.getX(),FReal(5.0)) );
            assert( FMath::LookEqual(posReal.getY(),FReal(15.0)) );
            assert( FMath::LookEqual(posReal.getZ(),FReal(105.0)) );
        }
    }

    // set test
    void SetTests(){
        AddTest(&TestConvert::ValideFMathEqual,"Test Equal");
        AddTest(&TestConvert::ConvertMorton,"Test Convert Morton");
    }
};


TestClass(TestConvert)


