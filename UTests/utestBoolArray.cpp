// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================
#include "FUTester.hpp"

#include "Containers/FBoolArray.hpp"


/**
* This file is a unit test for the FBoolArray class
*/


/** this class test the bool array container */
class TestArray : public FUTester<TestArray> {

    void TestGetSet(){
        FBoolArray array(500);
        for(int idx = 0 ; idx < 500 ; ++idx){
            uassert(!array.get(idx));
        }

        for(int idx = 0 ; idx < 500 ; ++idx){
            array.set(idx, true);
            uassert(array.get(idx));
            array.set(idx, false);
            uassert(!array.get(idx));
        }

        for(int idx = 0 ; idx < 500 ; ++idx){
            array.set(idx, true);
        }
        array.setToZeros();
        for(int idx = 0 ; idx < 500 ; ++idx){
            uassert(!array.get(idx));
        }
        array.setToOnes();
        for(int idx = 0 ; idx < 500 ; ++idx){
            uassert(array.get(idx));
        }
    }

    void TestGetSet2(){
        FBoolArray array(100);

        for(int idx = 0 ; idx < 100 ; ++idx){
            if(idx%3){
                array.set(idx, true);
                uassert(array.get(idx));
            }
            else{
                uassert(!array.get(idx));
            }
        }
    }

    void TestEqual(){
        FBoolArray array1(10);
        FBoolArray array2(10);


        uassert(array1 == array2);

        array1.set(1, true);
        uassert(array1 != array2);

        array2.set(1, true);
        uassert(array1 == array2);

        array1.set(5, true);
        array2 = array1;
        uassert(array1 == array2);
    }

    // set test
    void SetTests(){
        AddTest(&TestArray::TestGetSet,"Test Get & Set");
        AddTest(&TestArray::TestGetSet2,"Test Get & Set 2");
        AddTest(&TestArray::TestEqual,"Test Equal");
    }
};

// You must do this
TestClass(TestArray)


