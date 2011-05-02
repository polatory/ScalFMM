#include "FUTester.hpp"

#include "../Src/Containers/FBoolArray.hpp"


/**
* This file is a unit test for the FBoolArray class
*/


/** this class test the bool array container */
class TestArray : public FUTester<TestArray> {

        void TestGetSet(){
                FBoolArray array(10);
                for(int idx = 0 ; idx < 10 ; ++idx){
                    assert(!array.get(idx));
                }

                for(int idx = 0 ; idx < 10 ; ++idx){
                    array.set(idx, true);
                    assert(array.get(idx));
                }
	}

        void TestGetSet2(){
                FBoolArray array(10);

                for(int idx = 0 ; idx < 10 ; ++idx){
                    if(idx%3){
                        array.set(idx, true);
                        assert(array.get(idx));
                    }
                    else{
                        assert(!array.get(idx));
                    }
                }
        }

        void TestEqual(){
                FBoolArray array1(10);
                FBoolArray array2(10);


                assert(array1 == array2);

                array1.set(1, true);
                assert(array1 != array2);

                array2.set(1, true);
                assert(array1 == array2);

                array1.set(5, true);
                array2 = array1;
                assert(array1 == array2);
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


