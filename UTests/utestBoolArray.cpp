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

#include "../Src/Containers/FBoolArray.hpp"


/**
* This file is a unit test for the FBoolArray class
*/


/** this class test the bool array container */
class TestArray : public FUTester<TestArray> {

        void TestGetSet(){
                FBoolArray array(500);
                for(int idx = 0 ; idx < 500 ; ++idx){
                    assert(!array.get(idx));
                }

                for(int idx = 0 ; idx < 500 ; ++idx){
                    array.set(idx, true);
                    assert(array.get(idx));
                    array.set(idx, false);
                    assert(!array.get(idx));
                }

                for(int idx = 0 ; idx < 500 ; ++idx){
                    array.set(idx, true);
                }
                array.setToZeros();
                for(int idx = 0 ; idx < 500 ; ++idx){
                    assert(!array.get(idx));
                }
	}

        void TestGetSet2(){
                FBoolArray array(100);

                for(int idx = 0 ; idx < 100 ; ++idx){
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


