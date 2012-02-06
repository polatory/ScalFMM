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

// compile by g++ utestTest.cpp -o utestTest.exe

/** this class show a simple example of unit test */
class MyTest : public FUTester<MyTest> {
	void Before(){
		Print("Before running the test");
	}

	void TestOne(){
		assert(true);
                //or assert(false); make an error
		assert(1 == 1);
	}
	
	void TestTwo(){
		equal(1 , 1);
                different(1 , 2);
	}
	
	void After(){
		Print("After running the test");
	}
	
	void PreTest(){
		Print("Before each test");
	}
	
	void PostTest(){
		Print("After each test");
	}
	
	// You must implement it
	void SetTests(){
            AddTest(&MyTest::TestOne);
            AddTest(&MyTest::TestTwo,"My Second Test");
	}
};

// You must do this
TestClass(MyTest)


