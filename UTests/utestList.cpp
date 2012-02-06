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

#include "../Src/Containers/FList.hpp"

// compile by g++ utestList.cpp -o utestList.exe

/**
* This file is a unit test for the FList class
*/

/**
* This class is simply used to count alloc dealloc
*/
class TestObject{
public:
	static int counter;
	static int dealloced;

	TestObject(){
		++counter;
	}
	TestObject(const TestObject&){
		++counter;
	}
	~TestObject(){
		++dealloced;
	}
};

int TestObject::counter(0);
int TestObject::dealloced(0);


/** this class test the list container */
class TestList : public FUTester<TestList> {
	// Called before each test : simply set counter to 0
	void PreTest(){
		TestObject::counter = 0;
		TestObject::dealloced = 0;
	}

	// test size
	void TestSize(){
		FList<TestObject> list;
		list.push(TestObject());
		list.push(TestObject());
		list.push(TestObject());
		assert(list.getSize() == 3);
		
		assert((TestObject::counter - TestObject::dealloced) == list.getSize());

		list.clear();
		assert(list.getSize() == 0);

		assert(TestObject::counter == TestObject::dealloced);
	}
	
	// test copy
	void TestCopy(){
		FList<TestObject> list;
		list.push(TestObject());
		list.push(TestObject());
		list.push(TestObject());

		FList<TestObject> list2 = list;
		assert(list.getSize() == list2.getSize());
		
		assert((TestObject::counter - TestObject::dealloced) == (list.getSize() + list2.getSize()));
	}

	// test iter
	void TestIter(){		
		FList<TestObject> list;
		{
                        FList<TestObject>::ConstBasicIterator iter(list);
			assert(!iter.hasNotFinished());
		}
		{
			list.push(TestObject());
			list.push(TestObject());
			list.push(TestObject());

                        FList<TestObject>::ConstBasicIterator iter(list);
			assert(iter.hasNotFinished());

			int counter = 0;
                        while(iter.hasNotFinished()){
                            iter.gotoNext();
                            ++counter;
                        }

			assert(!iter.hasNotFinished());
			assert(counter == list.getSize());
		}
	}
		
	// set test
	void SetTests(){
            AddTest(&TestList::TestSize,"Test Size");
            AddTest(&TestList::TestCopy,"Test Copy");
            AddTest(&TestList::TestIter,"Test Iter");
	}
};

// You must do this
TestClass(TestList)


