#include "FUTester.hpp"

#include "../Sources/Containers/FList.hpp"

// compile by g++ flistUTest.cpp -o testFile.exe

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
		list.pushFront(TestObject());
		list.pushFront(TestObject());
		list.pushFront(TestObject());
		assert(list.getSize() == 3);
		
		assert((TestObject::counter - TestObject::dealloced) == list.getSize());

		list.clear();
		assert(list.getSize() == 0);

		assert(TestObject::counter == TestObject::dealloced);
	}
	
	// test copy
	void TestCopy(){
		FList<TestObject> list;
		list.pushFront(TestObject());
		list.pushFront(TestObject());
		list.pushFront(TestObject());

		FList<TestObject> list2 = list;
		assert(list.getSize() == list2.getSize());
		
		assert((TestObject::counter - TestObject::dealloced) == (list.getSize() + list2.getSize()));
	}

	// test iter
	void TestIter(){		
		FList<TestObject> list;
		{
			FList<TestObject>::BasicIterator iter(list);
			assert(!iter.isValide());
		}
		{
			list.pushFront(TestObject());
			list.pushFront(TestObject());
			list.pushFront(TestObject());

			FList<TestObject>::BasicIterator iter(list);
			assert(iter.isValide());

			int counter = 0;
			while(iter.isValide()){ iter.progress(); ++counter; }
			assert(!iter.isValide());
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


