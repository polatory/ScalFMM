#include "FUTester.hpp"

#include "../Src/Containers/FVector.hpp"


/**
* This file is a unit test for the FVector class
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


/** this class test the vector container */
class TestVector : public FUTester<TestVector> {
	// Called before each test : simply set counter to 0
	void PreTest(){
		TestObject::counter = 0;
		TestObject::dealloced = 0;
	}

	// test size
	void TestSize(){
                FVector<TestObject> vector;
                vector.push(TestObject());
                vector.push(TestObject());
                vector.push(TestObject());
                assert(vector.getSize() == 3);
		
                assert((TestObject::counter - TestObject::dealloced) == vector.getSize());

                vector.clear();
                assert(vector.getSize() == 0);

		assert(TestObject::counter == TestObject::dealloced);
	}
	
	// test copy
	void TestCopy(){
                FVector<TestObject> vector;
                vector.push(TestObject());
                vector.push(TestObject());
                vector.push(TestObject());

                {
                    FVector<TestObject> vector2(vector);
                    assert(vector.getSize() == vector2.getSize());
                    assert((TestObject::counter - TestObject::dealloced) == (vector.getSize() + vector2.getSize()));
                }
                {
                    FVector<TestObject> vector2(vector.getSize()/2);
                    vector2 = vector;
                    assert(vector.getSize() == vector2.getSize());
                    assert((TestObject::counter - TestObject::dealloced) == (vector.getSize() + vector2.getSize()));
                }
	}

	// test iter
	void TestIter(){		
                FVector<TestObject> vector;
		{
                        FVector<TestObject>::ConstBasicIterator iter(vector);
			assert(!iter.hasNotFinished());
		}
		{
                        vector.push(TestObject());
                        vector.push(TestObject());
                        vector.push(TestObject());

                        FVector<TestObject>::ConstBasicIterator iter(vector);
			assert(iter.hasNotFinished());

			int counter = 0;
			while(iter.hasNotFinished()){ iter.gotoNext(); ++counter; }
			assert(!iter.hasNotFinished());
                        assert(counter == vector.getSize());
		}
	}

        // test remove
        void TestRemove(){
                FVector<TestObject> vector;
                {
                        FVector<TestObject>::BasicIterator iter(vector);
                        assert(!iter.hasNotFinished());
                }
                {
                        vector.push(TestObject());
                        vector.push(TestObject());
                        vector.push(TestObject());

                        FVector<TestObject>::BasicIterator iter(vector);
                        assert(iter.hasNotFinished());

                        iter.gotoNext();
                        iter.remove();
                        assert(iter.hasNotFinished());

                        iter.gotoNext();
                        assert(!iter.hasNotFinished());

                        assert(2 == vector.getSize());
                }
        }
		
	// set test
	void SetTests(){
            AddTest(&TestVector::TestSize,"Test Size");
            AddTest(&TestVector::TestCopy,"Test Copy");
            AddTest(&TestVector::TestIter,"Test Iter");
            AddTest(&TestVector::TestRemove,"Test Remove");
	}
};

// You must do this
TestClass(TestVector)


