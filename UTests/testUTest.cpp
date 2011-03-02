#include "FUTester.hpp"

// compile by g++ testUTest.cpp -o test.exe

/** this class show a simple example of unit test */
class MyTest : public FUTester<MyTest> {
	void Before(){
		Print("Before running the test");
	}

	void TestOne(){
		assert(true);
		assert(false);
		assert(1 == 1);
	}
	
	void TestTwo(){
		equal(1 , 1);
		different(1 , 1);
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


