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
#ifndef UTESTER_HPP
#define UTESTER_HPP


#include <iostream>
#include <list>
#include <string>
#include <cstdio>


#define TestClass(X)\
int main(void){\
	X Controller;\
        return Controller.Run();\
}\


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* Please read the license
*
* This class is used to make simple unit test cases
*
* Please refer to testUTest.cpp to see an example
* @warning Create a derived class that implement SetTests() and use TestClass() macro
*/
template <class TestClass>
class FUTester{
	// Test function pointer
	typedef void (TestClass::*TestFunc)(void);
	
	/** Test descriptor */
	struct TestFuncDescriptor{
		TestFunc func;		//< Test adress
		std::string name;	//< Test name
	};
	

	std::list<TestFuncDescriptor> tests;	//< all tests
	
	int totalTests;				//< number of tests
	
	int currentTest;			//< current processing test in the run
	int currentStep;			//< current processing step in the run

	int failledSteps;			//< number of failled step in the current test
	int failledTests;			//< number of failled tests
	
protected:
	/** Constructor */
	FUTester(){
		totalTests = 0;
	}

	/** Callback before processing test */
	virtual void Before(){}
	
	/** Callback after processing test */
	virtual void After(){}
	
	/** Callback before each unit test */
	virtual void PreTest(){}
	
	/** Callback after each unit test */
	virtual void PostTest(){}
	
	/**
	* This function has to add tests
        * <code> AddTest(&MyTest::TestOne); </code>
	*/
        virtual void SetTests() = 0;

	/**
	* Add a test without giving a name
	* @param inFunc test function address
	*/
        void AddTest(TestFunc inFunc){
		char buff[256];
		sprintf(buff,"Unamed Test number %d",totalTests+1);
		AddTest(inFunc,buff);
	}
	
	/**
	* Add a test with a name
	* @param inFunc test function address
	* @param inFuncName function name
	*/
        void AddTest(TestFunc inFunc, const std::string& inFuncName){
		++totalTests;
		TestFuncDescriptor desc;
		desc.func = inFunc;
                desc.name = inFuncName;
		tests.push_back(desc);
	}
	
	/**
	* To print a message manually in the test
	* @param value a object that ostream can work on
	*/
        template <class Output>
	void Print(const Output& value){
		std::cout<< "--- Output from program : " << value << "\n";
	}
	
	/**
	* To test
	* @param result the test result
	* if result is false test failled
	*/
        void assert(const bool result){
		++currentStep;
		if(!result){
			std::cout << ">> Step " << currentStep << " Failled\n";
			++failledSteps;
		}
	}
	
	/**
	* To test equality
	* @param v1 value one
	* @param v2 value 2
	* if v1 is not equal v2 test failled
	*/
	template <class T>
        void equal(const T& v1, const T& v2){
                assert(v1 == v2);
	}
	
	/**
	* To test equality
	* @param v1 value one
	* @param v2 value 2
	* if v1 is equal v2 test failled
	*/
	template <class T>
        void different(const T& v1, const T& v2){
                assert(v1 != v2);
	}
	
public :
	/**
	* Processing the test
        * return application exit code (= nb of errors)
	*/
        int Run(){
		tests.clear();
		// register tests
		SetTests();
	
		TestClass* const toTest = static_cast<TestClass*>(this);
		currentTest = 0;
		failledTests = 0;

		Before();
		
		// for each tests
                const typename std::list<TestFuncDescriptor>::const_iterator end = tests.end();
                for(typename std::list<TestFuncDescriptor>::iterator iter = tests.begin() ; iter != end ; ++iter){
			currentStep = 0;
			failledSteps = 0;
			
			std::cout << "[Start] " << (*iter).name << "\n";
			
			PreTest();
			TestFunc ff = (*iter).func;
			(toTest->*ff)();
			PostTest();
			
			if(failledSteps){
				std::cout << "[Finished] FAILLED (" << failledSteps << "/" << currentStep<< " steps failled)\n";
				++failledTests;
			}
                        else{
				std::cout << "[Finished] PASSED (" << currentStep << " steps)\n";
			}
			
			++currentTest;
		}
		
		
		After();

		std::cout <<"Test is over, " << (totalTests-failledTests) << " Passed, " << failledTests << " Failled\n";

                return failledTests;
	}

};





#endif
