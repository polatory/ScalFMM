
// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#include "FUTester.hpp"

#include "../Src/Containers/FBlockAllocator.hpp"

#include <cstring>

/**
* This file is a unit test for the List block allocator
*/

/**
* This class is simply used to count alloc dealloc
*/
class TestObject{
public:
    static int counter;
    static int dealloced;

    int array[50];

    TestObject(){
        memset(array, 0, 50 * sizeof(int));
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
class TestBlock : public FUTester<TestBlock> {
    // Called before each test : simply set counter to 0
    void PreTest(){
        TestObject::counter = 0;
        TestObject::dealloced = 0;
    }

    // test copy
    void TestBlockFunction(){
        FListBlockAllocator<TestObject, 10> alloc;

        TestObject* ptr[100];
        for(int idx = 0 ; idx < 100 ; ++idx){
            TestObject* dl1 = alloc.newObject();
            TestObject* dl2 = alloc.newObject();
            alloc.deleteObject(dl1);
            ptr[idx] = alloc.newObject();
            alloc.deleteObject(dl2);
        }

        for(int idx = 0 ; idx < 100 ; ++idx){
            for(int idxval = 0 ; idxval < 50 ; ++idxval){
                ptr[idx]->array[idxval] += (idxval * idx);
            }
        }

        for(int idx = 0 ; idx < 100 ; ++idx){
            for(int idxval = 0 ; idxval < 50 ; ++idxval){
                uassert(ptr[idx]->array[idxval] == (idxval * idx));
            }
        }

        for(int idx = 0 ; idx < 100 ; ++idx){
            alloc.deleteObject(ptr[idx]);
        }

        uassert(TestObject::counter == 300);
        uassert(TestObject::counter == TestObject::dealloced);
    }

    // set test
    void SetTests(){
            AddTest(&TestBlock::TestBlockFunction,"Test Allocate Deallocate");
    }
};

// You must do this
TestClass(TestBlock)


