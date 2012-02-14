#ifndef FSMARTPOINTER_HPP
#define FSMARTPOINTER_HPP

template <class ClassType>
class FSmartPointer {
    ClassType* pointer;
    int* counter;

public:
    typedef ClassType ValueType;

    FSmartPointer() : pointer(0), counter(0) {
    }

    FSmartPointer(ClassType* const inPointer) : pointer(0), counter(0) {
        assign(inPointer);
    }

    FSmartPointer(const FSmartPointer& inPointer) : pointer(0), counter(0) {
        assign(inPointer);
    }

    ~FSmartPointer(){
        release();
    }

    void operator=(ClassType* const inPointer){
        assign(inPointer);
    }

    void operator=(const FSmartPointer& inPointer){
        assign(inPointer);
    }

    void assign(ClassType* const inPointer){
        release();
        pointer = inPointer;
        counter = new int;
        (*counter) = 1;
    }

    void assign(const FSmartPointer& inPointer){
        release();
        pointer = inPointer.pointer;
        counter = inPointer.counter;
        if(counter) (*counter) = (*counter) + 1;
    }

    void release(){
        if(counter){
            (*counter) = (*counter) - 1;
            if( (*counter) == 0 ){
                delete pointer;
            }
            pointer = 0;
            counter = 0;
        }
    }

    bool isAssigned() const{
        return pointer != 0;
    }

    bool isLast() const{
        return counter && (*counter) == 1;
    }

    ClassType* getPtr(){
        return pointer;
    }

    const ClassType* getPtr() const {
        return pointer;
    }

    ClassType& operator[](const int& index){
        return pointer[index];
    }

    const ClassType& operator[](const int& index) const {
        return pointer[index];
    }

    ClassType& operator*(){
        return (*pointer);
    }

    const ClassType& operator*() const {
        return (*pointer);
    }

    ClassType* operator->(){
        return pointer;
    }

    const ClassType* operator->() const {
        return pointer;
    }

    operator const ClassType*() const {
        return pointer;
    }

    operator ClassType*() {
        return pointer;
    }
};

#endif // FSMARTPOINTER_HPP
