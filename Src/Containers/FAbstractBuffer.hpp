#ifndef FABSTRACTBUFFER_HPP
#define FABSTRACTBUFFER_HPP

class FAbstractBufferReader {
public:
    virtual ~FAbstractBufferReader(){
    }

    virtual char* data() = 0;
    virtual const char* data() const  = 0;
    virtual int getSize() const = 0;
    virtual void seek(const int inIndex) = 0;
    virtual int tell() const  = 0;
    virtual void reserve(const int nbBytes) = 0;
    virtual void reset() = 0;
    template <class ClassType>
    ClassType getValue(){
        static_assert(sizeof(ClassType) == 0, "Your Buffer should implement getValue.");
        return ClassType();
    }
    template <class ClassType>
    void fillValue(ClassType* const){
        static_assert(sizeof(ClassType) == 0, "Your Buffer should implement fillValue.");
    }
    template <class ClassType>
    void fillArray(ClassType* const , const int ){
        static_assert(sizeof(ClassType) == 0, "Your Buffer should implement fillArray.");
    }
    template <class ClassType>
    FAbstractBufferReader& operator>>(ClassType& ){
        static_assert(sizeof(ClassType) == 0, "Your Buffer should implement operator>>.");
        return *this;
    }
};



class FAbstractBufferWriter {
public:
    virtual ~FAbstractBufferWriter(){
    }

    virtual char* data() = 0;
    virtual const char* data()  const = 0;
    virtual int getSize() const = 0;
    virtual void reset() = 0;

    template <class ClassType>
    void write(const ClassType& object){
        static_assert(sizeof(ClassType) == 0, "Your Buffer should implement write.");
    }
    template <class ClassType>
    void writeAt(const int position, const ClassType& object){
        static_assert(sizeof(ClassType) == 0, "Your Buffer should implement writeAt.");
    }
    template <class ClassType>
    void write(const ClassType* const objects, const int inSize){
        static_assert(sizeof(ClassType) == 0, "Your Buffer should implement write.");
    }
    template <class ClassType>
    FAbstractBufferWriter& operator<<(const ClassType& ){
        static_assert(sizeof(ClassType) == 0, "Your Buffer should implement operator<<.");
        return *this;
    }
};

#endif // FABSTRACTBUFFER_HPP
