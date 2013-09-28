#ifndef FBLOCKALLOCATOR_HPP
#define FBLOCKALLOCATOR_HPP

#include <list>
#include <cstring>

/**
 * What a cell allocator should implement
 */
template <class ObjectClass>
class FAbstractBlockAllocator{
public:
    FAbstractBlockAllocator(){}

    virtual ~FAbstractBlockAllocator(){
    }

    FAbstractBlockAllocator(const FAbstractBlockAllocator&) = delete;
    FAbstractBlockAllocator& operator=(const FAbstractBlockAllocator&) = delete;

    virtual ObjectClass* newObject() = 0;
    virtual void deleteObject(const ObjectClass*) = 0;
};

/**
 * Redirection to normal operators
 */
template <class ObjectClass>
class FBasicBlockAllocator : public FAbstractBlockAllocator<ObjectClass> {
public:
    ObjectClass* newObject(){
        return new ObjectClass;
    }

    void deleteObject(const ObjectClass* inObject){
        delete inObject;
    }
};


/**
 * Allocation per blocks
 */
template <class ObjectClass, int SizeOfBlock >
class FListBlockAllocator : public FAbstractBlockAllocator<ObjectClass>{
    class CellsBlock {
        int nbObjectInBlock;
        unsigned char usedBlocks[SizeOfBlock];
        unsigned char objectMemory[SizeOfBlock * sizeof(ObjectClass)];

    public:
        CellsBlock() : nbObjectInBlock(0){
            memset(usedBlocks, 0, sizeof(unsigned char) * SizeOfBlock);
        }

        bool isInsideBlock(const ObjectClass*const cellPtr) const{
            const unsigned char*const inPtr = reinterpret_cast<const unsigned char*>( cellPtr );
            return objectMemory <= inPtr && inPtr < (objectMemory+ SizeOfBlock * sizeof(ObjectClass));
        }

        int getPositionInBlock(const ObjectClass*const cellPtr) const{
            const unsigned char*const inPtr = reinterpret_cast<const unsigned char*>( cellPtr );
            return int((inPtr - objectMemory) / sizeof(ObjectClass));
        }

        void deleteObject(const int position){
            reinterpret_cast<ObjectClass*>(objectMemory)[position].~ObjectClass();
            nbObjectInBlock -= 1;
            usedBlocks[position] = 0x0;
        }

        ObjectClass* getNewObject(){
            for(int idx = 0 ; idx < SizeOfBlock ; ++idx){
                if(usedBlocks[idx] == 0){
                    nbObjectInBlock += 1;
                    usedBlocks[idx] = 0x1;
                    new (&reinterpret_cast<ObjectClass*>(objectMemory)[idx]) ObjectClass;
                    return &reinterpret_cast<ObjectClass*>(objectMemory)[idx];
                }
            }

            return 0;
        }

        int isEmpty() const{
            return nbObjectInBlock == 0;
        }

        int isFull() const{
            return nbObjectInBlock == SizeOfBlock;
        }
    };

    std::list<CellsBlock> blocks;

public:
    ObjectClass* newObject(){
        typename std::list<CellsBlock>::iterator iterator(blocks.begin());
        const typename std::list<CellsBlock>::iterator end(blocks.end());

        while( iterator != end){
            if( !(*iterator).isFull() ){
                return (*iterator).getNewObject();
            }
            ++iterator;
        }

        blocks.emplace_back();
        return blocks.back().getNewObject();
    }

    void deleteObject(const ObjectClass* inObject){
        typename std::list<CellsBlock>::iterator iterator(blocks.begin());
        const  typename std::list<CellsBlock>::iterator end(blocks.end());

        while( iterator != end ){
            if( (*iterator).isInsideBlock(inObject) ){
                (*iterator).deleteObject((*iterator).getPositionInBlock(inObject));
                if( (*iterator).isEmpty() ){
                    blocks.erase(iterator);
                }
                break;
            }
            ++iterator;
        }
    }
};

#endif // FBLOCKALLOCATOR_HPP
