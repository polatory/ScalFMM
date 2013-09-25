#ifndef FBLOCKALLOCATOR_HPP
#define FBLOCKALLOCATOR_HPP

#include <list>

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

    virtual ObjectClass* newCell() = 0;
    virtual void deleteCell(const ObjectClass*) = 0;
};

/**
 * Redirection to normal operators
 */
template <class ObjectClass>
class FBasicBlockAllocator : public FAbstractBlockAllocator<ObjectClass> {
public:
    ObjectClass* newCell(){
        return new ObjectClass;
    }

    void deleteCell(const ObjectClass* inCell){
        delete inCell;
    }
};


/**
 * Allocation per blocks
 */
template <class ObjectClass, int SizeOfBlock >
class FListBlockAllocator : public FAbstractBlockAllocator<ObjectClass>{
    class CellsBlock {
        int nbCellsInBlock;
        unsigned char usedBlocks[SizeOfBlock];
        unsigned char cellsMemory[SizeOfBlock * sizeof(ObjectClass)];

    public:
        CellsBlock() : nbCellsInBlock(0){
            memset(usedBlocks, 0, sizeof(unsigned char) * SizeOfBlock);
        }

        bool isInsideBlock(const ObjectClass*const cellPtr) const{
            const unsigned char*const inPtr = reinterpret_cast<const unsigned char*>( cellPtr );
            return cellsMemory <= inPtr && inPtr < (cellsMemory+ SizeOfBlock * sizeof(ObjectClass));
        }

        int getPositionInBlock(const ObjectClass*const cellPtr) const{
            const unsigned char*const inPtr = reinterpret_cast<const unsigned char*>( cellPtr );
            return int((inPtr - cellsMemory) / sizeof(ObjectClass));
        }

        void deleteCell(const int position){
            ((ObjectClass*)&cellsMemory)[position].~ObjectClass();
            nbCellsInBlock -= 1;
            usedBlocks[position] = 0x0;
        }

        ObjectClass* getNewCell(){
            for(int idx = 0 ; idx < SizeOfBlock ; ++idx){
                if(usedBlocks[idx] == 0){
                    nbCellsInBlock += 1;
                    usedBlocks[idx] = 0x1;
                    new (&((ObjectClass*)&cellsMemory)[idx]) ObjectClass;
                    return &((ObjectClass*)&cellsMemory)[idx];
                }
            }

            return 0;
        }

        int isEmpty() const{
            return nbCellsInBlock == 0;
        }

        int isFull() const{
            return nbCellsInBlock == SizeOfBlock;
        }
    };

    std::list<CellsBlock> blocks;

public:
    ObjectClass* newCell(){
        typename std::list<CellsBlock>::iterator iterator(blocks.begin());
        const typename std::list<CellsBlock>::iterator end(blocks.end());

        while( iterator != end){
            if( !(*iterator).isFull() ){
                return (*iterator).getNewCell();
            }
            ++iterator;
        }

        blocks.emplace_back();
        return blocks.back().getNewCell();
    }

    void deleteCell(const ObjectClass* inCell){
        typename std::list<CellsBlock>::iterator iterator(blocks.begin());
        const  typename std::list<CellsBlock>::iterator end(blocks.end());

        while( iterator != end ){
            if( (*iterator).isInsideBlock(inCell) ){
                (*iterator).deleteCell((*iterator).getPositionInBlock(inCell));
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
