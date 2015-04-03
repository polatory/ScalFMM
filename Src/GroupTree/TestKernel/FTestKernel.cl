/** This file contains the prototype for a kernel in opencl */
// @SCALFMM_PRIVATE


/***************************************************************************/
/***************************************************************************/
/************************CHANGE THINGS HERE*********************************/
/***************************************************************************/

typedef ___FSize___ FSize;
typedef ___FReal___ FReal;
typedef ___FParticleValueClass___ FParticleValueClass;
typedef long long int MortonIndex;

#define FOpenCLGroupOfCellsCellIsEmptyFlag  ((MortonIndex)-1)

#define NbAttributesPerParticle ___NbAttributesPerParticle___
#define NbSymbAttributes ___NbSymbAttributes___

#define FOpenCLGroupOfParticlesMemoryAlignementBytes  ___FP2PDefaultAlignement___
#define FOpenCLGroupOfParticlesMemoryAlignementParticles (FOpenCLGroupOfParticlesMemoryAlignementBytes/sizeof(FReal))
#define FOpenCLGroupOfParticlesLeafIsEmptyFlag ((MortonIndex)-1)

#define NULLPTR (0)

#define DefaultStructAlign ___DefaultStructAlign___

struct FSymboleCellClass {
    MortonIndex mortonIndex;
    int coordinates[3];
} __attribute__ ((aligned (DefaultStructAlign)));

typedef long long FTestCellPODData;
typedef FTestCellPODData FPoleCellClass;
typedef FTestCellPODData FLocalCellClass;

struct FWrappeCell{
    __global struct FSymboleCellClass* symb;
    __global FPoleCellClass* up;
    __global FLocalCellClass* down;        
};


/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

struct OutOfBlockInteraction{
    MortonIndex outIndex;
    MortonIndex insideIndex;
    int outPosition;
} __attribute__ ((aligned (DefaultStructAlign)));

#define Between(inValue, inMin, inMax)  ( (inMin) <= (inValue) && (inValue) < (inMax) )
#define pow2(power)  (1 << (power))
#define Abs(inV) (inV < 0 ? -inV : inV)

int3 GetPositionFromMorton(MortonIndex inIndex, const int inLevel){
    MortonIndex mask = 0x1LL;

    int3 coord;
    coord.x = 0;
    coord.y = 0;
    coord.z = 0;

    for(int indexLevel = 0; indexLevel < inLevel ; ++indexLevel){
        coord.z |= (int)(inIndex & mask);
        inIndex >>= 1;
        coord.y |= (int)(inIndex & mask);
        inIndex >>= 1;
        coord.x |= (int)(inIndex & mask);

        mask <<= 1;
    }

    return coord;
}

MortonIndex GetMortonIndex(const int3 coord, const int inLevel) {
    MortonIndex index = 0x0LL;
    MortonIndex mask = 0x1LL;
    // the ordre is xyz.xyz...
    MortonIndex mx = coord.x << 2;
    MortonIndex my = coord.y << 1;
    MortonIndex mz = coord.z;

    for(int indexLevel = 0; indexLevel < inLevel ; ++indexLevel){
        index |= (mz & mask);
        mask <<= 1;
        index |= (my & mask);
        mask <<= 1;
        index |= (mx & mask);
        mask <<= 1;

        mz <<= 2;
        my <<= 2;
        mx <<= 2;
    }

    return index;
}

int GetNeighborsIndexes(const int3 coord, const int OctreeHeight, MortonIndex indexes[26], int indexInArray[26]) {
    int idxNeig = 0;
    int limite = 1 << (OctreeHeight - 1);
    // We test all cells around
    for(int idxX = -1 ; idxX <= 1 ; ++idxX){
        if(!Between(coord.x + idxX,0, limite)) continue;

        for(int idxY = -1 ; idxY <= 1 ; ++idxY){
            if(!Between(coord.y + idxY,0, limite)) continue;

            for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                if(!Between(coord.z + idxZ,0, limite)) continue;

                // if we are not on the current cell
                if( idxX || idxY || idxZ ){
                    int3 other;

                    other.x = coord.x + idxX;
                    other.y = coord.y + idxY;
                    other.z = coord.z + idxZ;

                    indexes[ idxNeig ] = GetMortonIndex(other, OctreeHeight - 1);
                    indexInArray[ idxNeig ] = ((idxX+1)*3 + (idxY+1)) * 3 + (idxZ+1);
                    ++idxNeig;
                }
            }
        }
    }
    return idxNeig;
}

int GetInteractionNeighbors(const int3 coord, const int inLevel, MortonIndex inNeighbors[189], int inNeighborsPosition[189]) {
    // Then take each child of the parent's neighbors if not in directNeighbors
    // Father coordinate
    int3 parentCell;
    parentCell.x = coord.x>>1;
    parentCell.y = coord.y>>1;
    parentCell.z = coord.z>>1;

    // Limite at parent level number of box (split by 2 by level)
    const int limite = pow2(inLevel-1);

    int idxNeighbors = 0;
    // We test all cells around
    for(int idxX = -1 ; idxX <= 1 ; ++idxX){
        if(!Between(parentCell.x + idxX,0,limite)) continue;

        for(int idxY = -1 ; idxY <= 1 ; ++idxY){
            if(!Between(parentCell.y + idxY,0,limite)) continue;

            for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                if(!Between(parentCell.z + idxZ,0,limite)) continue;

                // if we are not on the current cell
                if( idxX || idxY || idxZ ){
                    int3 otherParent;

                    otherParent.x = parentCell.x + idxX;
                    otherParent.y = parentCell.y + idxY;
                    otherParent.z = parentCell.z + idxZ;

                    const MortonIndex mortonOther = GetMortonIndex(otherParent, inLevel-1);

                    // For each child
                    for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                        const int xdiff  = ((otherParent.x<<1) | ( (idxCousin>>2) & 1)) - coord.x;
                        const int ydiff  = ((otherParent.y<<1) | ( (idxCousin>>1) & 1)) - coord.y;
                        const int zdiff  = ((otherParent.z<<1) | (idxCousin&1)) - coord.z;

                        // Test if it is a direct neighbor
                        if(Abs(xdiff) > 1 || Abs(ydiff) > 1 || Abs(zdiff) > 1){
                            // add to neighbors
                            inNeighborsPosition[idxNeighbors] = ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3;
                            inNeighbors[idxNeighbors++] = (mortonOther << 3) | idxCousin;
                        }
                    }
                }
            }
        }
    }

    return idxNeighbors;
}


void FSetToNullptr343(struct FWrappeCell ptrs[343]){
    int idx;
    for( idx = 0 ; idx < 343 ; ++idx){
        ptrs[idx].symb = NULLPTR;
    }
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/


struct FOpenCLGroupAttachedLeaf {
    //< Nb of particles in the current leaf
    FSize nbParticles;
    //< Pointers to the positions of the particles
    __global FReal* positionsPointers[3];
    //< Pointers to the attributes of the particles
    __global FParticleValueClass* attributes[NbSymbAttributes+NbAttributesPerParticle];
};

struct FOpenCLGroupAttachedLeaf BuildFOpenCLGroupAttachedLeaf(const FSize inNbParticles, __global FReal* inPositionBuffer, const size_t inLeadingPosition,
                       __global FParticleValueClass* inAttributesBuffer, const size_t inLeadingAttributes){
    struct FOpenCLGroupAttachedLeaf leaf;
    leaf.nbParticles = (inNbParticles);
    // Redirect pointers to position
    leaf.positionsPointers[0] = inPositionBuffer;
    leaf.positionsPointers[1] = (__global FReal*)(((__global unsigned char*)inPositionBuffer) + inLeadingPosition);
    leaf.positionsPointers[2] = (__global FReal*)(((__global unsigned char*)inPositionBuffer) + inLeadingPosition*2);

    for(unsigned idxAttribute = 0 ; idxAttribute < NbSymbAttributes ; ++idxAttribute){
        leaf.attributes[idxAttribute] =  (__global FParticleValueClass*)(((__global unsigned char*)inPositionBuffer) + inLeadingPosition*(idxAttribute+3));
    }

    // Redirect pointers to data
    if(inAttributesBuffer){
        for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
            leaf.attributes[idxAttribute+NbSymbAttributes] = (__global FParticleValueClass*)(((__global unsigned char*)inAttributesBuffer) + idxAttribute*inLeadingAttributes);
        }
    }
    else{
        for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
            leaf.attributes[idxAttribute+NbSymbAttributes] = NULLPTR;
        }
    }
    return leaf;
}

struct FOpenCLGroupAttachedLeaf EmptyFOpenCLGroupAttachedLeaf(){
    struct FOpenCLGroupAttachedLeaf leaf;
    leaf.nbParticles = -1;
    // Redirect pointers to position
    leaf.positionsPointers[0] = NULLPTR;
    leaf.positionsPointers[1] = NULLPTR;
    leaf.positionsPointers[2] = NULLPTR;

    // Redirect pointers to data
    for(unsigned idxAttribute = 0 ; idxAttribute < NbSymbAttributes+NbAttributesPerParticle ; ++idxAttribute){
        leaf.attributes[idxAttribute] = NULLPTR;
    }
    return leaf;
}

bool FOpenCLGroupAttachedLeaf_isAttachedToSomething(const struct FOpenCLGroupAttachedLeaf* group){
    return (group->nbParticles != -1);
}
bool FOpenCLGroupAttachedLeaf_getNbParticles(const struct FOpenCLGroupAttachedLeaf* group){
    return (group->nbParticles);
}


/** One header is allocated at the beginning of each block */
struct FOpenCLGroupOfParticlesBlockHeader{
    MortonIndex startingIndex;
    MortonIndex endingIndex;
    int numberOfLeavesInBlock;
    int blockIndexesTableSize;

    //< The real number of particles allocated
    FSize nbParticlesAllocatedInGroup;
    //< Bytes difference/offset between position
    size_t positionOffset;
    //< Bytes difference/offset between attributes
    size_t attributeOffset;
    //< The total number of particles in the group
    FSize nbParticlesInGroup;
}__attribute__ ((aligned (DefaultStructAlign)));

/** Information about a leaf */
struct FOpenCLGroupOfParticlesLeafHeader {
    FSize nbParticles;
    size_t offSet;
}__attribute__ ((aligned (DefaultStructAlign)));


struct FOpenCLGroupOfParticles {
    //< The size of memoryBuffer in byte
    size_t allocatedMemoryInByte;
    //< Pointer to a block memory
    __global unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    __global struct FOpenCLGroupOfParticlesBlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    __global int*            blockIndexesTable;
    //< Pointer to leaves information
    __global struct FOpenCLGroupOfParticlesLeafHeader*     leafHeader;
    //< The total number of particles in the group
    const FSize nbParticlesInGroup;

    //< Pointers to particle position x, y, z
    __global FReal* particlePosition[3];

    //< Pointers to the particles data inside the block memory
    __global FParticleValueClass*      attributesBuffer;
    __global FParticleValueClass*      particleAttributes[NbSymbAttributes+NbAttributesPerParticle];
};

struct FOpenCLGroupOfParticles BuildFOpenCLGroupOfParticles(__global unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
                                                            __global unsigned char* inAttributeBuffer){
    struct FOpenCLGroupOfParticles group;
    group.allocatedMemoryInByte = (inAllocatedMemoryInByte);
    group.memoryBuffer = (inBuffer);

    // Move the pointers to the correct position
    group.blockHeader         = ((__global struct FOpenCLGroupOfParticlesBlockHeader*)group.memoryBuffer);
    group.blockIndexesTable   = ((__global int*)(group.memoryBuffer+sizeof(struct FOpenCLGroupOfParticlesBlockHeader)));
    group.leafHeader          = ((__global struct FOpenCLGroupOfParticlesLeafHeader*)(group.memoryBuffer+sizeof(struct FOpenCLGroupOfParticlesBlockHeader)+(group.blockHeader->blockIndexesTableSize*sizeof(int))));

    // Init particle pointers
    group.blockHeader->positionOffset = (sizeof(FReal) * group.blockHeader->nbParticlesAllocatedInGroup);
    group.particlePosition[0] = (__global FReal*) (( ((size_t)(group.leafHeader + group.blockHeader->numberOfLeavesInBlock))
                                                   +FOpenCLGroupOfParticlesMemoryAlignementBytes-1) & ~(FOpenCLGroupOfParticlesMemoryAlignementBytes-1));
    group.particlePosition[1] = (group.particlePosition[0] + group.blockHeader->nbParticlesAllocatedInGroup);
    group.particlePosition[2] = (group.particlePosition[1] + group.blockHeader->nbParticlesAllocatedInGroup);

    // Redirect pointer to data
    group.blockHeader->attributeOffset = (sizeof(FParticleValueClass) * group.blockHeader->nbParticlesAllocatedInGroup);
    __global unsigned char* previousPointer = ((__global unsigned char*)(group.particlePosition[2] + group.blockHeader->nbParticlesAllocatedInGroup));
    for(unsigned idxAttribute = 0 ; idxAttribute < NbSymbAttributes ; ++idxAttribute){
        group.particleAttributes[idxAttribute] = ((__global FParticleValueClass*)previousPointer);
        previousPointer += sizeof(FParticleValueClass)*group.blockHeader->nbParticlesAllocatedInGroup;
    }
    
    if(inAttributeBuffer){
        group.attributesBuffer = (__global FParticleValueClass*)inAttributeBuffer;    
        for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
            group.particleAttributes[idxAttribute+NbSymbAttributes] = ((__global FParticleValueClass*)inAttributeBuffer);
            inAttributeBuffer += sizeof(FParticleValueClass)*group.blockHeader->nbParticlesAllocatedInGroup;
        }
    }
    else{
        group.attributesBuffer = NULLPTR;
        for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
            group.particleAttributes[idxAttribute+NbSymbAttributes] = NULLPTR;
        }    
    }
    
    return group;
}
MortonIndex FOpenCLGroupOfParticles_getStartingIndex(const struct FOpenCLGroupOfParticles* group) {
    return group->blockHeader->startingIndex;
}
MortonIndex FOpenCLGroupOfParticles_getEndingIndex(const struct FOpenCLGroupOfParticles* group) {
    return group->blockHeader->endingIndex;
}
int FOpenCLGroupOfParticles_getNumberOfLeaves(const struct FOpenCLGroupOfParticles* group) {
    return group->blockHeader->numberOfLeavesInBlock;
}
bool FOpenCLGroupOfParticles_isInside(const struct FOpenCLGroupOfParticles* group, const MortonIndex inIndex) {
    return group->blockHeader->startingIndex <= inIndex && inIndex < group->blockHeader->endingIndex;
}
bool FOpenCLGroupOfParticles_exists(const struct FOpenCLGroupOfParticles* group, const MortonIndex inIndex) {
    return FOpenCLGroupOfParticles_isInside(group, inIndex) && (group->blockIndexesTable[inIndex-group->blockHeader->startingIndex] != FOpenCLGroupOfParticlesLeafIsEmptyFlag);
}
struct FOpenCLGroupAttachedLeaf FOpenCLGroupOfParticles_getLeaf(struct FOpenCLGroupOfParticles* group, const MortonIndex leafIndex){
    if(group->blockIndexesTable[leafIndex - group->blockHeader->startingIndex] != FOpenCLGroupOfParticlesLeafIsEmptyFlag){
        const int id = group->blockIndexesTable[leafIndex - group->blockHeader->startingIndex];
        return BuildFOpenCLGroupAttachedLeaf(group->leafHeader[id].nbParticles,
                                      group->particlePosition[0] + group->leafHeader[id].offSet,
                                        group->blockHeader->positionOffset,
                                      (group->attributesBuffer?group->particleAttributes[NbSymbAttributes] + group->leafHeader[id].offSet:NULLPTR),
                                        group->blockHeader->attributeOffset);
    }
    return EmptyFOpenCLGroupAttachedLeaf();
}


struct FOpenCLGroupOfCellsBlockHeader{
    MortonIndex startingIndex;
    MortonIndex endingIndex;
    int numberOfCellsInBlock;
    int blockIndexesTableSize;
} __attribute__ ((aligned (DefaultStructAlign)));


struct FOpenCLGroupOfCells {
    //< The size of the memoryBuffer
    size_t allocatedMemoryInByte;
    //< Pointer to a block memory
    __global unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    __global struct FOpenCLGroupOfCellsBlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    __global int*            blockIndexesTable;
    //< Pointer to the cells inside the block memory
    __global struct FSymboleCellClass*      blockCells;
    
    //< The multipole data
    __global FPoleCellClass* cellMultipoles;
    //< The local data
    __global FLocalCellClass* cellLocals;
};

struct FOpenCLGroupOfCells BuildFOpenCLGroupOfCells(__global unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
                __global unsigned char* inCellMultipoles, __global unsigned char* inCellLocals){
      struct FOpenCLGroupOfCells group;
      group.memoryBuffer = (inBuffer);
      group.allocatedMemoryInByte = (inAllocatedMemoryInByte);

    // Move the pointers to the correct position
    group.blockHeader         = (__global struct FOpenCLGroupOfCellsBlockHeader*)(inBuffer);
    inBuffer += sizeof(struct FOpenCLGroupOfCellsBlockHeader);
    group.blockIndexesTable   = (__global int*)(inBuffer);
    inBuffer += (group.blockHeader->blockIndexesTableSize*sizeof(int));
    group.blockCells          = (__global struct FSymboleCellClass*)(inBuffer);
    inBuffer += (group.blockHeader->numberOfCellsInBlock*sizeof(struct FSymboleCellClass));
    // Assert(((size_t)(inBuffer-group.memoryBuffer) == inAllocatedMemoryInByte);

    group.cellMultipoles = (__global FPoleCellClass*)inCellMultipoles;
    group.cellLocals = (__global FLocalCellClass*)inCellLocals;
    return group;
}
MortonIndex FOpenCLGroupOfCells_getStartingIndex(const struct FOpenCLGroupOfCells* group) {
    return group->blockHeader->startingIndex;
}
MortonIndex FOpenCLGroupOfCells_getEndingIndex(const struct FOpenCLGroupOfCells* group) {
    return group->blockHeader->endingIndex;
}
int FOpenCLGroupOfCells_getNumberOfCellsInBlock(const struct FOpenCLGroupOfCells* group) {
    return group->blockHeader->numberOfCellsInBlock;
}
int FOpenCLGroupOfCells_getSizeOfInterval(const struct FOpenCLGroupOfCells* group) {
    return group->blockHeader->blockIndexesTableSize;
}
bool FOpenCLGroupOfCells_isInside(const struct FOpenCLGroupOfCells* group, const MortonIndex inIndex){
    return group->blockHeader->startingIndex <= inIndex && inIndex < group->blockHeader->endingIndex;
}
bool FOpenCLGroupOfCells_exists(const struct FOpenCLGroupOfCells* group, const MortonIndex inIndex) {
    return FOpenCLGroupOfCells_isInside(group, inIndex) && (group->blockIndexesTable[inIndex-group->blockHeader->startingIndex] != FOpenCLGroupOfCellsCellIsEmptyFlag);
}
struct FWrappeCell FOpenCLGroupOfCells_getCompleteCell(struct FOpenCLGroupOfCells* group, const MortonIndex inIndex){
    size_t pos = group->blockIndexesTable[inIndex-group->blockHeader->startingIndex];
    struct FWrappeCell cell;
    cell.symb = &group->blockCells[pos];
    cell.up = &group->cellMultipoles[pos];
    cell.down = &group->cellLocals[pos];        
    return cell;
}

struct FWrappeCell FOpenCLGroupOfCells_getUpCell(struct FOpenCLGroupOfCells* group, const MortonIndex inIndex){
    size_t pos = group->blockIndexesTable[inIndex-group->blockHeader->startingIndex];
    struct FWrappeCell cell;
    cell.symb = &group->blockCells[pos];
    cell.up = &group->cellMultipoles[pos];
    cell.down = NULLPTR;       
    return cell;
}

struct FWrappeCell FOpenCLGroupOfCells_getDownCell(struct FOpenCLGroupOfCells* group, const MortonIndex inIndex){
    size_t pos = group->blockIndexesTable[inIndex-group->blockHeader->startingIndex];
        struct FWrappeCell cell;
    cell.symb = &group->blockCells[pos];
    cell.up = NULLPTR;
    cell.down =&group->cellLocals[pos];        
    return cell;
}

struct Uptr9{
    __global unsigned char* ptrs[9];
} __attribute__ ((aligned (DefaultStructAlign)));

struct size_t9{
    size_t v[9];
} __attribute__ ((aligned (DefaultStructAlign)));

struct Uptr343{
    __global unsigned char* ptrs[343];
};

/***************************************************************************/
/***************************************************************************/
/************************CHANGE THINGS HERE*********************************/
/***************************************************************************/


void P2M(struct FWrappeCell pole, const struct FOpenCLGroupAttachedLeaf particles, __global void* user_data) {
    *pole.up = particles.nbParticles;
}

void M2M(struct FWrappeCell  pole, struct FWrappeCell child[8], const int level, __global void* user_data) {
    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
        if(child[idxChild].symb){
            *pole.up += *child[idxChild].up;
        }
    }
}

void M2L(struct FWrappeCell const pole, const struct FWrappeCell distantNeighbors[343],
    const int size, const int level, __global void* user_data) {
    for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
        if(distantNeighbors[idxNeigh].symb){
            *pole.down += *distantNeighbors[idxNeigh].up;
        }
    }
}

void L2L(const struct FWrappeCell localCell, struct FWrappeCell child[8], const int level, __global void* user_data) {
    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
        if(child[idxChild].symb){
            *child[idxChild].down += *localCell.down;
        }
    }
}

void L2P(const struct FWrappeCell localCell, struct FOpenCLGroupAttachedLeaf particles, __global void* user_data){
    __global long long* partdown = particles.attributes[0];
    for(FSize idxPart = 0 ; idxPart < particles.nbParticles ; ++idxPart){
        partdown[idxPart] += *localCell.down;
    }
}

void P2P(const int3 pos,
             struct FOpenCLGroupAttachedLeaf  targets, const struct FOpenCLGroupAttachedLeaf sources,
             struct FOpenCLGroupAttachedLeaf directNeighborsParticles[27], int directNeighborsPositions[27], const int counter, __global void* user_data){
    long long cumul = sources.nbParticles-1;
    
    for(int idxNeigh = 0 ; idxNeigh < counter ; ++idxNeigh){
        if(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&directNeighborsParticles[idxNeigh])){
            cumul += directNeighborsParticles[idxNeigh].nbParticles;
        }
    }
    
    __global long long* partdown = targets.attributes[0];
    for(FSize idxPart = 0 ; idxPart < targets.nbParticles ; ++idxPart){
        partdown[idxPart] += cumul;
    }
}

void P2PRemote(const int3 pos,
             struct FOpenCLGroupAttachedLeaf  targets, const struct FOpenCLGroupAttachedLeaf  sources,
             struct FOpenCLGroupAttachedLeaf directNeighborsParticles, const int position, __global void* user_data){
    __global long long* partdown = targets.attributes[0];
    for(FSize idxPart = 0 ; idxPart < targets.nbParticles ; ++idxPart){
        partdown[idxPart] += directNeighborsParticles.nbParticles;
    }
}

int3 getCoordinate(const struct FWrappeCell cell) {
    int3 coord;
    coord.x = cell.symb->coordinates[0];
    coord.y = cell.symb->coordinates[1];
    coord.z = cell.symb->coordinates[2];
    return coord;    
}


/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

#define FOpenCLCheck( test ) { FOpenCLCheckCore((test), __FILE__, __LINE__); }
#define FOpenCLCheckAfterCall() { FOpenCLCheckCore((cudaGetLastError()), __FILE__, __LINE__); }
#define FOpenCLAssertLF(ARGS) if(!(ARGS)){ *((char*)0x09) = 'e'; }
//#define FOpenCLAssertLF(ARGS) ARGS;

#define FMGetOppositeNeighIndex(index) (27-(index)-1)
#define FMGetOppositeInterIndex(index) (343-(index)-1)

#define FOpenCLMax(x,y) ((x)<(y) ? (y) : (x))
#define FOpenCLMin(x,y) ((x)>(y) ? (y) : (x))


__kernel void FOpenCL__bottomPassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize,__global unsigned char* leafCellsUpPtr,
                                         __global unsigned char* containersPtr, size_t containersSize,
                                         __global void* userkernel ){
    struct FOpenCLGroupOfCells leafCells = BuildFOpenCLGroupOfCells(leafCellsPtr, leafCellsSize, leafCellsUpPtr, NULLPTR);
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr, containersSize, NULLPTR);

    const MortonIndex blockStartIdx = FOpenCLGroupOfCells_getStartingIndex(&leafCells);
    const MortonIndex blockEndIdx = FOpenCLGroupOfCells_getEndingIndex(&leafCells);
    
    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        if(FOpenCLGroupOfCells_exists(&leafCells, mindex)){
            struct FWrappeCell cell = FOpenCLGroupOfCells_getUpCell(&leafCells, mindex);
            FOpenCLAssertLF(cell.symb->mortonIndex == mindex);
            struct FOpenCLGroupAttachedLeaf particles = FOpenCLGroupOfParticles_getLeaf(&containers, mindex);
            FOpenCLAssertLF(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&particles));
            P2M(cell, particles, userkernel);
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////
/// Upward Pass
/////////////////////////////////////////////////////////////////////////////////////

__kernel void FOpenCL__upwardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize, __global unsigned char* currentCellsUpPtr,
                                  struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize, struct Uptr9 subCellGroupsUpPtr,
                                  int nbSubCellGroups, int idxLevel, __global void* userkernel){
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize, currentCellsUpPtr, NULLPTR);
    struct FOpenCLGroupOfCells subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx] = BuildFOpenCLGroupOfCells(subCellGroupsPtr.ptrs[idx], subCellGroupsSize.v[idx], subCellGroupsUpPtr.ptrs[idx], NULLPTR);
    }

    FOpenCLAssertLF(nbSubCellGroups != 0);
    const MortonIndex blockStartIdx = FOpenCLMax(FOpenCLGroupOfCells_getStartingIndex(&currentCells),
                                          FOpenCLGroupOfCells_getStartingIndex(&subCellGroups[0])>>3);
    const MortonIndex blockEndIdx   = FOpenCLMin(FOpenCLGroupOfCells_getEndingIndex(&currentCells),
                                          ((FOpenCLGroupOfCells_getEndingIndex(&subCellGroups[nbSubCellGroups-1])-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        if(FOpenCLGroupOfCells_exists(&currentCells, mindex)){
            struct FWrappeCell cell = FOpenCLGroupOfCells_getUpCell(&currentCells, mindex);
            FOpenCLAssertLF(cell.symb->mortonIndex == mindex);
            struct FWrappeCell child[8];

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( FOpenCLGroupOfCells_getEndingIndex(&subCellGroups[idxSubCellGroup]) <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                if(FOpenCLGroupOfCells_exists(&subCellGroups[idxSubCellGroup], (mindex<<3)+idxChild)){
                    child[idxChild] = FOpenCLGroupOfCells_getUpCell(&subCellGroups[idxSubCellGroup], (mindex<<3)+idxChild);
                }
                else{
                    child[idxChild].symb = NULLPTR;
                }
            }

            M2M(cell, child, idxLevel, userkernel);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass Mpi
/////////////////////////////////////////////////////////////////////////////////////


__kernel  void FOpenCL__transferInoutPassPerformMpi(__global unsigned char* currentCellsPtr, size_t currentCellsSize, __global unsigned char* currentCellsDownPtr,
                                             __global unsigned char* externalCellsPtr, size_t externalCellsSize, __global unsigned char* externalCellsUpPtr,
                                             int idxLevel, const __global struct OutOfBlockInteraction* outsideInteractions,
                                             size_t nbOutsideInteractions, __global void* userkernel){
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize, NULLPTR, currentCellsDownPtr);
    struct FOpenCLGroupOfCells cellsOther = BuildFOpenCLGroupOfCells(externalCellsPtr, externalCellsSize, externalCellsUpPtr, NULLPTR);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        if(FOpenCLGroupOfCells_exists(&cellsOther, outsideInteractions[outInterIdx].outIndex)){
            struct FWrappeCell interCell = FOpenCLGroupOfCells_getUpCell(&cellsOther, outsideInteractions[outInterIdx].outIndex);
            FOpenCLAssertLF(interCell.symb->mortonIndex == outsideInteractions[outInterIdx].outIndex);
            struct FWrappeCell cell = FOpenCLGroupOfCells_getDownCell(&currentCells, outsideInteractions[outInterIdx].insideIndex);
            FOpenCLAssertLF(cell.symb->mortonIndex == outsideInteractions[outInterIdx].insideIndex);

            struct FWrappeCell interactions[343];
            FSetToNullptr343(interactions);
            interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
            const int counter = 1;
            M2L( cell , interactions, counter, idxLevel, userkernel);
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass
/////////////////////////////////////////////////////////////////////////////////////



__kernel  void FOpenCL__transferInPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                        __global unsigned char* currentCellsUpPtr, __global unsigned char* currentCellsDownPtr,
                                       int idxLevel, __global void* userkernel){
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize, currentCellsUpPtr, currentCellsDownPtr);

    const MortonIndex blockStartIdx = FOpenCLGroupOfCells_getStartingIndex(&currentCells);
    const MortonIndex blockEndIdx = FOpenCLGroupOfCells_getEndingIndex(&currentCells);

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        if(FOpenCLGroupOfCells_exists(&currentCells, mindex)){
            struct FWrappeCell cell = FOpenCLGroupOfCells_getDownCell(&currentCells, mindex);
            FOpenCLAssertLF(cell.symb->mortonIndex == mindex);
            MortonIndex interactionsIndexes[189];
            int interactionsPosition[189];
            const int3 coord = (getCoordinate(cell));
            int counter = GetInteractionNeighbors(coord, idxLevel,interactionsIndexes,interactionsPosition);

            struct FWrappeCell interactions[343];
            FSetToNullptr343(interactions);
            int counterExistingCell = 0;

            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                    if(FOpenCLGroupOfCells_exists(&currentCells, interactionsIndexes[idxInter])){
                        struct FWrappeCell interCell = FOpenCLGroupOfCells_getUpCell(&currentCells, interactionsIndexes[idxInter]);
                        FOpenCLAssertLF(interCell.symb->mortonIndex == interactionsIndexes[idxInter]);
                        FOpenCLAssertLF(interactions[interactionsPosition[idxInter]].symb == NULLPTR);
                        interactions[interactionsPosition[idxInter]] = interCell;
                        counterExistingCell += 1;
                    }
                }
            }

            M2L( cell , interactions, counterExistingCell, idxLevel, userkernel);
        }
    }
}



__kernel void FOpenCL__transferInoutPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                         __global unsigned char*  currentCellsUpPtr, __global unsigned char*  currentCellsDownPtr,
                                         __global unsigned char* externalCellsPtr, size_t externalCellsSize,
                                         __global unsigned char* externalCellsUpPtr, __global unsigned char* externalCellsDownPtr,
                                         int idxLevel, const __global struct OutOfBlockInteraction* outsideInteractions,
                                         size_t nbOutsideInteractions, __global void* userkernel){
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize, currentCellsUpPtr, currentCellsDownPtr);
    struct FOpenCLGroupOfCells cellsOther = BuildFOpenCLGroupOfCells(externalCellsPtr, externalCellsSize, externalCellsUpPtr, externalCellsDownPtr);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        if(FOpenCLGroupOfCells_exists(&cellsOther, outsideInteractions[outInterIdx].outIndex)){
            struct FWrappeCell interCell = FOpenCLGroupOfCells_getCompleteCell(&cellsOther, outsideInteractions[outInterIdx].outIndex);
            FOpenCLAssertLF(interCell.symb->mortonIndex == outsideInteractions[outInterIdx].outIndex);
            struct FWrappeCell cell = FOpenCLGroupOfCells_getCompleteCell(&currentCells, outsideInteractions[outInterIdx].insideIndex);
            FOpenCLAssertLF(cell.symb->mortonIndex == outsideInteractions[outInterIdx].insideIndex);

            struct FWrappeCell interactions[343];
            FSetToNullptr343(interactions);
            interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
            const int counter = 1;
            M2L( cell , interactions, counter, idxLevel, userkernel);

            interactions[outsideInteractions[outInterIdx].outPosition].symb = NULLPTR;
            interactions[FMGetOppositeInterIndex(outsideInteractions[outInterIdx].outPosition)] = cell;
            M2L( interCell , interactions, counter, idxLevel, userkernel);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/// Downard Pass
/////////////////////////////////////////////////////////////////////////////////////


__kernel void FOpenCL__downardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize, __global unsigned char* currentCellsDownPtr,
                                   struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize, struct Uptr9 subCellGroupsDownPtr,
                                   int nbSubCellGroups, int idxLevel, __global void* userkernel){
    FOpenCLAssertLF(nbSubCellGroups != 0);
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize, NULLPTR, currentCellsDownPtr);
    struct FOpenCLGroupOfCells subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx] = BuildFOpenCLGroupOfCells(subCellGroupsPtr.ptrs[idx], subCellGroupsSize.v[idx], NULLPTR, subCellGroupsDownPtr.ptrs[idx]);
    }

    const MortonIndex blockStartIdx = FOpenCLMax(FOpenCLGroupOfCells_getStartingIndex(&currentCells),
                                          FOpenCLGroupOfCells_getStartingIndex(&subCellGroups[0])>>3);
    const MortonIndex blockEndIdx   = FOpenCLMin(FOpenCLGroupOfCells_getEndingIndex(&currentCells),
                                          ((FOpenCLGroupOfCells_getEndingIndex(&subCellGroups[nbSubCellGroups-1])-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        if(FOpenCLGroupOfCells_exists(&currentCells, mindex)){
            struct FWrappeCell cell = FOpenCLGroupOfCells_getDownCell(&currentCells, mindex);
            FOpenCLAssertLF(cell.symb->mortonIndex == mindex);
            struct FWrappeCell child[8];

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( FOpenCLGroupOfCells_getEndingIndex(&subCellGroups[idxSubCellGroup]) <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                if(FOpenCLGroupOfCells_exists(&subCellGroups[idxSubCellGroup], (mindex<<3)+idxChild)){
                    child[idxChild] = FOpenCLGroupOfCells_getDownCell(&subCellGroups[idxSubCellGroup], (mindex<<3)+idxChild);                
                }
                else{
                    child[idxChild].symb = NULLPTR;
                }
            }

            L2L(cell, child, idxLevel, userkernel);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass MPI
/////////////////////////////////////////////////////////////////////////////////////


__kernel void FOpenCL__directInoutPassPerformMpi(__global unsigned char* containersPtr, size_t containersSize, __global unsigned char* containersDownPtr,
                                          __global unsigned char* externalContainersPtr, size_t externalContainersSize, __global unsigned char* outsideInteractionsCl,
                                          const __global struct OutOfBlockInteraction* outsideInteractions,
                                          size_t nbOutsideInteractions, const int treeHeight, __global void* userkernel){
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr, containersSize, containersDownPtr);
    struct FOpenCLGroupOfParticles containersOther = BuildFOpenCLGroupOfParticles(externalContainersPtr, externalContainersSize, NULLPTR);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        struct FOpenCLGroupAttachedLeaf interParticles = FOpenCLGroupOfParticles_getLeaf(&containersOther, outsideInteractions[outInterIdx].outIndex);
        if(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&interParticles)){
            struct FOpenCLGroupAttachedLeaf particles = FOpenCLGroupOfParticles_getLeaf(&containers, outsideInteractions[outInterIdx].insideIndex);
            FOpenCLAssertLF(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&particles));

            P2PRemote( GetPositionFromMorton(outsideInteractions[outInterIdx].insideIndex, treeHeight-1), particles, particles , interParticles, outsideInteractions[outInterIdx].outPosition, userkernel);
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass
/////////////////////////////////////////////////////////////////////////////////////



__kernel void FOpenCL__directInPassPerform(__global unsigned char* containersPtr, size_t containersSize, __global unsigned char* containersDownPtr,
                                    const int treeHeight, __global void* userkernel){
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr, containersSize, containersDownPtr);

    const MortonIndex blockStartIdx = FOpenCLGroupOfParticles_getStartingIndex(&containers);
    const MortonIndex blockEndIdx = FOpenCLGroupOfParticles_getEndingIndex(&containers);

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        struct FOpenCLGroupAttachedLeaf particles = FOpenCLGroupOfParticles_getLeaf(&containers, mindex);
        if(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&particles)){
            MortonIndex interactionsIndexes[26];
            int interactionsPosition[26];
            const int3 coord = GetPositionFromMorton(mindex, treeHeight-1);
            int counter = GetNeighborsIndexes(coord, treeHeight,interactionsIndexes,interactionsPosition);

            struct FOpenCLGroupAttachedLeaf interactionsObjects[27];
            int neighPosition[26];
            int counterExistingCell = 0;

            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                    interactionsObjects[counterExistingCell] = FOpenCLGroupOfParticles_getLeaf(&containers, interactionsIndexes[idxInter]);
                    if(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&interactionsObjects[counterExistingCell])){
                        neighPosition[counterExistingCell] = interactionsPosition[idxInter];
                        counterExistingCell += 1;
                    }
                }
            }

            P2P( coord, particles, particles , interactionsObjects, neighPosition, counterExistingCell, userkernel);
        }
    }
}



__kernel void FOpenCL__directInoutPassPerform(__global unsigned char* containersPtr, size_t containersSize, __global unsigned char* containersDownPtr,
                                       __global unsigned char* externalContainersPtr, size_t externalContainersSize, __global unsigned char* externalContainersDownPtr,
                                       const __global struct OutOfBlockInteraction* outsideInteractions,
                                       size_t nbOutsideInteractions, const int treeHeight, __global void* userkernel){
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr, containersSize, containersDownPtr);
    struct FOpenCLGroupOfParticles containersOther = BuildFOpenCLGroupOfParticles(externalContainersPtr, externalContainersSize, externalContainersDownPtr);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        struct FOpenCLGroupAttachedLeaf interParticles = FOpenCLGroupOfParticles_getLeaf(&containersOther, outsideInteractions[outInterIdx].outIndex);
        if(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&interParticles)){
            struct FOpenCLGroupAttachedLeaf particles = FOpenCLGroupOfParticles_getLeaf(&containers, outsideInteractions[outInterIdx].insideIndex);
            FOpenCLAssertLF(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&particles));
            FOpenCLAssertLF(particles.nbParticles);
            FOpenCLAssertLF(interParticles.nbParticles);

            P2PRemote( GetPositionFromMorton(outsideInteractions[outInterIdx].insideIndex, treeHeight-1), particles, particles , interParticles, outsideInteractions[outInterIdx].outPosition, userkernel );

            P2PRemote( GetPositionFromMorton(outsideInteractions[outInterIdx].outIndex, treeHeight-1), interParticles, interParticles , particles, FMGetOppositeNeighIndex(outsideInteractions[outInterIdx].outPosition), userkernel);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/// Merge Pass
/////////////////////////////////////////////////////////////////////////////////////



__kernel void FOpenCL__mergePassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize, __global unsigned char* leafCellsDownPtr,
                                 __global unsigned char* containersPtr, size_t containersSize, __global unsigned char* containersDownPtr,
                                 __global void* userkernel){
    struct FOpenCLGroupOfCells leafCells = BuildFOpenCLGroupOfCells(leafCellsPtr,leafCellsSize, NULLPTR, leafCellsDownPtr);
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr,containersSize, containersDownPtr);

    const MortonIndex blockStartIdx = FOpenCLGroupOfCells_getStartingIndex(&leafCells);
    const MortonIndex blockEndIdx = FOpenCLGroupOfCells_getEndingIndex(&leafCells);

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        if(FOpenCLGroupOfCells_exists(&leafCells, mindex)){
            struct FWrappeCell cell = FOpenCLGroupOfCells_getDownCell(&leafCells, mindex);
            FOpenCLAssertLF(cell.symb->mortonIndex == mindex);
            struct FOpenCLGroupAttachedLeaf particles = FOpenCLGroupOfParticles_getLeaf(&containers, mindex);
            FOpenCLAssertLF(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&particles));
            L2P(cell, particles, userkernel);
        }
    }
}

