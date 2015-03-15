/** This file contains the prototype for a kernel in opencl */
// @SCALFMM_PRIVATE


/***************************************************************************/
/***************************************************************************/
/************************CHANGE THINGS HERE*********************************/
/***************************************************************************/

typedef ___FReal___ FReal;
typedef ___FParticleValueClass___ FParticleValueClass;
typedef long long int MortonIndex;

#define FCellClassSize ___FCellClassSize___
#define FCellUpOffset ___FCellUpOffset___
#define FCellDownOffset ___FCellDownOffset___
#define FCellMortonOffset ___FCellMortonOffset___
#define FCellCoordinateOffset ___FCellCoordinateOffset___
#define FOpenCLGroupOfCellsCellIsEmptyFlag  ((MortonIndex)-1)

#define NbAttributesPerParticle ___NbAttributesPerParticle___

#define FOpenCLGroupOfParticlesMemoryAlignementBytes  32
#define FOpenCLGroupOfParticlesMemoryAlignementParticles (FOpenCLGroupOfParticlesMemoryAlignementBytes/sizeof(FReal))
#define FOpenCLGroupOfParticlesLeafIsEmptyFlag ((MortonIndex)-1)

#define NULLPTR (0)

#define DefaultStructAlign ___DefaultStructAlign___

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

typedef struct OutOfBlockInteraction{
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


void FSetToNullptr343(__global const unsigned char* ptrs[343]){
    int idx;
    for( idx = 0 ; idx < 343 ; ++idx){
        ptrs[idx] = NULLPTR;
    }
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/


struct FOpenCLGroupAttachedLeaf {
    //< Nb of particles in the current leaf
    int nbParticles;
    //< Pointers to the positions of the particles
    __global FReal* positionsPointers[3];
    //< Pointers to the attributes of the particles
    __global FParticleValueClass* attributes[NbAttributesPerParticle];
};

struct FOpenCLGroupAttachedLeaf BuildFOpenCLGroupAttachedLeaf(const int inNbParticles, __global FReal* inPositionBuffer, const size_t inLeadingPosition,
                       __global FParticleValueClass* inAttributesBuffer, const size_t inLeadingAttributes){
    struct FOpenCLGroupAttachedLeaf leaf;
    leaf.nbParticles = (inNbParticles);
    // Redirect pointers to position
    leaf.positionsPointers[0] = inPositionBuffer;
    leaf.positionsPointers[1] = (__global FReal*)(((__global unsigned char*)inPositionBuffer) + inLeadingPosition);
    leaf.positionsPointers[2] = (__global FReal*)(((__global unsigned char*)inPositionBuffer) + inLeadingPosition*2);

    // Redirect pointers to data
    for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
        leaf.attributes[idxAttribute] = (__global FParticleValueClass*)(((__global unsigned char*)inAttributesBuffer) + idxAttribute*inLeadingAttributes);
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
    for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
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
    int nbParticlesAllocatedInGroup;
    //< Bytes difference/offset between position
    size_t positionOffset;
    //< Bytes difference/offset between attributes
    size_t attributeOffset;
    //< The total number of particles in the group
    int nbParticlesInGroup;
}__attribute__ ((aligned (DefaultStructAlign)));

/** Information about a leaf */
struct FOpenCLGroupOfParticlesLeafHeader {
    int nbParticles;
    size_t offSet;
}__attribute__ ((aligned (DefaultStructAlign)));


struct FOpenCLGroupOfParticles {
    //< The size of memoryBuffer in byte
    int allocatedMemoryInByte;
    //< Pointer to a block memory
    __global unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    __global struct FOpenCLGroupOfParticlesBlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    __global int*            blockIndexesTable;
    //< Pointer to leaves information
    __global struct FOpenCLGroupOfParticlesLeafHeader*     leafHeader;
    //< The total number of particles in the group
    const int nbParticlesInGroup;

    //< Pointers to particle position x, y, z
    __global FReal* particlePosition[3];

    //< Pointers to the particles data inside the block memory
    __global FParticleValueClass*      particleAttributes[NbAttributesPerParticle];
};

struct FOpenCLGroupOfParticles BuildFOpenCLGroupOfParticles(__global unsigned char* inBuffer, const size_t inAllocatedMemoryInByte){
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
    for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
        group.particleAttributes[idxAttribute] = ((__global FParticleValueClass*)previousPointer);
        previousPointer += sizeof(FParticleValueClass)*group.blockHeader->nbParticlesAllocatedInGroup;
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
                                        group->particleAttributes[0] + group->leafHeader[id].offSet,
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
    int allocatedMemoryInByte;
    //< Pointer to a block memory
    __global unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    __global struct FOpenCLGroupOfCellsBlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    __global int*            blockIndexesTable;
    //< Pointer to the cells inside the block memory
    __global unsigned char*      blockCells;
};

struct FOpenCLGroupOfCells BuildFOpenCLGroupOfCells(__global unsigned char* inBuffer, const size_t inAllocatedMemoryInByte){
      struct FOpenCLGroupOfCells group;
      group.memoryBuffer = (inBuffer);
      group.allocatedMemoryInByte = (inAllocatedMemoryInByte);

    // Move the pointers to the correct position
    group.blockHeader         = (__global struct FOpenCLGroupOfCellsBlockHeader*)(group.memoryBuffer);
    group.blockIndexesTable   = (__global int*)(group.memoryBuffer+sizeof(struct FOpenCLGroupOfCellsBlockHeader));
    group.blockCells          = (__global unsigned char*)(group.memoryBuffer+sizeof(struct FOpenCLGroupOfCellsBlockHeader)+(group.blockHeader->blockIndexesTableSize*sizeof(int)));
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
__global unsigned char* FOpenCLGroupOfCells_getCell(struct FOpenCLGroupOfCells* group, const MortonIndex inIndex){
    if( FOpenCLGroupOfCells_exists(group, inIndex) ) return &group->blockCells[FCellClassSize*group->blockIndexesTable[inIndex-group->blockHeader->startingIndex]];
    else return NULLPTR;
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


void P2M(__global unsigned char* pole, const struct FOpenCLGroupAttachedLeaf particles, __global void* user_data) {
    __global long long* up = (__global long long*)(pole+FCellUpOffset);
    (*up) = particles.nbParticles;
}

void M2M(__global unsigned char*  pole, __global unsigned char* child[8], const int level, __global void* user_data) {
    __global long long* up = (__global long long*)(pole+FCellUpOffset);

    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
        if(child[idxChild]){
            const __global long long* childup = (const __global long long*)(child[idxChild]+FCellUpOffset);
            (*up) += (*childup);
        }
    }
}

void M2L(__global unsigned char* const pole, const __global unsigned char* distantNeighbors[343],
    const int size, const int level, __global void* user_data) {
    __global long long* down = (__global long long*)(pole+FCellDownOffset);
    
    for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
        if(distantNeighbors[idxNeigh]){
            const __global long long* neighup = (const __global long long*)(distantNeighbors[idxNeigh]+FCellUpOffset);
            (*down) += (*neighup);
        }
    }
}

void L2L(__global const unsigned char* localCell, __global unsigned char* child[8], const int level, __global void* user_data) {
    const __global long long* down = (const __global long long*)(localCell+FCellDownOffset);
    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
        if(child[idxChild]){
            __global long long* childdown = (__global long long*)(child[idxChild]+FCellDownOffset);
            (*childdown) += (*down);
        }
    }
}

void L2P(__global const unsigned char* localCell, struct FOpenCLGroupAttachedLeaf particles, __global void* user_data){
    const __global long long* down = (const __global long long*)(localCell+FCellDownOffset);
    __global long long* partdown = particles.attributes[0];
    for(int idxPart = 0 ; idxPart < particles.nbParticles ; ++idxPart){
        partdown[idxPart] += (*down);
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
    for(int idxPart = 0 ; idxPart < targets.nbParticles ; ++idxPart){
        partdown[idxPart] += cumul;
    }
}

void P2PRemote(const int3 pos,
             struct FOpenCLGroupAttachedLeaf  targets, const struct FOpenCLGroupAttachedLeaf  sources,
             struct FOpenCLGroupAttachedLeaf directNeighborsParticles, const int position, __global void* user_data){
    __global long long* partdown = targets.attributes[0];
    for(int idxPart = 0 ; idxPart < targets.nbParticles ; ++idxPart){
        partdown[idxPart] += directNeighborsParticles.nbParticles;
    }
}

MortonIndex getMortonIndex(__global const unsigned char* cell, __global void* user_data) {
    __global MortonIndex* mindex = (__global MortonIndex*)(cell+FCellMortonOffset);
    return (*mindex);
}

int3 getCoordinate(__global const unsigned char* cell, __global void* user_data) {
    __global int* cellcoord = (__global int*)(cell+FCellCoordinateOffset);
    int3 coord;
    coord.x = cellcoord[0];
    coord.y = cellcoord[1];
    coord.z = cellcoord[2];
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


__kernel void FOpenCL__bottomPassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize,
                                         __global unsigned char* containersPtr, size_t containersSize,
                                         __global void* userkernel/*, __global int* output */){                                             
    struct FOpenCLGroupOfCells leafCells = BuildFOpenCLGroupOfCells(leafCellsPtr, leafCellsSize);
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr, containersSize);

    const MortonIndex blockStartIdx = FOpenCLGroupOfCells_getStartingIndex(&leafCells);
    const MortonIndex blockEndIdx = FOpenCLGroupOfCells_getEndingIndex(&leafCells);
    
    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        __global unsigned char* cell = FOpenCLGroupOfCells_getCell(&leafCells, mindex);
        if(cell){
            FOpenCLAssertLF(getMortonIndex(cell, userkernel) == mindex);
            struct FOpenCLGroupAttachedLeaf particles = FOpenCLGroupOfParticles_getLeaf(&containers, mindex);
            FOpenCLAssertLF(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&particles));
            P2M(cell, particles, userkernel);
            /*for(int idx = 0 ; idx < FCellClassSize ; ++idx){
                cell[idx] = 'z';            
            }*/
            /*output[0] = blockStartIdx;
            output[1] = blockEndIdx;
            output[2] = particles.nbParticles;
            output[3] = getMortonIndex(cell, userkernel);
            output[4] = mindex;
            output[5] = FOpenCLGroupOfCells_exists(&leafCells, mindex);
            output[6] = FOpenCLGroupOfParticles_getStartingIndex(&containers);
            output[7] = FOpenCLGroupOfParticles_getEndingIndex(&containers);
            output[8] = FOpenCLGroupOfParticles_getNumberOfLeaves(&containers);
            int count = 0;
            for(int idx = FOpenCLGroupOfParticles_getStartingIndex(&containers) ; idx < FOpenCLGroupOfParticles_getEndingIndex(&containers) ; ++idx){
                if((&containers)->blockIndexesTable[idx - (&containers)->blockHeader->startingIndex] != -1) count++;            
            }
            output[8] = FOpenCLGroupAttachedLeaf_isAttachedToSomething(&particles);
            output[9] = FOpenCLGroupAttachedLeaf_getNbParticles(&particles);//sizeof(struct FOpenCLGroupOfParticlesBlockHeader);//count;
            //return;//TODO
            //__global long long* up = (__global long long*)(((unsigned char*)output)+FCellUpOffset);
            //(*up) = particles.nbParticles;
            return;*/
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////
/// Upward Pass
/////////////////////////////////////////////////////////////////////////////////////

__kernel void FOpenCL__upwardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                  struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize,
                                  int nbSubCellGroups, int idxLevel, __global void* userkernel){
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize);
    struct FOpenCLGroupOfCells subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx] = BuildFOpenCLGroupOfCells(subCellGroupsPtr.ptrs[idx], subCellGroupsSize.v[idx]);
    }

    FOpenCLAssertLF(nbSubCellGroups != 0);
    const MortonIndex blockStartIdx = FOpenCLMax(FOpenCLGroupOfCells_getStartingIndex(&currentCells),
                                          FOpenCLGroupOfCells_getStartingIndex(&subCellGroups[0])>>3);
    const MortonIndex blockEndIdx   = FOpenCLMin(FOpenCLGroupOfCells_getEndingIndex(&currentCells),
                                          ((FOpenCLGroupOfCells_getEndingIndex(&subCellGroups[nbSubCellGroups-1])-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        __global unsigned char* cell = FOpenCLGroupOfCells_getCell(&currentCells, mindex);
        if(cell){
            FOpenCLAssertLF(getMortonIndex(cell, userkernel) == mindex);
            __global unsigned char* child[8] = {NULLPTR,NULLPTR,NULLPTR,NULLPTR,NULLPTR,NULLPTR,NULLPTR,NULLPTR};

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( FOpenCLGroupOfCells_getEndingIndex(&subCellGroups[idxSubCellGroup]) <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                child[idxChild] = FOpenCLGroupOfCells_getCell(&subCellGroups[idxSubCellGroup], (mindex<<3)+idxChild);
                FOpenCLAssertLF(child[idxChild] == NULLPTR || getMortonIndex(child[idxChild], userkernel) == ((mindex<<3)+idxChild));
            }

            M2M(cell, child, idxLevel, userkernel);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass Mpi
/////////////////////////////////////////////////////////////////////////////////////


__kernel  void FOpenCL__transferInoutPassPerformMpi(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                             __global unsigned char* externalCellsPtr, size_t externalCellsSize,
                                             int idxLevel, const __global struct OutOfBlockInteraction* outsideInteractions,
                                             size_t nbOutsideInteractions, __global void* userkernel){
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize);
    struct FOpenCLGroupOfCells cellsOther = BuildFOpenCLGroupOfCells(externalCellsPtr, externalCellsSize);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        __global unsigned char* interCell = FOpenCLGroupOfCells_getCell(&cellsOther, outsideInteractions[outInterIdx].outIndex);
        if(interCell){
            FOpenCLAssertLF(getMortonIndex(interCell, userkernel) == outsideInteractions[outInterIdx].outIndex);
            __global unsigned char* cell = FOpenCLGroupOfCells_getCell(&currentCells, outsideInteractions[outInterIdx].insideIndex);
            FOpenCLAssertLF(cell);
            FOpenCLAssertLF(getMortonIndex(cell, userkernel) == outsideInteractions[outInterIdx].insideIndex);

            const __global unsigned char* interactions[343];
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
                                       int idxLevel, __global void* userkernel){
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize);

    const MortonIndex blockStartIdx = FOpenCLGroupOfCells_getStartingIndex(&currentCells);
    const MortonIndex blockEndIdx = FOpenCLGroupOfCells_getEndingIndex(&currentCells);

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        __global unsigned char* cell = FOpenCLGroupOfCells_getCell(&currentCells, mindex);
        if(cell){
            FOpenCLAssertLF(getMortonIndex(cell, userkernel) == mindex);
            MortonIndex interactionsIndexes[189];
            int interactionsPosition[189];
            const int3 coord = (getCoordinate(cell, userkernel));
            int counter = GetInteractionNeighbors(coord, idxLevel,interactionsIndexes,interactionsPosition);

            const __global unsigned char* interactions[343];
            FSetToNullptr343(interactions);
            int counterExistingCell = 0;

            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                    __global unsigned char* interCell = FOpenCLGroupOfCells_getCell(&currentCells, interactionsIndexes[idxInter]);
                    if(interCell){
                        FOpenCLAssertLF(getMortonIndex(interCell, userkernel) == interactionsIndexes[idxInter]);
                        FOpenCLAssertLF(interactions[interactionsPosition[idxInter]] == NULLPTR);
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
                                         __global unsigned char* externalCellsPtr, size_t externalCellsSize,
                                         int idxLevel, const __global struct OutOfBlockInteraction* outsideInteractions,
                                         size_t nbOutsideInteractions, __global void* userkernel){
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize);
    struct FOpenCLGroupOfCells cellsOther = BuildFOpenCLGroupOfCells(externalCellsPtr, externalCellsSize);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        __global unsigned char* interCell = FOpenCLGroupOfCells_getCell(&cellsOther, outsideInteractions[outInterIdx].outIndex);
        if(interCell){
            FOpenCLAssertLF(getMortonIndex(interCell, userkernel) == outsideInteractions[outInterIdx].outIndex);
            __global unsigned char* cell = FOpenCLGroupOfCells_getCell(&currentCells, outsideInteractions[outInterIdx].insideIndex);
            FOpenCLAssertLF(cell);
            FOpenCLAssertLF(getMortonIndex(cell, userkernel) == outsideInteractions[outInterIdx].insideIndex);

            const __global unsigned char* interactions[343];
            FSetToNullptr343(interactions);
            interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
            const int counter = 1;
            M2L( cell , interactions, counter, idxLevel, userkernel);

            interactions[outsideInteractions[outInterIdx].outPosition] = NULLPTR;
            interactions[FMGetOppositeInterIndex(outsideInteractions[outInterIdx].outPosition)] = cell;
            M2L( interCell , interactions, counter, idxLevel, userkernel);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/// Downard Pass
/////////////////////////////////////////////////////////////////////////////////////


__kernel void FOpenCL__downardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                   struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize,
                                   int nbSubCellGroups, int idxLevel, __global void* userkernel){
    FOpenCLAssertLF(nbSubCellGroups != 0);
    struct FOpenCLGroupOfCells currentCells = BuildFOpenCLGroupOfCells(currentCellsPtr, currentCellsSize);
    struct FOpenCLGroupOfCells subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx] = BuildFOpenCLGroupOfCells(subCellGroupsPtr.ptrs[idx], subCellGroupsSize.v[idx]);
    }

    const MortonIndex blockStartIdx = FOpenCLMax(FOpenCLGroupOfCells_getStartingIndex(&currentCells),
                                          FOpenCLGroupOfCells_getStartingIndex(&subCellGroups[0])>>3);
    const MortonIndex blockEndIdx   = FOpenCLMin(FOpenCLGroupOfCells_getEndingIndex(&currentCells),
                                          ((FOpenCLGroupOfCells_getEndingIndex(&subCellGroups[nbSubCellGroups-1])-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        __global unsigned char* cell = FOpenCLGroupOfCells_getCell(&currentCells, mindex);
        if(cell){
            FOpenCLAssertLF(getMortonIndex(cell, userkernel) == mindex);
            __global unsigned char* child[8] = {NULLPTR,NULLPTR,NULLPTR,NULLPTR,NULLPTR,NULLPTR,NULLPTR,NULLPTR};

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( FOpenCLGroupOfCells_getEndingIndex(&subCellGroups[idxSubCellGroup]) <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                child[idxChild] = FOpenCLGroupOfCells_getCell(&subCellGroups[idxSubCellGroup], (mindex<<3)+idxChild);
                FOpenCLAssertLF(child[idxChild] == NULLPTR || getMortonIndex(child[idxChild], userkernel) == ((mindex<<3)+idxChild));
            }

            L2L(cell, child, idxLevel, userkernel);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass MPI
/////////////////////////////////////////////////////////////////////////////////////


__kernel void FOpenCL__directInoutPassPerformMpi(__global unsigned char* containersPtr, size_t containersSize,
                                          __global unsigned char* externalContainersPtr, size_t externalContainersSize,
                                          const __global struct OutOfBlockInteraction* outsideInteractions,
                                          size_t nbOutsideInteractions, const int treeHeight, __global void* userkernel){
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr, containersSize);
    struct FOpenCLGroupOfParticles containersOther = BuildFOpenCLGroupOfParticles(externalContainersPtr, externalContainersSize);

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



__kernel void FOpenCL__directInPassPerform(__global unsigned char* containersPtr, size_t containersSize,
                                    const int treeHeight, __global void* userkernel){
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr, containersSize);

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



__kernel void FOpenCL__directInoutPassPerform(__global unsigned char* containersPtr, size_t containersSize,
                                       __global unsigned char* externalContainersPtr, size_t externalContainersSize,
                                       const __global struct OutOfBlockInteraction* outsideInteractions,
                                       size_t nbOutsideInteractions, const int treeHeight, __global void* userkernel){
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr, containersSize);
    struct FOpenCLGroupOfParticles containersOther = BuildFOpenCLGroupOfParticles(externalContainersPtr, externalContainersSize);

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



__kernel void FOpenCL__mergePassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize,
                                 __global unsigned char* containersPtr, size_t containersSize,
                                 __global void* userkernel){
    struct FOpenCLGroupOfCells leafCells = BuildFOpenCLGroupOfCells(leafCellsPtr,leafCellsSize);
    struct FOpenCLGroupOfParticles containers = BuildFOpenCLGroupOfParticles(containersPtr,containersSize);

    const MortonIndex blockStartIdx = FOpenCLGroupOfCells_getStartingIndex(&leafCells);
    const MortonIndex blockEndIdx = FOpenCLGroupOfCells_getEndingIndex(&leafCells);

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        __global unsigned char* cell = FOpenCLGroupOfCells_getCell(&leafCells, mindex);
        if(cell){
            FOpenCLAssertLF(getMortonIndex(cell, userkernel) == mindex);
            struct FOpenCLGroupAttachedLeaf particles = FOpenCLGroupOfParticles_getLeaf(&containers, mindex);
            FOpenCLAssertLF(FOpenCLGroupAttachedLeaf_isAttachedToSomething(&particles));
            L2P(cell, particles, userkernel);
        }
    }
}

