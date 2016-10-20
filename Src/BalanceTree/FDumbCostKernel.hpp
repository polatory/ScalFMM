// ===================================================================================
// Copyright ScalFmm 2016 INRIA, Olivier Coulaud, Bérenger Bramas,
// Matthias Messner olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the
// FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
// An extension to the license is given to allow static linking of scalfmm
// inside a proprietary application (no matter its license).
// See the main license file for more details.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE


#ifndef DUMBCOSTKERNEL_HPP
#define DUMBCOSTKERNEL_HPP

#include <stdexcept>

#include "Utils/FGlobal.hpp"

#include "Utils/FSmartPointer.hpp"

#include "Components/FAbstractKernels.hpp"

class FTreeCoordinate;

/**
 * \author Quentin Khan, Matthias Messner (original file: FChebFlopsSymKernel)
 * \brief Dumb cost computation for FMM.
 *
 * Please read the license
 *
 * This kernel implements a dumb cost the cost computation for FMM operators
 * exploiting the symmetries in the far-field. It implements all interfaces
 * (P2P, P2M, M2M, M2L, L2L, L2P) which are required by the
 * FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * \tparam Freal Type of real number representation
 * \tparam CellClass Type of cell
 * \tparam ContainerClass Type of container to store particles
 * \tparam OctreeClass Class of the octree to work on.
 */
template < typename FReal, class CellClass, class ContainerClass, class OctreeClass>
class FDumbCostKernel : public FAbstractKernels<CellClass, ContainerClass>
{
public:
    // Class types available to the rest of the world
    /// Type of real number representation
    using _FReal = FReal;
    /// Type of cell
    using _CellClass = CellClass;
    /// Type of container to store particles
    using _ContainerClass = ContainerClass;
    /// Class of the octree to work on
    using _OctreeClass = OctreeClass;


    /// The tree the kernel is working on. Needed to attain cells in the P2P
    /// operator (we only get the particles containers otherwise)
    const OctreeClass* const _tree;

    /// The tree height
    const unsigned int _treeHeight;

    /// count permuted local and multipole expansions
    unsigned int _countExp[343];

    /// Flops count for each operator of the FMM. 
    unsigned long long flopsP2M = 0,
        flopsM2M = 0,
        flopsM2L = 0,
        flopsL2L = 0,
        flopsL2P = 0,
        flopsP2P = 0;

    /// Operators count.
    unsigned long long countP2M = 0,
        countM2M = 0,
        countM2L = 0,
        countL2L = 0,
        countL2P = 0,
        countP2P = 0;

public:
    /**
     * The constructor initializes all constant attributes and it reads the
     * precomputed and compressed M2L operators from a binary file (a
     * runtime_error is thrown if the required file is not valid).
     */
    FDumbCostKernel(OctreeClass* tree, const FReal Epsilon) :
          _tree(tree),
          _treeHeight(_tree->getHeight()) {
        //_countExp = new unsigned int [343];
    }

    /** Copy constructor */
    FDumbCostKernel(const FDumbCostKernel& other) :
        _tree(other._tree),
        _treeHeight(other._treeHeight) {
        // _countExp = new unsigned int [343];
    }

    /** Destructor */
    ~FDumbCostKernel() {
        // delete [] _countExp;
    }

    void P2M(CellClass* const cell, const ContainerClass* const SourceParticles) override {
        FSize tmpCost = SourceParticles->getNbParticles();
        flopsP2M += tmpCost;
        cell->addCost(tmpCost);
        countP2M++;
    }


    void M2M(CellClass* const FRestrict cell,
             const CellClass* const FRestrict* const FRestrict ChildrenCells,
             const int /*TreeLevel*/) override {
        FSize flops = 0;
        for ( unsigned int childIdx = 0; childIdx < 8; ++childIdx )
            if ( ChildrenCells[childIdx] )    
                flops += 1;
        flopsM2M += flops;
        cell->addCost(flops);
        countM2M++;
    }


    void M2L(CellClass* const FRestrict cell,
             const CellClass* /*SourceCells*/[],
             const int /*positions*/[],
             const int size,
             const int /* TreeLevel */) override {
        FSize flops = 0;

		memset(_countExp, 0, sizeof(int) * 343);
        for (int idx = 0; idx < size; ++idx) {
            flops += 1;
        }
            
        flopsM2L += flops;
        cell->addCost(flops);
        countM2L++;
    }


    void L2L(const CellClass* const FRestrict /* not needed */,
             CellClass* FRestrict *const FRestrict ChildCells,
             const int /* TreeLevel*/) override {
        FSize flops = 0;
        FSize tmpCost = 0;
        for (unsigned int childIdx = 0; childIdx < 8; ++childIdx ) {
            if (ChildCells[childIdx]) {
                tmpCost = 1;
                flops += tmpCost;
                ChildCells[childIdx]->addCost(tmpCost);
            }            
        }

        flopsL2L += flops;
        countL2L++;
    }



    void L2P(const CellClass* const cell,
             ContainerClass* const TargetParticles) override {
        FSize tmpCost = TargetParticles->getNbParticles();
        flopsL2P += tmpCost;
        cell->addCost(tmpCost);
        countL2P++;
    }

    

    void P2P(const FTreeCoordinate& LeafCellCoordinate, // needed for periodic boundary conditions
             ContainerClass* const FRestrict TargetParticles,
             const ContainerClass* const FRestrict SourceParticles,
             ContainerClass* const /*NeighborSourceParticles*/[],
             const int /*positions*/[],
             const int size) override {
        FSize srcPartCount = SourceParticles->getNbParticles();
        FSize tgtPartCount = TargetParticles->getNbParticles();

        FSize tmpCost = srcPartCount * tgtPartCount;


        CellClass* cell = _tree->getCell(
            LeafCellCoordinate.getMortonIndex(),
            _treeHeight - 1);

        flopsP2P += tmpCost;
        cell->addNearCost(tmpCost);
        countP2P++;
    }

    void P2POuter(const FTreeCoordinate& LeafCellCoordinate, // needed for periodic boundary conditions
             ContainerClass* const FRestrict TargetParticles,
             ContainerClass* const /*NeighborSourceParticles*/[],
             const int /*positions*/[],
             const int size) override {
        FSize tmpCost = 0;

        CellClass* cell = _tree->getCell(
            LeafCellCoordinate.getMortonIndex(),
            _treeHeight - 1);

        flopsP2P += tmpCost;
        cell->addNearCost(tmpCost);
        countP2P++;
    }
};





#endif 

// [--END--]
