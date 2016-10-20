// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================

#ifndef FABSTRACTBALANCEALGORITHM_H
#define FABSTRACTBALANCEALGORITHM_H


/**
 * @author Cyrille Piacibello
 * @class FAbstractBalanceAlgorithm
 *
 * @brief This class provide the methods that are used to balance a
 * tree FMpiTreeBuilder::EqualizeAndFillTree
 */
class FAbstractBalanceAlgorithm{
public:
  virtual ~FAbstractBalanceAlgorithm(){
  }

  /**
   * @brief Give the right leaves (ie the min) of the interval that
   * will be handle by idxOfProc
   * @param numberOfLeaves Total number of leaves that exist.
   * @param numberOfPartPerLeaf Array of lenght numberOfLeaves containing the number of particles in each leaf
   * @param numberOfPart Number of particles in the whole field
   * @param idxOfLeaves Array of lenght numberOfLeaves containing the Morton Index of each Leaf
   * @param numberOfProc Number of MPI processus that will handle the Octree
   * @param idxOfProc Idx of the proc calling.
   */
  virtual FSize getRight(const FSize numberOfLeaves,
                         const int numberOfProc, const int idxOfProc) = 0;

  /**
   * @brief Give the Leaft leaves (ie the max) of the interval that
   * will be handle by idxOfProc
   * @param numberOfLeaves Total number of leaves that exist.
   * @param numberOfPartPerLeaf Array of lenght numberOfLeaves containing the number of particles in each leaf
   * @param numberOfPart Number of particles in the whole field
   * @param idxOfLeaves Array of lenght numberOfLeaves containing the Morton Index of each Leaf
   * @param numberOfProc Number of MPI processus that will handle the Octree
   * @param idxOfProc Idx of the proc calling.
   */
  virtual FSize getLeft(const FSize numberOfLeaves,
                        const int numberOfProc, const int idxOfProc) = 0;

};

#endif //FABSTRACTBALANCEALGORITHM_H
