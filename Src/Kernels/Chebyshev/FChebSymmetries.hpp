#ifndef FCHEBSYMMETRIES_HPP
#define FCHEBSYMMETRIES_HPP

#include <climits>

#include "../../Utils/FBlas.hpp"

#include "./FChebTensor.hpp"

/**
 * @author Matthias Messner (matthias.matthias@inria.fr)
 * Please read the license
 */

/**
 * @class FChebSymmetries
 *
 * The class @p FChebSymmetries exploits all symmetries
 */
template <int ORDER>
class FChebSymmetries
{
	enum {nnodes = TensorTraits<ORDER>::nnodes};

	// index permutations (j<i)(k<i)(k<j)
	unsigned int perms[8][3];

	// 48 global permutations (combinations of 8 quadrants and 6 cones respectively)
	unsigned int permutations[64][nnodes];

	unsigned int getQuadIdx(const int i, const int j, const int k) const
	{
		// find right quadrant index (if < 0 then 0, else 1)
		const int si = ((unsigned int)i >> (sizeof(int)*CHAR_BIT-1));
		const int sj = ((unsigned int)j >> (sizeof(int)*CHAR_BIT-2)) & 2;
		const int sk = ((unsigned int)k >> (sizeof(int)*CHAR_BIT-3)) & 4;
		return  (sk | sj | si);
	}


	public:
	FChebSymmetries()
	{
		// permutations for 8 quadrants
		unsigned int quads[8][nnodes];
		// permutations for 6 cones in quadrant (+++), 2 and 5 do not exist
		unsigned int cones[8][nnodes];

		// set quads and cones permutations
		unsigned int evn[ORDER]; unsigned int odd[ORDER];
		for (unsigned int o=0; o<ORDER; ++o) {evn[o] = o;	odd[o] = ORDER-1 - o;}

		for (unsigned int i=0; i<ORDER; ++i) {
			for (unsigned int j=0; j<ORDER; ++j) {
				for (unsigned int k=0; k<ORDER; ++k) {
					const unsigned int index = k*ORDER*ORDER + j*ORDER + i;

					// global axis parallel symmetries (8 quads) ///////////
					quads[0][index] = evn[k]*ORDER*ORDER + evn[j]*ORDER + evn[i];	// - - -
					quads[1][index] = evn[k]*ORDER*ORDER + evn[j]*ORDER + odd[i]; // - - +
					quads[2][index] = evn[k]*ORDER*ORDER + odd[j]*ORDER + evn[i]; // - + -
					quads[3][index] = evn[k]*ORDER*ORDER + odd[j]*ORDER + odd[i]; // - + +
					quads[4][index] = odd[k]*ORDER*ORDER + evn[j]*ORDER + evn[i]; // + - -
					quads[5][index] = odd[k]*ORDER*ORDER + evn[j]*ORDER + odd[i]; // + - +
					quads[6][index] = odd[k]*ORDER*ORDER + odd[j]*ORDER + evn[i]; // + + -
					quads[7][index] = odd[k]*ORDER*ORDER + odd[j]*ORDER + odd[i]; // + + +

					// diagonal symmetries (j>i)(k>i)(k>j) /////////////////////
					cones[0][index] = evn[k]*ORDER*ORDER + evn[j]*ORDER + evn[i]; // (0) 000
					cones[1][index] = evn[j]*ORDER*ORDER + evn[k]*ORDER + evn[i];	// (1) 001
					// cones[2] does not exist
					cones[3][index] = evn[j]*ORDER*ORDER + evn[i]*ORDER + evn[k];	// (3) 011
					cones[4][index] = evn[k]*ORDER*ORDER + evn[i]*ORDER + evn[j];	// (4) 100
					// cones[5] does not exist
					cones[6][index] = evn[i]*ORDER*ORDER + evn[k]*ORDER + evn[j];	// (6) 110
					cones[7][index] = evn[i]*ORDER*ORDER + evn[j]*ORDER + evn[k];	// (7) 111

				}
			}
		}

		// set 48 global permutations (combinations of 8 quadrants and 6 cones respectively)
		for (unsigned int q=0; q<8; ++q)
			for (unsigned int c=0; c<8; ++c)
				if (c!=2 && c!=5) // cone 2 and 5 do not exist
					for (unsigned int n=0; n<nnodes; ++n)
						permutations[q*8 + c][n] = cones[c][quads[q][n]];
		
		// permutation of interaction indices (already absolute value)
		perms[0][0] = 0; perms[0][1] = 1; perms[0][2] = 2;
		perms[1][0] = 0; perms[1][1] = 2; perms[1][2] = 1;
		// perms[2] does not exist
		perms[3][0] = 2; perms[3][1] = 0; perms[3][2] = 1;
		perms[4][0] = 1; perms[4][1] = 0; perms[4][2] = 2;
		// perms[5] does not exist
		perms[6][0] = 1; perms[6][1] = 2; perms[6][2] = 0;
		perms[7][0] = 2; perms[7][1] = 1; perms[7][2] = 0;
	}



	const unsigned int getPermutationArrayAndIndex(const int i, const int j, const int k,
																								 unsigned int permutation[nnodes]) const
	{
		// find right quadrant index (if < 0 then 0, else 1)
		const unsigned int qidx = getQuadIdx(i,j,k);

		// store absolut values of (i,j,k) in (u[0],u[1],u[2]) 
		const int imask = i >> (sizeof(int)*CHAR_BIT-1);
		const int jmask = j >> (sizeof(int)*CHAR_BIT-1);
		const int kmask = k >> (sizeof(int)*CHAR_BIT-1);
		const unsigned int u[3] = {(i+imask)^imask, (j+jmask)^jmask, (k+kmask)^kmask};

		// find right cone index
		const int q0 = (u[1]>u[0]) << 2;
		const int q1 = (u[2]>u[0]) << 1;
		const int q2 = (u[2]>u[1]);
		const int cidx = (q2 | q1 | q0);

		// set permutation array /////////////////////////////////////////
		for (unsigned int n=0; n<nnodes; ++n)	permutation[n] = permutations[qidx*8 + cidx][n];
		
		// set permutation index /////////////////////////////////////////
		return (u[perms[cidx][0]]+3)*7*7 + (u[perms[cidx][1]]+3)*7 + (u[perms[cidx][2]]+3);
	}
	


	static void permuteMatrix(const unsigned int perm[nnodes], FReal *const Matrix)
	{
		// allocate temporary memory for matrix
		FReal temp[nnodes*nnodes];
		// permute rows
		for (unsigned int r=0; r<nnodes; ++r)
			FBlas::copy(nnodes, Matrix+r, nnodes, temp+perm[r], nnodes);
		// permute columns
		for (unsigned int c=0; c<nnodes; ++c)
			FBlas::copy(nnodes, temp + c * nnodes, Matrix + perm[c] * nnodes);
	}

};


#endif


////		if (uj>ui) FMemUtils::swap(uj,ui);
////		if (uk>ui) FMemUtils::swap(uk,ui);
////		if (uk>uj) FMemUtils::swap(uk,uj);

////	const unsigned int *const getPermQuadrant(const int i, const int j, const int k) const
////	{
////		// if < 0 then 0, else 1
////		const int si = ((unsigned int)i >> (sizeof(int)*CHAR_BIT - 1));
////		const int sj = ((unsigned int)j >> (sizeof(int)*CHAR_BIT - 2)) & 2;
////		const int sk = ((unsigned int)k >> (sizeof(int)*CHAR_BIT - 3)) & 4;
////		const int qidx = sk | sj | si;
////		return quads[qidx];
////	}

////	const void permuteMatrix(const int i, const int j, const int k, FReal *const Matrix) const
////	{	permuteMatrix(getPermQuadrant(i, j, k), Matrix); }


////for (unsigned int n=0; n<nnodes; ++n)	permutation[n] = cones[cidx][quads[qidx][n]];
