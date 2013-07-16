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
#ifndef FTAYLORKERNEL_HPP
#define FTAYLORKERNEL_HPP

#include "../../Components/FAbstractKernels.hpp"
#include "../../Utils/FMemUtils.hpp"
#include "../../Utils/FDebug.hpp"

#include "../P2P/FP2P.hpp"

/**
 * @author Cyrille Piacibello 
 * @class FTaylorKernel 
 * 
 * @brief This kernel is an implementation of the different operators
 * needed to compute the Fast Multipole Method using Taylor Expansion
 * for the Far fields interaction.
 */


//TODO spécifier les arguments.
template< class CellClass, class ContainerClass, int P, int order>
class FTaylorKernel : public FAbstractKernels<CellClass,ContainerClass> {
  
private:
  //Size of the multipole and local vectors
  static const int SizeVector = ((P+1)*(P+2)*(P+3))*order/6;

  
  ////////////////////////////////////////////////////
  // Object Attributes
  ////////////////////////////////////////////////////
  const FReal boxWidth;               //< the box width at leaf level
  const int   treeHeight;             //< The height of the tree
  const FReal widthAtLeafLevel;       //< width of box at leaf level
  const FReal widthAtLeafLevelDiv2;   //< width of box at leaf leve div 2
  const FPoint boxCorner;             //< position of the box corner
  
  FReal factorials[2*P+1];             //< This contains the factorial until P
  FReal arrayDX[P+2],arrayDY[P+2],arrayDZ[P+2] ; //< Working arrays


  // For debugging purpose
  FILE * out;

  ////////////////////////////////////////////////////
  // Private method
  ////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////
  // Precomputation
  ///////////////////////////////////////////////////////

  /** Compute the factorial from 0 to P
   * Then the data is accessible in factorials array:
   * factorials[n] = n! with n <= P
   */
  void precomputeFactorials(){
    factorials[0] = 1.0;
    FReal fidx = 1.0;
    for(int idx = 1 ; idx <= 2*P ; ++idx, ++fidx){
      factorials[idx] = fidx * factorials[idx-1];
    }
  }


  /** Return the position of a leaf from its tree coordinate */
  FPoint getLeafCenter(const FTreeCoordinate coordinate) const {
    return FPoint(
		  FReal(coordinate.getX()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX(),
		  FReal(coordinate.getY()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getY(),
		  FReal(coordinate.getZ()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getZ());
  }

  /** 
   * @brief Return the position of the center of a cell from its tree
   *  coordinate 
   * @param FTreeCoordinate
   * @param inLevel the current level of Cell
   */
  FPoint getCellCenter(const FTreeCoordinate coordinate, int inLevel)
  {

    //Set the boxes width needed
    FReal widthAtCurrentLevel = widthAtLeafLevel*FReal(1 << (treeHeight-(inLevel+1)));   
    FReal widthAtCurrentLevelDiv2 = widthAtCurrentLevel/FReal(2);

    //Get the coordinate
    int a = coordinate.getX();
    int b = coordinate.getY();
    int c = coordinate.getZ();
    
    //Set the center real coordinates from box corner and widths.
    FReal X = boxCorner.getX() + FReal(a)*widthAtCurrentLevel + widthAtCurrentLevelDiv2;
    FReal Y = boxCorner.getY() + FReal(b)*widthAtCurrentLevel + widthAtCurrentLevelDiv2;
    FReal Z = boxCorner.getZ() + FReal(c)*widthAtCurrentLevel + widthAtCurrentLevelDiv2;
    
    FPoint cCenter = FPoint(X,Y,Z);

    //For debug purpose
    //printf("%f,%f,%f\n",cCenter.getX(),cCenter.getY(),cCenter.getZ());
    return cCenter;
  }


  /** 
   * @brief Incrementation of powers in Taylor expansion 
   * Result : ...,[2,0,0],[1,1,0],[1,0,1],[0,2,0]...  3-tuple are sorted 
   * by size then alphabetical order.
   */
  void incPowers(int * const FRestrict a, int *const FRestrict b, int *const FRestrict c)
  {
    int t = (*a)+(*b)+(*c);
    if(t==0)
      {a[0]=1;}
    else{ if(t==a[0])
	{a[0]--;  b[0]++;}
      else{ if(t==c[0])
	  {a[0]=t+1;  c[0]=0;}
	else{ if(b[0]!=0)
	    {b[0]--; c[0]++;}
	  else{ 
	    b[0]=c[0]+1;
	    a[0]--;
	    c[0]=0;
	  }
	}
      }
    }
  }
  
  /**
   * @brief Give the index of array from the corresponding 3-tuple
   * powers.
   */
  int powerToIdx(const int a,const int b,const int c)
  {
    int t,p = a+b+c;
    if(p==0)
      {return 0;}
    else
      {
	int res = p*(p+1)*(p+2)/6;
	t=p-a;
	res+=t*(t+1)/2+c;
	return res;
      }
  }
  
  /* Return the factorial of a number
   */
  FReal fact(const int a){
    if(a<0) {
      printf("fact :: Error factorial negative!! a=%d\n",a);
      return FReal(0);
    }
    FReal result = 1;
    for(int i = 1 ; i <= a ; ++i){
      result *= FReal(i);
    }
    return result;
  }

  /* Return the product of factorial of 3 numbers
   */
  FReal fact3int(const int a,const int b,const int c)
  {
    return ( factorials[a]*factorials[b]* factorials[c]) ; 
  }

  /* Return the combine of a paire of number
   */
  FReal combin(const int& a, const int& b){
    if(a-b<0)  {printf("combin :: Error combin negative!! a=%d b=%d\n",a,b); exit(-1) ;  }
        return  factorials[a]/  (factorials[b]* factorials[a-b]) ; 
  }


  /** @brief Init the derivative array for using of following formula
   * from "A cartesian tree-code for screened coulomb interactions"
   *
   *  @todo METTRE les fonctions pour intialiser la recurrence. \f$x_i\f$ ?? \f$x_i\f$ ?? 
   *  @todo LA formule ci-dessous n'utilise pas k! 
   */

  void initDerivative(const FReal & dx ,const FReal & dy ,const FReal & dz  ,   FReal * tab)
  {
    FReal R2 = dx*dx+dy*dy+dz*dz;
    printf("dx : %f dy : %f dz : %f\n",dx,dy,dz);
    tab[0]=FReal(1)/FMath::Sqrt(R2);   
    FReal R3 = tab[0]/(R2);
    tab[1]= -dx*R3;                 //Derivative in (1,0,0) il doit y avoir un -
    tab[2]= -dy*R3;                 //Derivative in (0,1,0)
    tab[3]= -dz*R3;                 //Derivative in (0,0,1)
    FReal R5 = R3/R2;
    tab[4] = FReal(3)*dx*dx*R5-R3;  //Derivative in (2,0,0)
    tab[5] = FReal(3)*dx*dy*R5;     //Derivative in (1,1,0)
    tab[6] = FReal(3)*dx*dz*R5;     //Derivative in (1,0,1)
    tab[7] = FReal(3)*dy*dy*R5-R3;  //Derivative in (0,2,0)
    tab[8] = FReal(3)*dy*dz*R5;     //Derivative in (0,1,1)
    tab[9] = FReal(3)*dz*dz*R5-R3;  //Derivative in (0,0,2)
     for(int c=0 ; c<=9 ; ++c){
       printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,0,0,0,c,tab[c]);
       }
  }
 
  /** @brief Compute and store the derivative for a given tuple.
   *  Derivative are used for the M2L
   *
   *\f[
   * \Psi_{\mathbf{k}}^{c} \left [\left |\mathbf{k}\right |\times \left |
   * \mathbf{x}_i-\mathbf{x}_c\right |^2  \right ]\				
   * = (2\times \left |{\mathbf{k}}\right |-1)
   * \sum_{j=0}^{3}\left [ k_j (x_{i_j}-x_{c_j})
   * \Psi_{\mathbf{k}-e_j,i}^{c}\right ]\ 
   * -(\left |\mathbf{k}\right |-1)   \sum_{j=0}^{3}\left
   * [ k_j(k_j-1) \Psi_{\mathbf{k}-2 e_j,i}^{c} \right]
   * \f]
   *  where    \f$ \mathbf{k} = (k_1,k_2,k_3) \f$ 
   */
  void computeFullDerivative( FReal  dx,  FReal  dy,  FReal  dz, // Distance from distant center to local center
			      FReal * yetComputed)
  {
   
    initDerivative(dx,dy,dz,yetComputed);
    FReal dist2 =  dx*dx+dy*dy+dz*dz;
    int idxTarget;                      //Index of current yetComputed entry
    int idxSrc1, idxSrc2, idxSrc3,      //Indexes of needed yetComputed entries
      idxSrc4, idxSrc5, idxSrc6;        
    int a=0,b=0,c=0;                    //Powers of expansions
    
    for(c=3 ; c<=2*P ; ++c){
      //Computation of derivatives Psi_{0,0,c}
      // |x-y|^2 * Psi_{0,0,c} + (2*c-1) * dz *Psi_{0,0,c-1} + (c-1)^2 * Psi_{0,0,c-2} = 0
      idxTarget = powerToIdx(0,0,c);
      idxSrc1 = powerToIdx(0,0,c-1);
      idxSrc2 = powerToIdx(0,0,c-2);
      yetComputed[idxTarget] = -(FReal(2*c-1)*dz*yetComputed[idxSrc1] + FReal((c-1)*(c-1))*yetComputed[idxSrc2])/dist2;
      printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,c,idxTarget,yetComputed[idxTarget]);
    }
    //printf(" Psi_{0,0,c} computed \n");
    b=1;
    for(c=2 ; c<=2*P-1 ; ++c){
      //Computation of derivatives Psi_{0,1,c}
      // |x-y|^2 * Psi_{0,1,c} + (2*c) * dz *Psi_{0,1,c-1} + c*(c-1) * Psi_{0,1,c-2} + dy*Psi_{0,0,c} = 0
      idxTarget = powerToIdx(0,1,c);
      idxSrc1 = powerToIdx(0,1,c-1);
      idxSrc2 = powerToIdx(0,1,c-2);
      idxSrc3 = powerToIdx(0,0,c);
      yetComputed[idxTarget] = -(FReal(2*c)*dz*yetComputed[idxSrc1] + FReal(c*(c-1))*yetComputed[idxSrc2]+ dy*yetComputed[idxSrc3])/dist2;
      printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,c,idxTarget,yetComputed[idxTarget]);
    }
    //printf(" Psi_{0,1,c} computed \n");
    b=2;
    for(c=1 ; c<= 2*P-b ; ++c){
      //Computation of derivatives Psi_{0,2,c}
      //|x-y|^2 * Psi_{0,2,c} + (2*c) * dz *Psi_{0,2,c-1} + (c*(c-1)) * Psi_{0,2,c-2} + 3*dy * Psi_{0,1,c} + Psi_{0,0,c}  = 0
      idxTarget = powerToIdx(0,2,c);
      idxSrc1 = powerToIdx(0,2,c-1);
      idxSrc2 = powerToIdx(0,2,c-2);
      idxSrc3 = powerToIdx(0,1,c);
      idxSrc4 = powerToIdx(0,0,c);
      yetComputed[idxTarget] = -(FReal(2*c)*dz*yetComputed[idxSrc1] + FReal(c*(c-1))*yetComputed[idxSrc2]
				 + FReal(3)*dy*yetComputed[idxSrc3] + yetComputed[idxSrc4])/dist2;
      printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,c,idxTarget,yetComputed[idxTarget]);
    }
    //printf(" Psi_{0,2,c} computed \n");
  
    for(b=3 ; b<= 2*P ; ++b){
      //Computation of derivatives Psi_{0,b,0}
      // |x-y|^2 * Psi_{0,b,0} + (2*b-1) * dy *Psi_{0,b-1,0} + (b-1)^2 * Psi_{0,b-2,c} = 0
      idxTarget = powerToIdx(0,b,0);
      idxSrc1 = powerToIdx(0,b-1,0);
      idxSrc2 = powerToIdx(0,b-2,0);
      yetComputed[idxTarget] = -(FReal(2*b-1)*dy*yetComputed[idxSrc1] + FReal((b-1)*(b-1))*yetComputed[idxSrc2])/dist2;
      printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,0,idxTarget,yetComputed[idxTarget]);
   
      for(c=1 ; c<= 2*P-b ; ++c) {
	//Computation of derivatives Psi_{0,b,c}
	//|x-y|^2*Psi_{0,b,c} + (2*c)*dz*Psi_{0,b,c-1} + (c*(c-1))*Psi_{0,b,c-2} + (2*b-1)*dy*Psi_{0,b-1,c} + (b-1)^2 * Psi_{0,b-2,c}  = 0
	idxTarget = powerToIdx(0,b,c);
	idxSrc1 = powerToIdx(0,b,c-1);
	idxSrc2 = powerToIdx(0,b,c-2);
	idxSrc3 = powerToIdx(0,b-1,c);
	idxSrc4 = powerToIdx(0,b-2,c);
	yetComputed[idxTarget] = -(FReal(2*c)*dz*yetComputed[idxSrc1] + FReal(c*(c-1))*yetComputed[idxSrc2]
				   + FReal(2*b-1)*dy*yetComputed[idxSrc3] + FReal((b-1)*(b-1))*yetComputed[idxSrc4])/dist2;
	printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,c,idxTarget,yetComputed[idxTarget]);
      }
    }
    //printf(" Psi_{0,b,c} computed \n");
  
    a=1;
    b=0;
    for(c=2 ; c<= 2*P-1 ; ++c){
      //Computation of derivatives Psi_{1,0,c}
      //|x-y|^2 * Psi_{1,0,c} + (2*c)*dz*Psi_{1,0,c-1} + c*(c-1)*Psi_{1,0,c-2} + dx*Psi_{0,0,c}
      idxTarget = powerToIdx(1,0,c);
      idxSrc1 = powerToIdx(1,0,c-1);
      idxSrc2 = powerToIdx(1,0,c-2);
      idxSrc3 = powerToIdx(0,0,c);
      yetComputed[idxTarget] = -(FReal(2*c)*dz*yetComputed[idxSrc1] + FReal(c*(c-1))*yetComputed[idxSrc2] + dx*yetComputed[idxSrc3])/dist2;
      printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,c,idxTarget,yetComputed[idxTarget]);
    }
    //printf(" Psi_{1,0,c} computed \n");
    b=1;
    //Computation of derivatives Psi_{1,1,1}
    //|x-y|^2 * Psi_{1,1,1} + 2*dz*Psi_{1,1,0} + 2*dy*Psi_{1,0,1} + dx*Psi_{0,1,1}
    idxTarget = powerToIdx(1,1,1);
    idxSrc1 = powerToIdx(1,1,0);
    idxSrc2 = powerToIdx(1,0,1);
    idxSrc3 = powerToIdx(0,1,1);
    yetComputed[idxTarget] = -(FReal(2)*dz*yetComputed[idxSrc1] + FReal(2)*dy*yetComputed[idxSrc2] + dx*yetComputed[idxSrc3])/dist2;
    printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,1,1,1,idxTarget,yetComputed[idxTarget]);
    for(c=2 ; c<= 2*P-2 ; ++c){
      //Computation of derivatives Psi_{1,1,c}
      //|x-y|^2 * Psi_{1,1,c} + (2*c)*dz*Psi_{1,1,c-1} + c*(c-1)*Psi_{1,1,c-2} + 2*dy*Psi_{1,0,c} + dx*Psi_{0,1,c}
      idxTarget = powerToIdx(1,1,c);
      idxSrc1 = powerToIdx(1,1,c-1);
      idxSrc2 = powerToIdx(1,1,c-2);
      idxSrc3 = powerToIdx(1,0,c);
      idxSrc4 = powerToIdx(0,1,c);
      yetComputed[idxTarget] = -(FReal(2*c)*dz*yetComputed[idxSrc1] + FReal(c*(c-1))*yetComputed[idxSrc2] 
				 + FReal(2)*dy*yetComputed[idxSrc3]+ dx*yetComputed[idxSrc4])/dist2;
      printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,c,idxTarget,yetComputed[idxTarget]);
    }
    //printf(" Psi_{1,1,c} computed \n");
  
    for(b=2 ; b<= 2*P-a ; ++b){
      for(c=0 ; c<= 2*P-b-1 ; ++c){
	//Computation of derivatives Psi_{1,b,c}
	//|x-y|^2 * Psi_{1,b,c} + (2*b)*dy*Psi_{1,b-1,c} + b*(b-1)*Psi_{1,b-2,c} + (2*c)*dz*Psi_{1,b,c-1} + c*(c-1)*Psi_{1,b,c-2} + dx*Psi_{0,b,c}
	idxTarget = powerToIdx(1,b,c);
	idxSrc1 = powerToIdx(1,b-1,c);
	idxSrc2 = powerToIdx(1,b-2,c);
	idxSrc3 = powerToIdx(1,b,c-1);
	idxSrc4 = powerToIdx(1,b,c-2);
	idxSrc5 = powerToIdx(0,b,c);
	yetComputed[idxTarget] = -(FReal(2*b)*dy*yetComputed[idxSrc1] + FReal(b*(b-1))*yetComputed[idxSrc2] 
				   + FReal(2*c)*dz*yetComputed[idxSrc3]+ FReal(c*(c-1))*yetComputed[idxSrc4]
				   + dx*yetComputed[idxSrc5])/dist2;
	printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,c,idxTarget,yetComputed[idxTarget]);
      }
    }
    //printf(" Psi_{1,b,c} computed \n");
  
    for(a=2 ; a<=2*P ; ++a){
      //Computation of derivatives Psi_{a,0,0}
      // |x-y|^2 * Psi_{a,0,0} + (2*a-1) * dx *Psi_{a-1,0,0} + (a-1)^2 * Psi_{a-2,0,0} = 0
      idxTarget = powerToIdx(a,0,0);
      idxSrc1 = powerToIdx(a-1,0,0);
      idxSrc2 = powerToIdx(a-2,0,0);
      yetComputed[idxTarget] = -(FReal(2*a-1)*dx*yetComputed[idxSrc1] + FReal((a-1)*(a-1))*yetComputed[idxSrc2])/dist2;
      printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,0,0,idxTarget,yetComputed[idxTarget]);
      if(a <= 2*P-1){
	//Computation of derivatives Psi_{a,0,1}
	// |x-y|^2 * Psi_{a,0,1} + 2*dz*Psi_{a,0,0} + (2*a-1)*dx*Psi_{a-1,0,1} + (a-1)^2*Psi_{a-2,0,1} = 0
	idxSrc1 = idxTarget;
	idxTarget = powerToIdx(a,0,1);
	idxSrc2 = powerToIdx(a-1,0,1);
	idxSrc3 = powerToIdx(a-2,0,1);
	yetComputed[idxTarget] = -(FReal(2)*dz*yetComputed[idxSrc1] + FReal(2*a-1)*dx*yetComputed[idxSrc2] + FReal((a-1)*(a-1))*yetComputed[idxSrc3])/dist2;
	printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,0,1,idxTarget,yetComputed[idxTarget]);
	//Computation of derivatives Psi_{a,1,0}
	// |x-y|^2 * Psi_{a,1,0} + 2*dy*Psi_{a,0,0} + (2*a-1)*dx*Psi_{a-1,1,0} + (a-1)^2*Psi_{a-2,1,0} = 0
	idxTarget = powerToIdx(a,1,0);
	idxSrc2 = powerToIdx(a-1,1,0);
	idxSrc3 = powerToIdx(a-2,1,0);
	yetComputed[idxTarget] = -(FReal(2)*dy*yetComputed[idxSrc1] + FReal(2*a-1)*dx*yetComputed[idxSrc2] + FReal((a-1)*(a-1))*yetComputed[idxSrc3])/dist2;
	printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,1,0,idxTarget,yetComputed[idxTarget]);
	if(a <= 2*P-2){
	  b=0;
	  for(c=2 ; c <= 2*P-a ; ++c){
	    //Computation of derivatives Psi_{a,0,c}
	    // |x-y|^2 * Psi_{a,0,c} + 2*c*dz*Psi_{a,0,c-1} + c*(c-1)*Psi_{a,0,c-2} + (2*a-1)*dx*Psi_{a-1,0,c} + (a-1)^2*Psi_{a-2,0,c} = 0
	    idxTarget = powerToIdx(a,0,c);
	    idxSrc1 = powerToIdx(a,0,c-1);
	    idxSrc2 = powerToIdx(a,0,c-2);
	    idxSrc3 = powerToIdx(a-1,0,c);
	    idxSrc4 = powerToIdx(a-2,0,c);
	    yetComputed[idxTarget] = -(FReal(2*c)*dz*yetComputed[idxSrc1] + FReal(c*(c-1))*yetComputed[idxSrc2] 
				       + FReal(2*a-1)*dx*yetComputed[idxSrc3] + FReal((a-1)*(a-1))*yetComputed[idxSrc4])/dist2;
	    printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,0,c,idxTarget,yetComputed[idxTarget]);
	  }
	  b=1;
	  for(c=1 ; c <= 2*P-a-1 ; ++c){
	    //Computation of derivatives Psi_{a,1,c}
	    // |x-y|^2 * Psi_{a,1,c} + 2*c*dz*Psi_{a,1,c-1} + c*(c-1)*Psi_{a,1,c-2} + 2*a*dx*Psi_{a-1,1,c} + a*(a-1)*Psi_{a-2,1,c} + dy*Psi_{a,0,c}= 0
	    idxTarget = powerToIdx(a,1,c);
	    idxSrc1 = powerToIdx(a,1,c-1);
	    idxSrc2 = powerToIdx(a,1,c-2);
	    idxSrc3 = powerToIdx(a-1,1,c);
	    idxSrc4 = powerToIdx(a-2,1,c);
	    idxSrc5 = powerToIdx(a,0,c);
	    yetComputed[idxTarget] = -(FReal(2*c)*dz*yetComputed[idxSrc1] + FReal(c*(c-1))*yetComputed[idxSrc2] 
				       + FReal(2*a)*dx*yetComputed[idxSrc3] + FReal(a*(a-1))*yetComputed[idxSrc4]
				       + dy*yetComputed[idxSrc5])/dist2;
	    printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,1,c,idxTarget,yetComputed[idxTarget]);
	  }
	
	  for(b=2 ; b <= 2*P-a ; ++b){
	    //Computation of derivatives Psi_{a,b,0}
	    // |x-y|^2 * Psi_{a,b,0} + 2*b*dy*Psi_{a,b-1,0} + b*(b-1)*Psi_{a,b-2,0} + (2*a-1)*dx*Psi_{a-1,b,0} + (a-1)^2*Psi_{a-2,b,0} = 0
	    idxTarget = powerToIdx(a,b,0);
	    idxSrc1 = powerToIdx(a,b-1,0);
	    idxSrc2 = powerToIdx(a,b-2,0);
	    idxSrc3 = powerToIdx(a-1,b,0);
	    idxSrc4 = powerToIdx(a-2,b,0);
	    yetComputed[idxTarget] = -(FReal(2*b)*dy*yetComputed[idxSrc1] + FReal(b*(b-1))*yetComputed[idxSrc2]
				       + FReal(2*a-1)*dx*yetComputed[idxSrc3] + FReal((a-1)*(a-1))*yetComputed[idxSrc4])/dist2;
	    printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,0,idxTarget,yetComputed[idxTarget]);
	  
	    if(a+b < 2*P){
	      //Computation of derivatives Psi_{a,b,1}
	      // |x-y|^2 * Psi_{a,b,1} + 2*b*dy*Psi_{a,b-1,1} + b*(b-1)*Psi_{a,b-2,1} + 2*a*dx*Psi_{a-1,b,1} + a*(a-1)*Psi_{a-2,b,1} + dz*Psi_{a,b,0}= 0
	      idxTarget = powerToIdx(a,b,1);
	      idxSrc1 = powerToIdx(a,b-1,1);
	      idxSrc2 = powerToIdx(a,b-2,1);
	      idxSrc3 = powerToIdx(a-1,b,1);
	      idxSrc4 = powerToIdx(a-2,b,1);
	      idxSrc5 = powerToIdx(a,b,0);
	      yetComputed[idxTarget] = -(FReal(2*b)*dy*yetComputed[idxSrc1] + FReal(b*(b-1))*yetComputed[idxSrc2]
					 + FReal(2*a)*dx*yetComputed[idxSrc3] + FReal(a*(a-1))*yetComputed[idxSrc4]
					 + dz*yetComputed[idxSrc5])/dist2;
	      printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,1,idxTarget,yetComputed[idxTarget]);
	    }
	    for(c=2 ; c <= 2*P-b-a ; ++c){
	      //Computation of derivatives Psi_{a,b,c} with a >= 2
	      // |x-y|^2*Psi_{a,b,c} + (2*a-1)*dx*Psi_{a-1,b,c} + a*(a-2)*Psi_{a-2,b,c} + 2*b*dy*Psi_{a,b-1,c} + b*(b-1)*Psi_{a,b-2,c} + 2*c= 0
	      idxTarget = powerToIdx(a,b,c);
	      idxSrc1 = powerToIdx(a-1,b,c);
	      idxSrc2 = powerToIdx(a,b-1,c);
	      idxSrc3 = powerToIdx(a,b,c-1);
	      idxSrc4 = powerToIdx(a-2,b,c);
	      idxSrc5 = powerToIdx(a,b-2,c);
	      idxSrc6 = powerToIdx(a,b,c-2);
	      yetComputed[idxTarget] = -(FReal(2*a-1)*dx*yetComputed[idxSrc1] + FReal((a-1)*(a-1))*yetComputed[idxSrc4]
					 + FReal(2*b)*dy*yetComputed[idxSrc2] + FReal(b*(b-1))*yetComputed[idxSrc5]
					 + FReal(2*c)*dz*yetComputed[idxSrc3] + FReal(c*(c-1))*yetComputed[idxSrc6])/dist2;
	    
	      printf("just computed %f, a=%d, b=%d, c=%d target: %d %f\n",dx,a,b,c,idxTarget,yetComputed[idxTarget]);
	    }
	  }
	}
      }
    }
    //printf(" Psi_{a,b,c} computed \n");
  }


  /////////////////////////////////
  ///////// Public Methods ////////
  /////////////////////////////////

public:
  
  /*Constructor, need system information*/
  FTaylorKernel(const int inTreeHeight, const FReal inBoxWidth, const FPoint& inBoxCenter) :
    boxWidth(inBoxWidth),
    treeHeight(inTreeHeight),
    widthAtLeafLevel(inBoxWidth/FReal(1 << (inTreeHeight-1))),
    widthAtLeafLevelDiv2(widthAtLeafLevel/2),
    boxCorner(inBoxCenter.getX()-(inBoxWidth/2),inBoxCenter.getY()-(inBoxWidth/2),inBoxCenter.getZ()-(inBoxWidth/2))
  {
    this->precomputeFactorials() ;
  }
  
  /* Default destructor
   */
  virtual ~FTaylorKernel(){
    fclose(out);
  }

  /**P2M 
   * @brief Fill the Multipole with the field created by the cell
   * particles.
   *  
   * Formula :
   * \f[
   *   M_{k} = \sum_{j=0}^{N}{ q_j * \frac{|k|!}{k! k!} (x_c-x_j)^{k}}
   * \f]
   * where \f$x_c\f$ is the centre of the cell and \f$x_j\f$ the \f$j^{th}\f$ particles and \f$q_j\f$ its charge and  \f$N\f$ the particle number.
   */
  void P2M(CellClass* const pole, 
	   const ContainerClass* const particles)
  { 
    out = fopen("./res_3.data","a+");


    //Variables computed for each power of Multipole
    int a,b,c ;
    FReal facto, coeff; 
    //Copying cell center position once and for all
    const FPoint& cellCenter = getLeafCenter(pole->getCoordinate());
    printf("P2M :: pole : X: %f, Y: %f, Z:%f  \n",cellCenter.getX(),cellCenter.getY(),cellCenter.getZ());
    FReal * multipole = pole->getMultipole();
    FMemUtils::memset(multipole,0,SizeVector*FReal(0.0));
    //    
    // Iterator over Particles
    //
    int nbPart = particles->getNbParticles();
    const FReal* const * positions = particles->getPositions();
    const FReal* posX = positions[0];
    const FReal* posY = positions[1];
    const FReal* posZ = positions[2];
    
    const FReal* phyValue = particles->getPhysicalValues();
    //
    // Iterating over Particles
    //
    FReal xc = cellCenter.getX(), yc = cellCenter.getY(), zc = cellCenter.getZ() ;
    for(int idPart=0 ; idPart<nbPart ; ++idPart){
      printf("P2M :: part : X: %f, Y: %f, Z:%f   Q %f\n", posX[idPart] , posY[idPart] , posZ[idPart] ,phyValue[idPart]);

      // compute the distance to the centre	  
      FReal dx =  xc - posX[idPart] ;
      FReal dy =  yc - posY[idPart] ;
      FReal dz =  zc - posZ[idPart] ;

      // Precompute the an arrays of dx^i
      arrayDX[0] = 1.0 ;
      arrayDY[0] = 1.0 ;
      arrayDZ[0] = 1.0 ;
      for (int i = 1 ; i <= P ; ++i) 	{
	arrayDX[i] = dx * arrayDX[i-1] ;
	arrayDY[i] = dy * arrayDY[i-1] ; 
	arrayDZ[i] = dz * arrayDZ[i-1] ;
	printf("arrayD? ,i : %d, locForce : %f  %f  %f\n",i-1, arrayDX[i-1], arrayDY[i-1], arrayDZ[i-1] );
      }
      printf("arrayD? ,i : %d, locForce : %f  %f  %f\n",P, arrayDX[P], arrayDY[P], arrayDZ[P] );

      //
      //Iterating over MutlipoleVector
      //
      a = 0, b = 0, c = 0;
      for(int i=0 ; i<SizeVector ; ++i)
	{
	  //
	  // update needed values
	  //
	  facto = fact3int(a,b,c);
	  coeff = factorials[a+b+c]/(facto*facto);
	  //
	  //Computation
	  //
	  multipole[i] += coeff*phyValue[idPart]*arrayDX[a]*arrayDY[b]*arrayDZ[c];
	  
	  incPowers(&a,&b,&c);       //inc powers of expansion
	}  // end loop on multipoles
    }  // end loop on particles
    // Print the multipoles 
    //    if(xc == FReal(3)){
      a = b = c = 0; 
      for(int i=0 ; i<SizeVector ; ++i)
	{
	  fprintf(stdout," P2M :: cell : X = %f, Y = %f, Z = %f, %d = (%d,%d,%d) --> coeff %f   M= %f\n ",
		  cellCenter.getX(),cellCenter.getY(),cellCenter.getZ(),a+b+c,a,b,c,factorials[a+b+c]/fact3int(a,b,c),multipole[i]);
	  incPowers(&a,&b,&c);   
	} 	
      //    }
    std::cout << std::endl;
  //   for(int l=0 , idx = 0; l<= P ; ++l) // length of i + j + k = l
  //     {    
  //   	for( c=0 ; c <= l ; ++c)  
  //   	  {
  //   	    for( b = 0 ; b<= l-c ; ++b)
  //   	      {
  //   		for( a = l-c-b ;  a+b+c==l; --a, ++idx)
  //   		  {
  //   		    std::cout << "P2M>> "<< idx << " = (i,j,k) = ("<< a << " , " <<b << " , " << c << " ) " <<std::endl;
  //   		  } 
  //   	      } 
  //   	  } 
  //     } 
  // } 
    
  }
  

  /**
   * @brief Fill the parent multipole with the 8 values of child multipoles
   * 
   *
   */
  void M2M(CellClass* const FRestrict pole, 
	   const CellClass*const FRestrict *const FRestrict child, 
	   const int inLevel)
  {
    exit(0);
    printf("M2M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    //Powers of expansions
    int a=0,b=0,c=0;
    
    //Indexes of powers
    int idx_a,idx_b,idx_c;

    //Distance from current child to parent
    // FReal boxSize = (boxWidth/FReal(1 << inLevel));
    FReal dx = 0.0;
    FReal dy = 0.0;
    FReal dz = 0.0;
    //Center point of parent cell
    const FPoint& cellCenter = getLeafCenter(pole->getCoordinate());
    printf("M2M :: pole_target : X: %f, Y: %f, Z:%f\n",cellCenter.getX(),cellCenter.getY(),cellCenter.getZ());
    FReal * mult = pole->getMultipole();
    FMemUtils::memset(pole,FReal(0.0),SizeVector*FReal(0.0));
    
    //Iteration over the eight children
    int idxChild;
    FReal coeff;
    for(idxChild=0 ; idxChild<8 ; ++idxChild)
      {
	if(child[idxChild]){
	  
	  const FPoint& childCenter = getLeafCenter(child[idxChild]->getCoordinate());
	  const FReal * const multChild = child[idxChild]->getMultipole();

	  //Set the distance between centers of cells
	  dx = childCenter.getX()-cellCenter.getX();
	  dy = childCenter.getY()-cellCenter.getY();
	  dz = childCenter.getZ()-cellCenter.getZ();

	  // dz = ((FReal )(2*(1 & idxChild)-1))*boxSize;
	  // dy = ((FReal )(2*((1 << 1) & idxChild)-1))*boxSize;
	  // dx = ((FReal )(2*((1 << 2) & idxChild)-1))*boxSize;
	  // printf("Distances dans le M2M : %f %f %f boxSize : %f \n",dx,dy,dz,boxSize);
	  
	  a=0;
	  b=0;
	  c=0;
	  //Iteration over parent multipole array
	  for(int idxMult = 0 ; idxMult<SizeVector ; idxMult++)
	    {
	      FReal value = mult[idxMult];
	      
	      int idMultiChild;
	      //Iteration over the powers to find the cell multipole
	      //involved in the computation of the parent multipole
	      for(idx_a=0 ; idx_a <= a ; ++idx_a){
		for(idx_b=0 ; idx_b <= b ; ++idx_b){
		  for(idx_c=0 ; idx_c <= c ; ++idx_c){
		    //Computation
		    //Child multipole involved
		    idMultiChild=powerToIdx(a-idx_a,b-idx_b,c-idx_c);
		    coeff = fact3int(a-idx_a,b-idx_b,c-idx_c)/(fact(a-idx_a+b-idx_b+c-idx_c)*fact3int(idx_a,idx_b,idx_c));
		    
		    value+=multChild[idMultiChild]*FMath::pow(dx,idx_a)*FMath::pow(dy,idx_b)*FMath::pow(dz,idx_c)*coeff;
		  }
		}
	      }
	      //	      printf("M2M :: cell : X = %f, Y = %f, Z = %f,  %d,%d,%d --> %f\n",
	      //     cellCenter.getX(),cellCenter.getY(),cellCenter.getZ(),a,b,c,value);

	      mult[idxMult] += value;
	      incPowers(&a,&b,&c);
	    }
	}
      }
  }
  
  /**
   *@brief Convert the multipole expansion into local expansion The
   * operator do not use symmetries.
   *
   * Formula : \f[ L_{\mathbf{n}}^{c} = \frac{|n|!}{n! n!}
   * \sum_{\mathbf{k}=0}^{p} \left [ M_\mathbf{k}^c \,
   * \Psi_{\mathbf{,n+k}}( \mathbf{x}_c^{target})\right ] \f]
   * and \f[ \Psi_{\mathbf{,i}}^{c}(\mathbf{x}) =
   * \frac{\partial^i}{\partial x^i} \frac{1}{|x-x_c^{src}|} =  \frac{\partial^{i_1}}{\partial x_1^{i_1}} \frac{\partial^{i_2}}{\partial x_2^{i_2}} \frac{\partial^{i_3}}{\partial x_3^{i_3}} \frac{1}{|x-x_c^{src}|}\f] 
   *
   * Where \f$x_c^{src}\f$ is the centre of the cell where the
   * multiplole are considered,\f$x_c^{target}\f$ is the centre of the
   * current cell. The cell where we compute the local expansion.
   *
   */
  void M2L(CellClass* const FRestrict local,             // Target cell
	   const CellClass* distantNeighbors[343],       // Sources to be read     
	   const int /*size*/, const int inLevel)
  {
    printf("M2L\n");
    //Iteration over distantNeighbors
    int idxNeigh;
    int sizeDerivative = (2*P+1)*(P+1)*(2*P+3)/3; 
    

    FPoint locCenter = getCellCenter(local->getCoordinate(),inLevel);
    if(locCenter.getX() == FReal(-3)){
      fprintf(out,"M2l :: pole_target : X: %f, Y: %f, Z:%f\n",locCenter.getX(),locCenter.getY(),locCenter.getZ());
    }
    FReal * iterLocal = local->getLocal();
    FMemUtils::memset(iterLocal,0,SizeVector*sizeof(FReal(0.0)));
    FReal yetComputed[sizeDerivative];
    
    for(idxNeigh=0 ; idxNeigh<343 ; ++idxNeigh){

      //Need to test if current neighbor is one of the interaction list
      if(distantNeighbors[idxNeigh]){
	//Derivatives are computed iteratively
	
	FMemUtils::memset(yetComputed,0,sizeDerivative*sizeof(FReal(0.0)));
	// 
	// Compute derivatives on  locCenter - curDistCenter
	//                           target       source
	FPoint curDistCenter = getCellCenter(distantNeighbors[idxNeigh]->getCoordinate(),inLevel);
	FReal dx = locCenter.getX()-curDistCenter.getX();
	FReal dy = locCenter.getY()-curDistCenter.getY();
	FReal dz = locCenter.getZ()-curDistCenter.getZ();

	
	//Computation of all the derivatives needed
	computeFullDerivative(dx,dy,dz,yetComputed);

	//Iteration over Multipole / Local
	int al=0,bl=0,cl=0; // For local array
	int am,bm,cm;       // For distant array
	//	
	//Iterating over local array : n
	for(int i=0 ; i< SizeVector ; ++i){
	  FReal fctl   = fact3int(al,bl,cl);
	  FReal coeffL = factorials[al+bl+cl]/(fctl*fctl);
	  //
	  //Iterator over multipole array 
	  const FReal * multipole = distantNeighbors[idxNeigh]->getMultipole();
	  
	  //For debugging purposes
	  //FReal multipole[SizeVector];
	  //FMemUtils::memset(multipole,0,SizeVector*sizeof(FReal(0.0)));
	  //multipole[3]=FReal(1);

	  FReal tmp = 0.0 ;
	  //Iterating over multipole array : k
	  //  Loc(al,bl,cl) = N(al,bl,cl)/((al,bl,cl)!*(al,bl,cl)!) sum_(am,bm,cm) Psi[am+al,bm+bl,cm+cl] * M[am,bm,cm]  
	  //
	  am=0;	  bm=0;  cm=0;
	  printf("al= %d, bl=%d, cl=%d ==>    i =%d \n",al,bl,cl,i);
	  for(int j = 0 ; j < SizeVector ; ++j){ //corresponding powers am,bm,cm
	    int idxPsi = powerToIdx(al+am,bl+bm,cl+cm);
	    tmp += yetComputed[idxPsi]*multipole[j];
	    

	    printf(" j= %d, am=%d, bm=%d, cm=%d,, aml=%d, bml=%d, cml=%d, psi[%d]=%f\n",j,am,bm,cm,am+al,bm+bl,cm+cl,powerToIdx(al+am,bl+bm,cl+cm),yetComputed[powerToIdx(al+am,bl+bm,cl+cm)]);

	    //updating a,b,c
	    incPowers(&am,&bm,&cm);
	  }
	  iterLocal[i] = tmp*coeffL ;
	  incPowers(&al,&bl,&cl);
	}
	// For Debugging ..........................................................
	//if(locCenter.getX() == FReal(-3)){
	  int x=0,y=0,z=0;
	  FReal tot = FReal(0);
	  for(int dby=0 ; dby<SizeVector ; dby++)
	    {	
	      //tot += yetComputed[dby];
	      fprintf(stdout,"M2L :: source %f, %d,%d,%d ==> %f\n",curDistCenter.getX(),x,y,z,iterLocal[dby]);
	      incPowers(&x,&y,&z);
	    }
	  fprintf(out,"tot : %f\n\n\n",tot);
	  x = y = z = 0;
	  for(int dby=0 ; dby<sizeDerivative ; dby++)
	    {	
	      tot+=yetComputed[dby];
	      //fprintf(stdout,"M2L :: source %f, (%d,%d,%d) ==> derive : %f\n",curDistCenter.getX(),x,y,z,yetComputed[dby]);
	      incPowers(&x,&y,&z);
	    }
	  //}
      }
    }
  }

  
  /**
   *@brief Translate the local expansion of parent cell to child cell
   * \f[
   * L_{son} = \sum_{i=k}^{l} L_{i} (x_{son}-x_{father})^{i-k} \binom{i}{k}
   *\f]
   */
  void L2L(const CellClass* const FRestrict local, 
	   CellClass* FRestrict * const FRestrict child, 
	   const int /*inLevel*/)
  {
    exit(0);
    FPoint locCenter = getLeafCenter(local->getCoordinate());
    FReal dx;
    FReal dy;
    FReal dz;
    int a, b, c;
    // For all children
    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
      // if child exists
      if(child[idxChild]){
	a=0;
	b=0;
	c=0;
	FPoint childCenter = getLeafCenter(child[idxChild]->getCoordinate());
	//Set the distance between centers of cells
	dx = childCenter.getX()-locCenter.getX();
	dy = childCenter.getY()-locCenter.getY();
	dz = childCenter.getZ()-locCenter.getZ();
	//iterator over child's local expansion (to be filled)
	for(int k=0 ; k<SizeVector ; k++){
	  
	  //Iterator over parent's local array
	  for(int i=k ; i<SizeVector ; i++){
	    (child[idxChild]->getLocal())[k] = (local->getLocal())[i]*FMath::pow(dx,a)*FMath::pow(dy,b)*FMath::pow(dz,c)*combin(i,k);
	  }
	  incPowers(&a,&b,&c);
	}	
      }
    }
  }
 
  
  /**L2P
   *@brief Apply on the particles the force computed from the local expansion to all particles in the cell
   *
   *  
   * Formula :
   * \f[
   *   Potential = \sum_{j=0}^{nb_{particles}}{q_j \sum_{k=0}^{P}{ L_k * (x_j-x_c)^{k}}}
   * \f]
   *
   * where \f$x_c\f$ is the centre of the local cell and \f$x_j\f$ the
   * \f$j^{th}\f$ particles and \f$q_j\f$ its charge.
   */
    void L2P(const CellClass* const local, 
	     ContainerClass* const particles)
    {
      
      FPoint locCenter = getLeafCenter(local->getCoordinate());
      //Iterator over particles
      int nbPart = particles->getNbParticles();
      
      //	
      //Iteration over Local array
      //
      const FReal * iterLocal = local->getLocal();

      const FReal * const * positions = particles->getPositions();
      const FReal * posX = positions[0];
      const FReal * posY = positions[1];
      const FReal * posZ = positions[2];
      
      FReal * const forceX = particles->getForcesX();
      FReal * const forceY = particles->getForcesY();
      FReal * const forceZ = particles->getForcesZ();
      //
      FReal * const targetsPotentials = particles->getPotentials();
      
      printf("L2P : Cell : %f, fx = %f, fy = %f, fz = %f\n\n",locCenter.getX(),forceX[0],forceY[0],forceZ[0]);

      FReal * const phyValues = particles->getPhysicalValues();

      //Iteration over particles
      for(int i=0 ; i<nbPart ; ++i){
	//	
	FReal dx =  posX[i] - locCenter.getX();
	FReal dy =  posY[i] - locCenter.getY();
	FReal dz =  posZ[i] - locCenter.getZ();
	printf("L2P: dx = %f, dy = %f, dz = %f\n",dx,dy,dz);
	//
	// Precompute an arrays of Array[i] = dx^(i-1)
	arrayDX[0] = 0.0 ; 
	arrayDY[0] = 0.0 ;
	arrayDZ[0] = 0.0 ;
	
	arrayDX[1] = 1.0 ;
	arrayDY[1] = 1.0 ;
	arrayDZ[1] = 1.0 ;
	std::cout<< std::endl;
	printf("  ,(dx,dy,dz)  %f  %f  %f\n" ,dx, dy, dz);
	for (int d = 2 ; d <= P+1 ; ++d){ //Array is staggered : Array[i] = dx^(i-1)
	  arrayDX[d] = dx * arrayDX[d-1] ;
	  arrayDY[d] = dy * arrayDY[d-1] ; 
	  arrayDZ[d] = dz * arrayDZ[d-1] ;
	  printf("arrayD? ,j : %d, dx^j : %f  %f  %f\n",d-1, arrayDX[d-1], arrayDY[d-1], arrayDZ[d-1] );

	}
	std::cout<< std::endl;	
	FReal partPhyValue = phyValues[i]; 
	//
	FReal  locPot = 0.0, locForceX = 0.0, locForceY = 0.0, locForceZ = 0.0 ;
	int a=0,b=0,c=0;
	for(int j=0 ; j<SizeVector ; ++j){
	  FReal locForce     = iterLocal[j];
	  // compute the potential
	  locPot += iterLocal[j]*arrayDX[a+1]*arrayDY[b+1]*arrayDZ[c+1];
	  //Application of forces
	  locForceX += FReal(a)*locForce*arrayDX[a]*arrayDY[b+1]*arrayDZ[c+1];
	  locForceY += FReal(b)*locForce*arrayDX[a+1]*arrayDY[b]*arrayDZ[c+1];
	  locForceZ += FReal(c)*locForce*arrayDX[a+1]*arrayDY[b+1]*arrayDZ[c];
	  //
	  printf("  force X : %d,%d,%d,j : %d, locForceZ : %f  L_j  %f  %f  %f  %f \n",a,b,c,j,locForceZ,locForce,arrayDX[a],arrayDY[b+1],arrayDZ[c+1]);
	  incPowers(&a,&b,&c);
	}
	targetsPotentials[i]  = partPhyValue*locPot ;
	forceX[i]            = partPhyValue*locForceX ;
	forceY[i]            = partPhyValue*locForceY ;
	forceZ[i]            = partPhyValue*locForceZ ;
      printf("  Part %d,  Pot %f FX : %f  FY : %f   FZ : %f \n",i,	targetsPotentials[i],forceX[i] ,forceY[i] ,forceZ[i] );
      }

    }

  /**
   * P2P
   * Particles to particles
   * @param inLeafPosition tree coordinate of the leaf
   * @param targets current boxe targets particles
   * @param sources current boxe sources particles (can be == to targets)
   * @param directNeighborsParticles the particles from direct neighbors (this is an array of list)
   * @param size the number of direct neighbors
   */
  void P2P(const FTreeCoordinate& /*inLeafPosition*/,
	   ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict /*sources*/,
	   ContainerClass* const directNeighborsParticles[27], const int /*size*/)
  {
    FP2P::FullMutual(targets,directNeighborsParticles,14);
  }

};

#endif
