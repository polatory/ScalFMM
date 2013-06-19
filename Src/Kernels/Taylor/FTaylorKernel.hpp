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
  
  ////////////////////////////////////////////////////
  // Private method
  ////////////////////////////////////////////////////

  /** Return the position of a leaf from its tree coordinate */
  FPoint getLeafCenter(const FTreeCoordinate coordinate) const {
    return FPoint(
		  FReal(coordinate.getX()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX(),
		  FReal(coordinate.getY()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getY(),
		  FReal(coordinate.getZ()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getZ());
  }

  /** 
   * @brief Incrementation of powers in Taylor expansion 
   * Result : ...,[2,0,0],[1,1,0],[1,0,1]...  3-tuple are sorted 
   * by size and alphabetical order.
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
      printf("Error factorial negative!! a=%d\n",a);
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
    return (fact(a)*fact(b)*fact(c));
  }
  /* Return the combine of a paire of number
   */
  FReal combin(const int& a, const int& b){
    if(a-b<0) printf("Error combin negative!! a=%d b=%d\n",a,b);
    return fact(a) / (fact(b)*fact(a-b));
  }

  
  /** @brief Init the derivative array for using of following formula
   * from "A cartesian tree-code for screened coulomb interactions"
   *
   *\f[
   * \Psi_{\mathbf{k}}^{c}\times \left [\left |\mathbf{k}\right |\times \left |
   * \mathbf{x}_i-\mathbf{x}_c\right |^2 \times 
   * \frac{1}{\mathbf{k}!} \right ]\				
   * = (2\times \left |{\mathbf{k}}\right |-1)\times
   * \sum_{j=0}^{3}\left [(x_{i_j}-x_{c_j})\times
   * \Psi_{\mathbf{k}-e_j,i}^{c}\times \frac{1}{({\mathbf{k}}-e_j)!}\right ]\ 
   * -(\left |\mathbf{k}\right |-1)\times  \sum_{j=0}^{3}\left
   * [\Psi_{\mathbf{k}-2\times e_j,i}^{c}\times \frac{1}{(\mathbf{k}-2\times e_j)!}\right]
   * \f]
   */
  void initDerivative(const FPoint target, const FPoint src, FReal * tab)
  {
    FReal dx = target.getX()-src.getX();
    FReal dy = target.getY()-src.getY();
    FReal dz = target.getZ()-src.getZ();
    
    tab[0]=1/(FMath::Sqrt(dx*dx+dy*dy+dz*dz));
    FReal R3 = tab[0]/(dx*dx+dy*dy+dz*dz);
    tab[1]=dx*R3;
    tab[2]=dy*R3;
    tab[3]=dz*R3;
    FReal R5 = R3/(dx*dx+dy*dy+dz*dz);
    tab[4] = R3-dx*dx*R5;
    tab[5] = (-dx)*dy*R5;
    tab[6] = (-dx)*dz*R5;
    tab[7] = R3-dy*dy*R5;
    tab[8] = (-dy)*dz*R5;
    tab[9] = R3-dz*dz*R5;
  }
  
  /** @brief Compute and store the derivative for a given tuple.
   *
   */
  FReal computeDerivative(const int a, const int b, const int c, 
			  const FPoint target, 
			  const FPoint src, 
			  FReal * yetComputed)
  {
    int idx = powerToIdx(a,b,c);
    if((yetComputed[idx] =! 0))
      {return yetComputed[idx];}
    else
      {
	FReal dx = target.getX()-src.getX();
	FReal dy = target.getY()-src.getY();
	FReal dz = target.getZ()-src.getZ();
	FReal fct = fact3int(a,b,c);
	FReal dist = FMath::Sqrt(dx*dx+dy*dy+dz*dz);
	FReal temp_value = FReal(0);
	int idxt;
	if(a > 0){
	    idxt = powerToIdx(a-1,b,c);
	    temp_value += ((FReal)(2*(a+b+c)-1))*dx*yetComputed[idxt]*FReal(a)/fct;
	    if(a > 1){
	      idxt = powerToIdx(a-2,b,c);
	      temp_value -= FReal(a+b+c-1)*yetComputed[idxt]*FReal(a*(a-1))/fct;
	    }
	}
	if(b > 0){
	    idxt = powerToIdx(a,b-1,c);
	    temp_value += FReal(2*(a+b+c)-1)*dy*yetComputed[idxt]*FReal(b)/fct;
	    if(b > 1){
	      idxt = powerToIdx(a,b-2,c);
	      temp_value -= FReal(a+b+c-1)*yetComputed[idxt]*FReal(b*(b-1))/fct;
	    }
	}
	if(c > 0){
	    idxt = powerToIdx(a,b,c-1);
	    temp_value += FReal(2*(a+b+c)-1)*dz*yetComputed[idxt]*FReal(c)/fct;
	    if(c > 1){
	      idxt = powerToIdx(a,b,c-2);
	      temp_value -= FReal(a+b+c-1)*yetComputed[idxt]*FReal(c*(c-1))/fct;
	    }
	}
	FReal coeff = FReal(a+b+c)*dist/fct;
	temp_value = temp_value/coeff;
	yetComputed[idx] = temp_value;
	return temp_value;
      }
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
  }
  
  /* Default destructor
   */
  virtual ~FTaylorKernel(){
  }

  /**P2M 
   * @brief Fill the Multipole with the field created by the cell
   * particles.
   * 
   * Formula :
   * \f[
   *   M_{k} = \sum_{j=0}^{nb_particule} q_j * (x_c-x_j)^{k}*|k|!/k!
   * \f]
   */
  void P2M(CellClass* const pole, 
	   const ContainerClass* const particles)
  { 
    //Variables computed for each power of Multipole
    int a=0,b=0,c=0;
    FReal facto;
    FReal coeff; 
    //Copying cell center position once and for all
    const FPoint& cellCenter = getLeafCenter(pole->getCoordinate());
    
    FReal * multipole = pole->getMultipole();
    
    //Iterating over MutlipoleVector
    for(int i=0 ; i<SizeVector ; i++)
      {
	// Iterator over Particles
	int nbPart = particles->getNbParticles();
	const FReal* const * positions = particles->getPositions();
	const FReal* posX = positions[0];
	const FReal* posY = positions[1];
	const FReal* posZ = positions[2];

	const FReal* phyValue = particles->getPhysicalValues();

	//update needed values
	facto=fact3int(a,b,c);
	coeff=fact(a+b+c)/(facto*facto);
	
	//Iterating over Particles
	for(int idPart=0 ; idPart<nbPart ; idPart++){
	  
	  FReal dx = cellCenter.getX()-posX[idPart];
	  FReal dy = cellCenter.getY()-posY[idPart];
	  FReal dz = cellCenter.getZ()-posZ[idPart];
	  
	  const FReal potential = phyValue[idPart];
	  
	  //Computation
	  FReal value = FReal(0);
	  value+=potential*FMath::pow(dx,a)*FMath::pow(dy,b)*FMath::pow(dz,c)*coeff;
	  multipole[i] = value;
	  
	 
	}
	//inc powers of expansion
	incPowers(&a,&b,&c);
      }
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
    //Powers of expansions
    int a=0,b=0,c=0;
    int sum = 0; //a+b+c
    //Indexes of powers
    int idx_a,idx_b,idx_c;

    //Distance from current child to parent
    FReal boxSize = (boxWidth/FReal(1 << inLevel));
    FReal dx = 0;
    FReal dy = 0;
    FReal dz = 0;
    //Center point of parent cell
    //const FPoint cellCenter = getLeafCenter(pole->getCoordinate());
    
    
    //Iteration over the eight children
    int idxChild;
    for(idxChild=0 ; idxChild<8 ; ++idxChild)
      {
	if(child[idxChild]){
	  //Set the distance between centers of cells
	  dz = ((FReal )(2*(1 & idxChild)-1))*boxSize;
	  dy = ((FReal )(2*((1 << 1) & idxChild)-1))*boxSize;
	  dx = ((FReal )(2*((1 << 2) & idxChild)-1))*boxSize;
	  
	  //Iteration over parent multipole array
	  FReal * mult = pole->getMultipole();
	  a=0;
	  b=0;
	  c=0;
	  sum=0;
	  for(int idxMult = 0 ; idxMult<SizeVector ; idxMult++)
	    {
	      FReal value = mult[idxMult];
	      int idMultiChild;
	      //Iteration over the powers to find the cell multipole
	      //involved in the computation of the parent multipole
	      for(idx_a=0 ; idx_a<=sum ; ++idx_a){
		for(idx_b=0 ; idx_b<=sum-a ; ++idx_b){
		  for(idx_c=0 ; idx_c<sum-(a+b) ; ++idx_c){
		    //Computation
		    //Child multipole involved
		    idMultiChild=powerToIdx(idx_a,idx_b,idx_c);
		    value+=child[idxChild]->getMultipole()[idMultiChild]
		      *FMath::pow(dx,idx_a)*FMath::pow(dy,idx_b)*FMath::pow(dz,idx_c)/fact3int(idx_a,idx_b,idx_c);
		  }
		}
	      }
	      mult[idxMult] = value;
	      incPowers(&a,&b,&c);
	      sum=a+b+c;
	    }
	}
      }
  }
  
  /**
   *@brief Convert the multipole expansion into local expansion The
   * operator do not use symmetries.  
   * \f[ 
   * L_{\mathbf{n}}^{c} =
   * \frac{1}{\mathbf{n}!} \times \sum_{\mathbf{k}=0}^{p} \left [
   * \Psi_{\mathbf{n+k}}^{c}(\mathbf{x}\times M_{\mathbf{k}^c})\right
   * ] \f]
   */
  void M2L(CellClass* const FRestrict local, 
	   const CellClass* distantNeighbors[343],
	   const int /*size*/, const int /*inLevel*/)
  {
    //Iteration over distantNeighbors
    int idxNeigh;
    FPoint locCenter = getLeafCenter(local->getCoordinate());
    for(idxNeigh=0 ; idxNeigh<343 ; idxNeigh++){

      //Need to test if current neighbor is one of the interaction list
      if(distantNeighbors[idxNeigh]){
	//Derivatives are computed iteratively
	FReal yetComputed[(2*P+1)*(2*P+2)*(2*P+3)/6];
	FMemUtils::memset(yetComputed,0,SizeVector*sizeof(FReal(0)));
	FPoint curDistCenter = getLeafCenter(distantNeighbors[idxNeigh]->getCoordinate());
	initDerivative(locCenter, curDistCenter, yetComputed);

	//Iteration over Multipole / Local
	int al=0,bl=0,cl=0;
	int am=0,bm=0,cm=0;
	FReal * iterLocal = local->getLocal();
	
	
	//Iterating over local array : n
	for(int i=0 ; i< SizeVector ; i++){
	  FReal fctl = fact3int(al,bl,cl);
	  
	  //Iterator over multipole array 
	  const FReal * multipole = distantNeighbors[idxNeigh]->getMultipole(); 
	  
	  //Iterating over multipole array : k
	  for(int j = 0 ; j < SizeVector ; j++){
	    FReal data = computeDerivative(al+am,bl+bm,cl+cm,locCenter,getLeafCenter(distantNeighbors[idxNeigh]->getCoordinate()),yetComputed);
	    data *= multipole[j]/fctl;
	    iterLocal[i]+=data;
	    
	    //updating a,b,c and iterator
	    incPowers(&am,&bm,&cm);
	  }
	  am=0;
	  bm=0;
	  cm=0;
	  incPowers(&al,&bl,&cl);
	}
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
 
  
  /**
   *@brief Apply on the particles the force computed from the local expansion
   *
   */
    void L2P(const CellClass* const local, 
	   ContainerClass* const particles)
    {
      FPoint locCenter = getLeafCenter(local->getCoordinate());
      //Iterator over particles
      int nbPart = particles->getNbParticles();
      
      const FReal * const * positions = particles->getPositions();
      const FReal * posX = positions[0];
      const FReal * posY = positions[1];
      const FReal * posZ = positions[2];
      
      FReal * const forceX = particles->getForcesX();
      FReal * const forceY = particles->getForcesY();
      FReal * const forceZ = particles->getForcesZ();
      
      //Iteration over particles
      for(int i=0 ; i<nbPart ; i++){
	
	FReal dx = locCenter.getX() - posX[0];
	FReal dy = locCenter.getY() - posY[1];
	FReal dz = locCenter.getZ() - posZ[2];
	
	int a=0,b=0,c=0;
	
	//Iteration over Local array
	const FReal * iterLocal = local->getLocal();
	for(int j=0 ; j<SizeVector ; j++){
	  
	  	  
	  FReal locForce = iterLocal[j];
	  //Application of forces
	  forceX[i] = FMath::pow(dx,a)*locForce;
	  forceY[i] = FMath::pow(dy,b)*locForce;
	  forceZ[i] = FMath::pow(dz,c)*locForce;
	  
	  
	  incPowers(&a,&b,&c);
	}
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
