#ifndef FABSTRACTCELL_HPP
#define FABSTRACTCELL_HPP
// /!\ Please, you must read the license at the bottom of this page


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractCell
* @brief
* Please read the license
*
* This class define the method that every cell class
* has to implement.
* @warning Inherite from this class when implement a specific cell type
*/
class FAbstractCell{
public:	
	/** Default destructor */
        virtual ~FAbstractCell(){
	}

	/**
        * Must be implemented by each user Cell class
	* @return the position of the current cell
	*/
        virtual MortonIndex getMortonIndex() const = 0;


        /**
        * Must be implemented by each user Cell class
        * @param inIndex the position of the current cell
        */
        virtual void setMortonIndex(const MortonIndex inIndex) = 0;
};


#endif //FABSTRACTCELL_HPP

// [--LICENSE--]
