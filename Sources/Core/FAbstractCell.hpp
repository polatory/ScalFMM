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
*
* In fact FOctree & FFMMAlgorithm need this function to be implemented.
* But you cannot use this interface with the extension (as an example :
* because the compiler will faill to know if getMortonIndex is coming
* from this interface or from the extension)
*
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

        /**
        * Must be implemented by each user Cell class
        * @param inPosition the position of the current cell
        */
        virtual void setPosition(const F3DPosition& inPosition) = 0;

        /** Because the system can run in ToS mode
          * a cell has to express if it has sources
          * @return true if there are sources particules inside
          */
        virtual bool hasSourcesChild() const = 0;

        /** Because the system can run in ToS mode
          * a cell has to express if it has targets
          * @return true if there are targets particules inside
          */
        virtual bool hasTargetsChild() const = 0;

        /**
          * This function make the cell containing sources
          */
        virtual void setSourcesChildTrue() = 0;

        /**
          * This function make the cell containing targets
          */
        virtual void setTargetsChildTrue() = 0;
};


#endif //FABSTRACTCELL_HPP

// [--LICENSE--]
