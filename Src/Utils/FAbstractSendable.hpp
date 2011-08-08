#ifndef FABSTRACTSENDABLE_HPP
#define FABSTRACTSENDABLE_HPP
// /!\ Please, you must read the license at the bottom of this page


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractSendable
* Please read the license
*/
class FAbstractSendable {
protected:
    /** Empty Destructor */
    virtual ~FAbstractSendable(){}

    //static const int SerializedSizeUp = sizeof(?);
    virtual void serializeUp(void* const buffer) const  = 0;
    virtual void deserializeUp(const void* const buffer) = 0;

    //static const int SerializedSizeDown = sizeof(?);
    virtual void serializeDown(void* const buffer) const = 0;
    virtual void deserializeDown(const void* const buffer) = 0;
};

#endif //FABSTRACTSENDABLE_HPP

// [--LICENSE--]
