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

    /** To know the number of bytes needed */
    virtual int bytesToSend() const = 0;

    /** To put the object into the buffer */
    virtual int write(void* const buffer, const int limit) const = 0;

    /** To retreive data from a buffer */
    virtual int read(void* const buffer, const int limit) = 0;

};

#endif //FABSTRACTSENDABLE_HPP

// [--LICENSE--]
