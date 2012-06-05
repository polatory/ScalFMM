#ifndef FABSTRACTSERIALIZABLE_HPP
#define FABSTRACTSERIALIZABLE_HPP

class FBufferReader;
class FBufferWriter;

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractSerializable
* Please read the license
*
* This propose an interface to save and restore a class
*
* To make your particles are usable in the mpi fmm,
* they must provide this interface
*/
class FAbstractSerializable {
protected:
    virtual void save(FBufferWriter&) const  = 0;
    virtual void restore(FBufferReader&) = 0;
};

#endif // FABSTRACTSERIALIZABLE_HPP
