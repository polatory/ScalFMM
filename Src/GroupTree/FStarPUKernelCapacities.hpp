#ifndef FSTARPUKERNELCAPACITIES_HPP
#define FSTARPUKERNELCAPACITIES_HPP

/** A class used with the starpu system should
  * implement this interface in order to inform the algorithm about what the kernel
  * is doing.
  */
class FStarPUKernelCapacities {
public:
    virtual bool supportP2M() const = 0;
    virtual bool supportM2M() const = 0;
    virtual bool supportM2L() const = 0;
    virtual bool supportL2L() const = 0;
    virtual bool supportL2P() const = 0;
    virtual bool supportP2P() const = 0;
};

/**
 * This is for the kernels that implement all the methods.
 */
template <class BaseClass>
class FStarPUAllYesCapacities : public BaseClass, public FStarPUKernelCapacities {
public:
    using BaseClass::BaseClass;

    bool supportP2M() const override {
        return true;
    }
    bool supportM2M() const override {
        return true;
    }
    bool supportM2L() const override {
        return true;
    }
    bool supportL2L() const override {
        return true;
    }
    bool supportL2P() const override {
        return true;
    }
    bool supportP2P() const override {
        return true;
    }
};

#endif // FSTARPUKERNELCAPACITIES_HPP

