// @SCALFMM_PRIVATE
#ifndef FSTARPUKERNELCAPACITIES_HPP
#define FSTARPUKERNELCAPACITIES_HPP

#include "FStarPUUtils.hpp"

/** A class used with the starpu system should
  * implement this interface in order to inform the algorithm about what the kernel
  * is doing.
  */
class FStarPUKernelCapacities {
public:
    virtual bool supportP2M(const FStarPUTypes inPu) const = 0;
    virtual bool supportM2M(const FStarPUTypes inPu) const = 0;
    virtual bool supportM2L(const FStarPUTypes inPu) const = 0;
    virtual bool supportL2L(const FStarPUTypes inPu) const = 0;
    virtual bool supportL2P(const FStarPUTypes inPu) const = 0;
    virtual bool supportP2P(const FStarPUTypes inPu) const = 0;
};

/**
 * This is for the kernels that implement all the methods.
 */
template <class BaseClass>
class FStarPUAllYesCapacities : public BaseClass, public FStarPUKernelCapacities {
public:
    using BaseClass::BaseClass;

    bool supportP2M(const FStarPUTypes /*inPu*/) const override {
        return true;
    }
    bool supportM2M(const FStarPUTypes /*inPu*/) const override {
        return true;
    }
    bool supportM2L(const FStarPUTypes /*inPu*/) const override {
        return true;
    }
    bool supportL2L(const FStarPUTypes /*inPu*/) const override {
        return true;
    }
    bool supportL2P(const FStarPUTypes /*inPu*/) const override {
        return true;
    }
    bool supportP2P(const FStarPUTypes /*inPu*/) const override {
        return true;
    }
};

template <class BaseClass>
class FStarPUAllCpuCapacities : public BaseClass, public FStarPUKernelCapacities {
public:
    using BaseClass::BaseClass;

    bool supportP2M(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CPU_IDX;
    }
    bool supportM2M(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CPU_IDX;
    }
    bool supportM2L(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CPU_IDX;
    }
    bool supportL2L(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CPU_IDX;
    }
    bool supportL2P(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CPU_IDX;
    }
    bool supportP2P(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CPU_IDX;
    }
};

#endif // FSTARPUKERNELCAPACITIES_HPP

