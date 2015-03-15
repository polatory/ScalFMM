#ifndef FTESTOPENCLCODE_HPP
#define FTESTOPENCLCODE_HPP

#include "../../Utils/FGlobal.hpp"
#include "../FStarPUDefaultAlign.hpp"
#include "FTextReplacer.hpp"

// Initialize the types
class FTestOpenCLCode{
    FTextReplacer kernelfile;
    size_t dim;

public:
    //FTestOpenCLCode() : kernelfile("/home/berenger/Projets/ScalfmmGit/scalfmm/Src/GroupTree/OpenCl/FEmptyKernel.cl"){
    FTestOpenCLCode() : kernelfile("/home/berenger/Projets/ScalfmmGit/scalfmm/Src/GroupTree/OpenCl/FTestKernel.cl"){
        kernelfile.replaceAll("___FReal___", "double");
        kernelfile.replaceAll("___FParticleValueClass___", "long long");
        kernelfile.replaceAll("___NbAttributesPerParticle___", 2);
        const size_t structAlign = FStarPUDefaultAlign::StructAlign;
        kernelfile.replaceAll("___DefaultStructAlign___", structAlign);

        dim = 1;
    }

    const char* getKernelCode(const int /*inDevId*/){
        return kernelfile.getContent();
    }

    void releaseKernelCode(){
        kernelfile.clear();
    }

    size_t getNbDims() const {
        return 1;
    }

    const size_t* getNbGroups() const {
        // We return 1
        return &dim;
    }

    const size_t* getGroupSize() const {
        // We return 1
        return &dim;
    }
};


#endif // FTESTOPENCLCODE_HPP

