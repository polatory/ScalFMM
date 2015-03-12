#ifndef FTESTOPENCLCODE_HPP
#define FTESTOPENCLCODE_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Components/FTestCell.hpp"
#include "FTextReplacer.hpp"

struct FTestCell_Alignement{
    static const int dataUp;
    static const int dataDown;
    static const int mindex;
    static const int coord;
};

const int FTestCell_Alignement::dataUp = reinterpret_cast<std::size_t>(&((reinterpret_cast<FTestCell*>(0xF00))->dataUp)) - std::size_t(0xF00);
const int FTestCell_Alignement::dataDown = reinterpret_cast<std::size_t>(&((reinterpret_cast<FTestCell*>(0xF00))->dataDown)) - std::size_t(0xF00);
const int FTestCell_Alignement::mindex = reinterpret_cast<std::size_t>(&((reinterpret_cast<FTestCell*>(0xF00))->mortonIndex)) - std::size_t(0xF00);
const int FTestCell_Alignement::coord = reinterpret_cast<std::size_t>(&((reinterpret_cast<FTestCell*>(0xF00))->coordinate)) - std::size_t(0xF00);


// Initialize the types
class FTestOpenCLCode{
    FTextReplacer kernelfile;
    size_t dim;

public:
    //FTestOpenCLCode() : kernelfile("/home/berenger/Projets/ScalfmmGit/scalfmm/Src/GroupTree/OpenCl/FEmptyKernel.cl"){
    FTestOpenCLCode() : kernelfile("/home/berenger/Projets/ScalfmmGit/scalfmm/Src/GroupTree/OpenCl/FTestKernel.cl"){
        kernelfile.replaceAll("___FReal___", "double");
        kernelfile.replaceAll("___FParticleValueClass___", "long long");
        kernelfile.replaceAll("___FCellClassSize___", sizeof(FTestCell));
        kernelfile.replaceAll("___NbAttributesPerParticle___", 2);
        kernelfile.replaceAll("___FCellUpOffset___", FTestCell_Alignement::dataUp);
        kernelfile.replaceAll("___FCellDownOffset___", FTestCell_Alignement::dataDown);
        kernelfile.replaceAll("___FCellMortonOffset___", FTestCell_Alignement::mindex);
        kernelfile.replaceAll("___FCellCoordinateOffset___", FTestCell_Alignement::coord);

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

