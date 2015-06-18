#include <iostream>

#include "Files/FMPIParticleDivider.hpp"

using FReal = double;

#include "testFmmAlgorithmBalancedArgs.hpp"


int main(int argc, char** argv)
{

    loadFMAAndRunFMMArgs args(argc, argv);

    FMPIParticleDivider<FReal>
        divider(args.inFileName(),
                args.outFileName(),
                args.zoneCount(),
                args.treeHeight());

}
