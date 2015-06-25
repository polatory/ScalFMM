#include <iostream>

#include "Files/FMpiFmaDivider.hpp"

using FReal = double;

#include "testFmmAlgorithmBalancedArgs.hpp"


int main(int argc, char** argv)
{

    loadFMAAndRunFMMArgs args(argc, argv);

    FMpiFmaDivider<FReal>
        divider(args.inFileName(),
                args.zoneCount(),
                args.treeHeight());

}
