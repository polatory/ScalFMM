#include <iostream>

#include "Files/FMpiFmaDivider.hpp"

using FReal = double;

#include "testFMpiFmaDividerArgs.hpp"


int main(int argc, char** argv)
{

    testFMpiFmaDividerArgs args(argc, argv);


    FMpiFmaDivider<FReal>
        divider(args.inFileName(),
                args.outFileName(),
                args.outFileExt(),
                args.zoneCount(),
                args.treeHeight(),
                FMpiFmaDivider<FReal>::DispatchPolicy(args.dispatchPolicy()));

}
