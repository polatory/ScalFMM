// @FUSE_MPI

#include <iostream>
#include <string>

#include "Files/FMpiFmaDivider.hpp"

#include "testFMpiFmaDividerArgs.hpp"

using FReal = double;

int main(int argc, char** argv)
{

    testFMpiFmaDividerArgs args(argc, argv);


    FMpiFmaDivider<FReal>
        divider(args.inFileName(),
                args.outFileName() + "_" + args.dispatchPolicyString()
                                   + "_" + std::to_string(args.zoneCount())
                                   + "z_h" + std::to_string(args.treeHeight()),
                args.outFileExt(),
                args.zoneCount(),
                args.treeHeight(),
                FMpiFmaDivider<FReal>::DispatchPolicy(args.dispatchPolicy()));

}
