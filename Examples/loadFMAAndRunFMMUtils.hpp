#ifndef _LOADFMAANDRUNFMMUTILS_HPP_
#define _LOADFMAANDRUNFMMUTILS_HPP_

#include "CostZones.hpp"

template<class OctreeClass, class CellClass>
void writeZones(const loadFMAAndRunFMMArgs& args, const CostZones <OctreeClass,CellClass>& costzones)
{
    // GCC versions before 5.0 have not implemented move constructors to streams
    std::vector<std::unique_ptr<std::ofstream>> outfiles;
    for ( int zoneIdx = 0; zoneIdx < args.treeHeight(); zoneIdx++ ) {
        std::unique_ptr<std::ofstream> out(
            new std::ofstream( args.outFileName()
                               + "_" + std::to_string(args.zoneCount()) + "z" 
                               + "." + std::to_string(zoneIdx)
                               + args.outFileExt()));
        *out << "x,y,z,zone" << std::endl;
        outfiles.push_back(std::move(out));
    }
    
    auto zones = costzones.getZones();
    int zoneIdx = 0;
    for ( auto zone : zones) {
        for ( auto cell : zone) {
            *(outfiles[cell.first]) << cell.second->getCoordinate().getX() << ",";
            *(outfiles[cell.first]) << cell.second->getCoordinate().getY() << ",";
            *(outfiles[cell.first]) << cell.second->getCoordinate().getZ() << ",";
            *(outfiles[cell.first]) << zoneIdx << "," << cell.first << std::endl;
        }
        zoneIdx++;
    }

    if ( args.verboseLevel() > 0) {        
        auto& zonebounds = costzones.getZoneBounds();
        zoneIdx = 0;
        for ( auto zone : zonebounds ) {
            std::cout << std::endl << "Zone " << zoneIdx << std::endl;
            int level = 0;
            for ( auto levelbounds : zone ) {
                std::cout << "Level" << level << " : [" << levelbounds.first << ":" << levelbounds.second << "]\n";
                level++;
            }
            zoneIdx++;
        }
    }
}

template <class OctreeClass>
void loadTree(OctreeClass& tree, FFmaGenericLoader& loader)
{
    FReal  physicalValue;
    FPoint particlePosition;
    // insertion
    for ( int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart ) {
        loader.fillParticle(&particlePosition, &physicalValue);
        tree.insert(particlePosition);
    }
}

#endif
