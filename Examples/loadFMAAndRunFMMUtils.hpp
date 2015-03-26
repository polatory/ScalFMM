#ifndef _LOADFMAANDRUNFMMUTILS_HPP_
#define _LOADFMAANDRUNFMMUTILS_HPP_

#include "CostZones.hpp"
#include <assert.h>

/**
 * \brief Saves the costzones to files.
 *
 * One file is created per level, one particle is stored per line in the form :
 * x,y,z,zone.
 *
 * \param args The args that were given to the program
 * \param costzones The CostZones object that was used get the tree balance.
 */
template<class OctreeClass, class CellClass>
void writeZones(const loadFMAAndRunFMMArgs& args, const CostZones <OctreeClass,CellClass>& costzones)
{
    std::string outFileBaseName = args.outFileName();
    std::string outFileExt = args.outFileExt();
    int verboseLevel = args.verboseLevel();
    int treeHeight = args.treeHeight();

    auto zones = costzones.getZones();
    int zoneCount = zones.size();//args.zoneCount();

    std::cout << "Writing " << zoneCount << " zones." << std::endl;

    // GCC versions before 5.0 have not implemented move constructors to streams
    // we use unique pointers to get around this problem.
    std::vector<std::unique_ptr<std::ofstream>> outfiles;
    for ( int levelIdx = 0; levelIdx < treeHeight; levelIdx++ ) {
        std::unique_ptr<std::ofstream> out(
            new std::ofstream( outFileBaseName
                               + "_" + std::to_string(zoneCount) + "z" 
                               + "." + std::to_string(levelIdx)
                               + outFileExt));
        *out << "x,y,z,zone" << std::endl;
        outfiles.push_back(std::move(out));
    }
    
    int zoneIdx = 0;
    for ( auto zone : zones) {
        for ( auto cell : zone) {
            *(outfiles[cell.first]) << cell.second->getCoordinate().getX() << ",";
            *(outfiles[cell.first]) << cell.second->getCoordinate().getY() << ",";
            *(outfiles[cell.first]) << cell.second->getCoordinate().getZ() << ",";
            *(outfiles[cell.first]) << zoneIdx << std::endl;
        }
        zoneIdx++;
    }

    if ( verboseLevel > 0) {        
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

/**
 * \brief Loads a tree from a loader.
 * \param tree The the to load into.
 * \param loader The loader to load from.
 */
template <class OctreeClass>
void loadTree(OctreeClass& tree, FFmaGenericLoader& loader)
{
    FReal  physicalValue;
    FPoint particlePosition;
    // insertion
    for ( int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart ) {
        loader.fillParticle(&particlePosition, &physicalValue);
        tree.insert(particlePosition, idxPart);
    }
}

#endif
