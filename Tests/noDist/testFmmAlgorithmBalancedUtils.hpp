// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _LOADFMAANDRUNFMMUTILS_HPP_
#define _LOADFMAANDRUNFMMUTILS_HPP_

#include <fstream>
#include <memory>
#include <assert.h>

#include <array>
#include <vector>

/**
 * \brief Saves the far field costzones to files.
 *
 * One file is created per level, one particle is stored per line in the form :
 * x,y,z,zone.
 *
 * \param args The args that were given to the program
 * \param costzones The CostZones object that was used get the tree balance.
 */
template<class OctreeClass, class CellClass>
void writeZones(const loadFMAAndRunFMMArgs& args, const FCostZones <OctreeClass,CellClass>& costzones)
{
    const std::string outFileBaseName = args.outFileName();
    const std::string outFileExt = args.outFileExt();
    const int verboseLevel = args.verboseLevel();
    const int treeHeight = args.treeHeight();

    auto zones = costzones.getZones();
    long unsigned int zoneCount = zones.size();//args.zoneCount();

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
template <typename FReal, class OctreeClass>
void loadTree(OctreeClass& tree, FFmaGenericLoader<FReal>& loader)
{
    FReal  physicalValue;
    FPoint<FReal> particlePosition;
    // insertion
    for ( int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart ) {
        loader.fillParticle(&particlePosition, &physicalValue);
        tree.insert(particlePosition);
    }
}

/**
 * \brief Loads a tree from a loader.
 * \param tree The the to load into.
 * \param loader The loader to load from.
 */
template <typename FReal, class OctreeClass>
void loadTree(OctreeClass& tree, FRandomLoader<FReal>& loader)
{
    FPoint<FReal> particlePosition;
    // insertion
    for ( int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart ) {
        loader.fillParticle(&particlePosition);
        tree.insert(particlePosition);
    }
}


/**
 * \brief Prints the cost of each zone.
 *
 * This function prints the far-field and near-field costzones costs to stdout.
 * The costs are printed using the format :
 *     zone_number far_field_cost near_field_cost
 *
 * \param tree The tree.
 * \param costzones The FCostZones object that was used to compute the zones. 
 */
template <typename OctreeClass, typename CellClass>
void printZonesCosts(OctreeClass& tree, FCostZones<OctreeClass, CellClass>& costzones)
{
    using CostType = typename CellClass::costtype;
    using TreeIterator = typename OctreeClass::Iterator;
    //using BoundClass = typename FCostZones<OctreeClass, CellClass>::BoundClass;

    auto zones = costzones.getZones();


    // Zones costs (internal, leaf)
    std::vector<std::tuple<CostType,CostType>> zonecosts(zones.size());
    // Leaf colour zone costs
    std::vector<std::array<CostType, FCoordColour::range>> zonecolourcosts(zones.size());


    int zoneIdx = 0;
    for(auto z : zones) {
        for(auto cell : z) {
            std::get<0>(zonecosts.at(zoneIdx)) += cell.second->getCost();
        }
        zoneIdx++;
    }

    auto nearZones = costzones.getLeafZoneBounds();
    zoneIdx = 0;
    int colourIdx = 0;

    for(auto z : nearZones) {
        std::cerr << "zone : " << zoneIdx;

        colourIdx = 0;
        for(auto c : z) {
            std::cerr << "  " << colourIdx;
            
            TreeIterator it = c.first;
            int count = c.second;

            while(count > 0) {
                if( FCoordColour::coord2colour(
                        it.getCurrentCell()->getCoordinate())
                    == colourIdx) {
                    
                    std::get<1>(zonecosts.at(zoneIdx)) +=
                        it.getCurrentCell()->getNearCost();

                    zonecolourcosts.at(zoneIdx).at(colourIdx) +=
                        it.getCurrentCell()->getNearCost();

                    count--;
                }
                if (! it.moveRight() && count > 0) {
                    std::cerr << "Reached the end...";
                    break;
                }
            }

            colourIdx++;
        }
        std::cerr << std::endl;
        zoneIdx++;
    }

    zoneIdx = 0;
    for(auto z : zonecosts) {
        std::cout << "@@"
                  << zoneIdx        << " "
                  << std::get<0>(z) << " "
                  << std::get<1>(z) << " ";

        for(auto c : zonecolourcosts.at(zoneIdx))
            std::cout << c << " ";

        std::cout << std::endl;
        zoneIdx++;
    }


}


#endif
