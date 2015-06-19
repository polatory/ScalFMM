#ifndef _FMPIPARTICLEDIVIDER_HPP_
#define _FMPIPARTICLEDIVIDER_HPP_

#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "../Files/FFmaGenericLoader.hpp"

#include "../Containers/FCoordinateComputer.hpp"
#include "../Utils/FPoint.hpp"



/** \brief Loads a particle file and splits it for MPI processing.
 *
 * This is meant to be used to split a particle file and save each part in a
 * different file for MPI processes to pick up later. The aim is to distribute
 * the particles in a way that will allow to balance the computing load of all
 * processes during the FMM.
 *
 * The division is not done at the leaf level of the tree. The tree is cut at a
 * set height and the sub trees that result are distributed among the processes.
 *
 *                  *
 *               /     \
 *            *           *
 *      ____/___\_______/___\____ <-- division level
 *      |  *  |  *  |  *  |  *  |
 *      | / \ | / \ | / \ | / \ | <-- sub-trees
 *      |*   *|*   *|*   *|*   *|
 *
 * 
 */
template <typename FReal>
class FMPIParticleDivider {
private:

    /** \brief Particle structure.
     *
     * This structure is used to store particles read from an input file. It is
     * designed to ease sorting according to the Morton index of the tree cells
     * without building it.
     *
     * \note All particles are stored in the same box space which is divided in
     * 8^_divisionLevel cells.
     */
    class Particle {
    public:
        /// Center of the particles' box
        static FPoint<FReal> _boxCenter;
        /// Width of the particles' box
        static FReal _boxWidth;
        /// Division level of the particles' box
        static int   _divisionLevel;
        
        /// Particle position
        FPoint<FReal> position;
        /// Particle physical value
        FReal value;
        /// Contructor
        Particle(FPoint<FReal> pos, FReal val) :
            position(pos),
            value(val) {
        }
        
        /// Default constructor
        Particle() = default;
        /// Default copy constructor
        Particle(const Particle& other) = default;
        /// Default move constructor
        Particle(Particle&& other) = default;

        /** \brief Get containing cell Morton index.
         *
         * The index is computed from #position relative to #_divisionLevel,
         * #_boxWidth and #_boxCenter.
         *
         * \return The containing cell Morton index
         */
        MortonIndex getMortonIndex() const {
            return FCoordinateComputer::
                GetCoordinateFromPosition(_boxCenter, _boxWidth,
                                          _divisionLevel, position).
                getMortonIndex(_divisionLevel);
            
        }

        /// Equality operator
        bool operator==(const Particle& other) const {
            return position == other.position &&
                value == other.value;
        }

        /// Inequality operator
        bool operator!=(const Particle& other) const {
            return ! operator==(other);
        }
        
        /// Strictly less than operator
        bool operator<(const Particle& other) const {
            return getMortonIndex() < other.getMortonIndex();
        }

        /// Less than operator
        bool operator<=(const Particle& other) const {
            return ! operator>(other);
        }

        /// Strictly greater than operator
        bool operator>(const Particle& other) const {
            return getMortonIndex() > other.getMortonIndex();
        }

        /// Greater than operator
        bool operator>=(const Particle& other) const {
            return ! operator<(other);
        }

        /// Output stream operator
        friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
            return os << p.position << " " << p.value;
        }
    };

    /// The input file name
    std::string _filename;

    /// The output's basename
    std::string _outputBasename;

    /// The number of parts to split the file into
    int _splitCount; 

    /// The level of division of space
    int _divisionLevel;

    /// A set to sort the particles using their Morton index
    std::multiset<Particle> _particles;

    /// The particle loader
    FFmaGenericLoader<FReal> _loader;

    /// The loader box center
    FPoint<FReal> _boxCenter{0,0,0};
    /// The loader box width
    FReal _boxWidth = 1;

    /// Compares two particle based on their Morton index
    bool compareParticles(const Particle& lhs, const Particle& rhs) {
        return lhs.getMortonIndex() < rhs.getMortonIndex();
    }




public:
    /** \brief Constructor
     *
     * Build the particle divider and opens the loader.
     * \param filename[in] The particle input file.
     * \param outputBasename[in] The basename for the output files.
     * \param splitcount[in] The number of parts to split the input file into.
     * \param divisionLevel[in] The level in the to-be tree at which to do the division.
     */
    FMPIParticleDivider(const std::string filename,
                        const std::string outputBasename,
                        int splitcount,
                        int divisionLevel) :
        _filename(filename),
        _outputBasename(outputBasename),
        _splitCount(splitcount),
        _divisionLevel(divisionLevel),
        _loader(filename),
        _boxCenter(_loader.getCenterOfBox()),
        _boxWidth(_loader.getBoxWidth()) {
        // Set particles' space to enable comparisons.
        Particle::_divisionLevel = _divisionLevel;
        Particle::_boxCenter = _boxCenter;
        Particle::_boxWidth = _boxWidth;

        readFile();
        writeFiles();

    }

    /// Deleted copy constructor because of loader
    FMPIParticleDivider(const FMPIParticleDivider&) = delete;
    /// Deleted move constructor because of loader
    FMPIParticleDivider(const FMPIParticleDivider&&) = delete;


protected:
    /// Reads the particles from the loader and stores them sorted by Morton
    /// index.
    void readFile() {
        FReal physicalValue;
        FPoint<FReal> particlePosition;
        for ( int idxPart = 0 ; idxPart < _loader.getNumberOfParticles() ; ++idxPart ) {
            // Read one particle from file
            _loader.fillParticle(&particlePosition, &physicalValue);
            // And insert in the sorted set
            _particles.emplace(particlePosition, physicalValue);
        }
    }

    /// Writes the output file based on #_outputBasename.
    void writeFiles() {
        
        /*
          Creation of _splitCount files to write to.

          GCC versions before 5.0 have not implemented stream move
          constructors; we use unique pointers to get around this problem.
        */
        std::vector<std::unique_ptr<std::ofstream>> outfiles(_splitCount);

        for ( int fileIdx = 0; fileIdx < _splitCount; fileIdx++ ) {
            std::unique_ptr<std::ofstream> out(
                new std::ofstream( _outputBasename
                                   + "_" + std::to_string(_splitCount) + "z" 
                                   + "." + std::to_string(fileIdx),
                                   std::ofstream::trunc
                                   ));
            outfiles[fileIdx] = std::move(out);
        }


        // Temporary vector to store all particles with same Morton index
        std::vector<Particle> container;

        using size_type = typename std::vector<Particle>::size_type;
        size_type totalParticleCount = _particles.size();
        size_type currentGlobalParticleCount = 0;
        size_type currentFileParticleCount = 0;
        MortonIndex currentMortonIndex = 0;
        long long int currentFileIndex = 0;
        while(_particles.begin() != _particles.end()) {
            /* Particles are sorted, we pull out all those with the same Morton
             * index and put them in the temporary container. */
            while((*_particles.begin()).getMortonIndex() == currentMortonIndex) {
                container.push_back(*(_particles.begin()));
                _particles.erase(_particles.begin());
            }          
            // File number
            currentFileIndex = currentGlobalParticleCount * _splitCount
                / totalParticleCount;
            // Write particles to file.
            for(Particle p : container) {
                *outfiles[currentFileIndex] << p << std::endl;
            }

            currentGlobalParticleCount += container.size();
            currentMortonIndex += container.size();
            container.clear();
            currentMortonIndex = (*_particles.begin()).getMortonIndex();
        }

    }


};


template<typename FReal>
FPoint<FReal> FMPIParticleDivider<FReal>::Particle::_boxCenter = FPoint<FReal>{0,0,0};
template<typename FReal>
FReal FMPIParticleDivider<FReal>::Particle::_boxWidth = 1;
template<typename FReal>
int FMPIParticleDivider<FReal>::Particle::_divisionLevel = 4;



#endif
