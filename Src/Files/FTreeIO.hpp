#ifndef FTREEIO_HPP
#define FTREEIO_HPP

#include <iostream>
#include <fstream>

/** This class proposes static methods to save and load
  * a tree.
  * It used binary format (so FReal must be the same!)
  *
  * The format is :
  * [nb leaves]
  * [nb particles in leaf][particles.....] ...
  * [nb cells at level]
  * [morton index][cell] ....
  * ...
  */
class FTreeIO{
public:
    /** The serializer class call a method on the particles or cells
      * So they have to implement the write/read method
      */
    template <class CellClass, class ParticleClass>
    class Serializer {
    public:
        static void PutCell(std::ofstream*const stream, const CellClass*const cell) {
            cell->write(stream);
        }
        static void PutParticles(std::ofstream*const stream, const ParticleClass* const particles, const int nbParticles) {
            for( int idxParticle = 0 ; idxParticle < nbParticles ; ++idxParticle){
                particles[idxParticle]->write(stream);
            }
        }

        static void GetCell(std::ifstream*const stream, CellClass* const cell){
            cell->read(stream);
        }
        static void GetParticles(std::ifstream*const stream, const ParticleClass* const particles, const int nbParticles){
            for( int idxParticle = 0 ; idxParticle < nbParticles ; ++idxParticle){
                particles[idxParticle]->read(stream);
            }
        }
    };

    /** The copier method simple copy in memory using memcpy
      */
    template <class CellClass, class ParticleClass>
    class Copier {
    public:
        static void PutCell(std::ofstream*const stream, const CellClass* const cell) {
            stream->write((const char*)cell, sizeof(CellClass));
        }
        static void PutParticles(std::ofstream*const stream, const ParticleClass* const particles, const int nbParticles) {
            stream->write((const char*)particles, nbParticles * sizeof(ParticleClass));
        }

        static void GetCell(std::ifstream*const stream, CellClass* const cell){
            stream->read((char*)cell, sizeof(CellClass));
        }
        static void GetParticles(std::ifstream*const stream, const ParticleClass* const particles, const int nbParticles){
            stream->read((char*)particles, nbParticles * sizeof(ParticleClass));
        }
    };

    /** To save in memory */
    template <class OctreeClass, class CellClass, class ParticleClass, class ClassProptotype >
    static bool Save(const char filename[], OctreeClass& tree){
        std::ofstream file(filename, std::ofstream::binary | std::ofstream::out );

        if(!file.good()){
            return false;
        }

        const size_t sizeof_freal = sizeof(FReal);
        file.write((const char*)&sizeof_freal, sizeof(size_t));

        const int height = tree.getHeight();
        file.write((const char*)&height, sizeof(int));

        const int subHeight = tree.getSubHeight();
        file.write((const char*)&subHeight, sizeof(int));

        const FReal width = tree.getBoxWidth();
        file.write((const char*)&width, sizeof(FReal));

        file.write((const char*)&tree.getBoxCenter(), sizeof(F3DPosition));

        {
            typename OctreeClass::Iterator octreeIterator(&tree);

            int maxParticlesInLeaf = 0;
            int nbLeaf = 0;
            const std::ofstream::pos_type posNbLeaf = file.tellp();
            file.write((const char*)&nbLeaf,sizeof(int));
            file.write((const char*)&maxParticlesInLeaf,sizeof(int));

            octreeIterator.gotoBottomLeft();
            const bool useTargetSource = (octreeIterator.getCurrentListSrc() != octreeIterator.getCurrentListTargets());
            if( useTargetSource ){
                do{
                    const int nbParticlesInLeaf = (octreeIterator.getCurrentListSrc()->getSize() + octreeIterator.getCurrentListTargets()->getSize());
                    file.write((const char*)&nbParticlesInLeaf,sizeof(int));
                    ClassProptotype::PutParticles( &file, octreeIterator.getCurrentListSrc()->data(), octreeIterator.getCurrentListSrc()->getSize());
                    ClassProptotype::PutParticles( &file, octreeIterator.getCurrentListTargets()->data(), octreeIterator.getCurrentListTargets()->getSize());

                    ++nbLeaf;
                    if( maxParticlesInLeaf < nbParticlesInLeaf) maxParticlesInLeaf = nbParticlesInLeaf;
                } while(octreeIterator.moveRight());
            }
            else{
                do{
                    const int nbParticlesInLeaf = octreeIterator.getCurrentListSrc()->getSize();
                    file.write((const char*)&nbParticlesInLeaf,sizeof(int));
                    ClassProptotype::PutParticles( &file, octreeIterator.getCurrentListSrc()->data(), octreeIterator.getCurrentListSrc()->getSize());
                    ++nbLeaf;
                    if( maxParticlesInLeaf < nbParticlesInLeaf) maxParticlesInLeaf = nbParticlesInLeaf;
                } while(octreeIterator.moveRight());
            }

            const std::ofstream::pos_type currentPos = file.tellp();
            file.seekp(posNbLeaf);
            file.write((const char*)&nbLeaf,sizeof(int));
            file.write((const char*)&maxParticlesInLeaf,sizeof(int));
            file.seekp(currentPos);
        }

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = tree.getHeight() - 1 ; idxLevel > 1 ; --idxLevel ){
            int nbCells = 0;
            const std::ofstream::pos_type posNbCells = file.tellp();
            file.write((const char*)&nbCells,sizeof(int));

            do{
                const MortonIndex mindex = octreeIterator.getCurrentGlobalIndex();
                file.write((const char*)&mindex,sizeof(MortonIndex));;
                ClassProptotype::PutCell( &file, octreeIterator.getCurrentCell());
                ++nbCells;
            } while(octreeIterator.moveRight());

            const std::ofstream::pos_type currentPos = file.tellp();
            file.seekp(posNbCells);
            file.write((const char*)&nbCells,sizeof(int));
            file.seekp(currentPos);

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();
        }


        file.flush();
        file.close();

        return true;
    }


    /** To load from memory */
    template <class OctreeClass, class CellClass, class ParticleClass, class ClassProptotype >
    static bool Load(const char filename[], OctreeClass& tree){
        std::ifstream file(filename, std::ifstream::binary | std::ifstream::in );

        if(!file.good()){
            return false;
        }

        size_t sizeof_freal = 0;
        file.read((char*)&sizeof_freal, sizeof(size_t));
        if( sizeof_freal != sizeof(FReal)){
            std::cerr << "Error Freal do not coincide with file type:\n";
            std::cerr << "In file : " << sizeof_freal << " Real : " << sizeof(FReal) << "\n";
            return false;
        }

        int treeHeight = 0;
        file.read((char*)&treeHeight, sizeof(int));
        int treeSubHeight = 0;
        file.read((char*)&treeSubHeight, sizeof(int));
        FReal boxWidth = 0;
        file.read((char*)&boxWidth, sizeof(FReal));
        F3DPosition center;
        file.read((char*)&center, sizeof(F3DPosition));

        tree.~OctreeClass();
        new (&tree) OctreeClass(treeHeight,treeSubHeight,boxWidth,center);

        {
            int nbLeaf = 0;
            file.read((char*)&nbLeaf, sizeof(int));
            int maxParticlesInLeaf = 0;
            file.read((char*)&maxParticlesInLeaf, sizeof(int));

            ParticleClass* const particles = reinterpret_cast<ParticleClass*>(new char[maxParticlesInLeaf * sizeof(ParticleClass)]);

            for(int idxLeaf = 0 ; idxLeaf < nbLeaf ; ++idxLeaf){
                int particlesInLeaf = 0;
                file.read((char*)&particlesInLeaf, sizeof(int));
                ClassProptotype::GetParticles(&file, particles, particlesInLeaf);
                for(int idxParticle = 0 ; idxParticle < particlesInLeaf ; ++idxParticle){
                    tree.insert(particles[idxParticle]);
                }
            }

            delete[] reinterpret_cast<char*>(particles);
        }

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = tree.getHeight() - 1 ; idxLevel > 1 ; --idxLevel ){
            int nbCells = 0;
            file.read((char*)&nbCells, sizeof(int));

            do{
                MortonIndex mindex;
                file.read((char*)&mindex, sizeof(MortonIndex));
                if(mindex != octreeIterator.getCurrentGlobalIndex()){
                    std::cerr << "Error indexes are different\n";
                    return false;
                }

                ClassProptotype::GetCell(&file,octreeIterator.getCurrentCell());
                --nbCells;
            } while(octreeIterator.moveRight());

            if(nbCells != 0){
                std::cerr << "Wrong number of cells at level " << idxLevel << "\n";
                return false;
            }

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;
        }

        file.close();

        return true;
    }

};


#endif // FTREEIO_HPP
