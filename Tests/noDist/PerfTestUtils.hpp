#ifndef _PERFTESTUTILS_HPP_
#define _PERFTESTUTILS_HPP_

#include <string>

#include "Utils/FTic.hpp"
#include "Files/FFmaGenericLoader.hpp"

#include "Containers/FOctree.hpp"

/**
 * \brief Store the PerfTest program parameters.
 */
struct FPerfTestParams {
    int subTreeHeight = 2; ///< Subtree height.
    int treeHeight = 5;    ///< Tree height.
    int nbThreads = 1;     ///< Maximum number of threads (when used). 
    std::string filename = ""; ///< Particles file.
    std::string algo = "task"; ///< Algorithm to run.
    bool omp_static_schedule = false; ///< OpenMP static or dynamic schedule.
};


/**
 * \brief Base class for tree loaders.
 *
 * This class itself does not provide anything but a base on which to build tree
 * loaders. A tree loader should satisfy the following rules.
 *
 *    - Define the public typedefs : CellClass, ContainerClass, LeafClass,
 *      OctreeClass.
 *    - Provide public acces to a member of type OctreeClass _tree as the tree
 *      that is loaded.
 *    - Tree loading must happen at construction.
 *    - It may provide any other members or typdefs required by a special
 *      FKernelLoader or FAlgoLoader.
 *
 * For convenience, this class provides a timer and a basic loadTree method that
 * should be enough to load a tree from and FMA file.
 *
 * \note It is not mandatory that a loader inherit from this class. It must
 * however follow the aforementioned rules.
 */
class FTreeLoader {
public:
    /// A timer
    /** Is used to time the loadTree method.
     */
    FTic time;
protected:

    /**
     * \brief A template which type is always false.
     *
     * This template is only expanded by the compiler when it is requested
     * (ie. the compiler will not try to optimize out its value.). Must be used
     * to create false static_assert to catch unintended use of a template.
     */
    template <typename... Args>
    struct false_type {
        bool value = false;
    };  

    /**
     * \brief Failure method for unimplemented loadTree templates.
     *
     * This template will catch unspecialised call to the loadTree method and
     * will cause the compilation to fail with a (somewhat) meaningfull message.
     */
    template <typename...Args>
    void loadTree(Args...) {
        static_assert(false_type<Args...>::value,
                      "I don't know how to load this tree with this loader...");
    }


    /**
     * \brief Simple method to load a tree from a FMA file.
     *
     * The template parameters are usualy guessed by the compiler.
     *
     * \tparam OctreeClass The class of the tree to fill.
     * \tparam FReal The floating point type.
     *
     * \param loader The file loader to read from the file.
     * \param tree The tree to be filled.
     */
    template <class OctreeClass, typename FReal>
    void loadTree(FFmaGenericLoader<FReal>& loader, OctreeClass& tree) {
        std::cout << "Creating & inserting particles" << std::flush;

        time.tic();

        FPoint<FReal> position;
        FReal physicalValue = 0.0;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart) {
            // Read particle per particle from file
            loader.fillParticle(&position,&physicalValue);
            // put particle in octree
            tree.insert(position, idxPart, physicalValue);
        }

        time.tac();
        std::cout << " Done  (" << time.elapsed() << "s)." << std::endl;
    }

};

/**
 * \brief Base class for kernel loaders.
 *
 * This class itself does not provide anything but a base on which to build
 * kernel loaders. A kernel loader should satisfy the following rules.
 *
 *    - Define the public typedefs : TreeLoader, KernelClass.
 *    - Provide public acces to a member of type Kernelclass _kernel as the
 *      kernel that is loaded.
 *    - Kernel loading must happen at construction.
 *    - It may provide any other members or typdefs required by a special
 *      FAlgoLoader.
 *
 * For convenience, this class provides a timer.
 *
 * \tparam _TreeLoader The tree loader that was used.
 *
 * \note It is not mandatory that a loader inherit from this class. It must
 * however follow the aforementioned rules.
 */
template<class _TreeLoader>
class FKernelLoader {
    /// The tree loader that was used (see FTreeLoader).
    using TreeLoader = _TreeLoader;
public:
        FTic time;
};

/**
 * \brief Base class for algorithm loaders.
 *
 * This class itself does not provide anything but a base on which to build
 * algorithm loaders. A kernel loader should satisfy the following rules.
 *
 *    - Define the public typedefs : TreeLoader, KernelLoader.
 *    - Provide public acces to a member of type
 *      \link TreeLoader Treeloader::OctreeClass* \endlink` _algo`
 *      as the algorithm that is loaded. This pointer should be valid from the
 *      end of the ::run method to the destruction of the loader.
 *    - It may provide any other members or typdefs.
 *
 * For convenience, this class provides a timer.
 *
 * \tparam _TreeLoader The tree loader that was used.
 * \tparam _KernelLoader The kernel loader *template* that was used, the
 *         KernelLoader type will then be _KernelLoader<_TreeLoader>.
 *
 * \note It is not mandatory that a loader inherit from this class. It must
 * however follow the aforementioned rules.
 */
template <class _TreeLoader, template<typename> class _KernelLoader>
class FAlgoLoader {
    /// The tree loader that was used (see FTreeLoader).
    using TreeLoader = _TreeLoader;
    /// The kernel loader that was used (see FKernelLoader).
    using KernelLoader = _KernelLoader<TreeLoader>;
public:
    /// A timer.
    FTic time;
    /// Method that runs the algorithm.
    virtual void run() = 0;
};





#endif
