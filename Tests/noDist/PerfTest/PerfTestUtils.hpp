// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _PERFTESTUTILS_HPP_
#define _PERFTESTUTILS_HPP_

#include <string>

#ifdef SCALFMM_USE_MPI
#include "Utils/FMpi.hpp"
#endif

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
    int  omp_chunk_size = 0;   ///< OpenMP chunk size for basic algorithm (FFmmAlgorithmThread)
#ifdef SCALFMM_USE_MPI
    FMpi* mpiContext = nullptr;
#endif
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
    /// A timer used to time the loadTree method.
    FTic time;
protected:

    /**
     * \brief Load a tree from a file.
     *
     * \param loader The file loader to read from the file.
     * \param tree The tree to be filled.
     */
    virtual void loadTree() = 0;
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
    /// A timer
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

    /// Additionnal information for specific algorithm loader.
    /**  
     * The string should be formated as a key:value list separated by spaces.
     * For instance : "key1:value1 key2:value2 ". It may be a good idea to add a
     * space at the end of the string.
     */
    virtual std::string getRunInfoString() const {
        return "";
    }
};





#endif
