#ifndef FEQUALIZE_HPP
#define FEQUALIZE_HPP


#include <iostream>
#include <vector>

/**
 * This class proposes a method to distribute an interval own by a worker
 * to some others.
 * It returns a vector of Package which tell what interval to sent to which process.
 * The algorithm works only if the objectives (the intervals that the workers shoud obtain
 * are exclusive and given in ascendant order).
 * Also each of the element from the current interval must belong to someone!
 * Finally the current worker is included in the Package if its objective interval
 * is related to its starting interval.
 */
class FEqualize {
    // Just return the min
    static size_t Min(const size_t v1, const size_t v2){
        return v1 < v2 ? v1 : v2;
    }

public:
    /** To represent an interval to proceed */
    struct Package{
        int idProc;
        size_t elementFrom;
        size_t elementTo;
    };


    /**
     * To know what to send to who.
     * @param myCurrentInterval current process interval
     * @param allObjectives the intevals that each process should have (in ascendant order, exclusive)
     * @return the package that the current worker should sent to others
     */
    static std::vector<Package> GetPackToSend(const std::pair<size_t, size_t> myCurrentInterval,
                                              const std::vector< std::pair<size_t,size_t> >& allObjectives){
        std::vector<Package> packToSend;

        int idxProc = 0;

        // Find the first proc to send to
        while( idxProc != allObjectives.size()
               && allObjectives[idxProc].second < myCurrentInterval.first){
            idxProc += 1;
        }

        // We will from the first element for sure
        size_t currentElement = 0;

        // Check each proc to send to
        while( idxProc != allObjectives.size()
               && allObjectives[idxProc].first < myCurrentInterval.second){
            Package pack;
            pack.idProc = idxProc;
            // The current element must start where the previous one end
            // We can assert currentElement == allObjectives[idxProc].first - myCurrentInterval.first
            pack.elementFrom = currentElement;
            pack.elementTo   = Min(allObjectives[idxProc].second , myCurrentInterval.second) - myCurrentInterval.first;
            // Next time give from the previous end
            currentElement   = pack.elementTo;
            packToSend.push_back(pack);
            // Progress
            idxProc += 1;
        }
        // We can assert that currentElement == myCurrentInterval.second

        return packToSend;
    }
};

#endif // FEQUALIZE_HPP
