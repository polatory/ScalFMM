#ifndef FOMPPRIORITIES_HPP
#define FOMPPRIORITIES_HPP

#include "../../Utils/FGlobal.hpp"

class FOmpPriorities{
    int insertionPositionP2M;
    int insertionPositionM2M;

    int insertionPositionP2MSend;
    int insertionPositionM2MSend;

    int insertionPositionM2L;
    int insertionPositionM2LExtern;
    int insertionPositionM2LLastLevel;
    int insertionPositionL2L;
    int insertionPositionL2P;
    int insertionPositionP2P;
    int insertionPositionP2PExtern;

    int treeHeight;

    int maxprio;

public:
    FOmpPriorities(const int inTreeHeight) :
        insertionPositionP2M(0), insertionPositionM2M(0), insertionPositionP2MSend(0),
        insertionPositionM2MSend(0), insertionPositionM2L(0), insertionPositionM2LExtern(0),
        insertionPositionM2LLastLevel(0), insertionPositionL2L(0), insertionPositionL2P(0), insertionPositionP2P(0),
        insertionPositionP2PExtern(0), treeHeight(inTreeHeight) , maxprio(0){
        if(inTreeHeight > 2){
            int incPrio = 0;

            insertionPositionP2MSend = incPrio++;
            insertionPositionP2M     = incPrio++;

            insertionPositionM2MSend = incPrio++;
            insertionPositionM2M     = incPrio++;

            insertionPositionM2L     = incPrio++;
            insertionPositionM2LExtern = incPrio++;

            insertionPositionL2L     = incPrio++;

            incPrio += (treeHeight-3) - 1;   // M2L is done treeHeight-2 times
            incPrio += (treeHeight-3) - 1;   // M2L is done treeHeight-2 times
            incPrio += (treeHeight-3) - 1;   // L2L is done treeHeight-3 times

            insertionPositionP2P       = incPrio++;

            insertionPositionM2LLastLevel = incPrio++;

            insertionPositionL2P     = incPrio++;

            insertionPositionP2PExtern = incPrio++;

            assert(incPrio == 8 + (treeHeight-3) + (treeHeight-3) + (treeHeight-3));
            maxprio = incPrio;
        }
        else{
            int incPrio = 0;

            insertionPositionP2MSend = -1;
            insertionPositionP2M     = -1;

            insertionPositionM2MSend = -1;
            insertionPositionM2M     = -1;

            insertionPositionM2L     = -1;
            insertionPositionM2LExtern = -1;
            insertionPositionM2LLastLevel = -1;

            insertionPositionL2L     = -1;

            insertionPositionP2P     = incPrio++;
            insertionPositionP2PExtern = insertionPositionP2P;

            insertionPositionL2P     = -1;
            assert(incPrio == 1);
            maxprio = incPrio;
        }
    }

    int getMaxPrio() const{
        return maxprio;
    }

    int getInsertionPosP2M() const {
        return insertionPositionP2M;
    }
    int getInsertionPosM2M(const int /*inLevel*/) const {
        return insertionPositionM2M;
    }
    int getInsertionPosM2L(const int inLevel) const {
        return (inLevel==treeHeight-1? insertionPositionM2LLastLevel : insertionPositionM2L + (inLevel - 2)*3);
    }
    int getInsertionPosM2LExtern(const int inLevel) const {
        return (inLevel==treeHeight-1? insertionPositionM2LLastLevel : insertionPositionM2LExtern + (inLevel - 2)*3);
    }
    int getInsertionPosL2L(const int inLevel) const {
        return insertionPositionL2L + (inLevel - 2)*3;
    }
    int getInsertionPosL2P() const {
        return insertionPositionL2P;
    }
    int getInsertionPosP2P() const {
        return insertionPositionP2P;
    }
    int getInsertionPosP2PExtern() const {
        return insertionPositionP2PExtern;
    }
};

#endif // FOMPPRIORITIES_HPP

