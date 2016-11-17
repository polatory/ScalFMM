// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _FCOORDCOLOUR_HPP_
#define _FCOORDCOLOUR_HPP_

class FCoordColour {
    
public:
    enum {range = 3*3*3};

    static int coord2colour(const FTreeCoordinate& coord) {
        return (coord.getX() % 3) * 9 
            +  (coord.getY() % 3) * 3
            +  (coord.getZ() % 3);
    }
};

#endif
