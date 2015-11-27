// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FSVGRECT_HPP
#define FSVGRECT_HPP

// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"
#include "Utils/FAssert.hpp"

#include <vector>
#include <algorithm>


class FSvgRect {
protected:
    int dim;
    int svgSize;
    int margin;

    int rectCounter;
    double pixelperunknown;
    int descriptiony;

    FILE* fsvg;

    std::vector<std::tuple<int, int, int>> colors;

public:
    FSvgRect(const char inFilename[], const int inDim,
             const int inSvgSize = 2048, const int inMargin = 50)
        : dim(inDim), svgSize(inSvgSize), margin(inMargin), rectCounter(0),
          pixelperunknown(0), descriptiony(0), fsvg(NULL){
        fsvg = fopen(inFilename, "w");
        FAssertLF(fsvg);

        pixelperunknown = double(svgSize-margin-margin*2)/double(dim);

        fprintf(fsvg, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
        fprintf(fsvg, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n"
                "<svg xmlns=\"http://www.w3.org/2000/svg\"  xmlns:xlink=\"http://www.w3.org/1999/xlink\"  width=\"%d\" height=\"%d\">\n",
                svgSize, svgSize);
        fprintf(fsvg, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"fill:rgb(%d,%d,%d);stroke-width:1;stroke:rgb(0,0,0)\" />\n",
                            0, 0,svgSize, svgSize, 255, 255, 255);

        descriptiony = int(svgSize - margin);
        fprintf(fsvg, "<text x=\"%d\" y=\"%d\" fill=\"red\" id=\"res\" font-size=\"55\">Description: </text>\n", 0, descriptiony);

    }

    FSvgRect(const char inDirname[], const char inFilename[], const int inDim,
             const int inSvgSize = 2048, const int inMargin = 50)
        : FSvgRect((std::string(inDirname)+std::string(inFilename)).c_str(), inDim, inSvgSize, inMargin) {
    }

    ~FSvgRect(){
        fprintf(fsvg, "</svg>");
        fclose(fsvg);
    }

    void addRect(const int inX, const int inY, const int inWidth, const int inHeight, const int inLevel = 0){
        if(int(colors.size()) <= inLevel){
            for(int idxColor = int(colors.size()) ; idxColor <= inLevel ; ++idxColor){
                colors.push_back(std::tuple<int, int, int>(255*drand48(),255*drand48(),255*drand48()));
            }
        }
        fprintf(fsvg, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"fill:rgb(%d,%d,%d);stroke-width:1;stroke:rgb(0,0,0)\" />\n",
                            int(margin + inX*pixelperunknown),
                            int(margin + inY*pixelperunknown),
                            int(inWidth*pixelperunknown),
                            int(inHeight*pixelperunknown),
                            std::get<0>(colors[inLevel]),std::get<1>(colors[inLevel]),std::get<2>(colors[inLevel]));

        rectCounter += 1;
    }

    template <class... T>
    void addRectWithLegend(const int inX, const int inY, const int inWidth, const int inHeight, int inLevel = -1){
        if(int(colors.size()) <= inLevel){
            for(int idxColor = int(colors.size()) ; idxColor <= inLevel ; ++idxColor){
                colors.push_back(std::tuple<int, int, int>(255*drand48(),255*drand48(),255*drand48()));
            }
        }

        fprintf(fsvg, "<text x=\"%d\" y=\"%d\" fill=\"red\" id=\"res%d\" font-size=\"55\">Description: X = %d, Y = %d, W = %d, H = %d, L = % d </text>\n",
                0, descriptiony + svgSize, rectCounter,
                inX, inY, inWidth, inHeight, inLevel);

        fprintf(fsvg, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"fill:rgb(%d,%d,%d);stroke-width:1;stroke:rgb(0,0,0)\" ",
                int(margin + inX*pixelperunknown),
                int(margin + inY*pixelperunknown),
                int(inWidth*pixelperunknown),
                int(inHeight*pixelperunknown),
                std::get<0>(colors[inLevel]),std::get<1>(colors[inLevel]),std::get<2>(colors[inLevel]));

        fprintf(fsvg, " onmousedown=\"evt.target.parentNode.getElementById('res%d').setAttribute('y', '%d');\" ",
                rectCounter,descriptiony);
        fprintf(fsvg, " onmouseup=\"evt.target.parentNode.getElementById('res%d').setAttribute('y', '%d');\" ",
                rectCounter,descriptiony + svgSize);

        fprintf(fsvg, "/>\n");

        rectCounter += 1;
    }
};

#endif // FSVGRECT_HPP

