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
// 
// @SCALFMM_PRIVATE
// 
#ifndef FMATRIXIO_HPP
#define FMATRIXIO_HPP

// std includes
#include <numeric>
#include <stdexcept>
#include <string>
#include <sstream>
#include <fstream>
#include <typeinfo>

// ScalFMM includes
#include "Utils/FGlobal.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FMath.hpp"

/**
 * @brief The FMatrixIO class
 * Beware! Only handles sizes (nrows/ncols) representable by int32. 
 */
class FMatrixIO {

public:

    // Write matrix to binary file.
    template<class FReal>
    static void write(const FSize nrows, const FSize ncols, FReal* matrix, const std::string filename){

        // store into binary file
        std::ofstream stream(filename.c_str(),
                             std::ios::out | std::ios::binary | std::ios::trunc);
        if (stream.good()) {
            stream.seekp(0);
            // 1) write number of rows (int)
            //is_int(nrows);
            int _nrows = int(nrows);
            stream.write(reinterpret_cast<char*>(&_nrows), sizeof(int));
            // 2) write number of cols (int)
            //is_int(ncols);
            int _ncols = int(ncols);
            stream.write(reinterpret_cast<char*>(&_ncols), sizeof(int));
            // 1) write matrix
            stream.write(reinterpret_cast<char*>(matrix), sizeof(FReal)*nrows*ncols);
        }   
        else throw std::runtime_error("File could not be opened to write");
        stream.close();

    }

    // Read matrix from binary file. 
    // Also allocates memory! Please set matrix to nullptr before calling this function. 
    // nrows is just here to control the values that is read from file
    // ncols  is updated with the value read in the file (because it is not always known!)
    template<class FReal>
    static void read(const FSize nrows, FSize &ncols, FReal* &matrix, const std::string filename, const bool readSize = true){

        // start reading process
        if (matrix) throw std::runtime_error("Matrix is already set!");

        std::ifstream stream(filename.c_str(),
                             std::ios::in | std::ios::binary | std::ios::ate);
        const std::ifstream::pos_type size = stream.tellg();
        if (size<=0) {
            std::cout << "Info: The requested binary file " << filename
                      << " does not yet exist. Compute it now ... " << std::endl;
            return;
        } 
        if (stream.good()) {
            stream.seekg(0);
            if(readSize){
                // 1) read number of rows (int)
                int _nrows;
                stream.read(reinterpret_cast<char*>(&_nrows), sizeof(int));
                //is_int(nrows);
                if (_nrows!=int(nrows)) throw std::runtime_error("Read nrows and input nrows do not correspond");
                // 2) read number of cols (int)
                int _ncols;
                stream.read(reinterpret_cast<char*>(&_ncols), sizeof(int));
                //is_int(ncols);
                ncols=_ncols;
            }
            // 3) read matrix
            matrix = new FReal[nrows*ncols];
            stream.read(reinterpret_cast<char*>(matrix), sizeof(FReal)*nrows*ncols);
        }   else throw std::runtime_error("File could not be opened to read");
        stream.close();

    }

    template<class FReal>
    static bool read(const std::string filename, FReal** matrix, int* nbRows, int* nbCols){
        std::ifstream stream(filename.c_str(), std::ios::in | std::ios::binary | std::ios::ate);

        if (!stream.good()) {
            return false;
        }

        const std::ifstream::pos_type size = stream.tellg();
        if (size<=0) {
            std::cout << "Info: The requested binary file " << filename
                      << " does not yet exist. Compute it now ... " << std::endl;
            return false;
        }

        stream.seekg(0);

        stream.read(reinterpret_cast<char*>(nbRows), sizeof(int));
        stream.read(reinterpret_cast<char*>(nbCols), sizeof(int));

        delete[] (*matrix);
        (*matrix) = new FReal[(*nbRows)*(*nbCols)];
        stream.read(reinterpret_cast<char*>((*matrix)), sizeof(FReal)*(*nbRows)*(*nbCols));

        stream.close();
        return true;
    }
};


#endif // FMATRIXIO_HPP
