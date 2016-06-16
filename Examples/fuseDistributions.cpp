/**
 * \file
 * \brief Fuses FMA files to create a new distribution
 *
 * \author Quentin Khan
 * \copyright ScalFmm 2016 INRIA
 * \copyright [CeCILL-C licence](http://www.cecill.info)
 *
 *
 */

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


#include "Files/FFmaGenericLoader.hpp"

using FReal = double;

struct Particle {
    FPoint<FReal> pos;
    FReal val;
};

/// Distribution options
struct distribution {
    /// Distribution filename
    std::string filename = "";
    /// Distribution offset from center
    FPoint<FReal> offset = {0,0,0};
    /// Distribution rotation around its center
    FPoint<FReal> rot = {0,0,0};
    /// Distribution scale factor
    FReal scale = 1;
};



struct parameters {
    std::string output_filename;
    std::vector<distribution> distributions;
    FReal extra_length;
};

std::vector<distribution> subparse_file(const std::vector<std::string>& args, std::size_t& i) {
    std::stringstream sstr;
    // Grid layout
    unsigned int gx = 1, gy = 1, gz = 1;
    // Final distributions, one per grid part
    std::vector<distribution> distributions;
    // Distributions options
    distribution dist;

    while(i < args.size() && args[i] != "--file") {
        sstr.clear();
        if(false) {
        } else if(args[i] == "-s") {
            ++i;
            sstr.str(args.at(i));
            sstr >> dist.scale;
        } else if (args[i] == "-c") {
            ++i;
            char c; // Used to get the ':' from argument format
            sstr.str(args.at(i));
            sstr >> dist.offset[0] >> c >> dist.offset[1] >> c >> dist.offset[2];
        } else if(args[i] == "-g") {
            ++i;
            char c; // Used to get the ':' from argument format
            sstr.str(args.at(i));
            std::cerr << sstr.str() << '\n';
            sstr >> gx >> c >> gy >> c >> gz;
            std::cerr << gx << ' ' << gy << ' ' << gz << '\n';
        } else {
            if(dist.filename != "") {
                --i;
                break;
            }
            dist.filename = args[i];
        }
        ++i;
    }

    if(gx > 1 || gy > 1 || gz > 1) {

        // Compute offset of lowest left grid offset
        FFmaGenericLoader<FReal> loader(dist.filename);
        FReal box_width = loader.getBoxWidth() * dist.scale;
        dist.offset[0] -= (gx-1) * box_width / 2;
        dist.offset[1] -= (gy-1) * box_width / 2;
        dist.offset[2] -= (gz-1) * box_width / 2;

        // Create one distribution for each part of the grid layout
        for(unsigned int x = 0; x < gx; ++x) {
            for(unsigned int y = 0; y < gy; ++y) {
                for(unsigned int z = 0; z < gz; ++z) {
                    distribution tmp_dist = dist;
                    tmp_dist.offset[0] += x * box_width;
                    tmp_dist.offset[1] += y * box_width;
                    tmp_dist.offset[2] += z * box_width;
                    distributions.push_back(tmp_dist);
                }
            }
        }
    } else {
        distributions.push_back(dist);
    }

    return distributions;
}


parameters parse(const std::vector<std::string>& args) {
    parameters params;
    std::stringstream sstr;
    for(std::size_t i = 1; i < args.size(); ++i) {
        if(false) { // To reorder cases without breaking format
        } else if(args[i] == "--file") {
            ++i;
            auto ds = subparse_file(args, i);
            params.distributions.insert(params.distributions.end(),
                                        ds.begin(), ds.end());
        } else if(args[i] == "--extra-length") {
            ++i;
            sstr.str(args.at(i));
            sstr >> params.extra_length;
        } else if(args[i] == "-fout") {
            ++i;
            params.output_filename = args.at(i);
        } else {
            std::cerr << "Unknown or misplaced parameters: " << args[i] << '\n';
        }
    }
    return params;
}







int main(int argc, char** argv) {
    auto params = parse({argv,argv+argc});

    // Fail early if output file raises an error
    FFmaGenericWriter<FReal> writer(params.output_filename);

    // Fuse particle distributions
    std::vector<Particle> particles;
    FReal axis_max = 0;

    for(distribution& dist : params.distributions) {
        // Load particles into array
        FFmaGenericLoader<FReal> loader(dist.filename);
        const std::size_t count = loader.getParticleCount();
        // Particle array: x1, y1, z1, val1, x2, y2...
        particles.reserve(particles.size() + count);

        FPoint<FReal> center = loader.getBoxCenter();

        // Temp particle
        Particle p;
        for(std::size_t i = 0; i < count; ++i) {
            loader.fillParticle(&p.pos, &p.val);
            // Move distribution center to origin
            p.pos -= center;
            // Scale distribution
            p.pos *= dist.scale;
            // Rotate distribution
            // TODO

            // Move to new position
            p.pos += dist.offset;
            // Add particle to list
            particles.push_back(p);

            // Save particle x,y,z min/max to compute final box
            axis_max = std::max(std::abs(p.pos[0]), axis_max);
            axis_max = std::max(std::abs(p.pos[1]), axis_max);
            axis_max = std::max(std::abs(p.pos[2]), axis_max);
        }
    }


    // Write final distribution
    FPoint<FReal> center(0,0,0);
    // Compute final box width
    FReal box_width = 2 * (axis_max + params.extra_length);

    // Write header
    writer.writeHeader(center, box_width, particles.size(), 8, 4);

    // Write all particles

    // Buffer avoids duplicating particle vector
    std::vector<FReal> buffer;
    buffer.reserve(4*1024); // Avoid reallocations, size is a multiple of 4

    auto cur = particles.begin();
    auto sentinel = particles.end();

    // Fill and write buffer until we're done
    while(cur != sentinel) {
        buffer.clear();
        while(buffer.size() != buffer.capacity() && cur != sentinel) {
            buffer.push_back(cur->pos[0]);
            buffer.push_back(cur->pos[1]);
            buffer.push_back(cur->pos[2]);
            buffer.push_back(cur->val);
            ++cur;
        }
        writer.writeArrayOfReal(buffer.data(), 4, buffer.size()/4);
    }

}
