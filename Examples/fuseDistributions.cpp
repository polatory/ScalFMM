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

/**
 * \brief Print programe usage
 *
 * \param progname The name of the program.
 *
 * \note Use `argv[0]` as progname.
 */
void usage(const std::string& progname) {
    std::size_t start = progname.find_last_of('/');
    std::string name = progname.substr(start+1);
    std::cout <<
        "usage: " << name <<
        " --file [[-s scale] [-c cx:cy:cz] [-g gx:gy:gz]] filename"
        " -fout output_file"
        " [--extra-length length]"
        "\n"
        "\n"
        "Fuses multiple particle distributions into a bigger one."
        "\n"
        "\n"
        "Options:\n"
        "  -fout output_file\n"
        "    The output file name, must hase .bfma or .fma extension\n"
        "\n"
        "  --file [opts] filename [opts]\n"
        "    Add a .fma or .bfma distibution file. Multiple files may be specified by\n"
        "    adding more --file options. 'opts' is a combination of:\n"
        "      -s scale\n"
        "        Scale the distribution by 'scale' factor.\n"
        "      -c cx:cy:cz\n"
        "        Center the distribution at given coordinates. cx, cy and cz are\n"
        "        floating point numbers.\n"
        "      -g gx:gy:gz\n"
        "        Duplicate the distribution inside a grid of gx by gy by gz dimensions.\n"
        "        gx, gy and gz are integers. The grid center is governed by the -c \n"
        "        option.\n"
        "      -r rx:ry:rz\n"
        "        Rotate the distribution around its x, y and z axes. The rotation \n"
        "        center is the distribution center. rx, ry and rz are in radians.\n"
        "\n"
        "  --extra-length length\n"
        "    Length to be added to the final box width.\n"
        "\n"
        "  --help\n"
        "    Print this message."
        "\n"
        ;
}


using FReal = double;

/**
 * \brief Particle object
 */
struct Particle {
    FPoint<FReal> pos; ///< Particle position in space
    FReal val;         ///< Particle physical value
};

/**
 * \brief Distribution options
 *
 * Holds the tranformations that must be applied to a particle distribution.
 */
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

/**
 * \brief Main program parameters
 */
struct parameters {
    /// Output distribution file
    std::string output_filename;
    /// Particles distributions to combine
    std::vector<distribution> distributions;
    /// Length to add to the output distribution box
    FReal extra_length;
};

/**
 * \brief Helper method to parse the arguments that apply to one distribution.
 *
 * \param args A vector holding the command line arguments
 * \param i    Index of the current argument in args
 *
 * \return A vector of distribution objects
 */
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
        if(false) {// false for alignment of the other cases
        } else if(args[i] == "-s") {
            // Parse scale argument
            ++i;
            sstr.str(args.at(i));
            sstr >> dist.scale;
        } else if (args[i] == "-c") {
            // Parse offset from center argument
            ++i;
            char c; // Used to discard the ':' from argument format
            sstr.str(args.at(i));
            sstr >> dist.offset[0] >> c >> dist.offset[1] >> c >> dist.offset[2];
        } else if(args[i] == "-g") {
            // Parse grid layout
            ++i;
            char c; // Used to discard the ':' from argument format
            sstr.str(args.at(i));
            sstr >> gx >> c >> gy >> c >> gz;
        } else if(args[i] == "-r") {
            // Parse rotation angle
            ++i;
            char c; // Used to discard the ':' from argument format
            sstr.str(args.at(i));
            sstr >> dist.rot[0] >> c >> dist.rot[1] >> c >> dist.rot[2];
        } else {
            if(dist.filename != "") {
                --i;
                break;
            }
            dist.filename = args[i];
        }
        ++i;
    }

    if(i < args.size() && args[i] == "--file") {
        --i;
    }

    // If a grid layout was specified, duplicate the distribution for each cell
    // in the grid
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

/**
 * \brief Main argument parse function
 *
 * \param args A vector holding the command line arguments
 */
parameters parse(const std::vector<std::string>& args) {
    parameters params;
    std::stringstream sstr; // Use to convert strings to numbers
    for(std::size_t i = 1; i < args.size(); ++i) {
        if(args[i] == "--help") {
            usage(args[0]);
            exit(0);
        } else if(args[i] == "--file") {
            // Parse a new distribution specification
            ++i;
            auto ds = subparse_file(args, i);
            params.distributions.insert(params.distributions.end(),
                                        ds.begin(), ds.end());
        } else if(args[i] == "--extra-length") {
            // Parse extra length parameter
            ++i;
            sstr.str(args.at(i));
            sstr >> params.extra_length;
        } else if(args[i] == "-fout") {
            // Parse output file name
            ++i;
            params.output_filename = args.at(i);
        } else {
            std::cerr << "Unknown or misplaced parameters: " << args[i] << '\n';
        }
    }
    return params;
}

/**
 * \brief Rotate a particle around the origin
 *
 * \param p    particle to rotate
 * \param dist distribution transformation specification
 */
void rotate(Particle& p, const distribution& dist) {
    // Rotate around x axis
    if(dist.rot[0] > 1e-5 || dist.rot[0] < -1e-5) {
        FReal alpha = dist.rot[0];
        p.pos[1] = p.pos[1] * cos(alpha) - p.pos[2] * sin(alpha);
        p.pos[2] = p.pos[1] * sin(alpha) + p.pos[2] * cos(alpha);
    }
    // Rotate around y axis
    if(dist.rot[1] > 1e-5 || dist.rot[1] < -1e-5) {
        FReal alpha = dist.rot[1];
        p.pos[0] =  p.pos[0] * cos(alpha) + p.pos[2] * sin(alpha);
        p.pos[2] = -p.pos[0] * sin(alpha) + p.pos[2] * cos(alpha);
    }
    // Rotate around z axis
    if(dist.rot[2] > 1e-5 || dist.rot[2] < -1e-5) {
        FReal alpha = dist.rot[2];
        p.pos[0] = p.pos[0] * cos(alpha) - p.pos[1] * sin(alpha);
        p.pos[1] = p.pos[0] * sin(alpha) + p.pos[1] * cos(alpha);
    }
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

        // Load each particle and apply the distribution transformations
        Particle p;
        for(std::size_t i = 0; i < count; ++i) {
            loader.fillParticle(&p.pos, &p.val);
            // Move distribution center to origin
            p.pos -= center;
            // Scale distribution
            p.pos *= dist.scale;
            // Rotate distribution
            rotate(p, dist);
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
