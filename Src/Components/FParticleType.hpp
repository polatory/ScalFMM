// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================
#ifndef FPARTICLETYPE_HPP
#define FPARTICLETYPE_HPP

/**
 * @brief The FParticleType enum is to make a difference between Target and Source (Tsm)
 */
enum class FParticleType {
    FParticleTypeSource = 0,
    FParticleTypeTarget = 1
};

#endif // FPARTICLETYPE_HPP
