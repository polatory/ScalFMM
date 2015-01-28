// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas
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

#ifndef FALGORITHMTIMERS_HPP
#define FALGORITHMTIMERS_HPP

/**
 * @brief This class provide a way for the different algorithms to
 * store the time spent in each operator.
 *
 */
class FAlgorithmTimers{
protected:
    FTic * Timers;

public:

    enum FTimers{
        P2MTimer=0,
        M2MTimer=1,
        M2LTimer=2,
        L2LTimer=3,
        L2PTimer=4,
        P2PTimer=5,
        NearTimer=6,
    };

    const int nbTimers = 7;

    FAlgorithmTimers() : Timers(nullptr)
    {
        Timers = new FTic[nbTimers];
    }

    ~FAlgorithmTimers(){
        delete[] Timers;
    }

    const FTic * getAllTimers() const {
        return Timers;
    }

    int getNbOfTimerRecorded() const {
        return nbTimers;
    }

    double getTime(FTimers OpeTimer) const{
        //assert to verify size
        return Timers[OpeTimer].elapsed();
    }

    double getCumulatedTime(FTimers OpeTimer) const{
        //assert to verify size
        return Timers[OpeTimer].cumulated();
    }

};

#endif // FALGORITHMTIMERS_HPP
