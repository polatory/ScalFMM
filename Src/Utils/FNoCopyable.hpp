// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FNOCOPYABLE_HPP
#define FNOCOPYABLE_HPP

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* This class has to be inherited to forbid copy
* @todo use C++0x ?
*/
class FNoCopyable {
private:
        /** Forbiden copy constructor */
        FNoCopyable(const FNoCopyable&);
        /** Forbiden copy operator */
        FNoCopyable& operator=(const FNoCopyable&);
protected:
        /** Empty constructor */
        FNoCopyable(){}
};

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* This class has to be inherited to forbid assignement
*/
class FNoAssignement {
private:
        /** Forbiden copy operator */
        FNoAssignement& operator=(const FNoAssignement&);
protected:
        /** Empty constructor */
        FNoAssignement(){}
};

#endif // FNOCOPYABLE_HPP
