#ifndef FNOCOPYABLE_HPP
#define FNOCOPYABLE_HPP

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* This class has to be inherited to forbid copy
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
