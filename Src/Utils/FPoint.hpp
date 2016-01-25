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
//
#ifndef FPOINT_HPP
#define FPOINT_HPP

#include <array>
#include <iterator>
#include <ostream>

#include "FMath.hpp"
#include "FConstFuncs.hpp"

/** N-dimentional coordinates
 *
 * \author Berenger Bramas <berenger.bramas@inria.fr>, Quentin Khan <quentin.khan@inria.fr>
 *
 * A fixed size array that represents coordinates in space. This class adds
 * a few convenience operations such as addition, scalar multiplication and
 * division and formated stream output.
 *
 * \tparam Real The floating number type
 * \tparam Dim The space dimension
 **/
template<typename _Real, std::size_t _Dim = 3>
class FPoint : public std::array<_Real, _Dim> {
public:
    /// Floating number type
    using Real = _Real;
    /// Space dimension count
    constexpr static const std::size_t Dim = _Dim;

private:

    /// Type used in SFINAE to authorize arithmetic types only in template parameters
    template<class T>
    using must_be_arithmetic = typename std::enable_if<std::is_arithmetic<T>::value>::type*;
    /// Type used in SFINAE to authorize floating point types only in template parameters
    template<class T>
    using must_be_floating = typename std::enable_if<std::is_floating_point<T>::value>::type*;
    /// Type used in SFINAE to authorize integral types only in template parameters
    template<class T>
    using must_be_integral = typename std::enable_if<std::is_integral<T>::value>::type*;


    /** Recursive template for constructor */
    template<typename A = Real, typename... Args>
    void _fill_data(const A& arg, const Args... args) {
        this->data()[Dim-sizeof...(args)-1] = arg;
        _fill_data(args...);
    }

    /** Recursive template end condition for constructor */
    template<typename A = Real>
    void _fill_data(const A& arg) {
        this->data()[Dim-1] = arg;
    }

public:

    /** Default constructor */
    FPoint() : std::array<Real, Dim>{{0}} {};
    /** Copy constructor */
    FPoint(const FPoint&) = default;

    /** Copy constructor from other point type */
    template<typename A, must_be_arithmetic<A> = nullptr>
    FPoint(const FPoint<A, Dim>& other) {
        for(std::size_t i = 0; i < Dim; ++i) {
            this->data()[i] = other.data()[i];
        }
    }

    /** Constructor from array */
    FPoint(const Real array[Dim]) {
        for(std::size_t i = 0; i < Dim; ++i) {
            this->data()[i] = array[i];
        }
    }

    /** Constructor from args */
    template<typename A = Real, typename... Args>
    FPoint(const Real& arg, const Args... args) {
        static_assert(sizeof...(args)+1 == Dim, "FPoint argument list isn't the right length.");
        _fill_data(arg, args...);
    }

    /** Additive contructor, same as FPoint(other + add_value) */
    FPoint(const FPoint& other, const Real add_value) {
        for(std::size_t i = 0; i < Dim; ++i) {
            this->data()[i] = other.data()[i] + add_value;
        }
    }

    /** Copies #Dim values from an iterable container
     *
     * \param other Container that defines begin() and end() methods.
     */
    template<class T>
    void copy(const T& other) {

        auto other_it = other.begin();
        auto this_it  = this->begin();
        for(std::size_t i = 0; i < Dim; ++i, ++this_it, ++other_it) {
            *this_it = *other_it;
        }
    }

    /** Assignment operator
     *
     * \param other A FPoint object.
     */
    template<typename T, must_be_arithmetic<T> = nullptr>
    FPoint<Real, Dim>& operator=(const FPoint<T, Dim>& other) {
        this->copy(other);
        return *this;
    }


    /** Sets the point value */
    template<typename A = Real, typename... Args>
    void setPosition(const Real& X, const Args... YandZ) {
        static_assert(sizeof...(YandZ)+1 == Dim, "FPoint argument list isn't the right length.");
        _fill_data(X, YandZ...);
    }

    /** \brief Get x
     * \return this->data()[0]
     */
    Real getX() const{
        return this->data()[0];
    }


    /** \brief Get y
     * \return this->data()[1]
     */
    Real getY() const{
        return this->data()[1];
    }


    /** \brief Get z
     * \return this->data()[2]
     */
    Real getZ() const{
        return this->data()[2];
    }


    /** \brief Set x
     * \param the new x
     */
    void setX(const Real inX){
        this->data()[0] = inX;
    }


    /** \brief Set y
     * \param the new y
     */
    void setY(const Real inY){
        this->data()[1] = inY;
    }


    /** \brief Set z
     * \param the new z
     */
    void setZ(const Real inZ){
        this->data()[2] = inZ;
    }


    /** \brief Add to the x-dimension the inX value
     * \param  inXthe increment in x
     */
    void incX(const Real inX){
        this->data()[0] += inX;
    }


    /** \brief Add to the y-dimension the inY value
     * \param  in<<<<<<y the increment in y
     */
    void incY(const Real inY){
        this->data()[1] += inY;
    }


    /** \brief Add to z-dimension the inZ value
     * \param inZ the increment in z
     */
    void incZ(const Real inZ){
        this->data()[2] += inZ;
    }

    /** \brief Get a pointer on the coordinate of FPoint<Real>
     * \return the data value array
     */
    Real * getDataValue(){
        return this->data() ;
    }

    /** \brief Get a pointer on the coordinate of FPoint<Real>
     * \return the data value array
     */
    const Real *  getDataValue()  const{
        return this->data() ;
    }

    /** \brief Compute the distance to the origin
     * \return the norm of the FPoint
     */
    Real norm() const {
        return FMath::Sqrt(norm2()) ;
    }

    /** \brief Compute the distance to the origin
     * \return the square norm of the FPoint
     */
    Real norm2() const {
        Real square_sum = 0;
        for(std::size_t i = 0; i < Dim; ++i) {
            square_sum += this->data()[i]*this->data()[i];
        }
        return square_sum;
    }


    /** Addition assignment operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint& operator +=(const FPoint<T, Dim>& other) {
        auto other_it = other.begin();
        auto this_it  = this->begin();
        for(std::size_t i = 0; i < Dim; ++i, ++this_it, ++other_it) {
            *this_it += *other_it;
        }

        return *this;
    }

    /** Addition operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<Real, Dim> operator+(FPoint<Real, Dim> lhs, const FPoint<T, Dim>& rhs) {
        lhs += rhs;
        return lhs;
    }

    /** Scalar assignment addition */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint& operator+=(const T& val) {
        for(std::size_t i = 0; i < Dim; ++i) {
            this->data()[i] += val;
        }
        return *this;
    }

    /** Scalar addition */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<Real, Dim> operator+(FPoint<Real, Dim> lhs, const T& val) {
        lhs += val;
        return lhs;
    }

    /** Subtraction assignment operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint& operator -=(const FPoint<T, Dim>& other) {
        auto other_it = other.begin();
        auto this_it  = this->begin();
        for(std::size_t i = 0; i < Dim; i++, ++this_it, ++other_it) {
            *this_it -= *other_it;
        }
        return *this;
    }

    /** Subtraction operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<Real, Dim> operator-(FPoint<Real, Dim> lhs, const FPoint<T, Dim>& rhs) {
        lhs -= rhs;
        return lhs;
    }

    /** Scalar subtraction assignment */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint& operator-=(const T& val) {
        for(std::size_t i = 0; i < Dim; ++i) {
            this->data()[i] -= val;
        }
        return *this;
    }

    /** Scalar subtraction */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<Real, Dim> operator-(FPoint<Real, Dim> lhs, const T& rhs) {
        lhs -= rhs;
        return lhs;
    }


    /** Right scalar multiplication assignment operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint<Real, Dim>& operator *=(const T& val) {
        for(std::size_t i = 0; i < Dim; ++i) {
            this->data()[i] *= val;
        }
        return *this;
    }

    /** Right scalar multiplication operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<Real, Dim> operator*(FPoint<Real, Dim> lhs, const T& val) {
        lhs *= val;
        return lhs;
    }

    /** Left scalar multiplication operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<Real, Dim> operator*(const T& val, FPoint<Real, Dim> rhs) {
        rhs *= val;
        return rhs;
    }

    /** Data to data division assignment */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint<Real, Dim>& operator /=(const FPoint<T, Dim>& other) {
        for(std::size_t i = 0; i < Dim; ++i) {
            this->data()[i] *= other.data()[i];
        }
        return *this;
    }

    /** Data to data division operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<Real, Dim> operator/(FPoint<Real, Dim> lhs, const FPoint<Real, Dim>& rhs) {
        lhs /= rhs;
        return lhs;
    }

    /** Right scalar division assignment operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint<Real, Dim>& operator /=(const T& val) {
        auto this_it  = this->begin();
        for(std::size_t i = 0; i < Dim; i++, ++this_it) {
            *this_it /= val;
        }

        return *this;
    }

    /** Right scalar division operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<Real, Dim> operator/(FPoint<Real, Dim> lhs, const T& val) {
        lhs /= val;
        return lhs;
    }

    /** Equality test operator */
    template<class T, must_be_integral<T> = nullptr>
    friend bool operator==(const FPoint<Real, Dim>& lhs, const FPoint<T, Dim>& rhs) {
        auto lhs_it = lhs.begin(), rhs_it  = rhs.begin();
        for(std::size_t i = 0; i < Dim; i++, ++lhs_it, ++rhs_it) {
            if( *lhs_it != *rhs_it) {
                return false;
            }
        }
        return true;
    }

    /** Equality test operator */
    template<class T, must_be_floating<T> = nullptr>
    friend bool operator==(const FPoint<Real, Dim>& lhs, const FPoint<T, Dim>& rhs) {
        auto lhs_it = lhs.begin(), rhs_it  = rhs.begin();
        for(std::size_t i = 0; i < Dim; i++, ++lhs_it, ++rhs_it) {
            if(! Ffeq(*lhs_it, *rhs_it)) {
                return false;
            }
        }
        return true;
    }

    /** Non equality test operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend bool operator!=(const FPoint<Real, Dim>& lhs, const FPoint<T, Dim>& rhs) {
        return ! (lhs == rhs);
    }

    /** Formated output stream operator */
    friend std::ostream& operator<<(std::ostream& os, const FPoint<Real, Dim>& pos) {
        os << "[";
        for(auto it = pos.begin(); it != pos.end()-1; it++)
            os << *it << ", ";
        os << pos.back() << "]";
        return os;
    }

    /** Formated input stream operator */
    friend std::istream& operator>>(std::istream& is, const FPoint<Real, Dim>& pos) {
        for(std::size_t i = 0; i < Dim; ++i) {
            is >> pos.data()[i];
        }
        return is;
    }

    /** \brief Save current object */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const {
        for(std::size_t i = 0; i < Dim; ++i) {
            buffer << this->data()[i];
        }
    }

    /** \brief Retrieve current object */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer) {
        for(std::size_t i = 0; i < Dim; ++i) {
            buffer >> this->data()[i];
        }
    }

};


#endif
