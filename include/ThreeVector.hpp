//------------------------------------------------------------------------------
/// \file ThreeVector.hpp
//
// Author(s):
//    Frank Tackmann
//
// Copyright:
//    Copyright (C) 2012 MIT, Frank Tackmann
//
//    This file is part of the Geneva MC framework. Geneva is distributed under
//    the terms of the GNU General Public License version 3 (GPLv3), see the
//    COPYING file that comes with this distribution for details.
//    Please respect the academic usage guidelines in the GUIDELINES file.
//
// Description:
//    Definition of class ThreeVector
//------------------------------------------------------------------------------

#ifndef THREE_VECTOR_HPP
#define THREE_VECTOR_HPP

#include <iosfwd>
#include <cmath>

using namespace std;

//------------------------------------------------------------------------------
/**
 * \class ThreeVector
 * \ingroup Core
 *
 * \brief A threevector with Euclidean metric (+,+,+) using cartesian coordinates
 *
 * The class implements a momentum threevector \f$ \vec{p} \f$ or position
 * threevector \f$ \vec{x} \f$ in terms of its cartesian coordinates
 * \f$ \vec{p} = (p_x, p_y, p_z) \f$ or \f$ \vec{x} = (x, y, z)) \f$.
 *
 * Most of the member functions and operators are implemented inline to make the
 * class lightweight. The members p1(), p2(), p3() give direct access to
 * the components for both reading and writing.
 *
 * The following algebraic operations are defined with the obvious meanings
 * (where p, q, r are ThreeVectors, and c is a double):
 *
 * - p = q + r, p = q - r (component-wise addition/subtraction)
 * - p = c * q, p = q * c (multiplication by real scalar)
 * - c = q * p            (scalar product of p and q)
 * - p += q, p -= q       (computes p = p +/- q)
 * - p *= c               (computes p = c * p)
 * - p == q, p != q       (component-wise comparison)
 *
 * There are various physics related queries which return the square and
 * magnitude and manipulations like rotations and normalization.
 */
//------------------------------------------------------------------------------
class ThreeVector
{
   protected:
      double _p1;                      ///< 1-component (px or x)
      double _p2;                      ///< 2-component (py or y)
      double _p3;                      ///< 3-component (pz or z)

   public:
      /// \name Constructors etc.
      ///@{

      /// Creates a zero %ThreeVector.
      ThreeVector();

      /// Creates a %ThreeVector (p1, p2, p3).
      ThreeVector(double p1, double p2, double p3);

      /// %ThreeVector copy constructor
      ThreeVector(const ThreeVector& other);

      /// %ThreeVector assignment operator
      ThreeVector& operator=(const ThreeVector& rhs);

      ~ThreeVector() {}
      ///@}

      /// \name Accessing and Querying Components
      ///@{
      double p1() const;               ///< Returns 1-component (px or x).
      double p2() const;               ///< Returns 2-component (py or y).
      double p3() const;               ///< Returns 3-component (pz or z).

      double& p1();                    ///< Access 1-component (px or x)
      double& p2();                    ///< Access 2-component (py or y)
      double& p3();                    ///< Access 3-component (pz or z)

      // Returns the component indexed by mu.
      double operator[](int mu) const;

      // Gives access to the component indexed by mu.
      double& operator[](int mu);

      // Returns all components.
      void getComponents(double& p1, double& p2, double& p3) const;

      // Sets all components.
      void setComponents(double p1, double p2, double p3);

      /// Sets all components to zero.
      void setToZero();
      ///@}

      /// \name Algebra
      ///@{
      // ThreeVector equality comparison operator
      bool operator==(const ThreeVector& rhs) const;

      // ThreeVector nonequality comparison operator
      bool operator!=(const ThreeVector& rhs) const;

      // ThreeVector comparison operator, based on magnitudes
      bool operator<(const ThreeVector& rhs) const;

      // Test if the ThreeVector is nonzero.
      bool isNonzero() const;

      // ThreeVector addition assignment operator
      ThreeVector& operator+=(const ThreeVector& rhs);

      // ThreeVector subtraction assignment operator
      ThreeVector& operator-=(const ThreeVector& rhs);

      // ThreeVector multiplication with real number assignment operator
      ThreeVector& operator*=(double rhs);

      // ThreeVector division by real number assignment operator
      ThreeVector& operator/=(double rhs);
      ///@}

      /// \name Physics
      /// Physics related queries and manipulations
      ///@{
      // Returns the invariant square of the ThreeVector.
      double square() const;

      // Returns the magnitude of the ThreeVector.
      double magnitude() const;

      // Returns the %ThreeVector's transverse magnitude relative to the z-direction.
      double perp() const;

      // Returns the ThreeVector's transverse magnitude squared relative to the z-direction.
      double perp2() const;

      // Returns the ThreeVector's magnitude perpendicular to a given direction.
      double perp(const ThreeVector& n) const;

      // Returns the %ThreeVector's magnitude squared perpendicular to a given direction.
      double perp2(const ThreeVector& n) const;

      // Applies a transformation to the ThreeVector.
      template <class Transformation>
      ThreeVector& apply(const Transformation& trans);
      ///@}

      /// \name Output
      ///@{
      // ThreeVector stream insertion operator
      friend ostream& operator<<(ostream& out, const ThreeVector& p);
      ///@}

   protected:
      /// Streams a string in the format (p1, p2, p3) to \a out.
      ostream& print(ostream& out) const;
};

// ThreeVector addition operator
const ThreeVector operator+(const ThreeVector& lhs, const ThreeVector& rhs);

// ThreeVector subtraction operator
const ThreeVector operator-(const ThreeVector& lhs, const ThreeVector& rhs);

// ThreeVector negation operator
const ThreeVector operator-(const ThreeVector& lhs);

// ThreeVector multiplication operator with real number from the left
const ThreeVector operator*(double lhs, const ThreeVector& rhs);

// ThreeVector multiplication operator with real number from the right
const ThreeVector operator*(const ThreeVector& lhs, double rhs);

// ThreeVector division operator with real number
const ThreeVector operator/(const ThreeVector& lhs, double rhs);

// ThreeVector multiplication operator corresponding to inner product
double operator*(const ThreeVector& lhs, const ThreeVector& rhs);

// ThreeVector cross product
const ThreeVector crossProduct(const ThreeVector& lhs, const ThreeVector& rhs);

// Computes the cosine of the angle between two ThreeVectors.
double cosAngleBetween(const ThreeVector& lhs, const ThreeVector& rhs);


////////////////////////////////////////////////////////////////////////////////
// Inline Definitions
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
/**
 * \ingroup Core
 *
 * \brief Approximate comparison of two doubles to given precision and accuracy
 * \param x1 The 1st double to compare.
 * \param x2 The 2nd double to compare.
 * \param delta The relative precision \f$ \delta \f$ to be used. The internal
 *        default precision is used as default value.
 * \param epsilon The absolute accuracy \f$ \epsilon \f$ to be used. The
 *        internal default accuracy is used as default value.
 * \return An integer 1, 0, or -1 corresponding to the sign of \f$ x_1 - x_2 \f$
 *
 * Checks if \a x1 and \a x2 are approximately equal in terms of both absolute
 * accuracy and relative precision (for the latter using the algorithm by
 * D.E. Knuth in Seminumerical and used in the floating-point comparison from
 * the GSL). That is, \a x1 and \a x2 are considered equal (return value 0) if
 * their absolute difference is less or equal either than \a epsilon or roughly
 * \a delta times the larger of them.
 *
 * Specifically, the return value is computed as follows:
 * \f[
 * \begin{aligned}
 * -1 &\qquad\text{for}\quad& x_1 - x_2 &< -\varepsilon \\
 * 0  &\qquad\text{for}\quad& -\varepsilon \leq x_1 - x_2 &\leq \varepsilon \\
 * 1  &\qquad\text{for}\quad& \varepsilon < x_1 - x_2 &
 * \end{aligned}
 * \f]
 *
 * where \f$ \varepsilon \f$ is given as
 * \f$ \varepsilon = \mathrm{max} \{ \epsilon, \delta \times 2^n \}\f$
 * where \f$ n \f$ is the maximum base-2 exponent of \a x1 and \a x2.
 *
 * \note Since \a x1 and \a x2 are also compared to absolute accuracy, the value
 * of \a epsilon should be roughly chosen as \a delta times the typical size of
 * \a x1 and \a x2. To compare using relative precision only (i.e. to arbitrary
 * accuracy) use compareRelative().
 */
inline int compare(const double x1, const double x2,
                   const double delta = 1e-10,
                   const double epsilon = 1e-10)
{
   // compute sign and absolute value of x1 - x2
   int sign;
   double difference = x1 - x2;
   if (difference >= 0.)
      sign = 1;
   else {
      sign = -1;
      difference = -difference;
   }
   // if |x1 - x2| <= epsilon (absolute precision) we are done
   if (difference <= epsilon)
      return 0;

   // check the relative precision
   int exponent;
   frexp(abs(x1) > abs(x2) ? x1 : x2, &exponent);
   const double eps = ldexp(delta, exponent);
   if (difference <= eps)
      return 0;
   return sign;
}

//------------------------------------------------------------------------------
/**
 * \ingroup Core
 *
 * \brief Approximate comparison of a double with zero to a given accuracy
 * \param x The double to compare with zero.
 * \param epsilon The absolute accuracy \f$ \epsilon \f$ to be used. If
 *        \a epsilon is omitted, the internal default accuracy is used.
 * \return An integer 1, 0, or -1 corresponding to the sign of \a x
 *
 * Checks if \a x is approximately equal to zero.
 *
 * The return value is given by
 * \f[
 * \begin{aligned}
 * -1 &\qquad\text{for}\quad& x &< -\epsilon \\
 * 0  &\qquad\text{for}\quad& -\epsilon \leq x &\leq \epsilon \\
 * 1  &\qquad\text{for}\quad& \epsilon < x &
 * \end{aligned}
 * \f]
 *
 * \warning Comparing to an absolute accuracy is often not very meaningful
 * without further information where that zero came from. Hence, direct use of
 * this function should be avoided in favor of compare(). For example, instead of
 * \code compareToZero(x1 - x2 - a*x3) \endcode
 * one should rather use one of
 * \code
 * compare(x1 - x2, a*x3)
 * compare(x1, x2 + a*x3)
 * \endcode
 */
inline int compareToZero(const double x, const double epsilon = 1e-10)
{
   if (x > epsilon)
      return 1;
   if (x < -epsilon)
      return -1;
   return 0;
}

//------------------------------------------------------------------------------
/**
 * \ingroup Core
 *
 * \brief Approximate comparison of two doubles to given relative precision only
 * \param x1 The 1st double to compare.
 * \param x2 The 2nd double to compare.
 * \param delta The relative precision \f$ \delta \f$ to be used. If \a delta
 *        is omitted, the internal default precision is used.
 * \return An integer 1, 0, or -1 corresponding to the sign of \f$ x_1 - x_2 \f$
 *
 * Checks if \a x1 and \a x2 are approximately equal using relative precision
 * only as given by \a delta. This is equivalent to using compare() with
 * \a epsilon set to 0.0 (i.e. using arbitrary accuracy).
 *
 * \note Since \a x1 and \a x2 are only compared to relative precision, this
 * function should not be used to compare a double to 0.0.
 */
inline int compareRelative(const double x1, const double x2,
                           const double delta = 1e-10)
{
   int exponent;
   frexp(abs(x1) > abs(x2) ? x1 : x2, &exponent);
   const double eps = ldexp(delta, exponent);
   return compareToZero(x1 - x2, eps);
}

//------------------------------------------------------------------------------
// Constructors etc.
//------------------------------------------------------------------------------
// Creates a zero ThreeVector.
inline ThreeVector::ThreeVector()
   : _p1(0.), _p2(0.), _p3(0.)
{}

//------------------------------------------------------------------------------
// Creates a ThreeVector (p1, p2, p3).
inline ThreeVector::ThreeVector(double p1, double p2, double p3)
   : _p1(p1), _p2(p2), _p3(p3)
{}

//------------------------------------------------------------------------------
// ThreeVector copy constructor
inline ThreeVector::ThreeVector(const ThreeVector& other)
   : _p1(other._p1), _p2(other._p2), _p3(other._p3)
{}

//------------------------------------------------------------------------------
// ThreeVector assignment operator
inline ThreeVector& ThreeVector::operator=(const ThreeVector& rhs)
{
   _p1 = rhs._p1;
   _p2 = rhs._p2;
   _p3 = rhs._p3;
   return *this;
}

//------------------------------------------------------------------------------
// Accessing and Querying Components
//------------------------------------------------------------------------------
inline double ThreeVector::p1() const { return _p1; }
inline double ThreeVector::p2() const { return _p2; }
inline double ThreeVector::p3() const { return _p3; }

inline double& ThreeVector::p1()      { return _p1; }
inline double& ThreeVector::p2()      { return _p2; }
inline double& ThreeVector::p3()      { return _p3; }

//------------------------------------------------------------------------------
/**
 * \brief Returns all components.
 * \param[out] p1 px or x component
 * \param[out] p2 py or y component
 * \param[out] p3 pz or z component
 *
 * The components of the %ThreeVector are returned in the provided parameters.
 */
inline void ThreeVector::getComponents(double& p1, double& p2, double& p3) const
{
   p1 = _p1;
   p2 = _p2;
   p3 = _p3;
}

//------------------------------------------------------------------------------
/**
 * \brief Sets all components.
 * \param p1 px or x component
 * \param p2 py or y component
 * \param p3 pz or z component
 *
 * Sets all components of the %ThreeVector to the given values.
 */
inline void ThreeVector::setComponents(double p1, double p2, double p3)
{
   _p1 = p1;
   _p2 = p2;
   _p3 = p3;
}

//------------------------------------------------------------------------------
// Sets all components to zero.
inline void ThreeVector::setToZero()
{
   setComponents(0., 0., 0.);
}

//------------------------------------------------------------------------------
// Algebra
//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector equality comparison operator
 * \param rhs The %ThreeVector to compare to.
 * \return \c true if all components are equal, \c false otherwise.
 *
 * The comparison is done with the default accuracy and precision separately for
 * each of the components.
 */
inline bool ThreeVector::operator==(const ThreeVector& rhs) const
{
   return !operator!=(rhs);
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector inequality comparison operator
 * \param rhs The %ThreeVector to compare to.
 * \return \c true if at least one of the components are different,
 *         \c false otherwise.
 *
 * The comparison is done with the default accuracy and precision separately for
 * each of the components.
 */
inline bool ThreeVector::operator!=(const ThreeVector& rhs) const
{
   return compare(_p1, rhs._p1) || compare(_p2, rhs._p2) || compare(_p3, rhs._p3);
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector less than comparison operator, to generate ordering
 * \param rhs The %ThreeVector to compare to.
 * \return \c true if magnitude is smaller, or 
 *              checks if p1, p2, then p3 are smaller recursively.
 *         \c false otherwise.
 *
 * The comparison is done with the default accuracy and precision separately for
 * each of the components.
 */
inline bool ThreeVector::operator<(const ThreeVector& rhs) const
{
   int compsq = compare(square(), rhs.square());
   if (compsq != 0) {
      return (compsq == -1);
   } else { // if the squares are equal, start comparing p1
      int comp1 = compare(_p1, rhs._p1);
      if (comp1 != 0) {
         return (comp1 == -1);
      } else { // if p1 values equal, move on to p2
         int comp2 = compare(_p2, rhs._p2);
         if (comp2 != 0) {
            return (comp2 == -1);
         } else { // if p2 values equal, move on to p3
            int comp3 = compare(_p3, rhs._p3);
            if (comp3 != 0) {
               return (comp3 == -1);
            } else { // if p3 values equal, vectors are equal and just return false
               return false;
            }
         }
      }
   }
}

//------------------------------------------------------------------------------
/**
 * \brief Test if the %ThreeVector is nonzero.
 * \return \c true if at least one component is nonzero, \c false otherwise.
 *
 * Compares each components to 0.0 using the default absolute accuracy.
 *
 * \warning This function is here for completeness. Directly comparing to zero
 * is usually not very meaningfull for doubles since the overall scale is
 * unknown, so use of this function should be avoided. For example, one should
 * compare two %ThreeVectors with each other rather than their difference to
 * zero.
 */
inline bool ThreeVector::isNonzero() const
{
   return compareToZero(_p1) || compareToZero(_p2) || compareToZero(_p3);
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector addition assignment operator
 * \param rhs The %ThreeVector to be added.
 * \return A reference to the %ThreeVector.
 *
 * \a rhs is added to the %ThreeVector.
 */
inline ThreeVector& ThreeVector::operator+=(const ThreeVector& rhs)
{
   _p1 += rhs._p1;
   _p2 += rhs._p2;
   _p3 += rhs._p3;
   return *this;
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector subtraction assignment operator
 * \param rhs The %ThreeVector to be subtracted.
 * \return A reference to the %ThreeVector.
 *
 * \a rhs is subtracted from the %ThreeVector.
 */
inline ThreeVector& ThreeVector::operator-=(const ThreeVector& rhs)
{
   _p1 -= rhs._p1;
   _p2 -= rhs._p2;
   _p3 -= rhs._p3;
   return *this;
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector multiplication with real number assignment operator
 * \param rhs The number to be multiplied by.
 * \return A reference to the %ThreeVector.
 *
 * Each component of the %ThreeVector is multiplied by \a rhs.
 */
inline ThreeVector& ThreeVector::operator*=(double rhs)
{
   _p1 *= rhs;
   _p2 *= rhs;
   _p3 *= rhs;
   return *this;
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector division by real number assignment operator
 * \param rhs The number to be divided by.
 * \return A reference to the %ThreeVector.
 *
 * Each component of the %ThreeVector is divided by \a rhs.
 * \note Division by zero is not checked.
 */
inline ThreeVector& ThreeVector::operator/=(double rhs)
{
   return operator*=(1. / rhs);
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector addition operator
 * \param lhs The 1st %ThreeVector to add.
 * \param rhs The 2nd %ThreeVector to add.
 * \return A %ThreeVector given by the sum of \a lhs and \a rhs.
 */
inline const ThreeVector operator+(const ThreeVector& lhs, const ThreeVector& rhs)
{
   return ThreeVector(lhs.p1() + rhs.p1(), lhs.p2() + rhs.p2(), lhs.p3() + rhs.p3());
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector subtraction operator
 * \param lhs The 1st %ThreeVector
 * \param rhs The 2nd %ThreeVector to be subtracted from \a lhs.
 * \return A %ThreeVector given by the difference of \a lhs and \a rhs.
 */
inline const ThreeVector operator-(const ThreeVector& lhs, const ThreeVector& rhs)
{
   return ThreeVector(lhs.p1() - rhs.p1(), lhs.p2() - rhs.p2(), lhs.p3() - rhs.p3());
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector negation operator
 * \param lhs The %ThreeVector to be negated.
 * \return A %ThreeVector given by the negative of \a lhs.
 */
inline const ThreeVector operator-(const ThreeVector& lhs)
{
   return ThreeVector(-lhs.p1(), -lhs.p2(), -lhs.p3());
}

//------------------------------------------------------------------------------
/**
 * \brief  %ThreeVector multiplication operator with real number from the left
 * \param lhs The number to multiply with.
 * \param rhs The %ThreeVector to be multiplied by \a lhs.
 * \return A %ThreeVector given by the product of \a lhs times \a rhs.
 */
inline const ThreeVector operator*(double lhs, const ThreeVector& rhs)
{
   return ThreeVector(lhs * rhs.p1(), lhs * rhs.p2(), lhs * rhs.p3());
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector multiplication operator with real number from the right
 * \param lhs The %ThreeVector to be multiplied by \a rhs.
 * \param rhs The number to multiply with.
 * \return A %ThreeVector given by the product of \a lhs times \a rhs.
 */
inline const ThreeVector operator*(const ThreeVector& lhs, double rhs)
{
   return rhs * lhs;
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector division operator with real
 * \param lhs The %ThreeVector to be divided by \a rhs.
 * \param rhs The number to divide by.
 * \return A %ThreeVector given by \a lhs divided by \a rhs.
 * \note Division by zero is not checked.
 */
inline const ThreeVector operator/(const ThreeVector& lhs, double rhs)
{
   return (1./rhs) * lhs;
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector multiplication operator corresponding to inner product
 * \param lhs The 1st %ThreeVector.
 * \param rhs The 2nd %ThreeVector.
 * \return The result from taking the inner product of \a lhs and \a rhs.
 *
 * The inner product is computed as
 * \f$ \vec{p} \cdot \vec{q} = p_x q_x + p_y q_y + p_z q_z \f$.
 */
inline double operator*(const ThreeVector& lhs, const ThreeVector& rhs)
{
   return lhs.p1()*rhs.p1() + lhs.p2()*rhs.p2() + lhs.p3()*rhs.p3();
}

//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector cross product.
 * \param lhs The 1st %ThreeVector.
 * \param rhs The 2nd %ThreeVector.
 * \return The result from taking the cross product of \a lhs and \a rhs.
 */
inline const ThreeVector crossProduct(const ThreeVector& lhs, const ThreeVector& rhs)
{
   ThreeVector tmp(lhs.p2() * rhs.p3() - lhs.p3() * rhs.p2(),
                   lhs.p3() * rhs.p1() - lhs.p1() * rhs.p3(),
                   lhs.p1() * rhs.p2() - lhs.p2() * rhs.p1());
   return tmp;
}

//------------------------------------------------------------------------------
// Physics
//------------------------------------------------------------------------------
/**
 * \brief Returns the square of the %ThreeVector.
 * \return The value of \f$ \vec{p}^2 \f$.
 */
inline double ThreeVector::square() const
{
   return _p1*_p1 + _p2*_p2 + _p3*_p3;
}

//------------------------------------------------------------------------------
/**
 * \brief Returns the magnitude of the %ThreeVector.
 * \return The value of \f$ \lvert \vec{p} \rvert = \sqrt{\vec{p}^2} \f$.
 */
inline double ThreeVector::magnitude() const
{
   return sqrt(square());
}

//------------------------------------------------------------------------------
/**
 * \brief Returns the %ThreeVector's transverse magnitude relative to the z-direction.
 * \return The value of \f$ \lvert \vec{p}_T \rvert = \sqrt{p_x^2 + p_y^2} \f$.
 */
inline double ThreeVector::perp() const
{
   return sqrt(perp2());
}

//------------------------------------------------------------------------------
/**
 * \brief Returns the %ThreeVector's transverse magnitude squared relative to the z-direction.
 * \return The value of \f$ \lvert \vec{p}_T \rvert^2 = p_x^2 + p_y^2 \f$.
 */
inline double ThreeVector::perp2() const
{
   return _p1*_p1 + _p2*_p2;
}

//------------------------------------------------------------------------------
/**
 * \brief Returns the %ThreeVector's magnitude perpendicular to a given direction.
 * \param n The unit %ThreeVector \f$ \vec{n} \f$ to use as reference direction.
 * \return The value of \f$ \lvert \vec{p}_{\perp n} \rvert = \sqrt{\vec{p}^2 - (\vec{n}\cdot\vec{p})^2} \f$.
 *
 * \note The reference vector \a n is assumed to be normalized
 * \f$ \vec{n}^2 = 1 \f$.
 */
inline double ThreeVector::perp(const ThreeVector& n) const
{
   return sqrt(perp2(n));
}

//------------------------------------------------------------------------------
/**
 * \brief Returns the %ThreeVector's magnitude squared perpendicular to a given direction.
 * \param n The unit %ThreeVector \f$ \vec{n} \f$ to use as reference direction.
 * \return The value of \f$ \lvert \vec{p}_{\perp n} \rvert^2 = p^2 - (\vec{n}\cdot \vec{p})^2 \f$.
 *
 * \note The reference vector \a n is assumed to be normalized
 * \f$ \vec{n}^2 = 1 \f$.
 */
inline double ThreeVector::perp2(const ThreeVector& n) const
{
   return square() - pow(n * (*this), 2);
}

//------------------------------------------------------------------------------
/**
 * \brief Applies a transformation to the %ThreeVector.
 * \param trans The transformation to be applied.
 * \return A reference to the %ThreeVector.
 *
 * \a trans can be any function object which implements the function evaluation
 * operator with the signature
 * \code
 * const ThreeVector operator()(const ThreeVector& operand)
 * \endcode
 *
 * The result of the transformation is assigned to the %ThreeVector.
 */
template <class Transformation>
ThreeVector& ThreeVector::apply(const Transformation& trans)
{
   return *this = trans(*this);
}

//------------------------------------------------------------------------------
/**
 * \brief Computes the cosine of the angle between two %ThreeVectors.
 * \param lhs The 1st %ThreeVector.
 * \param rhs The 2nd %ThreeVector.
 * \return The value of \f$ \cos\theta \f$ where \f$ \theta \f$ is the angle
 *         between \a lhs and \a rhs.
 *
 * The angle is computed as
 * \f$ \cos\theta = \vec{p}_1\cdot\vec{p}_2/(\lvert\vec{p}_1\rvert \lvert\vec{p}_2\rvert)\f$.
 * The result is guaranteed to be in the range [-1,1].
 */
inline double cosAngleBetween(const ThreeVector& lhs, const ThreeVector& rhs)
{
   double result = (lhs * rhs) / sqrt(lhs.square() * rhs.square());
   if (result > 1.)
      return 1.;
   if (result < -1.)
      return -1.;
   return result;
}

//------------------------------------------------------------------------------
// Output
//------------------------------------------------------------------------------
/**
 * \brief %ThreeVector stream insertion operator
 * \param out The output stream to be used.
 * \param p The %ThreeVector to be streamed.
 * \return A reference to \a out.
 *
 * Streams a string in the format "(p1, p2, p3)" to \a out.
 */
inline ostream& operator<<(ostream& out, const ThreeVector& p)
{
   return p.print(out);
}

#endif // THREE_VECTOR_HPP
