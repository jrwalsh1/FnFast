//------------------------------------------------------------------------------
/// \file ThreeVector.cpp
//
// Author(s):
//    Jon Walsh
//
// Copyright:
//    Copyright (C) 2015  LBL
//
//    This file is part of the FnFast library. FnFast is distributed under the
//    terms of the GNU General Public License version 3 (GPLv3), see the COPYING
//    file that comes with this distribution for details.
//    Please respect the academic usage guidelines in the GUIDELINES file.
//
// Description:
//    Definition of class ThreeVector
//    Modified version of the ThreeVector class from the Geneva MC framework.
//------------------------------------------------------------------------------

#include "ThreeVector.hpp"

#include <iostream>
#include <ostream>

namespace fnfast {

//------------------------------------------------------------------------------
// Accessing and Querying Components
//------------------------------------------------------------------------------
/**
 * \brief Returns the component indexed by mu.
 * \param mu The index of the component 1 <= mu <= 3.
 * \return The value of the component.
 *
 * If mu is out of range an out_of_range error is thrown.
 */
double ThreeVector::operator[](int mu) const
{
   if (mu == 1)
      return _p1;
   if (mu == 2)
      return _p2;
   if (mu != 3)
      std::cout << "ThreeVector::operator[] : mu index " << mu << " is out of range" << std::endl;
   return _p3;
}

//------------------------------------------------------------------------------
/**
 * \brief Gives access to the component indexed by mu.
 * \param mu The index of the component 1 <= mu <= 3.
 * \return A reference to the component.
 *
 * If mu is out of range an out_of_range error is thrown.
 */
double& ThreeVector::operator[](int mu)
{
   if (mu == 1)
      return _p1;
   if (mu == 2)
      return _p2;
   if (mu != 3)
      std::cout << "ThreeVector::operator[] : mu index " << mu << " is out of range" << std::endl;
   return _p3;
}

//------------------------------------------------------------------------------
// Private member functions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Streams a string in the format (p1, p2, p3) to out.
std::ostream& ThreeVector::print(std::ostream& out) const
{
   return out << "(" << _p1 << ", " << _p2 << ", " << _p3 << ")";
}

} // namespace fnfast
