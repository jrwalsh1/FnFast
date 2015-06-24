//------------------------------------------------------------------------------
/// \file Random.hpp
//
// Author(s):
//    Frank Tackmann
//
// Copyright:
//    Copyright (C) 2012 MIT
//
//    This file is part of the Geneva MC framework. Geneva is distributed under
//    the terms of the GNU General Public License version 3 (GPLv3), see the
//    COPYING file that comes with this distribution for details.
//    Please respect the academic usage guidelines in the GUIDELINES file.
//
// Description:
//    Random number generator definitions
//------------------------------------------------------------------------------

#ifndef RANDOM_HPP
#define RANDOM_HPP

/// Random number between 0 and 1
double randomInterval();

/// Default random number generator; called by randomInterval()
double Ran2(long*);

/// Set global random number seed
void setRandomSeed(long newSeed);

#endif // RANDOM_HPP
