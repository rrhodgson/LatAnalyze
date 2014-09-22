/*
 * Minimizer.hpp, part of LatAnalyze 3
 *
 * Copyright (C) 2013 - 2014 Antonin Portelli
 *
 * LatAnalyze 3 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LatAnalyze 3 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LatAnalyze 3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef Latan_Minimizer_hpp_
#define Latan_Minimizer_hpp_

#include <LatAnalyze/Global.hpp>
#include <LatAnalyze/Function.hpp>
#include <LatAnalyze/Mat.hpp>
#include <LatAnalyze/Solver.hpp>

BEGIN_NAMESPACE

/******************************************************************************
 *                        Abstract minimizer class                            *
 ******************************************************************************/

class Minimizer: public Solver
{
public:
    // constructor
    Minimizer(void) = default;
    // destructor
    virtual ~Minimizer(void) = default;
    // minimization
    virtual const DVec & operator()(const DoubleFunction &f) = 0;
};

END_NAMESPACE

#endif // Latan_Minimizer_hpp_
