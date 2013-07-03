/*
 * Run.h
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of RD_3D.
 *
 * RD_3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RD_3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with RD_3D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 30 Oct 2012
 *      Author: robinsonm
 */

#ifndef RUN_H_
#define RUN_H_

#include "Operator.h"
#include <initializer_list>

//#include <boost/preprocessor/iteration/iterate.hpp>
//#define BOOST_PP_ITERATION_LIMITS (1, 20)
//#define BOOST_PP_FILENAME_1       "Run_spec.h"
//#include BOOST_PP_ITERATE()

namespace Tyche {
void run(double for_time, double dt, Operator& operators);
void run(double for_time, double dt, const std::initializer_list<Operator*>& operators);

void reset_then_run(double for_time, double dt, Operator& operators);
void reset_then_run(double for_time, double dt, const std::initializer_list<Operator*>& operators);
}


#endif /* RUN_H_ */
