/*
 * Control.h
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
 *  Created on: 4 Dec 2012
 *      Author: robinsonm
 */

#ifndef CONTROL_H_
#define CONTROL_H_

#include "Operator.h"
#include "NextSubvolumeMethod.h"
#include "Output.h"

namespace Tyche {

template<typename T>
class Control: public Operator {
public:
	Control(T& geometry):geometry(geometry) {
	}

protected:
	void integrate(const double dt) {
	}


	T& geometry;
};

template<typename T>
class GrowingInterface: public Control<T> {
public:
	struct my_accumulate {
		my_accumulate(std::vector<int>& vect_to_sum):vect_to_sum(vect_to_sum) {}
		int operator()(int x, int y) {
			return x+vect_to_sum[y];
		}
		std::vector<int>& vect_to_sum;
	};

	GrowingInterface(T& geometry, NextSubvolumeMethod& nsm,
					   const double move_by, const T& shrink_to, const T& grow_to,
					   const double check_dt,
					   const int grow_threshold, const int shrink_threshold):
		Control<T>(geometry),nsm(nsm),
		check_dt(check_dt),move_by(move_by),
		shrink_to(shrink_to),grow_to(grow_to),
		grow_threshold(grow_threshold),shrink_threshold(shrink_threshold) {
		last_check = 0;
	}
protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tGrowing/Shrinking Interface at "<< this->geometry;
	}


	NextSubvolumeMethod& nsm;
	const double check_dt;
	const double move_by;
	const T& shrink_to,grow_to;
	const int grow_threshold,shrink_threshold;
	double last_check;
};


}

#include "Control.impl.h"

#endif /* CONTROL_H_ */
