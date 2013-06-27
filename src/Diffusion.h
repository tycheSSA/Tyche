/*
 * BrownianDynamics.h
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
 *  Created on: 11 Oct 2012
 *      Author: robinsonm
 */

#ifndef DIFFUSION_H_
#define DIFFUSION_H_


#include "Species.h"
#include "Operator.h"
#include <vector>
#include "MyRandom.h"

namespace Tyche {


class Diffusion: public Operator {
public:
	Diffusion():norm(generator,boost::normal_distribution<>(0,1)) {}
	void operator()(const double dt);
	void add_species(Species &s);
	virtual void print(std::ostream& out) {
		out << "\tDiffusion";
	}
protected:
	double calc_step_length(Species &s, const double dt) {
		return sqrt(2.0*s.D*dt);
	}
	boost::variate_generator<base_generator_type&, boost::normal_distribution<> > norm;
	std::vector<int> compartment_rate_indicies;
};


Diffusion create_diffusion(Species &s);
Diffusion create_diffusion();

class Diffusion1D: public Diffusion {
public:
	Diffusion1D(const int dim):Diffusion(),dim(dim) {}
	void operator()(const double dt);
	const int dim;
	virtual void print(std::ostream& out) {
		out << "\t1D Diffusion along axis "<<dim;
	}
};

}

#endif /* DIFFUSION_H_ */
