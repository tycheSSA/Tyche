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
#include <set>
#include "MyRandom.h"
#include "NextSubvolumeMethod.h"

namespace Tyche {


class Diffusion: public Operator {
public:
	Diffusion():norm(generator,boost::normal_distribution<>(0,1)) {}
	static std::auto_ptr<Operator> New() {
		return std::auto_ptr<Operator>(new Diffusion());
	}


protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tDiffusion";
	}
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
	const int dim;

protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) {
		out << "\t1D Diffusion along axis "<<dim;
	}
};

template<typename T>
class DiffusionWithTracking : public Diffusion {
public:
  DiffusionWithTracking(const T& geometry, NextSubvolumeMethod *nsm) : nsm(nsm) {
    std::vector<int> slice_indices;
    std::set<int> indices;
    const StructuredGrid &subvolumes = nsm->get_grid();
    subvolumes.get_slice(geometry, slice_indices);
    const int n = slice_indices.size();
    for (int i = 0; i < n; ++i) {
      const std::vector<int>& neighbrs = subvolumes.get_neighbour_indicies(slice_indices[i]);
      const int nn = neighbrs.size();
      for (int j = 0; j < nn; ++j)
	if (!subvolumes.is_in(geometry, neighbrs[j]))
	  indices.insert(neighbrs[j]);
    }
    indices_to_consider = std::vector<int>(indices.begin(), indices.end());
  };
  static std::auto_ptr<Operator> New(const T& geometry, NextSubvolumeMethod *_nsm) {
    return std::auto_ptr<Operator>(new DiffusionWithTracking(geometry,_nsm));
  }

protected:
  virtual void integrate(const double dt);

private:
  NextSubvolumeMethod *nsm;

  std::vector<int> indices_to_consider;
};

#include "Diffusion.impl.h"

}

#endif /* DIFFUSION_H_ */
