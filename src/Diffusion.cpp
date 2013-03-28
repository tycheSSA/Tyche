/*
 * BrownianDynamics.cpp
 *
 *  Created on: 11 Oct 2012
 *      Author: robinsonm
 */

#include "Diffusion.h"
#include "Log.h"
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <math.h>

namespace Tyche {

void Diffusion::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	const int n = all_species.size();
	for (int i = 0; i < n; ++i) {
		Species &s = *(all_species[i]);
		const double step_length = calc_step_length(s, dt);
		BOOST_FOREACH(Vect3d &r,s.mols.r) {
			r += step_length * Vect3d(norm(),norm(),norm());
		}
	}
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}


std::ostream& operator<< (std::ostream& out, Diffusion &b) {
	return out << "\tDiffusion";
}

void Diffusion::add_species(Species& s) {
	Operator::add_species(s);
}

void Diffusion1D::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	const int n = all_species.size();
	for (int i = 0; i < n; ++i) {
		Species &s = *(all_species[i]);
		const double step_length = calc_step_length(s, dt);
		BOOST_FOREACH(Vect3d &r,s.mols.r) {
			r[dim] += step_length * norm();
		}
	}
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

std::ostream& operator <<(std::ostream& out, Diffusion1D& b) {
	return out << "\t1D Diffusion along axis "<<b.dim;
}

}
