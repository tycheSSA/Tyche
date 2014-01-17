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

void Diffusion::integrate(const double dt) {

	const int n = get_species().size();
	for (int i = 0; i < n; ++i) {
		Species &s = *(get_species()[i]);
		const double step_length = calc_step_length(s, dt);
		const int n = s.mols.size();
		for (int j = 0; j < n; ++j) {
			Vect3d& r = s.mols.get_position_nonconst(i);
			Vect3d& r0 = s.mols.get_old_position_nonconst(i);
			r0 = r;
			r += step_length * Vect3d(norm(),norm(),norm());
		}
	}

}



Diffusion create_diffusion() {
	return Diffusion();
}

Diffusion create_diffusion(Species &s) {
	Diffusion to_return; to_return.add_species(s);
	return to_return;
}

void Diffusion1D::integrate(const double dt) {

	const int n = get_species().size();
	for (int i = 0; i < n; ++i) {
		Species &s = *(get_species()[i]);
		const double step_length = calc_step_length(s, dt);
		const int nm = s.mols.size();
		for (int j = 0; j < nm; ++j) {
			s.mols.get_position_nonconst(j)[dim] += step_length * norm();
		}
	}

}

}
