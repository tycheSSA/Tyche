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
		const Vect3d step_length = calc_step_length(s, dt);
		const int n = s.mols.size();
		for (int j = 0; j < n; ++j) {
			s.mols.r0[j] = s.mols.r[j];
			s.mols.r[j] += step_length.cwiseProduct(Vect3d(norm(),norm(),norm()));
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


}
