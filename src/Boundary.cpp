/*
 * Boundary.cpp
 *
 *  Created on: 18 Oct 2012
 *      Author: robinsonm
 */


#include "Boundary.h"

namespace Tyche {

void FluxBoundary::integrate(const double dt) {

	BOOST_FOREACH(Species *s, get_species()) {
		Particles& mols = s->mols;
		boost::poisson_distribution<> p_dist(dt*rate);
		boost::variate_generator<base_generator_type&, boost::poisson_distribution<> > poisson(generator, p_dist);
		const unsigned int n = poisson();
		for (int i = 0; i < n; ++i) {
			mols.add_molecule(p + uni1()*t1 + uni2()*t2);
		}
	}
}



}

