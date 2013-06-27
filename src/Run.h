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

#include <boost/preprocessor/iteration/iterate.hpp>
#define BOOST_PP_ITERATION_LIMITS (1, 20)
#define BOOST_PP_FILENAME_1       "Run_spec.h"
#include BOOST_PP_ITERATE()

namespace Tyche {

void run(double for_time, double dt, OperatorList& operators) {
	const int iterations = int(for_time/dt) + 1;
	const double num_out = 100;
	const int it_per_out = int(iterations/num_out+0.5);

	boost::timer timer;
	double time = 0;
	LOG(1, "Running simulation for "<< for_time << " seconds with dt = " << dt << " seconds.");
	LOG(1, "Will perform "<<iterations<<" iterations for a total simulation time of "<<iterations*dt<<" seconds.");
	LOG(1, "Operators to be run (in order) are:"<<std::endl<<operators);
	for (int i = 0; i < iterations; ++i) {
		if ((it_per_out==0)||(i%it_per_out==0)) {
			int nmol = 0;
			int ncompart = 0;
			int ncomparts = 0;
			for (int j = 0; j < operators.get_species().size(); ++j) {
				nmol += operators.get_species()[j]->mols.size();
				ncompart += std::accumulate(operators.get_species()[j].copy_numbers.begin(),operators.get_species()[j].copy_numbers.end(),0);
				ncomparts += (unsigned int)operators.get_species()[j].copy_numbers.size();
			}
			LOG(1, "Simulation " << double(i)/iterations*100 << "% complete. There are " << nmol <<
					" molecules in free space and " << ncompart <<
					" molecules in " << ncomparts << " compartments");
		}
		timer.restart();
		for (int j = 0; j < operators.get_species().size(); ++j) {
			operators.get_species()[j]->mols.save_indicies();
		}

		operators(dt);

		time += timer.elapsed();
	}

	LOG(1,"Total time = " << time << " s");
	LOG(1, operators.get_global_time());
	LOG(1, "Operator times:");
	LOG(1, BOOST_PP_REPEAT(n, RUN_print_times, none) " ");
}


void reset_then_run(double for_time, double dt, OperatorList& operators) {
	LOG(1, "Reseting all operators...");
	operators.reset();
	LOG(1, "done.");
	run(for_time, dt, operators);
}
}


#endif /* RUN_H_ */
