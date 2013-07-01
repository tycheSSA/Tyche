/*
 * Run_spec.h
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
#include <boost/preprocessor.hpp>
#include <numeric>
#include <boost/timer.hpp>
#include "Log.h"

#define n BOOST_PP_ITERATION()

#define RUN_print(z, n, notused) BOOST_PP_CAT(arg,n) << std::endl <<
#define RUN_exec_ops_on_nsm(z, n, notused) BOOST_PP_CAT(arg,n).init_reactions(nsm);
#define RUN_arguements(z, n, notused) BOOST_PP_CAT(T,n)& BOOST_PP_CAT(arg,n)
#define RUN_(z, n, notused) BOOST_PP_CAT(T,n)& BOOST_PP_CAT(arg,n)
#define RUN_reset(z, n, notused) BOOST_PP_CAT(arg,n).reset();
#define RUN_execute(z, n, notused) BOOST_PP_CAT(arg,n)(dt);
#define RUN_print_times(z, n, notused) BOOST_PP_CAT(arg,n) << "\t" << BOOST_PP_CAT(arg,n).get_time_string() << std::endl <<
#define RUN_template(n) template <BOOST_PP_ENUM_PARAMS(n, typename T)>
#define RUN_empty(n)

namespace Tyche {

BOOST_PP_IF(n, RUN_template , RUN_empty)(n)
void run(Species& s, double for_time, double dt BOOST_PP_COMMA_IF(n) BOOST_PP_ENUM(n, RUN_arguements, none)) {
	const int iterations = int(for_time/dt) + 1;
	const double num_out = 100;
	const int it_per_out = int(iterations/num_out+0.5);

	boost::timer timer;
	double time = 0;
	LOG(1, "Running simulation for "<< for_time << " seconds with dt = " << dt << " seconds.");
	LOG(1, "Will perform "<<iterations<<" iterations for a total simulation time of "<<iterations*dt<<" seconds.");
	LOG(1, "Operators to be run (in order) are:");
	LOG(1, BOOST_PP_REPEAT(n, RUN_print, none) " ");
	//BOOST_PP_REPEAT(n, RUN_start_timers, none)
	for (int i = 0; i < iterations; ++i) {
		const int nn = s.copy_numbers.size();
		for (int ii = 0; ii < nn; ++ii) {
			if (s.copy_numbers[ii] < 0) {
				std::cout << "cell "<<ii<<" has "<<s.copy_numbers[ii]<<" molecules" <<std::endl;
				exit(-1);
			}
		}
		if ((it_per_out==0)||(i%it_per_out==0)) {
			int nmol = 0;
			for (int j = 0; j < arg0.get_species().size(); ++j) {
				nmol += arg0.get_species()[j]->mols.size();
			}
			LOG(1, "Simulation " << double(i)/iterations*100 << "% complete. There are " << nmol <<
					" molecules in free space and " << std::accumulate(s.copy_numbers.begin(),s.copy_numbers.end(),0) <<
					" molecules in " << (unsigned int)s.copy_numbers.size() << " compartments");
		}
		timer.restart();
		for (int j = 0; j < arg0.get_species().size(); ++j) {
			arg0.get_species()[j]->mols.save_indicies();
		}
		BOOST_PP_REPEAT(n, RUN_execute, none)

		time += timer.elapsed();
	}

	LOG(1,"Total time = " << time << " s");
	LOG(1, arg0.get_global_time());
	LOG(1, "Operator times:");
	LOG(1, BOOST_PP_REPEAT(n, RUN_print_times, none) " ");
}

BOOST_PP_IF(n, RUN_template , RUN_empty)(n)
void reset_then_run(Species& s, double for_time, double dt BOOST_PP_COMMA_IF(n) BOOST_PP_ENUM(n, RUN_arguements, none)) {
	LOG(1, "Reseting all operators...");
	BOOST_PP_REPEAT(n, RUN_reset, none)
	LOG(1, "done.");
	run(s, for_time, dt BOOST_PP_COMMA_IF(n) BOOST_PP_ENUM_PARAMS(n, arg));
}
}
#undef n
