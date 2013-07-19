/*
 * run.cpp
 * 
 * Copyright 2013 Martin Robinson
 *
 * This file is part of Tyche.
 *
 * Tyche is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Tyche is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Tyche.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 1 Jul 2013
 *      Author: robinsonm
 */

#include "Run.h"
namespace Tyche {

void run(double for_time, double dt, const std::initializer_list<Operator*>& operators) {
	OperatorList list = operators;
	run(for_time, dt, list);
}

void run(double for_time, double dt, Operator& operators) {
	const int iterations = int(for_time/dt) + 1;
	const double num_out = 100;
	const int it_per_out = int(iterations/num_out+0.5);

	boost::timer timer;
	double time = 0;
	LOG(1, "Running simulation for "<< for_time << " seconds with dt = " << dt << " seconds.");
	LOG(1, "Will perform "<<iterations<<" iterations for a total simulation time of "<<iterations*dt<<" seconds.");
	LOG(1, "Operators to be run (in order) are:"<<std::endl<<operators);
	LOG(1, "Number of species:"<<operators.get_species().size());
	for (int i = 0; i < iterations; ++i) {
		if ((it_per_out==0)||(i%it_per_out==0)) {
			int nmol = 0;
			int ncompart = 0;
			int ncomparts = 0;
			for (int j = 0; j < operators.get_species().size(); ++j) {
				nmol += operators.get_species()[j]->mols.size();
				ncompart += std::accumulate(operators.get_species()[j]->copy_numbers.begin(),operators.get_species()[j]->copy_numbers.end(),0);
			}
			ncomparts += (unsigned int)operators.get_species()[0]->copy_numbers.size();

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
	LOG(1, operators << "\t" << operators.get_time_string() << std::endl);

}

void reset_then_run(double for_time, double dt, const std::initializer_list<Operator*>& operators) {
	OperatorList list = operators;
	reset_then_run(for_time, dt, list);
}

void reset_then_run(double for_time, double dt, Operator& operators) {
	LOG(1, "Reseting all operators...");
	operators.reset();
	LOG(1, "done.");
	run(for_time, dt, operators);
}
}
