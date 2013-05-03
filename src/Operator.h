/*
 * Operator.h
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

#ifndef OPERATOR_H_
#define OPERATOR_H_

#include "Species.h"
#include <boost/timer.hpp>
#include <initializer_list>

namespace Tyche {

class Operator {
public:
	Operator();
	void add_species(Species &s);
	void add_species(std::initializer_list<Species> s);

	void operator()(const double dt);
	std::string get_time();
	std::string get_global_time();
	std::string get_time_percentage();
	void reset();
	int get_species_index(Species& s);
protected:
	void resume_timer();
	void stop_timer();

	boost::timer timer;
	double total_time;
	static double total_global_time;
	double time;
	//boost::timer global_timer;
	std::vector<Species*> all_species;
};

class CountMolsOnGrid: public Operator {
public:
	CountMolsOnGrid() {}
	void operator()(const double dt);
};

std::ostream& operator<< (std::ostream& out, CountMolsOnGrid& b);
}

#endif /* OPERATOR_H_ */
