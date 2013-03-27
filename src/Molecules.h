/*
 * MoleculesSimple.h
 * 
 * Copyright 2013 Martin Robinson
 *
 * This file is part of PDE_BD.
 *
 * PDE_BD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PDE_BD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PDE_BD.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 21 Feb 2013
 *      Author: robinsonm
 */

#ifndef MOLECULES_H_
#define MOLECULES_H_

#include <vector>
#include "zip.h"
#include <algorithm>

namespace Tyche {

class Molecules {
public:
	typedef double ST;
	typedef std::vector<ST> vector;
	Molecules();
	void add_particle(const ST x, const ST Y, const ST z);
	void remove_particle(const int i);
	void remove_particles(std::vector<int>& to_delete);

	void diffuse(const double dt, const double D);
	void reflective_boundaries(const double xmin,const double xmax,
			const double ymin, const double ymax,
			const double zmin, const double zmax);
	vector& get_x() {return x;}
	vector& get_y() {return y;}
	vector& get_z() {return z;}
	const int size() {return x.size();}
private:
	vector x,y,z;
};

}

#endif /* MOLECULESSIMPLE_H_ */
