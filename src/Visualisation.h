/*
 * Visualisation.h
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
 *  Created on: 22 Nov 2012
 *      Author: robinsonm
 */

#ifndef VISUALISATION_H_
#define VISUALISATION_H_

#include "Operator.h"

namespace Tyche {
class Visualisation: public Operator {
public:
	Visualisation(const double dt);
	void operator()(const double dt);
	virtual void print(std::ostream& out);
	void add_geometry(const xplane& geometry);
	void add_geometry(const yplane& geometry);
	void add_geometry(const zplane& geometry);
	void add_compartment_slice(const xplane& geometry);
	void add_compartment_slice(const yplane& geometry);
	void add_compartment_slice(const zplane& geometry);
	void reset() {
		Operator::reset();
		next_vis = time + vis_dt;
	}

private:
	void add_planes_to_vis();
	void add_plots_to_vis();
	void add_molecules_to_vis(Species& s);
	void add_compartments_to_vis(Species& s);
	void show_vis();
	std::vector<xplane> xplanes;
	std::vector<yplane> yplanes;
	std::vector<zplane> zplanes;
	std::vector<xplane> compartment_slices_x;
	std::vector<yplane> compartment_slices_y;
	std::vector<zplane> compartment_slices_z;
	Vect3d low, high;
	double vis_dt;
	double next_vis;
	struct MyVTKdata;
	MyVTKdata* my_vtk_data;
};


class Plot2d: public Operator {
public:
	Plot2d(const double dt, const std::vector<double>& x, const std::vector<double>& y,
			const char* x_label, const char* y_label, const char* title);
	void operator()(const double dt);
	virtual void print(std::ostream& out) {
		out << "\t2d plot showing \"" << x_label << "\" versus \"" << y_label << "\" entitled \"" << title << "\"";
	}
	void reset() {
		Operator::reset();
		next_vis = time + vis_dt;
	}

private:
	const std::vector<double>& x;
	const std::vector<double>& y;
	std::string x_label;
	std::string y_label;
	std::string title;
	double xlow,xhigh,ylow,yhigh;
	double vis_dt;
	double next_vis;
	struct MyVTKdata;
	MyVTKdata* my_vtk_data;
};

}
#endif /* VISUALISATION_H_ */
