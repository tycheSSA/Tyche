/*
 * Output.cpp
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
 *  Created on: Nov 18, 2012
 *      Author: mrobins
 */

#include "Output.h"
#include <numeric>

namespace Tyche {

void OutputMolecularConcentrations::operator ()(const double dt) {
	Operator::resume_timer();
	OutputMultiFilename::operator ()(dt);
	if (f.is_open()) {
		LOG(2, "Starting Operator: " << *this);
		const int n = all_species.size();
		ASSERT(n==1,"OutputConcentrations only implemented for one species");
		Species &s = *(all_species[0]);

		bins.assign(num_bins,0);
		const double scale = num_bins/(high-low);
		BOOST_FOREACH(Vect3d r,s.mols.r) {
			const double scaled_position = (r[0]-low)*scale;
			if ((scaled_position < 0) || (scaled_position > num_bins)) {
				//printf("outside area: position = (%f %f %f)\n",r[0],r[1],r[2]);
				continue;
			}
			const int index = int(scaled_position);
			bins[index]++;
		}

		for (int i = 0; i < num_bins; ++i) {
			f << (i+0.5)*(high-low)/num_bins + low << " 0 0 " << bins[i] << std::endl;
		}

		LOG(2, "Stopping Operator: " << *this);
	}
	Operator::stop_timer();
}

std::ostream& operator <<(std::ostream& out, OutputMolecularConcentrations& b) {
	return out << "\tOutput molecular concentrations";
}

void OutputCompartmentConcentrations::operator ()(const double dt) {
	Operator::resume_timer();
	OutputMultiFilename::operator ()(dt);
	if (f.is_open()) {
		LOG(2, "Starting Operator: " << *this);
		const int n = all_species.size();
		ASSERT(n==1,"OutputConcentrations only implemented for one species");
		Species &s = *(all_species[0]);
		const int num_bins = s.grid.get_cells_along_axes()[0];
		bins.assign(num_bins,0);
		const int nc = s.grid.size();
		for(int i = 0; i < nc; i++) {
			Vect3i cell_indices = s.grid.get_cell_indicies(i);
			bins[cell_indices[0]] += s.copy_numbers[i];
			//if ((cell_indices[0]==5)||(cell_indices[0]==4))  std::cout <<  "copy number for "<<cell_indices[0]<<" = " << s.copy_numbers[i] << std::endl;
		}

		const double low = s.grid.get_low()[0];
		const double high = s.grid.get_high()[0];
		for (int i = 0; i < num_bins; ++i) {
			f << (i+0.5)*(high-low)/num_bins + low << " 0 0 " << bins[i] << std::endl;
		}

		LOG(2, "Stopping Operator: " << *this);
	}
	Operator::stop_timer();
}

std::ostream& operator <<(std::ostream& out,
		OutputCompartmentConcentrations& b) {
	return out << "\tOutput compartments concentrations";
}

void OutputConcentrations::operator ()(const double dt) {
	Operator::resume_timer();
	Output::operator ()(dt);
	if (is_execute_time()) {
		LOG(2, "Starting Operator: " << *this);
		const int n = all_species.size();
		ASSERT(n==1,"OutputConcentrations only implemented for one species");
		Species &s = *(all_species[0]);
		std::vector<double> mol_con, comp_con;
		s.get_concentrations(grid,mol_con,comp_con);
		ASSERT(grid.size()==mol_con.size(), "copy numbers size is not the same as grid!");
		const int ngrid = mol_con.size();
		for (DataType::iterator i = data.begin(); i != data.end(); i++) {
			i->second.clear();
		}
		for (int i = 0; i < ngrid; ++i) {
			const Vect3d cc = grid.get_cell_centre(i);
			data["x"].push_back(cc[0]);
			data["y"].push_back(cc[1]);
			data["z"].push_back(cc[2]);
			data["Concentration(M)"].push_back(mol_con[i]);
			data["Concentration(C)"].push_back(comp_con[i]);
			data["Concentration"].push_back(mol_con[i] + comp_con[i]);
		}
		write(time);
		LOG(2, "Stopping Operator: " << *this);
	}
	Operator::stop_timer();
}

std::ostream& operator <<(std::ostream& out, OutputConcentrations& b) {
	return out << "\tOutput concentrations";
}



void Output::write() {
	write(filename + ".dat");
}

void Output::write(double time) {
	std::ostringstream s;
	s << filename << "_time_" << time  << ".dat";
	write(s.str());
}

void Output::write(std::string filename) {
	f.open(filename.c_str());
	f << "# ";
	const int n = data.begin()->second.size();
	for (DataType::const_iterator i = data.begin(); i != data.end(); i++) {
		ASSERT(n == i->second.size(), "all data assumed to be the same length");
		f << i->first << ' ';
	}
	f << std::endl;
	for (int i = 0; i < n; ++i) {
		for (DataType::const_iterator j = data.begin(); j != data.end(); j++) {
			f << j->second[i] << ' ';
		}
		f << std::endl;
	}
	f.close();
}

}









