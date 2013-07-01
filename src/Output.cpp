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


void OutputConcentrations::integrate(const double dt) {

	if (is_execute_time()) {
		const int n = get_species().size();
		ASSERT(n==1,"OutputConcentrations only implemented for one species");
		Species &s = *(get_species()[0]);
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
		write(get_time());
	}
}


void OutputSumConcentrations::integrate(const double dt) {

	if (is_execute_time()) {
		double sumsum_m = 0;
		double sumsum_c = 0;
		const int n = get_species().size();
		for (int i = 0; i < n; ++i) {
			Species &s = *(get_species()[i]);

			std::vector<double> mol_con, comp_con;
			s.get_concentrations(grid,mol_con,comp_con);
			ASSERT(grid.size()==mol_con.size(), "copy numbers size is not the same as grid!");
			const int ngrid = mol_con.size();
			double sum_m = 0;
			double sum_c = 0;
			for (int i = 0; i < ngrid; ++i) {
				sum_m += mol_con[i]*grid.get_cell_volume(i);
				sum_c += comp_con[i]*grid.get_cell_volume(i);
			}

			std::ostringstream ss;
			ss << s.id;
			data["Concentration(M)("+ss.str()+")"].push_back(sum_m);
			data["Concentration(C)("+ss.str()+")"].push_back(sum_c);
			data["Concentration("+ss.str()+")"].push_back(sum_m+sum_c);

			sumsum_m += sum_m;
			sumsum_c += sum_c;
		}

		data["Time"].push_back(get_time());
		data["Concentration(M)"].push_back(sumsum_m);
		data["Concentration(C)"].push_back(sumsum_c);
		data["Concentration"].push_back(sumsum_m+sumsum_c);
		write();
	}

}

void OutputSumConcentrations::add_species_execute(Species &s) {
	std::ostringstream ss;
	ss << s.id;
	data.insert(std::pair<std::string,std::vector<double> >("Concentration(M)("+ss.str()+")",std::vector<double>()));
	data.insert(std::pair<std::string,std::vector<double> >("Concentration(C)("+ss.str()+")",std::vector<double>()));
	data.insert(std::pair<std::string,std::vector<double> >("Concentration("+ss.str()+")",std::vector<double>()));

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









