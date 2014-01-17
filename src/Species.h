/*
 * species.h
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
 *  Created on: 11 Oct 2012
 *      Author: robinsonm
 */

#ifndef SPECIES_H_
#define SPECIES_H_

#include <vector>
#include <boost/foreach.hpp>
#include "MyRandom.h"
#include "Vector.h"
#include "StructuredGrid.h"

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>



namespace Tyche {

#define DATA_typename ParticleInfo
#define DATA_names   (r)(r0)(alive)(id)(saved_index)
#define DATA_types   (Vect3d)(Vect3d)(bool)(int)(int)
#include "Data.h"

template<int DataSize = 1>
class Particles {
public:
	static const int SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE = -1;

	Particles() {
		next_id = 0;
	}
	void fill_uniform(const Vect3d low, const Vect3d high, const unsigned int N) {
		//TODO: assumes a 3d rectangular region
		boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, boost::uniform_real<>(0,1));
		const Vect3d dist = high-low;
		for(int i=0;i<N;i++) {
			add_particle(Vect3d(uni()*dist[0],uni()*dist[1],uni()*dist[2])+low);
		}
	}
	void delete_particle(const unsigned int i) {
		const int last_index = info.size()-1;
		if (i != last_index) {
			info[i] = info[last_index];
			data[i] = data[last_index];
		}
		info.pop_back();
		data.pop_back();
	}
	void delete_particles() {
		int i = 0;
		while (i < info.size()) {
			if (!info.alive[i]) {
				delete_particle(i);
			} else {
				i++;
			}
		}
	}
	void clear() {
		info.clear();
		data.clear();
	}
	size_t size() const {
		return info.size();
	}
	void add_particle(const Vect3d& position) {
		info.push_back(position, position, true, next_id++, SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE);
		data.resize(info.size());
	}
	void add_particle(const Vect3d& position, const Vect3d& old_position) {
		info.push_back(position, old_position, true, next_id++, SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE);
		data.resize(info.size());
	}
	const std::vector<Vect3d>& get_position() const {
		return info.r;
	}
	const Vect3d& get_position(const unsigned int i) const {
		return info.r[i];
	}
	Vect3d& get_position_nonconst(const unsigned int i) {
		return info.r[i];
	}
	const Vect3d& get_old_position(const unsigned int i) const {
		return info.r0[i];
	}
	Vect3d& get_old_position_nonconst(const unsigned int i) {
		return info.r0[i];
	}
	const double* get_data(const unsigned int i) const {
		return data[i];
	}
	double* get_data_noncost(const unsigned int i) {
		return data[i];
	}
	bool is_alive(const unsigned int i) const {
		return info.alive[i];
	}
	const unsigned int& get_id(const unsigned int i) const {
		return info.id[i];
	}

	void mark_for_deletion(const unsigned int i) {
		info.alive[i] = false;
	}

	void save_indicies() {
		const int n = info.size();
		for (int i = 0; i < n; ++i) {
			info.saved_index[i] = i;
		}
	}

	vtkSmartPointer<vtkUnstructuredGrid> get_vtk_grid() {
		vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkIntArray> newInt = vtkSmartPointer<vtkIntArray>::New();
		newInt->SetName("id");
		const vtkIdType n = size();
		newPts->SetNumberOfPoints(n);
		newInt->SetNumberOfValues(n);
		for (int i = 0; i < n; ++i) {
			//std::cout << "adding mol to vtk at position "<<mols.r[i]<<std::endl;
			newPts->SetPoint(i,get_position(i)[0],get_position(i)[1],get_position(i)[2]);
			newInt->SetValue(n,info.id[i]);
		}
		newPts->ComputeBounds();

		grid->SetPoints(newPts);
		grid->GetPointData()->AddArray(newInt);

		return grid;
	}
private:
	ParticleInfo info;
	std::vector<boost::array<double, DataSize> > data;
	int next_id;
};





static const StructuredGrid empty_grid;

class Species {
public:
	Species(double D):D(D),grid(NULL) {
		id = species_count++;
		clear();
	}
	Species(double D, const StructuredGrid* grid):D(D),grid(grid)  {
		id = species_count++;
		clear();
	}
	static std::auto_ptr<Species> New(double D) {
		return std::auto_ptr<Species>(new Species(D));
	}
	void clear() {
		mols.clear();
		if (grid!=NULL) copy_numbers.assign(grid->size(),0);
	}
	void set_grid(const StructuredGrid* new_grid) {
		grid = new_grid;
		if (grid!=NULL) copy_numbers.assign(grid->size(),0);
	}
	//void fill_uniform(const Vect3d low, const Vect3d high, const unsigned int N);
	template<typename T>
	void fill_uniform(const T& geometry, const unsigned int N);
	template<typename T>
	unsigned int count_mols_in(const T& geometry);
	void get_concentrations(const StructuredGrid& calc_grid, std::vector<double>& mol_concentrations, std::vector<double>& compartment_concentrations) const;
	void get_concentration(const StructuredGrid& calc_grid, std::vector<double>& concentration) const;
	void get_concentration(const Vect3d low, const Vect3d high, const Vect3i n, std::vector<double>& concentration) const;
	vtkSmartPointer<vtkUnstructuredGrid> get_vtk();
	std::string get_status_string() const;
	friend std::ostream& operator<<( std::ostream& out, const Species& b ) {
		return out << b.get_status_string();
	}

	double D;
	Particles<1> mols;
	std::vector<int> copy_numbers;
	std::vector<int> mol_copy_numbers;
	const StructuredGrid* grid;
	int id;
	std::vector<double> tmpx,tmpy,tmpz;
private:
	static int species_count;
};

template<typename T>
void Species::fill_uniform(const T& geometry, const unsigned int N) {
	LOG(2,"Adding "<<N<<" molecules of Species ("<<id<<") within "<<geometry);
	for(int i=0;i<N;i++) {
		mols.add_particle(geometry.get_random_point_in());
	}
}

template<typename T>
unsigned int Species::count_mols_in(const T& geometry) {
	unsigned int count;
	const int n = mols.size();
	for (int i = 0; i < n; ++i) {
		if (geometry.is_in(mols.get_position(i))) count++;
	}
	return count;
}


extern Species null_species;
}


#endif /* SPECIES_H_ */
