/*
 * species.cpp
 *
 *  Created on: 11 Oct 2012
 *      Author: robinsonm
 */

#include "Species.h"
#include <boost/random.hpp>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

#include <vtkPointData.h>


namespace Tyche {
int Species::species_count = 0;
Species null_species(0);

vtkSmartPointer<vtkUnstructuredGrid> Molecules::get_vtk_grid() {
	/*
	 * setup points
	 */
	vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
	const int num_points = r.size();
	for (int i = 0; i < num_points; i++) {
		newPts->InsertNextPoint(r[i][0],r[i][1],r[i][2]);
	}

	/*
	 * setup grid
	 */
	vtkSmartPointer<vtkUnstructuredGrid> vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtk_grid->SetPoints(newPts);

	return vtk_grid;
}


int Molecules::delete_molecule(const unsigned int i) {
	const int last_index = this->size()-1;
	if (i != last_index) {
		(*this)[i] = (*this)[last_index];
	}
	this->pop_back();
}


void Molecules::fill_uniform(const Vect3d low, const Vect3d high,
		const unsigned int n) {
	//TODO: assumes a 3d rectangular region
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, boost::uniform_real<>(0,1));
	const Vect3d dist = high-low;
	for(int i=0;i<n;i++) {
		add_molecule(Vect3d(uni()*dist[0],uni()*dist[1],uni()*dist[2])+low);
	}
}

void Species::fill_uniform(const Vect3d low, const Vect3d high, const unsigned int N) {
	LOG(2,"Adding "<<N<<" molecules of Species ("<<id<<") within the rectangle defined by "<<low<<" and "<<high);
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, boost::uniform_real<>(0,1));
	const Vect3d dist = high-low;
	for(int i=0;i<N;i++) {
		const Vect3d pos = Vect3d(uni()*dist[0],uni()*dist[1],uni()*dist[2])+low;
		mols.add_molecule(pos);
	}
}


vtkSmartPointer<vtkUnstructuredGrid> Species::get_vtk() {
	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkIntArray> newInt = vtkSmartPointer<vtkIntArray>::New();
	newInt->SetName("id");
	const vtkIdType n = mols.size();
	newPts->SetNumberOfPoints(n);
	newInt->SetNumberOfValues(n);
	for (int i = 0; i < n; ++i) {
		//std::cout << "adding mol to vtk at position "<<mols.r[i]<<std::endl;
		newPts->SetPoint(i,mols.r[i][0],mols.r[i][1],mols.r[i][2]);
		newInt->SetValue(n,mols.id[i]);
	}
	newPts->ComputeBounds();

	grid->SetPoints(newPts);
	grid->GetPointData()->AddArray(newInt);

	return grid;
}

int Molecules::add_molecule(const Vect3d& position) {
	this->push_back(position, position, true, next_id++, SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE);
}
int Molecules::add_molecule(const Vect3d& position, const Vect3d& old_position) {
	this->push_back(position, old_position, true, next_id++, SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE);
}

int Molecules::delete_molecules() {
   int i = 0;
   while (i < this->size()) {
		if (!alive[i]) {
			delete_molecule(i);
		} else {
		   i++;
		}
	}
}

int Molecules::mark_for_deletion(const unsigned int i) {
	alive[i] = false;
}

void Molecules::save_indicies() {
	const int n = this->size();
	for (int i = 0; i < n; ++i) {
		saved_index[i] = i;
	}
}


void Species::get_concentrations(const StructuredGrid& calc_grid,
		std::vector<double>& mol_concentrations,
		std::vector<double>& compartment_concentrations) const {

	mol_concentrations.assign(calc_grid.size(),0);
	BOOST_FOREACH(Vect3d r,mols.r) {
		if (calc_grid.is_in(r)) {
			mol_concentrations[calc_grid.get_cell_index(r)]++;
		}
	}

	const int n = calc_grid.size();
	compartment_concentrations.assign(calc_grid.size(),0);
	if ((grid!=NULL)&&(copy_numbers.size() != 0)) {
		for (int i = 0; i < n; ++i) {
			std::vector<int> indicies;
			std::vector<double> volume_ratio;
			grid->get_overlap(calc_grid.get_low_point(i),calc_grid.get_high_point(i),indicies,volume_ratio);
			const int noverlap = indicies.size();
			//double sum_of_volume_ratios = 0;
			for (int j = 0; j < noverlap; ++j) {
				//sum_of_volume_ratios += volume_ratio[j];
				compartment_concentrations[i] += copy_numbers[indicies[j]]*volume_ratio[j];
				//std::cout << " compartment "<<i<<" overlap with compartment "<<indicies[j]<<" with volume ratio"<<volume_ratio[j]<<std::endl;
			}
			//std::cout <<" sum of vol ratio = "<<sum_of_volume_ratios<<std::endl;
		}
	}

	for (int i = 0; i < n; ++i) {
		const double inv_vol = 1.0/calc_grid.get_cell_volume(i);
		compartment_concentrations[i] *= inv_vol;
		mol_concentrations[i] *= inv_vol;
	}
}

void Species::get_concentration(const StructuredGrid& calc_grid,
		std::vector<double>& concentration) const {
	std::vector<double> compartment_concentrations;
	get_concentrations(calc_grid,concentration,compartment_concentrations);
	const int n = concentration.size();
//	double totalm = 0;
//	double totalc = 0;
//	double totalv = 0;
	for (int i = 0; i < n; ++i) {
//		if (i >= n/2) {
//			totalm += concentration[i]*calc_grid.get_cell_volume(i);
//			totalc += compartment_concentrations[i]*calc_grid.get_cell_volume(i);
//			totalv += calc_grid.get_cell_volume(i);
//		}
		concentration[i] += compartment_concentrations[i];
	}
	//std::cout <<" there are "<<totalm<<" particles and "<<totalc<<" compartment mols in "<<totalv<<" volume"<<std::endl;
}

void Species::get_concentration(const Vect3d low, const Vect3d high, const Vect3i n,
		std::vector<double>& concentration) const {

	const Vect3d spacing = (high-low).cwiseQuotient(n.cast<double>());
	StructuredGrid calc_grid(low,high,spacing);
	get_concentration(calc_grid,concentration);

//	const Vect3d inv_spacing = Vect3d(1,1,1).cwiseQuotient(spacing);
//
//	const double volume = spacing.prod();
//	const int nnn = n.prod();
//	const int num_cells_along_yz = n[2]*n[1];
//
//
//	concentration.assign(nnn,0);
//	for (Vect3d r: mols.r) {
//		if (((r.array() >= low.array()).all()) && ((r.array() < high.array()).all())) {
//			const Vect3i celli = ((r-low).cwiseProduct(inv_spacing)).cast<int>();
//			const int index = celli[0] * num_cells_along_yz + celli[1] * n[1] + celli[2];
//			concentration[index]++;
//		}
//	}
}

std::string Species::get_status_string() const {
	const int ncomp = std::accumulate(copy_numbers.begin(),copy_numbers.end(),0);
	std::ostringstream ss;
	ss << "Species "<<id<<":\t" << mols.size() + ncomp << " particles.";
	return ss.str();
}

}




