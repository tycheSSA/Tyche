/*
 * Boundary.impl.h
 *
 *  Created on: 19 Oct 2012
 *      Author: robinsonm
 */

#ifndef BOUNDARY_IMPL_H_
#define BOUNDARY_IMPL_H_

#include <boost/foreach.hpp>
#include <set>
#include "Vector.h"
#include "Boundary.h"
#include "Log.h"

namespace Tyche {

template<typename T>
void JumpBoundary<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	BOOST_FOREACH(Species *s, this->all_species) {
		Molecules& mols = s->mols;
		const int n = mols.size();
		for (int i = 0; i < n; ++i) {
			while (this->geometry.at_boundary(mols.r[i])) {
				mols.r[i] += jump_by;
				//mols.saved_index[i] = SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE;
			}
		}
	}
	LOG(2,"Stopping Operator: " << *this);
	Operator::stop_timer();
}

template<typename T>
void DiffusionCorrectedBoundary<T>::timestep_initialise(const double dt) {
   const int n = this->all_species.size();
   for (int i = 0; i < n; ++i) {
      Species& s = *(this->all_species[i]);
      std::vector<double>& prev_distance = *(all_prev_distance[i]);
      std::vector<double>& curr_distance = *(all_curr_distance[i]);
      Molecules& mols = s.mols;

      recalc_constants(s, dt);
      if (prev_distance.size() == 0) init_prev_distance(mols, prev_distance);
      const int nmol = mols.size();
      curr_distance.resize(nmol);
      for (int ii = 0; ii < nmol; ++ii) {
         const double dist_to_wall = this->geometry.distance_to_boundary(mols.r[ii]);
         curr_distance[ii] = dist_to_wall;
      }
   }
}

template<typename T>
void DiffusionCorrectedBoundary<T>::timestep_finalise() {
   const int n = this->all_species.size();
   for (int i = 0; i < n; ++i) {
      std::vector<double> *tmp = all_prev_distance[i];
      all_prev_distance[i] = all_curr_distance[i];
      all_curr_distance[i] = tmp;
   }
}

template<typename T>
bool DiffusionCorrectedBoundary<T>::particle_crossed_boundary(const int p_i, const int s_i) {

	const double dist_to_wall = (*all_curr_distance[s_i])[p_i];
	if (dist_to_wall <= 0) {
		return true;
	} else if (dist_to_wall < test_this_distance_from_wall) {
	   const int saved_index = this->all_species[s_i]->mols.saved_index[p_i];
	   if (saved_index == SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE) return false;
	   const double old_dist_to_wall = (*all_prev_distance[s_i])[saved_index];
	   if (old_dist_to_wall <= 0) return false;
		const double P = exp(-dist_to_wall*old_dist_to_wall/D_dt);
		if (uni() < P) {
			return true;
		}
	}
	return false;
}

template<typename T>
void DiffusionCorrectedBoundary<T>::add_species(Species &s) {
   Boundary<T>::add_species(s);
   std::vector<double>* prev_distance = new std::vector<double>;
   std::vector<double>* curr_distance = new std::vector<double>;
   init_prev_distance(s.mols, *prev_distance);
   all_curr_distance.push_back(curr_distance);
   all_prev_distance.push_back(prev_distance);
}

template<typename T>
void JumpBoundaryWithCorrection<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	DiffusionCorrectedBoundary<T>::timestep_initialise(dt);

	const int s_n = this->all_species.size();
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->all_species[s_i]);
		const int p_n = s.mols.size();
		for (int p_i = 0; p_i < p_n; ++p_i) {
			if (DiffusionCorrectedBoundary<T>::particle_crossed_boundary(p_i, s_i)) {
				s.mols.r[p_i] += jump_by;
				while (this->geometry.at_boundary(s.mols.r[p_i])) {
					s.mols.r[p_i] += jump_by;
				}
				//s.mols.saved_index[p_i] = SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE;
			}
		}
	}

	DiffusionCorrectedBoundary<T>::timestep_finalise();
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

template<typename T>
void RemoveBoundaryWithCorrection<T>::add_species(Species& s) {
	DiffusionCorrectedBoundary<T>::add_species(s);
	removed_molecules.push_back(Molecules());
}

template<typename T>
Molecules& RemoveBoundaryWithCorrection<T>::get_removed(Species& s) {
	const int s_i = Operator::get_species_index(s);
	return removed_molecules[s_i];
}

template<typename T>
void RemoveBoundaryWithCorrection<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	DiffusionCorrectedBoundary<T>::timestep_initialise(dt);

	const int s_n = this->all_species.size();
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->all_species[s_i]);
		const int p_n = s.mols.size();
		for (int p_i = 0; p_i < p_n; ++p_i) {
			if (DiffusionCorrectedBoundary<T>::particle_crossed_boundary(p_i, s_i)) {
				s.mols.mark_for_deletion(p_i);
				removed_molecules[s_i].add_molecule(s.mols.r[p_i]);
			}
		}
		s.mols.delete_molecules();
	}
	DiffusionCorrectedBoundary<T>::timestep_finalise();
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

template<typename T>
void RemoveBoundary<T>::add_species(Species& s) {
	Boundary<T>::add_species(s);
	removed_molecules.push_back(Molecules());
}

template<typename T>
Molecules& RemoveBoundary<T>::get_removed(Species& s) {
	const int s_i = Operator::get_species_index(s);
	return removed_molecules[s_i];
}

template<typename T>
void RemoveBoundary<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);

	const int s_n = this->all_species.size();
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->all_species[s_i]);
		const int p_n = s.mols.size();
		for (int p_i = 0; p_i < p_n; ++p_i) {
			if (this->geometry.at_boundary(s.mols.r[p_i])) {
				s.mols.mark_for_deletion(p_i);
				removed_molecules[s_i].add_molecule(s.mols.r[p_i]);
			}
		}
		s.mols.delete_molecules();
	}
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

template<typename T>
void ReflectiveBoundary<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	BOOST_FOREACH(Species *s, this->all_species) {
		Molecules& mols = s->mols;
		const int n = mols.size();
		for (int i = 0; i < n; ++i) {
			Vect3d nv;
			if (this->geometry.at_boundary(mols.r[i],nv)) {
				mols.r[i] += 2.0*nv;
				//mols.saved_index[i] = SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE;
			}
		}
	}
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}


template<typename T>
void DestroyBoundary<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	BOOST_FOREACH(Species *s, this->all_species) {
		Molecules& mols = s->mols;
		for (int i = 0; i < mols.size(); ++i) {
			if (Boundary<T>::geometry.at_boundary(mols.r[i])) {
				mols.delete_molecule(i);
			}
		}
	}
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

template<typename T>
void CouplingBoundary_M_to_C<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	DiffusionCorrectedBoundary<T>::timestep_initialise(dt);
	const int s_n = this->all_species.size();
	std::set<int> dirty_indicies;
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->all_species[s_i]);
		const int p_n = s.mols.size();
		for (int p_i = 0; p_i < p_n; ++p_i) {
			if (DiffusionCorrectedBoundary<T>::particle_crossed_boundary(p_i, s_i)) {
				const double dist_to_wall = (*(this->all_curr_distance[s_i]))[p_i];
				//if (dist_to_wall < 0) {
				//if (dist_to_wall < 0) {
				const Vect3d r = s.mols.r[p_i];
				Vect3d comp_r = r;
				if (dist_to_wall > 0) {
					comp_r += (1.0 + s.grid.get_tolerance())*this->geometry.shortest_vector_to_boundary(r);
				}
//				else {
//					comp_r += (1.0 - s.grid.get_tolerance())*this->geometry.shortest_vector_to_boundary(r);
//				}
				const int i = s.grid.get_cell_index(comp_r);
				dirty_indicies.insert(i);
				//dirty_indicies.push_back(i);
				s.copy_numbers[i]++;
				s.mols.mark_for_deletion(p_i);
			}
		}

		s.mols.delete_molecules();
		//std::cout << count << "particles moved to compartments. Free space = "<<s.mols.size() <<" compart = "<< std::accumulate(s.copy_numbers.begin(),s.copy_numbers.end(),0) << std::endl;
	}
	BOOST_FOREACH(int i, dirty_indicies) {
		nsm.recalc_priority(i);
	}
	DiffusionCorrectedBoundary<T>::timestep_finalise();
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

template<typename T>
void CouplingBoundary_C_to_M<T>::add_species(Species& s, const double dt) {
	Boundary<T>::add_species(s);
	/*
	 * setup interface reactions in the boundary grid cells
	 */
	nsm.set_interface(s,this->geometry,dt);
//	boundary_compartment_indicies.push_back(std::vector<int>());
//	boundary_intersections.push_back(std::vector<typename T::SurfaceElementType>());
//	std::vector<int>& boundary_compartment_indicies_ref = *(boundary_compartment_indicies.end()-1);
//	std::vector<typename T::SurfaceElementType>& boundary_intersections_ref = *(boundary_intersections.end()-1);
//	const int nc = s.copy_numbers.size();
//	for (int i = 0; i < nc; ++i) {
//		if (s.grid.geometry_intersects_cell(i, this->geometry)) {
//			boundary_compartment_indicies_ref.push_back(i);
//			boundary_intersections_ref.push_back(s.grid.get_intersection_of_cell(i,this->geometry));
//		}
//	}
}


template<typename T>
void CouplingBoundary_C_to_M<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	const int s_n = this->all_species.size();
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->all_species[s_i]);
		/*
		 * find where boundaries intersect
		 */
		std::vector<typename T::SurfaceElementType> boundary_intersections;
		std::vector<int> boundary_compartment_indicies;
		s.grid.get_slice(this->geometry,boundary_compartment_indicies);
		BOOST_FOREACH(int i,boundary_compartment_indicies) {
			boundary_intersections.push_back(s.grid.get_intersection_of_cell(i,this->geometry));
		}

		/*
		 * setup interface reactions in the boundary grid cells
		 */
		if (old_dt != dt) {
			nsm.set_interface_reactions(s,boundary_compartment_indicies,dt);
			old_dt = dt;
		}

		/*
		 * move particles in compartments to free space
		 */
		const double step_length = sqrt(2.0*s.D*dt);
		const int c_n = boundary_compartment_indicies.size();
		for (int c_i = 0; c_i < c_n; ++c_i) {
			const int num_particles_to_add = s.copy_numbers[boundary_compartment_indicies[c_i]];
			typename T::SurfaceElementType& intersect =  boundary_intersections[c_i];
			for (int i = 0; i < num_particles_to_add; ++i) {
				Vect3d newr,newn;
				intersect.get_random_point_and_normal_triangle(newr, newn);
				const double P = uni();
				const double P2 = pow(P,2);
				const double dist_from_intersect = step_length*(0.729614*P - 0.70252*P2)/(1.0 - 1.47494*P + 0.484371*P2);
				newr += newn*dist_from_intersect;
				s.mols.add_molecule(newr);
			}
			s.copy_numbers[boundary_compartment_indicies[c_i]] = 0;
		}
		//std::cout << count << "particles moved to free-space. Free space = "<<s.mols.size() <<" compart = "<< std::accumulate(s.copy_numbers.begin(),s.copy_numbers.end(),0) << std::endl;
	}
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

}

#endif /* BOUNDARY_IMPL_H_ */
