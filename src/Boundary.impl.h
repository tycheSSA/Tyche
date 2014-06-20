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
void JumpBoundary<T>::integrate(const double dt) {
	BOOST_FOREACH(Species *s, this->get_species()) {
		Molecules& mols = s->mols;
		const int n = mols.size();
		for (int i = 0; i < n; ++i) {
			if (this->geometry.lineXsurface(mols.r0[i],mols.r[i])) {
				mols.r[i] += jump_by;
				//mols.r0[i] += jump_by;
				//mols.saved_index[i] = SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE;
			}
//			while (this->geometry.distance_to_boundary(mols.r[i]) < 0) {
//				mols.r[i] += jump_by;
//				mols.r0[i] += jump_by;
//				//mols.saved_index[i] = SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE;
//			}
		}
	}

}

template<typename T>
void DiffusionCorrectedBoundary<T>::timestep_initialise(const double dt) {
   const int n = this->get_species().size();
   for (int i = 0; i < n; ++i) {
      Species& s = *(this->get_species()[i]);
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
   const int n = this->get_species().size();
   for (int i = 0; i < n; ++i) {
      std::vector<double> *tmp = all_prev_distance[i];
      all_prev_distance[i] = all_curr_distance[i];
      all_curr_distance[i] = tmp;
   }
}

template<typename T>
bool DiffusionCorrectedBoundary<T>::particle_crossed_boundary(const int p_i, const int s_i) {
	Molecules& mols = this->get_species()[s_i]->mols;
	const double dist_to_wall = (*all_curr_distance[s_i])[p_i];
	if (this->geometry.lineXsurface(mols.r0[p_i],mols.r[p_i])) {
		return true;
	} else if (dist_to_wall < test_this_distance_from_wall) {
	   const int saved_index = mols.saved_index[p_i];
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
void DiffusionCorrectedBoundary<T>::add_species_execute(Species &s) {
   std::vector<double>* prev_distance = new std::vector<double>;
   std::vector<double>* curr_distance = new std::vector<double>;
   init_prev_distance(s.mols, *prev_distance);
   all_curr_distance.push_back(curr_distance);
   all_prev_distance.push_back(prev_distance);
}

template<typename T>
void JumpBoundaryWithCorrection<T>::integrate(const double dt) {

	DiffusionCorrectedBoundary<T>::timestep_initialise(dt);

	const int s_n = this->get_species().size();
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->get_species()[s_i]);
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

}

template<typename T>
void RemoveBoundaryWithCorrection<T>::add_species_execute(Species& s) {

	removed_molecules.push_back(Molecules());
}

template<typename T>
Molecules& RemoveBoundaryWithCorrection<T>::get_removed(Species& s) {
	const int s_i = Operator::get_species_index(s);
	return removed_molecules[s_i];
}

template<typename T>
void RemoveBoundaryWithCorrection<T>::integrate(const double dt) {

	DiffusionCorrectedBoundary<T>::timestep_initialise(dt);

	const int s_n = this->get_species().size();
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->get_species()[s_i]);
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

}

template<typename T>
void RemoveBoundary<T>::add_species_execute(Species& s) {

	removed_molecules.push_back(Molecules());
}

template<typename T>
Molecules& RemoveBoundary<T>::get_removed(Species& s) {
	const int s_i = Operator::get_species_index(s);
	return removed_molecules[s_i];
}

template<typename T>
void RemoveBoundary<T>::integrate(const double dt) {


	const int s_n = this->get_species().size();
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->get_species()[s_i]);
		Molecules& mols = s.mols;
		const int p_n = s.mols.size();
		for (int p_i = 0; p_i < p_n; ++p_i) {
			if (this->geometry.lineXsurface(mols.r0[p_i],mols.r[p_i])) {
				s.mols.mark_for_deletion(p_i);
				removed_molecules[s_i].add_molecule(s.mols.r[p_i],mols.r0[p_i]);
			}
		}
		s.mols.delete_molecules();
	}

}

template<typename T>
void ReflectiveBoundary<T>::integrate(const double dt) {

	BOOST_FOREACH(Species *s, this->get_species()) {
		Molecules& mols = s->mols;
		const int n = mols.size();
		for (int i = 0; i < n; ++i) {
			Vect3d nv,ip;
			int iteration = 0;
			while (this->geometry.lineXsurface(mols.r0[i],mols.r[i],&ip,&nv)&&(iteration++<3)) {
				Vect3d v = ip-mols.r0[i];
				v = v-2.0*(v.dot(nv)*nv);
				v /= v.norm();
				//std::cout << "intersection from "<<mols.r0[i]<<" to "<<mols.r[i]<<" collision at "<<ip<<" reflecting to "<<ip+v<<std::endl;
				//mols.r0[i] = ip + GEOMETRY_TOLERANCE*v;
				mols.r[i] = ip + (mols.r[i] - ip).norm()*v;
				//mols.saved_index[i] = SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE;
			}
//			Vect3d nv;
//			const double d = this->geometry.distance_to_boundary(mols.r[i]);
//			if (d<0) {
//				const Vect3d vect_to_wall = this->geometry.shortest_vector_to_boundary(mols.r[i]);
//				mols.r0[i] = mols.r[i] + vect_to_wall;
//				mols.r[i] += 2.0*vect_to_wall;
//				//mols.saved_index[i] = SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE;
//			}
		}
	}

}


template<typename T>
void DestroyBoundary<T>::integrate(const double dt) {

	BOOST_FOREACH(Species *s, this->get_species()) {
		Molecules& mols = s->mols;
		for (int i = 0; i < mols.size(); ++i) {
			if (this->geometry.lineXsurface(mols.r0[i],mols.r[i])) {
				mols.delete_molecule(i);
			}
//			if (Boundary<T>::geometry.distance_to_boundary(mols.r[i])<0) {
//				mols.delete_molecule(i);
//			}
		}
	}

}

template<typename T>
void CouplingBoundary_M_to_C<T>::integrate(const double dt) {

	DiffusionCorrectedBoundary<T>::timestep_initialise(dt);
	const int s_n = this->get_species().size();
	std::set<int> dirty_indicies;
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->get_species()[s_i]);
		const int p_n = s.mols.size();
		for (int p_i = 0; p_i < p_n; ++p_i) {
			if (DiffusionCorrectedBoundary<T>::particle_crossed_boundary(p_i, s_i)) {
				const double dist_to_wall = (*(this->all_curr_distance[s_i]))[p_i];
				//if (dist_to_wall < 0) {
				//if (dist_to_wall < 0) {
				const Vect3d r = s.mols.r[p_i];
				Vect3d comp_r = r;
				if (dist_to_wall > 0) {
					comp_r += (1.0 + s.grid->get_tolerance())*this->geometry.shortest_vector_to_boundary(r);
				}
//				else {
//					comp_r += (1.0 - s.grid.get_tolerance())*this->geometry.shortest_vector_to_boundary(r);
//				}
				const int i = s.grid->get_cell_index(comp_r);
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

}

template<typename T>
void CouplingBoundary<T>::integrate(const double dt) {
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator,boost::uniform_real<>(0,1));

	const int s_n = this->get_species().size();
	std::set<int> dirty_indicies;
	for (int s_i = 0; s_i < s_n; ++s_i) {
		Species &s = *(this->get_species()[s_i]);
		const int p_n = s.mols.size();
		for (int p_i = 0; p_i < p_n; ++p_i) {
			const Vect3d r = s.mols.r[p_i];
			const Vect3d rold = s.mols.r0[p_i];
			if (this->geometry.is_in(r)) {
				const int i = s.grid->get_cell_index(r);
				ASSERT(i>=0, "Invalid negative compartment index!");
				dirty_indicies.insert(i);
				s.copy_numbers[i]++;
				s.mols.mark_for_deletion(p_i);
			} else if (corrected) {
				const double old_dist_to_wall = this->geometry.distance_to_boundary(rold);
				if (old_dist_to_wall==0.0) continue;
				const Vect3d vect_to_wall = this->geometry.shortest_vector_to_boundary(r);
				const double dist_to_wall = vect_to_wall.norm();
				//std::cout << "rold = "<<rold<<" old_dist_to_wall = "<<old_dist_to_wall<<" r = "<<r<<" dist_to_wall = "<<dist_to_wall<<std::endl;
				//TODO: assumes isotropic diffusion
				const double P = exp(-dist_to_wall*old_dist_to_wall/(s.D.maxCoeff()*dt));
				if (uni() < P) {
					const int i = s.grid->get_cell_index(r + 1.000001*vect_to_wall);
					ASSERT(i>=0, "Invalid negative compartment index!");
					dirty_indicies.insert(i);
					s.copy_numbers[i]++;
					s.mols.mark_for_deletion(p_i);
				}
			}
		}
		s.mols.delete_molecules();
		//std::cout << count << "particles moved to compartments. Free space = "<<s.mols.size() <<" compart = "<< std::accumulate(s.copy_numbers.begin(),s.copy_numbers.end(),0) << std::endl;
	}
	BOOST_FOREACH(int i, dirty_indicies) {
		nsm->recalc_priority(i);
	}
}



}

#endif /* BOUNDARY_IMPL_H_ */
