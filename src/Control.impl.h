/*
 * Control.h
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
 *  Created on: 4 Dec 2012
 *      Author: robinsonm
 */

#ifndef CONTROL_IMPL_H_
#define CONTROL_IMPL_H_

namespace Tyche {

template<typename T>
void GrowingInterface<T>::integrate(const double dt) {

	//Operator::operator()();
	if (last_check > check_dt) {
		last_check = 0;
		ASSERT(this->get_species().size() > 0,"no species!!");
		//TODO: assumes only one species;
		Species &s = *(this->get_species()[0]);

		std::vector<int> grid_indices_current, grid_indices_shrink;
		s.grid.get_slice(this->geometry,grid_indices_current);
		this->geometry -= move_by;
		if (this->geometry.is_between(shrink_to,grow_to)) {
			s.grid.get_slice(this->geometry,grid_indices_shrink);

			const int total_copy_number_shrink = std::accumulate(grid_indices_shrink.begin(),grid_indices_shrink.end(),0,my_accumulate(s.copy_numbers));
			const bool is_meaningful_shrink = (grid_indices_shrink.size() > 0) || (grid_indices_current.size() > 0);
			if ((is_meaningful_shrink)&&(total_copy_number_shrink < shrink_threshold)) {
				LOG(1.5,"Shrinking interface at " << this->geometry <<". found " << total_copy_number_shrink <<" molecules behind interface");
				nsm.unset_interface_reactions(s,grid_indices_current);
				/*
				 * If there are any molecules in the compartments, randomly
				 * pick a position for them and either convert them to particles
				 * (if outside interface) or keep them in compartments
				 */
				BOOST_FOREACH(int i,grid_indices_shrink) {
					const int n = s.copy_numbers[i];
					for (int ii = 0; ii < n; ++ii) {
						const Vect3d r = s.grid.get_random_point(i);
						if (this->geometry.distance_to_boundary(r) >= 0) {
							s.mols.add_molecule(r);
							s.copy_numbers[i]--;
							ASSERT(s.copy_numbers[i]>=0,"can't have negative copy numbers");
						}
					}
				}
				nsm.set_interface_reactions(s,grid_indices_shrink,dt);
				//					OutputMolecularConcentrations out_mol("simpleReact_mi_after_shrink_mols_",0,0.0,1.0,100);out_mol.add_species(s);
				//					OutputCompartmentConcentrations out_compart("simpleReact_mi_after_shrink_compart_",0);out_compart.add_species(s);
				//					out_mol(this->time);out_compart(this->time);
				return;
			}
		}

		this->geometry += 2.0*move_by;

		if (this->geometry.is_between(shrink_to,grow_to)) {
			//std::fill(s.mol_copy_numbers.begin(),s.mol_copy_numbers.end(),0);
			/*
			 * count and store indicies of all molecules behind new geometry
			 */
			std::vector<int> mol_indices;
			const int n = s.mols.size();
			for (int i = 0; i < n; ++i) {
				const double dist = this->geometry.distance_to_boundary(s.mols.r[i]);
				if ((dist < 0)&&(dist > -move_by)) {
					mol_indices.push_back(i);
				}
			}


			std::vector<int> grid_indices_grow;
			s.grid.get_slice(this->geometry,grid_indices_grow);
			const int total_copy_number_grow = mol_indices.size();
			//const int total_copy_number_grow = std::accumulate(grid_indices_grow.begin(),grid_indices_grow.end(),0,my_accumulate(s.mol_copy_numbers));
			const bool is_meaningful_grow = (grid_indices_grow.size() > 0) || (grid_indices_current.size() > 0);
			if ((is_meaningful_grow)&&(total_copy_number_grow > grow_threshold)) {
				LOG(1.5,"Growing interface at " << this->geometry <<". found " << total_copy_number_grow <<" molecules in front of interface");

				/*
				 * delete particles behind interface and add them to compartments
				 */
				BOOST_FOREACH(int i, mol_indices) {
					s.copy_numbers[s.grid.get_cell_index(s.mols.r[i])]++;
					s.mols.delete_molecule(i);
				}
				nsm.unset_interface_reactions(s, grid_indices_current);
				nsm.set_interface_reactions(s,grid_indices_grow,dt);
				//				nsm.clear_reactions(grid_indices_grow);
				//				nsm.copy_reactions(grid_indices_shrink[0],)

				//					OutputMolecularConcentrations out_mol("simpleReact_mi_after_growth_mols_",0,0.0,1.0,100);out_mol.add_species(s);
				//					OutputCompartmentConcentrations out_compart("simpleReact_mi_after_growth_compart_",0);out_compart.add_species(s);
				//					out_mol(this->time);out_compart(this->time);

				return;
			}
		}
		//return to normal
		this->geometry -= move_by;
	} else {
		last_check += dt;
	}

}

}

#endif /* CONTROL_IMPL_H_ */
