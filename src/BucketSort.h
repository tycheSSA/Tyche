/*
 * BucketSort.h
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
 *  Created on: 22 Oct 2012
 *      Author: robinsonm
 */

#ifndef BUCKETSORT_H_
#define BUCKETSORT_H_

#include "Vector.h"
#include "Constants.h"
#include "Log.h"
#include <vector>
#include <iostream>

namespace Tyche {

const int CELL_EMPTY = -1;

class BucketSort {
public:
	BucketSort(Vect3d low, Vect3d high, Vect3b periodic):
		low(low),high(high),domain_size(high-low),periodic(periodic) {
		LOG(2,"Creating bucketsort data structure with lower corner = "<<low<<" and upper corner = "<<high);
		const double dx = (high-low).maxCoeff()/10.0;
		reset(low, high, dx);
	}

	void reset(const Vect3d& low, const Vect3d& high, double _max_interaction_radius);
	inline const Vect3d& get_low() {return low;}
	inline const Vect3d& get_high() {return high;}

	void embed_points(std::vector<Vect3d> &positions);
	std::vector<int>& find_broadphase_neighbours(const Vect3d& r, const int my_index, const bool self);

	Vect3d correct_position_for_periodicity(const Vect3d& source_r, const Vect3d& to_correct_r);
	Vect3d correct_position_for_periodicity(const Vect3d& to_correct_r);

private:
	inline int vect_to_index(const Vect3i& vect) {
		return vect[0] * num_cells_along_yz + vect[1] * num_cells_along_axes[1] + vect[2];
	}
	inline int find_cell_index(const Vect3d &r) {
		const Vect3i celli = ((r-low).cwiseProduct(inv_cell_size) + Vect3d(1.0,1.0,1.0)).cast<int>();
		ASSERT((celli[0] >= 0) && (celli[0] < num_cells_along_axes[0]), "position is outside of x-range");
		ASSERT((celli[1] >= 0) && (celli[1] < num_cells_along_axes[1]), "position is outside of y-range");
		ASSERT((celli[2] >= 0) && (celli[2] < num_cells_along_axes[2]), "position is outside of z-range");
		//std::cout << " looking in cell " << celli <<" out of total cells " << num_cells_along_axes << " at position " << r<< std::endl;
		return vect_to_index(celli);
	}

    std::vector<int> cells;
    std::vector<std::vector<int> > ghosting_indices_pb;
    std::vector<std::pair<int,int> > ghosting_indices_cb;
    std::vector<int> dirty_cells;
	std::vector<int> linked_list;
	std::vector<int> neighbr_list;
	Vect3d low,high,domain_size;
	const Vect3b periodic;
	Vect3d cell_size,inv_cell_size;
	Vect3i num_cells_along_axes;
	int num_cells_along_yz;
	double max_interaction_radius;
	std::vector<int> surrounding_cell_offsets;
};

}

#endif /* BUCKETSORT_H_ */
