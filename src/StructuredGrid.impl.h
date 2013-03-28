#ifndef STRUCTUREDGRID_IMPL_H_
#define STRUCTUREDGRID_IMPL_H_

namespace Tyche {
template<typename T>
bool StructuredGrid::geometry_intersects_cell(const int i, const T geometry) const {

	const Vect3d low_point = index_to_vect(i).cast<double>().cwiseProduct(cell_size) + low;

	bool at_least_one_zero = false;
	bool at_least_one_negative = false;
	bool at_least_one_positive = false;

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				const Vect3d test_point = low_point + Vect3d(i,j,k).cwiseProduct(cell_size);
				const double dist = geometry.distance_to_boundary(test_point);
				if (std::abs(dist) < tolerance) {
					at_least_one_zero = true;
					continue;
				}
				const bool sign = std::signbit(dist);
				if (sign == true) {
					at_least_one_negative = true;
				} else {
					at_least_one_positive = true;
				}
			}
		}
	}

	return (at_least_one_negative && at_least_one_positive) || (at_least_one_zero && at_least_one_positive);
}



template<unsigned int DIM>
AxisAlignedRectangle<DIM> StructuredGrid::get_intersection_of_cell(const int i,
		const Intersection<AxisAlignedPlane<DIM>,AxisAlignedPlane<DIM> >& geometry) const {
	const double dim_map[3][2] = {{1,2},{0,2},{0,1}};
	Vect3d low_point = index_to_vect(i).cast<double>().cwiseProduct(cell_size)+low;
	Vect3d high_point = low_point + cell_size;
	const double dist = geometry.distance_to_boundary(low_point);
	low_point[DIM] += std::abs(dist);
	high_point[DIM] = low_point[DIM];
	int norm;
	if (dist < 0) {
		norm = +1;
	} else {
		norm = -1;
	}
	return AxisAlignedRectangle<DIM>(low_point, high_point, norm);
}

template<unsigned int DIM>
AxisAlignedRectangle<DIM> StructuredGrid::get_intersection_of_cell(const int i,
		const AxisAlignedPlane<DIM>& geometry) const {
	const double dim_map[3][2] = {{1,2},{0,2},{0,1}};
	Vect3d low_point = index_to_vect(i).cast<double>().cwiseProduct(cell_size)+low;
	Vect3d high_point = low_point + cell_size;
	low_point[DIM] = geometry.get_coord();
	high_point[DIM] = low_point[DIM];
	const int norm = geometry.get_normal();
	return AxisAlignedRectangle<DIM>(low_point, high_point, norm);
}

template<unsigned int DIM>
void StructuredGrid::get_slice(const Intersection<AxisAlignedPlane<DIM>,AxisAlignedPlane<DIM> >& surface, std::vector<int>& indices) const {
	const int coord_index1 = (surface.geometry1.get_coord()-low[DIM])*inv_cell_size[DIM];
	const int coord_index2 = (surface.geometry2.get_coord()-low[DIM])*inv_cell_size[DIM];
	get_slice(surface.geometry1,indices);
	if (coord_index2 != coord_index1) {
		get_slice(surface.geometry2,indices);
	}
}
template<unsigned int DIM>
void StructuredGrid::get_slice(const AxisAlignedPlane<DIM>& surface, std::vector<int>& indices) const {
	const double coord_index_double = (surface.get_coord()-low[DIM])*inv_cell_size[DIM];
	const int coord_index_int = int(floor(coord_index_double + surface.get_normal()*tolerance));

	if ((coord_index_int < 0) || (coord_index_int > num_cells_along_axes[DIM]-1)) return;
	static const int dim_map[][2] = {{1,2}, {0,2}, {0,1}};
	int ret_index = indices.size();
	const int size = ret_index + num_cells_along_axes[dim_map[DIM][0]]*num_cells_along_axes[dim_map[DIM][1]];

	indices.resize(size);
	Vect3i vect_index(0,0,0);
	vect_index[DIM] = coord_index_int;

	for (int i = 0; i < num_cells_along_axes[dim_map[DIM][0]]; ++i) {
		vect_index[dim_map[DIM][0]] = i;
		for (int j = 0; j < num_cells_along_axes[dim_map[DIM][1]]; ++j) {
			vect_index[dim_map[DIM][1]] = j;
			const int index = vect_to_index(vect_index);
			ASSERT(ret_index < size,"return index is out of bounds");
			indices[ret_index] = index;
			ret_index++;
		}
	}
}
}

#endif /* STRUCTUREDGRID_IMPL_H_ */
