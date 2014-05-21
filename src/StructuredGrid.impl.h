#ifndef STRUCTUREDGRID_IMPL_H_
#define STRUCTUREDGRID_IMPL_H_

namespace Tyche {

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
AxisAlignedRectangle<DIM> StructuredGrid::get_intersection_of_cell(const int i,
		const AxisAlignedRectangle<DIM>& geometry) const {
	const double dim_map[3][2] = {{1,2},{0,2},{0,1}};
	Vect3d low_point = index_to_vect(i).cast<double>().cwiseProduct(cell_size)+low;
	Vect3d high_point = low_point + cell_size;
	low_point[DIM] = geometry.get_coord();
	high_point[DIM] = low_point[DIM];
	const int norm = geometry.get_normal();
	return AxisAlignedRectangle<DIM>(low_point, high_point, norm);
}

template<unsigned int DIM>
void StructuredGrid::get_slice(const AxisAlignedRectangle<DIM>& surface, std::vector<int>& indices) const {
	const double coord_index_double = (surface.get_coord()-low[DIM])*inv_cell_size[DIM];
	const int coord_index_int = int(floor(coord_index_double + surface.get_normal()*tolerance));
	if ((coord_index_int < 0) || (coord_index_int > num_cells_along_axes[DIM]-1)) return;
	static const int dim_map[][2] = {{1,2}, {0,2}, {0,1}};

	Vect3d low = surface.get_low();
	low[dim_map[DIM][0]] += tolerance;
	low[dim_map[DIM][1]] += tolerance;
	Vect3i min_index = get_cell_index_vector(low);

	Vect3d high = surface.get_high();
	high[dim_map[DIM][0]] -= tolerance;
	high[dim_map[DIM][1]] -= tolerance;
	Vect3i max_index = get_cell_index_vector(high);

	ASSERT(max_index[DIM]==min_index[DIM],"surface not an axis aligned rectangle");

	int num_indicies = (max_index[dim_map[DIM][0]] - min_index[dim_map[DIM][0]] + 1) *
			(max_index[dim_map[DIM][1]] - min_index[dim_map[DIM][1]] + 1);

	const int size = indices.size() + num_indicies;
	int ret_index = indices.size();

	indices.resize(size);
	Vect3i vect_index(0,0,0);
	vect_index[DIM] = coord_index_int;
	for (int i = min_index[dim_map[DIM][0]]; i <= max_index[dim_map[DIM][0]]; ++i) {
		vect_index[dim_map[DIM][0]] = i;
		for (int j = min_index[dim_map[DIM][1]]; j <= max_index[dim_map[DIM][1]]; ++j) {
			vect_index[dim_map[DIM][1]] = j;
			const int index = vect_to_index(vect_index);
			ASSERT(ret_index < size,"return index is out of bounds");
			indices[ret_index] = index;
			ret_index++;
		}
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
