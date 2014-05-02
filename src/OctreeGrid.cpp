#include "OctreeGrid.h"

namespace Tyche {

void Octree::refine(const Vect3d& point, const double target_edge) {
  if(!is_in(point))
    return;
  if (target_edge>=edge_length)
    return;
  if (is_leaf()) {
    double h = edge_length/2.;
    for (int i = 0; i < 8; i++) {
      Vect3i direction(2*((i>>0)&1)-1,2*((i>>1)&1)-1,2*((i>>2)&1)-1);
      Vect3d new_point = center+h/2.*direction.cast<double>();
      Octree *o = new Octree(new_point, h, level+1, this);
      subtrees[i] = o;
    }
  }
  for (int i = 0; i < 8; i++)
    subtrees[i]->refine(point, target_edge);
}

std::ostream& operator<< (std::ostream& out, const Octree& o) {
  out << "Octree Leaf " << o.cell_index << " at (" << o.center[0] << "," << o.center[1] << "," << o.center[2] << ") with edge length " << o.edge_length << " and neighbours ";
  for (auto ni : o.get_neighbour_indices())
    out << ni << ",";
  return out << std::endl;
}

Vect3d OctreeGrid::get_random_point(const int i) const {
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator,boost::uniform_real<>(0,1));
  return cells[i]->get_low() + cells[i]->get_edge_length()*Vect3d(uni(),uni(),uni());
}

void OctreeGrid::get_slice(const Geometry& geometry, std::vector<int>& indices) const {
  indices.clear();
  const int nedges = 14;
  const double edges[nedges][2][3] = {{{0,0,0},{0,0,1}},
				      {{0,0,0},{0,1,0}},
				      {{0,0,0},{1,0,0}},
				      {{0,0,1},{0,1,1}},
				      {{0,0,1},{1,0,0}},
				      {{0,1,0},{1,1,0}},
				      {{0,1,0},{0,1,1}},
				      {{1,0,0},{1,1,0}},
				      {{1,0,0},{1,0,1}},
				      {{0,1,1},{1,1,1}},
				      {{1,1,0},{1,1,1}},
				      {{1,0,1},{1,1,1}},
				      {{0,0,0},{0.5,0.5,0.5}},
				      {{0.5,0.5,0.5},{1,1,1}}};
  
  for (int i = 0; i < cells.size(); ++i) {
    const Vect3d low_point = cells[i]->get_low();
    const double h = cells[i]->get_edge_length();
    for (int j = 0; j < nedges; ++j) {
      const Vect3d p1 = low_point + h*Vect3d(edges[j][0][0],edges[j][0][1],edges[j][0][2]);
      const Vect3d p2 = low_point + h*Vect3d(edges[j][1][0],edges[j][1][1],edges[j][1][2]);
      if (geometry.lineXsurface(p1,p2)) {
	indices.push_back(i);
	break;
      }
    }
  }
}

void OctreeGrid::get_region(const Geometry& geometry, std::vector<int>& indices) const {
  indices.clear();
  for (int i = 0; i < cells.size(); ++i) {
    if (is_in(geometry,i)) {
      indices.push_back(i);
    }
  }
}

bool OctreeGrid::is_in(const Geometry& geometry, const int index) const {
  if (geometry.is_in(cells[index]->get_center())) return true;
  const Vect3d low_point = cells[index]->get_low();
  const double h = cells[index]->get_edge_length();
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
	const Vect3d test_point = low_point + h*Vect3d(i,j,k);
	if (geometry.is_in(test_point)) {
	  return true;
	}
      }
    }
  }
  return false;
}

void OctreeGrid::get_overlap(const Vect3d& overlap_low, const Vect3d& overlap_high, std::vector<int>& indicies, std::vector<double>& volume) const {
  indicies.clear();
  volume.clear();
  if ((overlap_low.array() >= high.array()).any()) return;
  if ((overlap_high.array() <= low.array()).any()) return;
  Vect3d snap_low = overlap_low + Vect3d(tolerance,tolerance,tolerance);
  Vect3d snap_high = overlap_high - Vect3d(tolerance,tolerance,tolerance);
  for (int i = 0; i < 3; ++i) {
    if (snap_low[i] < low[i]) {
      snap_low[i] = low[i];
    }
    if (snap_high[i] > high[i]) {
      snap_high[i] = high[i]-tolerance;
    }
  }
  
  Box overlap(snap_low, snap_high, true);
  get_region(overlap, indicies);
  for (auto j : indicies) {
    Vect3d low_point = cells[j]->get_low();
    Vect3d high_point = cells[j]->get_high();
    const double h = cells[j]->get_edge_length();
    const double inv_volume = 1./(h*h*h);
    for (int i = 0; i < 3; ++i) {
      if (high_point[i] > overlap_high[i]) {
	high_point[i] = overlap_high[i];
      }
      if (low_point[i] < overlap_low[i]) {
	low_point[i] = overlap_low[i];
      }
    }
    volume.push_back((high_point-low_point).prod()*inv_volume);
  }
}

Rectangle OctreeGrid::get_face_between(const int i, const int j) const {
  Vect3d p1 = cells[i]->get_center();
  Vect3d p2 = cells[j]->get_center();
  Vect3d low1 = cells[i]->get_low();
  Vect3d high1 = cells[i]->get_high();
  Vect3d low2 = cells[j]->get_low();
  Vect3d high2 = cells[j]->get_high();
  double h1 = cells[i]->get_edge_length();
  double h2 = cells[j]->get_edge_length();
  
  Vect3d normal = Vect3d::Zero();
  for (int d=0; d<3; d++) {
    if (p1[d]+h1/2==p2[d]-h2/2) {
      low1[d] += h1;
      high2[d] -= h2;
      normal[d] = 1.;
    } else if (p1[d]-h1/2==p2[d]+h2/2) {
      low2[d] += h2;
      high1[d] -= h1;
      normal[d] = -1.;
    } else {
      continue;
    }
    break;
  }
  
  Vect3d low,high;
  if ((low1.array() < low2.array()).any()) {
    low = low2;
  } else {
    low = low1;
  }
  if ((high1.array() < high2.array()).any()) {
    high = high1;
  } else {
    high = high2;
  }
  
  const Vect3d mid = (low+high)/2.;
  const Vect3d in_plane = (low-mid).cross(normal);
  const Vect3d low_right = mid+in_plane;
  const Vect3d up_left = mid-in_plane;
  return Rectangle(low, low_right, up_left);
}

double OctreeGrid::get_laplace_coefficient(const int i, const int j) const {
  const double hi = cells[i]->get_edge_length();
  const double hj = cells[j]->get_edge_length();
  const double hm = fmin(hi,hj);
  const Rectangle face = get_face_between(i,j);
  const Vect3d rij = cells[j]->get_center()-cells[i]->get_center();
  const double dij = rij.norm();
  return hm*hm*face.get_normal().dot(rij)/(dij*dij*hi*hi*hi);
}

void OctreeGrid::calculate_neighbours() {
  std::vector<Octree*> leafs;
  for (auto o : octrees)
    o->get_leafs(leafs);
  for (int i=0; i<leafs.size(); i++)
    leafs[i]->clear_neighbours();
  for (int i=0; i<leafs.size(); i++) {
    for (int j=i+1; j<leafs.size(); j++) {
      Octree *a = leafs[i];
      Octree *b = leafs[j];
      if (a->is_adjacent(*b)) {
	a->insert_neighbour(b);
	b->insert_neighbour(a);
      }
    }
  }
}


}
