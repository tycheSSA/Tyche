/*
 * OctreeGrid.h
 *
 *  Created on: 1 May 2014
 *      Author: Ulrich Dobramysl
 */

#ifndef OCTREEGRID_H_
#define OCTREEGRID_H_

#include "Grid.h"
#include <set>

namespace Tyche {

class Octree {
public:
  Octree(Vect3d center, double edge_length, unsigned int level, Octree *parent) :
    center(center),
    level(level),
    edge_length(edge_length),
    parent(parent)
  {
    low = center.array()-edge_length/2.;
    high = center.array()+edge_length/2.;
    cell_index = -1;

    for (int i = 0; i < 8; i++) {
      subtrees[i] = NULL;
    }
  }
  Octree(Vect3d center, double edge_length) : 
    center(center),
    level(0),
    edge_length(edge_length),
    parent(NULL)
  {
    low = center.array()-edge_length/2.;
    high = center.array()+edge_length/2.;
    cell_index = -1;

    for (int i = 0; i < 8; i++) {
      subtrees[i] = NULL;
    }
  }

  ~Octree() {
    for (int i = 0; i < 8; i++) {
      delete subtrees[i];
    }
  }
  
  inline bool is_leaf() const {
    return subtrees[0]==NULL;
  }

  inline bool is_in(const Vect3d& point) const {
    return ((point.array() > low.array()).all() && (point.array() < high.array()).all());
  }

  inline bool is_adjacent(const Octree& other) const {
    const Vect3d olow = other.get_low();
    const Vect3d ohigh = other.get_high();
    if ((high.array() < olow.array()).any())
      return false;
    if ((low.array() > ohigh.array()).any())
      return false;
    int cond_sum = (high.array() == olow.array()).cast<int>().sum() + (low.array() == ohigh.array()).cast<int>().sum();
    if (cond_sum >= 2)
      return false;
    return true;
  }

  void refine(const Vect3d& point, const double target_edge);
  
  int get_cell_index(const Vect3d& point) const {
    if (!is_in(point))
      return -1;
    if (is_leaf())
      return cell_index;
    for (int i = 0; i < 8; i++) {
      if (subtrees[i]->is_in(point)) {
	return subtrees[i]->get_cell_index(point);
      }
    }
    ERROR("Found no matching cell for point " << point << "!");
    return -1;
  }

  int get_cell_index() const {
    return cell_index;
  }

  void set_cell_index(const int index) {
    cell_index = index;
  }

  void get_leafs(std::vector<Octree*> &leafs) {
    if (is_leaf()) {
      leafs.push_back(this);
      return;
    }
    for(int i=0; i<8; i++) {
      subtrees[i]->get_leafs(leafs);
    }
  }

  inline Vect3d get_center() const {
    return center;
  }
  inline Vect3d get_low() const {
    return low;
  }
  inline Vect3d get_high() const {
    return high;
  }

  inline double get_edge_length() const {
    return edge_length;
  }

  void insert_neighbour(Octree *n) {
    neighbours.insert(n);
  }

  void clear_neighbours() {
    neighbours.clear();
  }

  std::vector<int> get_neighbour_indices() const {
    std::vector<int> ret;
    for (auto n : neighbours) {
      ret.push_back(n->get_cell_index());
    }
    return ret;
  }

private:
  friend std::ostream& operator<< (std::ostream& out, const Octree& o);
  Vect3d center,low,high;
  double edge_length;
  unsigned int level;
  unsigned int cell_index;
  Octree *subtrees[8];
  Octree *parent;
  std::set<Octree*> neighbours;
};

std::ostream& operator<< (std::ostream& out, const Octree& o);

class OctreeGrid : public Grid {
public:
  OctreeGrid(const Vect3d& center, const double domain_size, const unsigned int num_trees_along_axes) : 
    low(center.array()-domain_size/2.),
    high(center.array()+domain_size/2.),
    tree_size(domain_size/num_trees_along_axes),
    num_trees_along_axes(num_trees_along_axes)
  {
    for (int i = 0; i < num_trees_along_axes; i++) {
      for (int j = 0; j < num_trees_along_axes; j++) {
	for (int k = 0; k < num_trees_along_axes; k++) {
	  Vect3d position = low + tree_size*Vect3d(i+0.5,j+0.5,k+0.5);
	  octrees.push_back(new Octree(position, tree_size));
	}
      }
    }
  }

  ~OctreeGrid() {
    for (auto o : octrees) {
      delete o;
    }
  }

  virtual Vect3d get_random_point(const int i) const;
  
  //virtual void get_slice(const xrect& surface, std::vector<int>& indices) const = 0;
  //virtual void get_slice(const yrect& surface, std::vector<int>& indices) const = 0;
  //virtual void get_slice(const zrect& surface, std::vector<int>& indices) const = 0;
  
  //virtual void get_slice(const xplane& surface, std::vector<int>& indices) const = 0;
  //virtual void get_slice(const yplane& surface, std::vector<int>& indices) const = 0;
  //virtual void get_slice(const zplane& surface, std::vector<int>& indices) const = 0;
  
  virtual void get_slice(const Geometry& geometry, std::vector<int>& indices) const;
  
  virtual void get_region(const Geometry& geometry, std::vector<int>& indices) const;
  
  virtual inline bool is_in(const Vect3d& r) const {
    return ((r.array() >= get_low().array()).all()) && ((r.array() < get_high().array()).all());
  }
  virtual bool is_in(const Geometry& geometry, const int index) const;

  virtual void get_overlap(const Vect3d& overlap_low, const Vect3d& overlap_high, std::vector<int>& indicies, std::vector<double>& volume) const;
  
  virtual inline const Vect3d& get_low() const {
    return low;
  }
  virtual inline const Vect3d& get_high() const {
    return high;
  }
  virtual inline int size() const {
    return cells.size();
  }
  virtual inline Vect3d get_cell_size() const {
    return tree_size*Vect3d::Ones();
  }
  virtual inline double get_cell_volume(const int i) const {
    double h = cells[i]->get_edge_length();
    return h*h*h;
  }
  virtual double get_tolerance() const {
    return tolerance;
  }
  virtual Vect3d get_cell_centre(const int i) const {
    return cells[i]->get_center();
  }
  virtual Vect3d get_low_point(const int i) const {
    return cells[i]->get_low();
  }
  virtual Vect3d get_high_point(const int i) const {
    return cells[i]->get_high();
  }
  Vect3i get_cells_along_axes() const {
    return num_trees_along_axes*Vect3i::Ones();
  }
  virtual const std::vector<int> get_neighbour_indicies(const int i) const {
    return cells[i]->get_neighbour_indices();
  }
  virtual const std::vector<double> get_neighbour_distances(const int i) const {
    std::vector<double> ret;
    Vect3d p1 = cells[i]->get_center();
    for (auto n : cells[i]->get_neighbour_indices())
      ret.push_back(get_distance_between(i,n));
    return ret;
  }
  
  virtual inline int get_cell_index(const Vect3d &point) const {
    if (!is_in(point))
      return -1;
    return octrees[get_tree_index(point)]->get_cell_index(point);
  }

  inline std::vector<int> get_leaf_indices(const Vect3i tree_vec) const {
    std::vector<int> ret;
    std::vector<Octree*> leafs;
    octrees[get_tree_index(tree_vec)]->get_leafs(leafs);
    for(auto l : leafs)
      ret.push_back(l->get_cell_index());
    return ret;
  }

  virtual Rectangle get_face_between(const int i, const int j) const;

  virtual double get_laplace_coefficient(const int i, const int j) const;

  virtual double get_distance_between(const int i, const int j) const {
    return (cells[j]->get_center()-cells[i]->get_center()).norm();
  }

  void refine(const Vect3d& point, const double target_edge) {
    if (!is_in(point))
      return;
    Octree *o = octrees[get_tree_index(point)];
    o->refine(point, target_edge);
  }

  void finalise() {
    cells.clear();
    for (auto o : octrees) {
      o->get_leafs(cells);
    }
    double min_edge = 1e10;
    for (int i = 0; i<cells.size(); i++) {
      cells[i]->set_cell_index(i);
      const double h = cells[i]->get_edge_length();
      if (min_edge>h)
	min_edge = h;
    }
    tolerance = min_edge/100000.0;
    calculate_neighbours();
  }

  std::vector<Octree*> get_cells() const {
    return cells;
  }

private:
  inline int get_tree_index(const Vect3d &point) const {
    Vect3i indices = ((point-low)/tree_size).cast<int>();
    return get_tree_index(indices);
  }
  inline int get_tree_index(const Vect3i &indices) const {
    return ((indices[0]*num_trees_along_axes+indices[1])*num_trees_along_axes)+indices[2];
  }

  void calculate_neighbours();

  Vect3d low,high;
  int num_trees_along_axes;
  double tree_size;
  std::vector<Octree*> octrees;
  std::vector<Octree*> cells;
  double tolerance;
};
}

#endif /* OCTREEGRID_H_ */
