/*
 * Grid.h
 *
 *  Created on: 1 May 2014
 *      Author: Ulrich Dobramysl
 */

#ifndef GRID_H_
#define GRID_H_

#include "Vector.h"
#include "Geometry.h"
#include <vector>

namespace Tyche {
class Grid {
public:
  virtual Vect3d get_random_point(const int i) const = 0;
  
  //virtual void get_slice(const xrect& surface, std::vector<int>& indices) const = 0;
  //virtual void get_slice(const yrect& surface, std::vector<int>& indices) const = 0;
  //virtual void get_slice(const zrect& surface, std::vector<int>& indices) const = 0;
  
  //virtual void get_slice(const xplane& surface, std::vector<int>& indices) const = 0;
  //virtual void get_slice(const yplane& surface, std::vector<int>& indices) const = 0;
  //virtual void get_slice(const zplane& surface, std::vector<int>& indices) const = 0;
  
  virtual void get_slice(const Geometry& geometry, std::vector<int>& indices) const = 0;
  
  virtual void get_region(const Geometry& geometry, std::vector<int>& indices) const = 0;
  
  virtual bool is_in(const Geometry& geometry, const int index) const = 0;
  virtual inline bool is_in(const Vect3d& r) const {
    return ((r.array() >= get_low().array()).all()) && ((r.array() < get_high().array()).all());
  }

  virtual void get_overlap(const Vect3d& low, const Vect3d& high, std::vector<int>& indicies, std::vector<double>& volume) const = 0;
  
  virtual const Vect3d& get_low() const = 0;
  virtual const Vect3d& get_high() const = 0;
  virtual int size() const = 0;
  virtual Vect3d get_cell_size() const = 0;
  virtual Vect3i get_cells_along_axes() const = 0;
  virtual double get_cell_volume(const int i) const = 0;
  virtual double get_tolerance() const = 0;
  virtual Vect3d get_cell_centre(const int i) const = 0;
  virtual Vect3d get_low_point(const int i) const = 0;
  virtual Vect3d get_high_point(const int i) const = 0;
  virtual const std::vector<int> get_neighbour_indicies(const int i) const = 0;
  virtual const std::vector<double> get_neighbour_distances(const int i) const = 0;
  
  virtual int get_cell_index(const Vect3d &r) const = 0;
  virtual Rectangle get_face_between(const int i, const int j) const = 0;
  virtual double get_laplace_coefficient(const int i, const int j) const = 0;
  virtual double get_distance_between(const int i, const int j) const = 0;
};
}

#endif /* GRID_H_ */
