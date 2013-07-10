/*
 * Geometry.h
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
 *  Created on: 19 Oct 2012
 *      Author: robinsonm
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "Vector.h"
#include "MyRandom.h"

namespace Tyche {

class NullGeometry {
public:
	NullGeometry() {}
};

class Cuboid {
public:
	Cuboid(const Vect3d min, const Vect3d max);
	bool at_boundary(const Vect3d r) const;
private:
	Vect3d min,max;

};

template<int DIM>
class AxisAlignedRectangle;

static const int dim_map[][2] = {{1,2}, {0,2}, {0,1}};

template<unsigned int DIM>
class AxisAlignedPlane: public NullGeometry {
public:
	typedef AxisAlignedRectangle<DIM> SurfaceElementType;
	static const int dim = DIM;

	AxisAlignedPlane(const double coord, const int normal):
		coord(coord),
		normal(normal)
	{}
	AxisAlignedPlane(const double coord):
			coord(coord),
			normal(1.0)
		{}
	AxisAlignedPlane(const AxisAlignedPlane& arg):
			coord(arg.coord),
			normal(arg.normal)
	{}

	bool is_between(const AxisAlignedPlane<DIM>& plane1, const AxisAlignedPlane<DIM>& plane2) const {
		if (plane1.coord < plane2.coord) {
			return ((coord > plane1.coord) && (coord < plane2.coord));
		} else {
			return ((coord < plane1.coord) && (coord > plane2.coord));
		}
	}
	inline bool at_boundary(const Vect3d& r) const {
		return normal*(r[DIM]-coord) < 0;
	}
	inline bool at_boundary(const Vect3d& r, Vect3d& shortest) const {
		//shortest = Vect3d::Zero();
		shortest[dim_map[DIM][0]] = 0;
		shortest[dim_map[DIM][1]] = 0;
		shortest[DIM] = coord-r[DIM];
		return normal*(-shortest[DIM]) < 0;
	}

	const Vect3d shortest_vector_to_boundary(const Vect3d& r) const {
	   Vect3d shortest = Vect3d::Zero();
		shortest[DIM] = coord-r[DIM];
		return shortest;
	}
	inline double distance_to_boundary(const Vect3d& r) const {
		return normal*(r[DIM]-coord);
	}
	const Vect3d operator-(const AxisAlignedPlane& arg) const {
	   Vect3d diff = Vect3d::Zero();
		diff[DIM] = coord - arg.coord;
		return diff;
	}
	void operator +=(const double move_by) {
		coord += normal*move_by;
	}
	void operator -=(const double move_by) {
		coord -= normal*move_by;
	}

	void set_coord(const double arg) {
		coord = arg;
	}

	const double& get_coord() const {
		return coord;
	}

	const int& get_normal() const {
		return normal;
	}

protected:
	double coord;
	int normal;

};


template<unsigned int DIM>
std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<DIM>& p) {
	return out << "Plane aligned with dimension " << DIM;
}



std::ostream& operator<< (std::ostream& out, const NullGeometry& b);
std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<0>& p);
std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<1>& p);
std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<2>& p);

typedef AxisAlignedPlane<0> xplane;
typedef AxisAlignedPlane<1> yplane;
typedef AxisAlignedPlane<2> zplane;



template<int DIM>
class AxisAlignedRectangle: public AxisAlignedPlane<DIM> {
public:
	AxisAlignedRectangle(const Vect3d _low, const Vect3d _high, const int _normal):
	   AxisAlignedPlane<DIM>(_low[DIM], _normal),
	   low(_low),high(_high),
	   normal_vector(Vect3d::Zero()),
	   uni1(generator,boost::uniform_real<>(_low[dim_map[DIM][0]],_high[dim_map[DIM][0]])),
	   uni2(generator,boost::uniform_real<>(_low[dim_map[DIM][1]],_high[dim_map[DIM][1]])),
	   tri1(generator,boost::triangle_distribution<>(_low[dim_map[DIM][0]],
			   	   	   	   	   0.5*(_low[dim_map[DIM][0]] + _high[dim_map[DIM][0]]),
			   	   	   	   	   	   	   	             _high[dim_map[DIM][0]] )),
	   tri2(generator,boost::triangle_distribution<>(_low[dim_map[DIM][1]],
			   0.5*(_low[dim_map[DIM][1]] + _high[dim_map[DIM][1]]),
			   	   	   	   	   				         _high[dim_map[DIM][1]] )) {
	   high[DIM] = low[DIM];
	   normal_vector[DIM] = this->normal;
	}
	AxisAlignedRectangle(const AxisAlignedRectangle<DIM>& arg):
		AxisAlignedPlane<DIM>(arg),
		low(arg.low),
		high(arg.high),
		normal_vector(arg.normal_vector),
		uni1(generator,boost::uniform_real<>(arg.low[dim_map[DIM][0]],arg.high[dim_map[DIM][0]])),
		uni2(generator,boost::uniform_real<>(arg.low[dim_map[DIM][1]],arg.high[dim_map[DIM][1]])),
		tri1(generator,boost::triangle_distribution<>(arg.low[dim_map[DIM][0]],
				0.5*(arg.low[dim_map[DIM][0]] + arg.high[dim_map[DIM][0]]),
				arg.high[dim_map[DIM][0]] )),
				tri2(generator,boost::triangle_distribution<>(arg.low[dim_map[DIM][1]],
						0.5*(arg.low[dim_map[DIM][1]] + arg.high[dim_map[DIM][1]]),
						arg.high[dim_map[DIM][1]] ))
						{}
	AxisAlignedRectangle<DIM>& operator=(const AxisAlignedRectangle<DIM>& arg) {
		AxisAlignedPlane<DIM>::operator=(arg);
		low = arg.low;
		high = arg.high;
		normal_vector = arg.normal_vector;
		boost::uniform_real<double>& dist1 = uni1.distribution();
		boost::uniform_real<double>& dist2 = uni2.distribution();
		dist1.param(boost::uniform_real<double>::param_type(arg.low[dim_map[DIM][0]],arg.high[dim_map[DIM][0]]));
		dist2.param(boost::uniform_real<double>::param_type(arg.low[dim_map[DIM][1]],arg.high[dim_map[DIM][1]]));
		boost::triangle_distribution<double>& dist3 = tri1.distribution();
		boost::triangle_distribution<double>& dist4 = tri2.distribution();
		dist3.param(boost::triangle_distribution<double>::param_type(arg.low[dim_map[DIM][0]],
				0.5*(arg.low[dim_map[DIM][0]] + arg.high[dim_map[DIM][0]]),
				arg.high[dim_map[DIM][0]] ));
		dist4.param(boost::triangle_distribution<double>::param_type(arg.low[dim_map[DIM][1]],
				0.5*(arg.low[dim_map[DIM][1]] + arg.high[dim_map[DIM][1]]),
				arg.high[dim_map[DIM][1]] ));
	}
	void get_random_point_and_normal(Vect3d& p, Vect3d& n) {
	   p = get_random_point();
	   n = normal_vector;
	}
	Vect3d get_random_point() {
		Vect3d ret;
		ret[DIM] = this->coord;
		ret[dim_map[DIM][0]] = uni1();
		ret[dim_map[DIM][1]] = uni2();
		return ret;
	}

	void get_random_point_and_normal_triangle(Vect3d& p, Vect3d& n) {
		p = get_random_point_triangle();
		n = normal_vector;
	}
	Vect3d get_random_point_triangle() {
		Vect3d ret;
		ret[DIM] = this->coord;
		ret[dim_map[DIM][0]] = tri1();
		ret[dim_map[DIM][1]] = tri2();
		return ret;
	}

private:
	Vect3d low,high,normal_vector;

	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni1, uni2;
	boost::variate_generator<base_generator_type&, boost::triangle_distribution<> > tri1, tri2;
};

template<unsigned int DIM>
std::ostream& operator<< (std::ostream& out, const AxisAlignedRectangle<DIM>& p) {
	return out << "Rectangle aligned with dimension " << DIM << ". Lower point in other dimensions is "<<p.low<<". Upper point in other dimensions is "<<p.high<<".";
}

}

#endif /* GEOMETRY_H_ */
