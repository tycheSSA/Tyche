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


template<unsigned int DIM>
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

	static std::auto_ptr<AxisAlignedPlane<DIM> > New(const double coord, const int normal) {
		return std::auto_ptr<AxisAlignedPlane<DIM> >(new AxisAlignedPlane<DIM>(coord,normal));
	}

	bool lineXsurface(const Vect3d& p1, const Vect3d& p2, Vect3d *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
		if (((p2[DIM]>=coord)&&(p1[coord]<coord))||((p2[DIM]<coord)&&(p1[coord]>=coord))) {
			if (intersect_point != NULL) {
				(*intersect_point)[DIM] = coord;
				(*intersect_point)[dim_map[DIM][0]] = 0.5*(p1[dim_map[DIM][0]] + p2[dim_map[DIM][0]]);
				(*intersect_point)[dim_map[DIM][1]] = 0.5*(p1[dim_map[DIM][1]] + p2[dim_map[DIM][1]]);
			}
			if (intersect_normal != NULL) {
				(*intersect_normal)[DIM] = 1.0;
				(*intersect_normal)[dim_map[DIM][0]] = 0.0;
				(*intersect_normal)[dim_map[DIM][1]] = 0.0;
			}
			return true;
		} else {
			return false;
		}
	}

	bool is_between(const AxisAlignedPlane<DIM>& plane1, const AxisAlignedPlane<DIM>& plane2) const {
		if (plane1.coord < plane2.coord) {
			return ((coord > plane1.coord) && (coord < plane2.coord));
		} else {
			return ((coord < plane1.coord) && (coord > plane2.coord));
		}
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

	const double& get_coord() const {
		return coord;
	}

	const int& get_normal() const {
		return normal;
	}

	void swap_normal() {
		normal = -normal;
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



template<unsigned int DIM>
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

	static std::auto_ptr<AxisAlignedRectangle<DIM> > New(const Vect3d _low, const Vect3d _high, const int _normal) {
		return std::auto_ptr<AxisAlignedRectangle<DIM> >(new AxisAlignedRectangle<DIM>(_low,_high,_normal));
	}

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
		return *this;
	}

	bool lineXsurface(const Vect3d& p1, const Vect3d& p2, Vect3d *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
		if (((p2[DIM]>=low[DIM])&&(p1[DIM]<low[DIM]))||((p2[DIM]<low[DIM])&&(p1[DIM]>=low[DIM]))) {
			const double intersect0 = 0.5*(p1[dim_map[DIM][0]] + p2[dim_map[DIM][0]]);
			if ((intersect0 >= low[dim_map[DIM][0]]) && (intersect0 < high[dim_map[DIM][0]])) {
				const double intersect1 = 0.5*(p1[dim_map[DIM][1]] + p2[dim_map[DIM][1]]);
				if ((intersect1 >= low[dim_map[DIM][1]]) && (intersect1 < high[dim_map[DIM][1]])) {
					if (intersect_point != NULL) {
						(*intersect_point)[DIM] = low[DIM];
						(*intersect_point)[dim_map[DIM][0]] = intersect0;
						(*intersect_point)[dim_map[DIM][1]] = intersect1;
					}
					if (intersect_normal != NULL) {
						(*intersect_normal)[DIM] = 1.0;
						(*intersect_normal)[dim_map[DIM][0]] = 0.0;
						(*intersect_normal)[dim_map[DIM][1]] = 0.0;
					}
					return true;
				}
			}
		}
		return false;
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

	const Vect3d& get_low() const {return low;}
	const Vect3d& get_high() const {return high;}

private:
	Vect3d low,high,normal_vector;

	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni1, uni2;
	boost::variate_generator<base_generator_type&, boost::triangle_distribution<> > tri1, tri2;
};

typedef AxisAlignedRectangle<0> xrect;
typedef AxisAlignedRectangle<1> yrect;
typedef AxisAlignedRectangle<2> zrect;

template<unsigned int DIM>
std::ostream& operator<< (std::ostream& out, const AxisAlignedRectangle<DIM>& p) {
	return out << "Rectangle aligned with dimension " << DIM << ". Lower point in other dimensions is "<<p.get_low()<<". Upper point in other dimensions is "<<p.get_high()<<".";
}

class Rectangle {
public:
	Rectangle(const Vect3d& lower_corner,
			const Vect3d& upper_left_corner,
			const Vect3d& lower_right_corner):
				low(lower_corner),high(upper_left_corner) {
		l = upper_left_corner-lower_corner;
		r = lower_right_corner-lower_corner;
		normal = l.cross(r);
		normal.normalize();
	}

	bool lineXsurface(const Vect3d& p1, const Vect3d& p2, Vect3d *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
		for (int d = 0; d < 3; ++d) {
			if (((p1[d]<low[d]) && (p2[d]<low[d])) || ((p1[d]>=high[d]) && (p2[d]>=high[d]))) {
				return false;
			}
		}
		const double denominator = (p2-p1).dot(normal);
		if (denominator==0) return false;
		const double numerator = (low-p1).dot(normal);
		if (intersect_point==NULL) {
			intersect_point = new Vect3d();
		}
		*intersect_point = (numerator/denominator) * (p2-p1) + p1;
		for (int d = 0; d < 3; ++d) {
			if (((*intersect_point)[d]>=high[d]) && ((*intersect_point)[d]<low[d])) {
				return false;
			}
		}
		if (intersect_normal != NULL) {
			*intersect_normal = normal;
		}
		return true;
	}


	void get_random_point_and_normal(Vect3d& p, Vect3d& n) {
		p = get_random_point();
		n = normal;
	}
	Vect3d get_random_point() {
		boost::variate_generator<base_generator_type&, boost::uniform_real<> >
			uni(generator,boost::uniform_real<>(0,1));
		return low + l*uni() + r*uni();
	}

	void get_random_point_and_normal_triangle(Vect3d& p, Vect3d& n) {
		p = get_random_point_triangle();
		n = normal;
	}
	Vect3d get_random_point_triangle() {
		boost::variate_generator<base_generator_type&, boost::triangle_distribution<> >
					tri(generator,boost::triangle_distribution<>(0,0.5,1));
		return low + l*tri() + r*tri();
	}
	const Vect3d& get_low() const {return low;}
	const Vect3d& get_l() const {return l;}
	const Vect3d& get_r() const {return r;}
	const Vect3d& get_normal() const {return normal;}
private:
	Vect3d low,high;
	Vect3d l,r;
	Vect3d normal;
};

std::ostream& operator<< (std::ostream& out, const Rectangle& p);

class Box {
public:
	Box(const Vect3d& lower_corner,
			const Vect3d& upper_corner,
			const bool in):low(lower_corner),high(upper_corner),in(in) {
	}
	bool is_in(const Vect3d& point) {
		const bool inside = ((point.array() >= low.array()).all() && (point.array() < high.array()).all());
		return inside == in;
	}
private:
	Vect3d low,high;
    bool in;
};

std::ostream& operator<< (std::ostream& out, const Box& p);

}

#endif /* GEOMETRY_H_ */
