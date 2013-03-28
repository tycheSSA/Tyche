/*
 * Union.h
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
 *  Created on: 16 Nov 2012
 *      Author: robinsonm
 */

#ifndef UNION_H_
#define UNION_H_

namespace Tyche {
template<typename T1, typename T2>
class Intersection {
public:
	Intersection(const T1& geometry1, const T2& geometry2):
	   geometry1(geometry1), geometry2(geometry2) {};

	typedef T1 GeometryType1;
	typedef T2 GeometryType2;

	typedef AxisAlignedRectangle<T1::dim> SurfaceElementType;


	bool at_boundary(const Vect3d& r) const {
		return geometry1.at_boundary(r) && geometry2.at_boundary(r);
	}
	bool at_boundary(const Vect3d& r, Vect3d& shortest) const {
		bool ret = geometry1.at_boundary(r,shortest);
		Vect3d shortest2;
		ret &= geometry2.at_boundary(r,shortest2);
		if (shortest.squaredNorm() > shortest2.squaredNorm()) {
			shortest = shortest2;
		}
		return ret;
	}

	const Vect3d shortest_vector_to_boundary(const Vect3d& r) const {
		Vect3d ret1 = geometry1.shortest_vector_to_boundary(r);
		Vect3d ret2 = geometry2.shortest_vector_to_boundary(r);
		if (ret1.squaredNorm() < ret2.squaredNorm()) {
			return ret1;
		} else {
			return ret2;
		}
	}
	double distance_to_boundary(const Vect3d& r) const {
		const double dist1 = geometry1.distance_to_boundary(r);
		const double dist2 = geometry2.distance_to_boundary(r);

		if (dist1 > 0) {
		   if (dist2 > 0) {
		      return std::min(dist1,dist2);
		   } else {
		      return dist1;
		   }
		} else {
		   if (dist2 > 0) {
		      return dist2;
		   } else {
		      return std::max(dist1,dist2);
		   }
		}
	}
	void operator +=(const double move_by) {
		geometry1 += move_by;
		geometry2 += move_by;
	}
	void operator -=(const double move_by) {
		geometry1 -= move_by;
		geometry2 -= move_by;
	}


	const T1& geometry1;
	const T2& geometry2;
};

template<typename T1, typename T2>
std::ostream& operator<< (std::ostream& out, const Intersection<T1,T2> &b) {
   return out << "intersection of " << b.geometry1 << " and " << b.geometry2;
}
}
//template<typename T1, typename T2>
//Intersection<T1,T2> make_union(const T1& geometry1, const T2& geometry2) {
//	return Intersection<T1, T2>(geometry1, geometry2);
//}

#endif /* VOLUME_H_ */
