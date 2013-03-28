/*
 * Geometry.cpp
 *
 *  Created on: 19 Oct 2012
 *      Author: robinsonm
 */

#include "Geometry.h"
#include "Constants.h"

namespace Tyche {


Cuboid::Cuboid(const Vect3d min, const Vect3d max):min(min),max(max) {

}

bool Cuboid::at_boundary(const Vect3d r) const {
	bool in = true;
	for (int d = 0; d < NDIM; ++d) {
		in &= ((r[d] >= min[d]) && (r[d] <= max[d]));
	}
	return in;
}


std::ostream& operator<< (std::ostream& out, const NullGeometry& b) {
	return out << "Null Geometry";
}

std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<0>& p) {
	return out << "x = " << p.get_coord() << " with normal " << p.get_normal();
}

std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<1>& p) {
	return out << "y = " << p.get_coord() << " with normal " << p.get_normal();
}

std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<2>& p) {
	return out << "z = " << p.get_coord() << " with normal " << p.get_normal();
}
}
