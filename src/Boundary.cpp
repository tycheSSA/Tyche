/*
 * Boundary.cpp
 *
 *  Created on: 18 Oct 2012
 *      Author: robinsonm
 */


#include "Boundary.h"

namespace Tyche {

void FluxBoundary::integrate(const double dt) {

	BOOST_FOREACH(Species *s, get_species()) {
		Molecules& mols = s->mols;
		boost::poisson_distribution<> p_dist(dt*rate);
		boost::variate_generator<base_generator_type&, boost::poisson_distribution<> > poisson(generator, p_dist);
		const unsigned int n = poisson();
		for (int i = 0; i < n; ++i) {
			mols.add_molecule(p + uni1()*t1 + uni2()*t2);
		}
	}
}


void BoundaryController::add_jump_boundary(Geometry &geometry,  const Vect3d jump_by) {
	param_type params;
	params[0] = JUMP;
	params[1] = jump_by[0];
	params[2] = jump_by[1];
	params[3] = jump_by[2];
	geometries.push_back(data_type(&geometry,params));
}

void BoundaryController::add_reflect_boundary(Geometry &geometry) {
	param_type params;
	params[0] = REFLECT;
	geometries.push_back(data_type(&geometry,params));
}

void BoundaryController::add_absorb_boundary(Geometry &geometry) {
	param_type params;
	params[0] = ABSORB;
	geometries.push_back(data_type(&geometry,params));
}


void BoundaryController::integrate(const double dt) {
	const int n = get_species().size();
	for (int s_i = 0; s_i < n; ++s_i) {
		Species& s = *(get_species()[s_i]);
		bool delete_molecules = false;
		Molecules& mols = s.mols;
		const int n = mols.size();
		for (int i = 0; i < n; ++i) {
			Vect3d nv,nv_min;
			double ip,ip_min;
			int iteration = 0;
			while (iteration++<3) {
				int min = -1;
				int index = 0;
				for (auto g: geometries) {
					if (g.first->lineXsurface(mols.r0[i],mols.r[i],&ip,&nv) && ((min == -1) || (ip < ip_min))) {
						ip_min = ip;
						nv_min = nv;
						min = index;
					}
					index++;
				}
				if (min != -1) {
					Vect3d v = ip_min*(mols.r[i]-mols.r0[i]);
					const Vect3d impact_point = mols.r0[i] + v;

					switch((int)geometries[min].second[0]) {
					case REFLECT:
						v = v-2.0*(v.dot(nv)*nv);
						v /= v.norm();
						//std::cout << "intersection from "<<mols.r0[i]<<" to "<<mols.r[i]<<" collision at "<<ip<<" reflecting to "<<ip+v<<std::endl;
						mols.r0[i] = impact_point + GEOMETRY_TOLERANCE*v;
						mols.r[i] = impact_point + (mols.r[i] - impact_point).norm()*v;
						break;
					case ABSORB:
						mols.mark_for_deletion(i);
						delete_molecules = true;
						break;
					case JUMP:
						const Vect3d jump_by(geometries[min].second[1],geometries[min].second[2],geometries[min].second[3]);
						mols.r[i] += jump_by;
						mols.r0[i] = impact_point + jump_by;
						break;
					}
				} else {
					break;
				}
			}

		}
		if (delete_molecules) mols.delete_molecules();

	}
}
}
