/*
 * Boundary.h
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
 *  Created on: 12 Oct 2012
 *      Author: robinsonm
 */

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "Species.h"
#include "Vector.h"
#include "Geometry.h"
#include "Operator.h"
#include "Union.h"
#include "NextSubvolumeMethod.h"
#include <iostream>
#include <array>


namespace Tyche {


class BoundaryController: public Operator {
public:
	enum BOUNDARY_TYPE {REFLECT,JUMP,ABSORB};
	BoundaryController() {
		 bt_names[REFLECT] = "Reflection";
		 bt_names[JUMP] = "Jumping";
		 bt_names[ABSORB] = "Absorption";
	}
	static std::auto_ptr<BoundaryController> New() {
		return std::auto_ptr<BoundaryController> (new BoundaryController());
	}
	void add_jump_boundary(Geometry &geometry, const Vect3d jump_by);
	void add_reflect_boundary(Geometry &geometry);
	void add_absorb_boundary(Geometry &geometry);
protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const  {
		out << "\tBoundary Controller at";
		for (auto i: geometries) {
			out << " "<<bt_names[(int)i.second[0]]<< " "<< *(i.first);
		}
	}
private:
	typedef std::array<double,10> param_type;
	typedef std::pair<Geometry*,param_type> data_type;
	std::array<std::string,3> bt_names;
	std::vector<data_type> geometries;

};

template<typename T>
class Boundary: public Operator {
public:
	Boundary<T>(const T& geometry):geometry(geometry) {
	}
	const T& geometry;

protected:
	virtual void print(std::ostream& out) const  {
		out << " at "<< this->geometry;
		Operator::print(out);
	}
};


template<typename T>
class DestroyBoundary: public Boundary<T> {
public:
	DestroyBoundary(const T& geometry):
		Boundary<T>(geometry) {}
protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const  {
		out << "\tDestroy Boundary at "<< this->geometry;
	}

};


template<typename T>
class JumpBoundary: public Boundary<T> {
public:
	JumpBoundary(const T& geometry, const Vect3d jump_by):
		Boundary<T>(geometry),jump_by(jump_by) {}
	static std::auto_ptr<Operator> New(const T& geometry, const Vect3d jump_by) {
			return std::auto_ptr<Operator> (new JumpBoundary(geometry,jump_by));
	}
protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tJump Boundary at "<< this->geometry;
	}

	const Vect3d jump_by;

};


template<typename T>
JumpBoundary<T> create_jump_boundary(T& geometry, const Vect3d jump_by) {
	return JumpBoundary<T>(geometry,jump_by);
}


template<typename T>
class DiffusionCorrectedBoundary: public Boundary<T> {
public:
   DiffusionCorrectedBoundary(const T& geometry):
      Boundary<T>(geometry),
      uni(generator,boost::uniform_real<>(0,1)) {

   }
   virtual ~DiffusionCorrectedBoundary() {
	   BOOST_FOREACH(std::vector<double>* i, all_prev_distance) {
		   delete i;
	   }
	   BOOST_FOREACH(std::vector<double>* i, all_curr_distance) {
		   delete i;
	   }
   }

protected:
   virtual void add_species_execute(Species &s);



   void timestep_initialise(const double dt);
   void timestep_finalise();
   bool particle_crossed_boundary(const int p_i, const int s_i);
   void recalc_constants(Species &s, const double new_dt) {
      D_dt = s.D*new_dt;
      test_this_distance_from_wall = 5.0*sqrt(2.0*s.D*new_dt);
      //test_this_distance_from_wall = 0;
   }
   void init_prev_distance(Molecules& mols, std::vector<double>& prev_distance) {
      const int n = mols.size();
      prev_distance.resize(n);
      for (int i = 0; i < n; ++i) {
         prev_distance[i] = this->geometry.distance_to_boundary(mols.r[i]);
      }
   }
   std::vector<std::vector<double>* > all_prev_distance, all_curr_distance;
   double D_dt;
   double test_this_distance_from_wall;
   boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni;
};

template<typename T>
class RemoveBoundaryWithCorrection: public DiffusionCorrectedBoundary<T> {
public:
	RemoveBoundaryWithCorrection(const T& geometry):
		DiffusionCorrectedBoundary<T>(geometry) {}
	Molecules& get_removed(Species& s);
protected:
	virtual void integrate(const double dt);
	virtual void add_species_execute(Species& s);
	virtual void print(std::ostream& out) const {
		out << "\tRemove Boundary With Correction at "<< this->geometry;
	}

private:
	std::vector<Molecules> removed_molecules;
};


template<typename T>
RemoveBoundaryWithCorrection<T> create_remove_boundary_corrected(T& geometry) {
	return RemoveBoundaryWithCorrection<T>(geometry);
}

template<typename T>
RemoveBoundaryWithCorrection<T> create_remove_boundary_corrected(T& geometry, Species& s) {
	RemoveBoundaryWithCorrection<T> to_return(geometry); to_return.add_species(s);
	return to_return;
}

template<typename T>
class RemoveBoundary: public Boundary<T> {
public:
	RemoveBoundary(const T& geometry):
		Boundary<T>(geometry) {}
	Molecules& get_removed(Species& s);

protected:

	virtual void integrate(const double dt);
	virtual void add_species_execute(Species& s);
	virtual void print(std::ostream& out) const {
		out << "\tRemove Boundary at "<< this->geometry;
	}

private:
	std::vector<Molecules> removed_molecules;
};


template<typename T>
RemoveBoundary<T> create_remove_boundary(T& geometry) {
	return RemoveBoundary<T>(geometry);
}

template<typename T>
RemoveBoundary<T> create_remove_boundary(T& geometry, Species& s) {
	RemoveBoundary<T> to_return(geometry); to_return.add_species(s);
	return to_return;
}


template<typename T>
class JumpBoundaryWithCorrection: public DiffusionCorrectedBoundary<T> {
public:
	JumpBoundaryWithCorrection(const T& geometry, const Vect3d jump_by):
		DiffusionCorrectedBoundary<T>(geometry),jump_by(jump_by) {}
	static std::auto_ptr<Operator> New(const T& geometry, const Vect3d jump_by) {
		return std::auto_ptr<Operator>(new JumpBoundaryWithCorrection(geometry,jump_by));
	}
protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tJump Boundary With Correction at "<< this->geometry;
	}


   const Vect3d jump_by;

};


template<typename T>
JumpBoundaryWithCorrection<T> create_jump_boundary_corrected(T& geometry, const Vect3d jump_by) {
	return JumpBoundaryWithCorrection<T>(geometry,jump_by);
}

template<typename T>
JumpBoundaryWithCorrection<T> create_jump_boundary_corrected(T& geometry, const Vect3d jump_by, Species& s) {
	JumpBoundaryWithCorrection<T> to_return(geometry,jump_by); to_return.add_species(s);
	return to_return;
}


template<typename T>
class ReflectiveBoundary: public Boundary<T> {
public:
	ReflectiveBoundary(const T& geometry):
		Boundary<T>(geometry) {}
	static std::auto_ptr<Operator> New(const T& geometry) {
		return std::auto_ptr<Operator>(new ReflectiveBoundary(geometry));
	}
protected:

	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tReflective Boundary";
		Boundary<T>::print(out);
	}
};

template<typename T>
ReflectiveBoundary<T> create_reflective_boundary(T& geometry) {
	return ReflectiveBoundary<T>(geometry);
}

template<typename T>
ReflectiveBoundary<T> create_reflective_boundary(T& geometry, Species& s) {
	ReflectiveBoundary<T> to_return(geometry); to_return.add_species(s);
	return to_return;
}

class FluxBoundary: public Boundary<NullGeometry> {
public:
	FluxBoundary(const Vect3d p, const Vect3d t1, const Vect3d t2, const double rate):
		Boundary<NullGeometry>(NullGeometry()),
		p(p),
		t1(t1),
		t2(t2),
		rate(rate),
		uni1(generator,boost::uniform_real<>(0,t1.norm())),
		uni2(generator,boost::uniform_real<>(0,t2.norm()))
		{}
	static std::auto_ptr<Operator> New(const Vect3d p, const Vect3d t1, const Vect3d t2, const double rate) {
		return std::auto_ptr<Operator>(new FluxBoundary(p,t1,t2,rate));
	}

	const double rate;
	const Vect3d p,t1,t2;

protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tFlux Boundary at (x,y,z) = " << p << " + s*" << t1 << " + t*" << t2 << " with s = (0 -> 1) and t = (0 -> 1)";
	}
private:
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni1, uni2;
};



template<typename T>
class CouplingBoundary_M_to_C: public DiffusionCorrectedBoundary<T> {
public:
	CouplingBoundary_M_to_C(const T& geometry, NextSubvolumeMethod& nsm):
		DiffusionCorrectedBoundary<T>(geometry),nsm(nsm) {}
	static std::auto_ptr<Operator> New(const T& geometry, NextSubvolumeMethod& nsm) {
		return std::auto_ptr<Operator>(new CouplingBoundary_M_to_C(geometry,nsm));
	}
protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const  {
		out << "\tCoupling Boundary from Molecules to Compartments at "<< std::endl << "\t\t"<< this->geometry;
	}

private:
	NextSubvolumeMethod& nsm;

};
template<typename T>
class CouplingBoundary: public Boundary<T> {
public:
	CouplingBoundary(const T& geometry, NextSubvolumeMethod* nsm, const bool corrected):
		Boundary<T>(geometry),nsm(nsm),corrected(corrected) {
	}
	static std::auto_ptr<Operator> New(const T& geometry, NextSubvolumeMethod* nsm, const bool corrected) {
		return std::auto_ptr<Operator>(new CouplingBoundary(geometry,nsm,corrected));
	}
protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const  {
		out << "\tCoupling Boundary from Molecules to Compartments at "<< std::endl << "\t\t"<< this->geometry;
	}

private:
	NextSubvolumeMethod* nsm;
	bool corrected;

};




}





#include "Boundary.impl.h"



#endif /* BOUNDARY_H_ */
