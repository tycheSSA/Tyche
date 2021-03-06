/*
 * Reaction.h
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
 *  Created on: 17 Oct 2012
 *      Author: robinsonm
 */

#ifndef REACTION_H_
#define REACTION_H_


#include "Species.h"
#include "Constants.h"
#include "BucketSort.h"
#include "Operator.h"
#include "ReactionEquation.h"
#include "Log.h"

#include <vector>
#include <list>
#include <utility>
#include <boost/random.hpp>
#include <boost/function.hpp>

namespace Tyche {

class Reaction: public Operator {
public:
	Reaction(const double rate):rate(rate) {};
protected:
	double rate;
};

class ZeroOrderMolecularReaction: public Reaction {
public:
	ZeroOrderMolecularReaction(const double rate, const Vect3d min, const Vect3d max):
		Reaction(rate),min(min),max(max) {
	};

	static std::auto_ptr<Operator> New(const double rate, const Vect3d min, const Vect3d max) {
		return std::auto_ptr<Operator>(new ZeroOrderMolecularReaction(rate,min,max));
	}

	void add_species(const double rate, Species& s);

protected:
	virtual void integrate(const double dt);
	virtual void add_species_execute(Species& s);
	virtual void print(std::ostream& out) const;
private:
	std::vector<double> rates;
	Vect3d min,max;
};


class UniMolecularReaction: public Reaction {
public:
	UniMolecularReaction(const double rate,const ReactionEquation& eq);
//	UniMolecularReaction(const double rate, const std::vector<Species*>& reactants,const std::vector<Species*>& products, const double init_radius=0.0):Reaction(rate) {
//		UniMolecularReaction(rate,ReactionSide(reactants)>>ReactionSide(products),init_radius);
//	}

	static std::auto_ptr<Operator> New(const double rate, const ReactionEquation& eq) {
		return std::auto_ptr<Operator>(new UniMolecularReaction(rate,eq));
	}
	size_t add_reaction(const double rate, const ReactionEquation& eq);
	void report_dt_suitability(const Vect3d pos, const double dt);
	void set_geometry(size_t index, Geometry *g) {
		geometries[index] = g;
	}
	std::auto_ptr<Geometry> get_geometry(size_t index) {
		return std::auto_ptr<Geometry>(geometries[index]);
	}
	void set_init_radius(size_t index, double r) {
		init_radii[index] = r;
	}
	double get_init_radius(size_t index) {
		return init_radii[index];
	}
protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tUnimolecular Reaction with reactions:";
		BOOST_FOREACH(ReactionEquation eq, eqs) {
			out << eq<<" (rate = "<<rate<<")";
		}
	}
private:
	double calc_mass_action(ReactionEquation& eq, int i) {
		double action = 1.0;
		for (std::vector<ReactionComponent>::iterator rc=eq.lhs.begin();rc!=eq.lhs.end();rc++) {
			if (rc->compartment_index == SpeciesType::LATTICE) {
				action *= pow(rc->species->copy_numbers[i]/rc->species->grid->get_cell_volume(i)
								,rc->multiplier);
			} else if (rc->compartment_index == SpeciesType::PDE) {
				action *= pow(rc->species->concentrations[i]
												,rc->multiplier);
			}
		}
		ASSERT(action >= 0, "calculated action is less than zero!!");
		return action;
	}
	void calculate_probabilities(const int cell_index, const Vect3d pos, const double dt);
	ReactionEquation& get_random_reaction(const double rand);
	std::vector<ReactionEquation> eqs;
	std::vector<double> probabilities;
	std::vector<double> rates;
	std::vector<double> init_radii;
	std::vector<Geometry *> geometries;
	double total_probability;
	double total_rate;
	std::vector<bool> do_mass_action;
};

class ZeroOrderMolecularReactionLattice: public Reaction {
public:
	ZeroOrderMolecularReactionLattice(const double rate,const ReactionEquation& eq);
	static std::auto_ptr<Operator> New(const double rate, const ReactionEquation& eq) {
		return std::auto_ptr<Operator>(new ZeroOrderMolecularReactionLattice(rate,eq));
	}
	void set_geometry(Geometry *g) {
		geometry = g;
	}
	std::auto_ptr<Geometry> get_geometry() {
		return std::auto_ptr<Geometry>(geometry);
	}
protected:
	virtual void integrate(const double dt);
	double calc_mass_action(int i) {
		double action = 1.0;
		//for (auto& rc : rs.lhs) {
		for (std::vector<ReactionComponent>::iterator rc=eq.lhs.begin();rc!=eq.lhs.end();rc++) {
			if (rc->compartment_index == SpeciesType::LATTICE) {
				action *= pow(rc->species->copy_numbers[i]/rc->species->grid->get_cell_volume(i)
												,rc->multiplier);
			} else {
				ASSERT(rc->compartment_index == SpeciesType::PDE,"SpeciesType cannot be OFF_LATTICE");
				action *= pow(rc->species->concentrations[i]
												,rc->multiplier);
			}
		}
		ASSERT(action >= 0, "calculated action is less than zero!!");
		return action;
	}
	virtual void print(std::ostream& out) const {
		out << "\tReaction Lattice with reaction:";
		out << "\t1"<<eq<<" (rate = "<<rate<<")";
	}
private:
	ReactionEquation eq;
	Geometry *geometry;
};



class BindingReaction: public Reaction {
public:
		BindingReaction(const double rate,
				const double diss_rate,
				Species& species,
				const double binding,
				const double unbinding,
				const double dt,
				Vect3d pos,
				const int binding_sites,
				const int initial_state);

	static std::auto_ptr<Operator> New(const double rate,
					   const double diss_rate,
					   Species& species,
					   const double binding,
					   const double unbinding,
					   const double dt,
					   Vect3d pos,
					   const int binding_sites,
					   const int initial_state) {
	  return std::auto_ptr<Operator>(new BindingReaction(rate,diss_rate,species,binding,unbinding,dt,pos,binding_sites,initial_state));
	}

	double get_rate() const {return this->rate;}
	double get_P_lambda() const {return this->P_lambda;}
	double get_binding_radius() const {return binding_radius;}
	double get_unbinding_radius() const {return unbinding_radius;}
        int get_site_state() const {return site_state;}
	std::list< std::pair< int, double > > get_state_sequence(bool clear=false) {
		std::list< std::pair< int, double > > retlist = std::list< std::pair< int, double > >(state_sequence);
		if (clear) state_sequence.clear();
		return retlist;
	}
	void report_dt_suitability(const double dt);

	void set_state_changed_cb(const boost::function< void(double, int) > callback)
	{
	  have_state_changed_cb = true;
	  state_changed_cb = callback;
	}

	void set_position(const Vect3d& pos) {
	  position = pos;
	}

protected:

	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tBinding Reaction: ("<<get_species()[0]->id<<") (rate = "<<rate<<", binding radius = "<<get_binding_radius()<<")";
	}

	double calculate_lambda(const double dt);

	void suggest_binding_unbinding(const double dt);
	void suggest_binding(const double dt);

	double binding_radius_dt;
	Vect3d position;
	double binding_radius,unbinding_radius;
	double P_lambda;
	double P_diss;
	int binding_sites;
	int site_state;
	bool have_state_changed_cb;
	std::list< std::pair< int, double > > state_sequence;
	boost::function< void(double, int) > state_changed_cb;
};

template<typename T>
class BiMolecularReaction: public Reaction {
public:
	BiMolecularReaction(const double rate, const ReactionEquation& eq,
			const double binding,
			const double unbinding,
			const double dt,
			Vect3d low, Vect3d high, Vect3b periodic,
			const bool reversible=false);
//	BiMolecularReaction(const double rate, const std::vector<Species*>& reactants,const std::vector<Species*>& products,
//				const double binding,
//				const double unbinding,
//				const double dt,
//				Vect3d low, Vect3d high, Vect3b periodic,
//				const bool reversible=false);
	BiMolecularReaction(const double rate, const ReactionEquation& eq, const double dt,
				Vect3d low, Vect3d high, Vect3b periodic, const bool reversible);

	static std::auto_ptr<Operator> New(const double rate,  const ReactionEquation& eq,
			const double binding,
			const double unbinding,
			const double dt,
			Vect3d low, Vect3d high, Vect3b periodic,
			const bool reversible=false) {
		return std::auto_ptr<Operator>(new BiMolecularReaction(rate,eq,binding,unbinding,dt,low,high,periodic,reversible));
	}

	static std::auto_ptr<Operator> New(const double rate,  const ReactionEquation& eq,
				const double dt,
				Vect3d low, Vect3d high, Vect3b periodic,
				const bool reversible=false) {
			return std::auto_ptr<Operator>(new BiMolecularReaction(rate,eq,dt,low,high,periodic,reversible));
		}

	double get_rate() const {return this->rate;}
	double get_P_lambda() const {return this->P_lambda;}
	double get_binding_radius() const {return binding_radius;}
	double get_unbinding_radius() const {return unbinding_radius;}
	void report_dt_suitability(const double dt);

protected:

	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		if (self_reaction) {
			out << "\tBimolecular Reaction: 1("<<get_species()[0]->id<<") + 1("<<get_species()[0]->id<<") >> ";
		} else {
			out << "\tBimolecular Reaction: 1("<<get_species()[0]->id<<") + 1("<<get_species()[1]->id<<") >> ";
		}
		for (auto component : products) {
			out << "1("<<component.species->id<<") ";
		}
		out << "(rate = "<<rate<<", binding radius = "<<get_binding_radius()<<")";
	}

	double calculate_lambda_reversible(const double dt);
	double calculate_lambda_irreversible(const double dt);

	void suggest_binding_unbinding(const double dt);
	void suggest_binding(const double dt);

	double binding_radius_dt;
	ReactionSide products;
	bool reversible;
	double binding_radius,unbinding_radius;
	double P_lambda;
	T neighbourhood_search;
	bool self_reaction;
};


class TriMolecularReaction: public Reaction {
public:
	TriMolecularReaction(const double rate, const ReactionEquation& eq,
			const double dt,
			Vect3d low, Vect3d high, Vect3b periodic);

	static std::auto_ptr<Operator> New(const double rate,  const ReactionEquation& eq,
			const double dt,
			Vect3d low, Vect3d high, Vect3b periodic) {
		return std::auto_ptr<Operator>(new TriMolecularReaction(rate,eq,dt,low,high,periodic));
	}

protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tTriMolecularReaction Reaction: 1("<<get_species()[0]->id<<") + 1("<<get_species()[1]->id<<") + 1("<<get_species()[2]->id<<") >> ";
		for (auto component : products) {
			out << "1("<<component.species->id<<") ";
		}
		out << "(rate = "<<rate<<")";
	}

	ReactionSide products;
	double radius_check;
	BucketSort neighbourhood_search2,neighbourhood_search3;
	double D1,D2,D3;
	double Dbar1,Dbar2;
	double invDbar1,invDbar2;
};


}
#endif /* REACTION_H_ */
