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
#include <boost/random.hpp>

namespace Tyche {

class Reaction: public Operator {
public:
	Reaction(const double rate):rate(rate) {};
protected:
	const double rate;
};

class ZeroOrderMolecularReaction: public Reaction {
public:
	ZeroOrderMolecularReaction(const double rate, const Vect3d min, const Vect3d max):
		Reaction(rate),min(min),max(max) {
	};

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
	UniMolecularReaction(const double rate,const ReactionEquation& eq, const double init_radius=0.0);
	void add_reaction(const double rate, const ReactionEquation& eq, const double init_radius=0.0);
	void report_dt_suitability(const double dt);
protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const {
		out << "\tUnimolecular Reaction with reactions:";
		BOOST_FOREACH(ReactionSide side, product_list) {
			out << "\t1("<<get_species()[0]->id<<") >> "<<side<<" (rate = "<<rate<<")";
		}
	}
private:
	void calculate_probabilities(const double dt);
	ReactionSide& get_random_reaction(const double rand);
	std::vector<ReactionSide> product_list;
	std::vector<double> probabilities;
	std::vector<double> rates;
	std::vector<double> init_radii;
	double total_probability;
	double total_rate;
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
	BiMolecularReaction(const double rate, const ReactionEquation& eq, const double dt,
				Vect3d low, Vect3d high, Vect3b periodic, const bool reversible);

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

	double binding_radius_dt;
	ReactionSide products;
	bool reversible;
	double binding_radius,unbinding_radius;
	double P_lambda;
	T neighbourhood_search;
	bool self_reaction;
};

}
#endif /* REACTION_H_ */
