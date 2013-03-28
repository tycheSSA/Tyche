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

class UniMolecularReaction: public Reaction {
public:
	UniMolecularReaction(const double rate,const ReactionEquation& eq):Reaction(rate) {
		CHECK((eq.lhs.size()==1) && (eq.lhs[0].multiplier == 1), "Reaction equation is not unimolecular!");
		this->add_species(*(eq.lhs[0].species));
		product_list.push_back(eq.rhs);
		probabilities.push_back(0);
		total_rate = rate;
		rates.push_back(rate);
	};
	void add_reaction(const double rate, const ReactionEquation& eq) {
		CHECK((eq.lhs.size()==1) && (eq.lhs[0].multiplier == 1), "Reaction equation is not unimolecular!");
		CHECK(eq.lhs[0].species == all_species[0], "Reactant is different from previous equation");
		product_list.push_back(eq.rhs);
		probabilities.push_back(0);
		rates.push_back(rate);
		total_rate += rate;
	}
	void operator()(const double dt);
	friend std::ostream& operator<< (std::ostream& out, UniMolecularReaction &r);
private:
	void calculate_probabilities(const double dt);
	ReactionSide& get_random_reaction(const double rand);
	std::vector<ReactionSide> product_list;
	std::vector<double> probabilities;
	std::vector<double> rates;
	double total_probability;
	double total_rate;
};

std::ostream& operator<< (std::ostream& out, UniMolecularReaction &r);


template<typename T>
class BiMolecularReaction: public Reaction {
public:
	BiMolecularReaction(const double rate, const ReactionEquation& eq, const double dt, Vect3d low, Vect3d high, Vect3b periodic);
	void operator()(const double dt);
	double get_rate() {return this->rate;}
	double get_binding_radius() {return binding_radius;}
protected:
	double binding_radius_dt;
	ReactionSide products;
	double binding_radius, binding_radius2;
	T neighbourhood_search;
	bool self_reaction;
};

template<typename T>
std::ostream& operator<< (std::ostream& out, BiMolecularReaction<T> &r) {
	return out << "\tBimolecular Reaction with rate = "<<r.get_rate()<<" and binding radius = "<<r.get_binding_radius();
}
}
#endif /* REACTION_H_ */
