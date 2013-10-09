/*
 * NextSubvolumeMethod.h
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
 *  Created on: 15 Oct 2012
 *      Author: robinsonm
 */

#ifndef NEXTSUBVOLUMEMETHOD_H_
#define NEXTSUBVOLUMEMETHOD_H_

#include <boost/heap/fibonacci_heap.hpp>
#include <vector>
#include "MyRandom.h"
#include "Species.h"
#include "Operator.h"
#include "StructuredGrid.h"
#include "ReactionEquation.h"

namespace Tyche {

#define LOWEST_PRIORITY 0

struct HeapNode {
	HeapNode(double time_at_next_reaction, unsigned int subvolume_index):time_at_next_reaction(time_at_next_reaction),subvolume_index(subvolume_index) {}
	HeapNode(const HeapNode& hn) {
		time_at_next_reaction = hn.time_at_next_reaction;
		subvolume_index = hn.subvolume_index;
	}

	double time_at_next_reaction;
	unsigned int subvolume_index;

	bool operator<(HeapNode const & rhs) const {
		return time_at_next_reaction > rhs.time_at_next_reaction;
	}
};
typedef boost::heap::fibonacci_heap<HeapNode> PriorityHeap;
typedef boost::heap::fibonacci_heap<HeapNode>::handle_type HeapHandle;



struct ReactionsWithSameRateAndLHS {
	ReactionsWithSameRateAndLHS(const double rate, const ReactionSide& lhs, const ReactionSide& rhs):
		lhs(lhs),
		rate(rate) {
		all_rhs.push_back(rhs);
		//std::cout <<"added new reaction. size of lhs = "<<lhs.size()<<std::endl;
	}
//	~ReactionsWithSameRateAndLHS() {
//		all_rhs.clear();
//	}
	bool add_if_same_lhs(const double rate_to_add, const ReactionSide& lhs_to_add, const ReactionSide& rhs_to_add);
	ReactionSide& pick_random_rhs(const double rand);

	int size() {
		return all_rhs.size();
	}

	ReactionSide lhs;
	double rate;
	std::vector<ReactionSide> all_rhs;
};



class ReactionList {
public:
	ReactionList():
		total_propensity(0),my_size(0),
		inv_total_propensity(0) {}
//	~ReactionList() {
//		reactions.clear();
//		propensities.clear();
//	}
	void operator=(const ReactionList& arg) {
		reactions = arg.reactions;
		propensities.assign(arg.propensities.size(),0);
		total_propensity = 0;
		inv_total_propensity = 0;
		my_size = arg.my_size;
	}
	void list_reactions();
	void add_reaction(const double rate, const ReactionEquation& eq);
	double delete_reaction(const ReactionEquation& eq);
	void clear();
	ReactionEquation pick_random_reaction(const double rand);
	double recalculate_propensities();
	double get_propensity() {
		return total_propensity;
	}
	int size() {
		return my_size;
	}
private:
	double total_propensity;
	double my_size;
	std::vector<ReactionsWithSameRateAndLHS> reactions;
	std::vector<double> propensities;
	double inv_total_propensity;
};

//template<typename T>
//void set_interface(Operator *op, const T& geometry, const double dt, const bool corrected) {
//	NextSubvolumeMethod *nsm = dynamic_cast<NextSubvolumeMethod *>(op);
//	if (!nsm) ERROR("Must pass NextSubvolumeMethod operator to set_interface");
//	nsm->set_interface(geometry,dt,corrected);
//}
	

class NextSubvolumeMethod: public Operator {
public:
	NextSubvolumeMethod(StructuredGrid& subvolumes);
	static std::auto_ptr<NextSubvolumeMethod> New(StructuredGrid& subvolumes) {
		return std::auto_ptr<NextSubvolumeMethod>(new NextSubvolumeMethod(subvolumes));
	}

	void list_reactions();
	template<typename T>
	void set_interface(const T& geometry, const double dt, const bool corrected) {
		std::vector<int> to_indicies,from_indicies;
		T to_geometry = geometry;
		to_geometry += 0.5*subvolumes.get_cell_size()[T::dim];
		T from_geometry = geometry;
		from_geometry.swap_normal();
		from_geometry += 0.5*subvolumes.get_cell_size()[T::dim];
		subvolumes.get_slice(from_geometry,from_indicies);
		subvolumes.get_slice(to_geometry,to_indicies);
		set_interface_reactions(from_indicies,to_indicies,dt,corrected);
	}

	template<typename T>
	void unset_interface(const T& geometry) {
		std::vector<int> to_indicies,from_indicies;
		subvolumes.get_slice(geometry,to_indicies);
		T opposite_normal = geometry;
		opposite_normal.swap_normal();
		subvolumes.get_slice(geometry,from_indicies);
		unset_interface_reactions(from_indicies,to_indicies);
	}
	void set_interface_reactions(std::vector<int>& from_indicies, std::vector<int>& to_indicies, const double dt, const bool corrected);
	void unset_interface_reactions(std::vector<int>& from_indicies, std::vector<int>& to_indicies);
	void add_diffusion(Species &s);
	void add_reaction(const double rate, ReactionEquation eq);
	template<typename T>
	void add_reaction_on(const double rate, ReactionEquation eq, const T& geometry) {
		std::vector<int> indicies;
		subvolumes.get_slice(geometry,indicies);
		for (unsigned int i = 0; i < indicies.size(); ++i) {
			add_reaction_to_compartment(rate,eq,i);
		}
	}

	template<typename T>
	void add_reaction_in(const double rate, ReactionEquation eq, const T& geometry) {
		std::vector<int> indicies;
		subvolumes.get_region(geometry,indicies);
		for (unsigned int i = 0; i < indicies.size(); ++i) {
			add_reaction_to_compartment(rate,eq,i);
		}
	}

	template<typename T>
	void add_diffusion_between(Species &s, const double rate, T& geometry_from, T& geometry_to) {
		std::vector<int> from,to;
		subvolumes.get_slice(geometry_from,from);
		subvolumes.get_slice(geometry_to,to);
		add_diffusion_between(s,rate,from,to);
	}
	void add_diffusion_between(Species &s, const double rate, std::vector<int>& from, std::vector<int>& to);

	void clear_reactions(std::vector<int>& cell_indicies) {
		//for (int i : cell_indicies) {
		for (std::vector<int>::iterator i=cell_indicies.begin();i!=cell_indicies.end();i++) {
			//subvolume_reactions[i].clear();
			subvolume_reactions[*i].clear();
		}
	}

	void reset_all_priorities();
	void reset_priority(const int i);
	void recalc_priority(const int i);
	double get_next_event_time() {
		return heap.top().time_at_next_reaction;
	}
	double get_time() {return time;}
	const StructuredGrid& get_grid() const { return subvolumes; }
	void add_reaction_to_compartment(const double rate, ReactionEquation eq, int i);

protected:
	virtual void add_species_execute(Species &s);

	virtual void reset_execute();
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const;
private:
	void react(ReactionEquation& r);
	StructuredGrid& subvolumes;
	PriorityHeap heap;
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni;
	double time;
	std::vector<ReactionList> subvolume_reactions;
	std::vector<ReactionList> saved_subvolume_reactions;
	std::vector<HeapHandle> subvolume_heap_handles;
};


}

#endif /* NEXTSUBVOLUMEMETHD_H */
