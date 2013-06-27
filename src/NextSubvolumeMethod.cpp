/*
 * NextSubvolumeMethod.cpp
 *
 *  Created on: 15 Oct 2012
 *      Author: robinsonm
 */

#include "NextSubvolumeMethod.h"
#include "Log.h"
#include "Constants.h"
#include <boost/foreach.hpp>

namespace Tyche {
bool ReactionsWithSameRateAndLHS::add_if_same_lhs(const double rate_to_add, const ReactionSide& lhs_to_add, const ReactionSide& rhs_to_add) {
	//std::sort(lhs_to_add);
	if ((lhs_to_add == lhs)&&(rate_to_add == rate)) {
		all_rhs.push_back(rhs_to_add);
		//			if (lhs[0].compartment_index==0) {
		//				std::cout<<"found duplicate lhs, adding new rhs with ci = "<<rhs_to_add[0].compartment_index<<". new all_rhs.size() = "<<all_rhs.size()<<std::endl;
		//			}
		return true;
	} else {
		return false;
	}
}

ReactionSide& ReactionsWithSameRateAndLHS::pick_random_rhs(const double rand) {
	const int n = all_rhs.size();
	if (n==1) return all_rhs[0];
	return all_rhs[int(floor(rand*n))];
}


void ReactionList::list_reactions() {
	BOOST_FOREACH(ReactionsWithSameRateAndLHS& r, reactions) {
		std::cout <<"With rate = " << r.rate << ":" << std::endl;
		BOOST_FOREACH(ReactionSide& rhs, r.all_rhs) {
			BOOST_FOREACH(ReactionComponent& c, r.lhs) {
				std::cout << "(" << c.multiplier << "*" << c.species->id << "<" << c.compartment_index << ">) ";
			}
			std::cout << "-> ";
			BOOST_FOREACH(ReactionComponent& c, rhs) {
				std::cout << "(" << c.multiplier << "*" << c.species->id << "<" << c.compartment_index << ">) ";
			}
			std::cout << std::endl;
		}
	}
}

void ReactionList::add_reaction(const double rate, const ReactionEquation& eq) {

	//		if (eq.rhs.size() > 0) {
	//			if (eq.lhs.size() > 0) {
	//				std::cout << "adding reaction with lhs = (" << eq.lhs[0].compartment_index << "," << eq.lhs[0].multiplier << ") rhs = (" << eq.rhs[0].compartment_index << "," << eq.rhs[0].multiplier << ")" << std::endl;
	//			} else {
	//				std::cout << "adding reaction with lhs = (null) rhs = (" << eq.rhs[0].compartment_index << "," << eq.rhs[0].multiplier << ")" << std::endl;
	//			}
	//		} else {
	//			std::cout << "adding reaction with lhs = (" << eq.lhs[0].compartment_index << "," << eq.lhs[0].multiplier << ") rhs = (null)" << std::endl;
	//		}
	bool added = false;
	//std::cout <<"added unsorted reaction. size of lhs = "<<lhs.size()<<std::endl;
	ReactionSide sorted_lhs(eq.lhs);
	std::sort(sorted_lhs.begin(),sorted_lhs.end());
	//std::cout <<"added sorted reaction. size of lhs = "<<sorted_lhs.size()<<std::endl;
	BOOST_FOREACH(ReactionsWithSameRateAndLHS& r, reactions) {
		added |= r.add_if_same_lhs(rate, sorted_lhs, eq.rhs);
		if (added) break;
	}
	if (!added) {
		reactions.push_back(ReactionsWithSameRateAndLHS(rate, sorted_lhs, eq.rhs));
		propensities.push_back(0);
	}
	my_size++;
}

void ReactionList::clear() {
	reactions.clear();
	propensities.clear();
	my_size = 0;
}

bool ReactionList::delete_reaction(const ReactionEquation& eq) {
	const int n = reactions.size();
	bool found = false;
	int i,j;
	for (i = 0; i < n; ++i) {
		if (eq.lhs == reactions[i].lhs) {
			const int n_r = reactions[i].all_rhs.size();
			for (j = 0; j < n_r; ++j) {
				found |= (eq.rhs == reactions[i].all_rhs[j]);
				if (found) break;
			}
			if (found) break;
		}
	}
//	if (!found) {
//		std::cerr << "ERROR: could not find reaction to erase!!!" << std::endl;
//		return;
//	}
	if (found) {
		reactions[i].all_rhs.erase(reactions[i].all_rhs.begin() + j);
		if (reactions[i].all_rhs.size() == 0) {
			reactions.erase(reactions.begin() + i);
			propensities.erase(propensities.begin()+i);
		}
		my_size--;
	}
	return found;
}

ReactionEquation ReactionList::pick_random_reaction(const double rand) {
	double last_sum_propensities = 0;
	double sum_propensities = 0;
	const int n = reactions.size();
	const double rand_times_total_propensity = rand*total_propensity;
	for (int i = 0; i < n; i++) {
		sum_propensities += propensities[i];
		if (rand_times_total_propensity < sum_propensities) {
			ASSERT(propensities[i] > 0, "chosen reaction with propensity less than or equal to zero");
			const double scaled_rand = (rand_times_total_propensity-last_sum_propensities)/(sum_propensities-last_sum_propensities);
			return ReactionEquation(reactions[i].lhs,reactions[i].pick_random_rhs(scaled_rand));
		}
		last_sum_propensities = sum_propensities;
	}

	std::cerr << "ERROR: should have picked a reaction. rand is either not 0->1 or total_propensity != sum of propensities!!!!!!" << std::endl;
	exit(-1);
	//std::cout << "returning reaction with lhs.size() = " << reactions[n].size()<<std::endl;
	//		const double scaled_rand = (rand_times_total_propensity - sum_propensities)/(total_propensity-sum_propensities);
	//		return ReactionEquation(reactions[n].lhs,reactions[n].pick_random_rhs(rand));
}

double ReactionList::recalculate_propensities() {
	total_propensity = 0;
	inv_total_propensity = 0;
	const int n = reactions.size();
	for (int i = 0; i < n; i++) {
		ReactionsWithSameRateAndLHS& rs = reactions[i];
		double& propensity = propensities[i];
		const int n = rs.lhs.size();
		propensity = 1.0;
		int beta = 0;
		BOOST_FOREACH(ReactionComponent rc, rs.lhs) {
			int copy_number = rc.species->copy_numbers[rc.compartment_index];
			beta += rc.multiplier;
			ASSERT(copy_number >= 0, "copy number is less than zero!!");
			if (copy_number < rc.multiplier) {
				propensity = 0.0;
				break;
			}
			for (int k = 1; k < rc.multiplier; ++k) {
				copy_number *= copy_number-k;
			}
			propensity *= copy_number;
		}
		propensity *= rs.size()*rs.rate;
		ASSERT(propensity >= 0, "calculated propensity is less than zero!!");
		total_propensity += propensity;
		//			if (reactions[i].lhs[0].compartment_index==0) {
		//				std::cout << "reaction with neighbour "<<reactions[i].all_rhs[0][0].compartment_index<<" and num particles = " <<
		//						reactions[i].lhs[0].species->copy_numbers[reactions[i].lhs[0].compartment_index] << " has propensity = "<<propensities[i]<<std::endl;
		//			}
	}
	if (total_propensity != 0) inv_total_propensity = 1.0/total_propensity;
	return inv_total_propensity;
}

static const double LONGEST_TIME = 100000;


NextSubvolumeMethod::NextSubvolumeMethod(const StructuredGrid& subvolumes):
		subvolumes(subvolumes),
		uni(generator,boost::uniform_real<>(0,1)),
		time(0) {
	const int n = subvolumes.size();
	//std::cout << "created "<<n<<" subvolumes"<<std::endl;
	heap.clear();
	for (int i = 0; i < n; ++i) {
		//subvolume_heap_handles.push_back(heap.push(HeapNode(LONGEST_TIME,i)));
		subvolume_heap_handles.push_back(HeapHandle());
		subvolume_reactions.push_back(ReactionList());
		saved_subvolume_reactions.push_back(ReactionList());
	}
}

void NextSubvolumeMethod::reset() {
	Operator::reset();
	reset_all_priorities();
}

void  NextSubvolumeMethod::add_reaction(const double rate, ReactionEquation eq) {
	const int n = subvolumes.size();
	for (int i = 0; i < n; ++i) {
		add_reaction_to_compartment(rate,eq,i);
	}
}

void  NextSubvolumeMethod::add_reaction_to_compartment(const double rate, ReactionEquation eq, const int i) {
	eq.lhs.set_compartment_index(i);
	eq.rhs.set_compartment_index(i);
	const int beta = eq.lhs.get_num_reactants();
	if (beta == 0) {
		subvolume_reactions[i].add_reaction(rate*subvolumes.get_cell_volume(i),eq);

	} else if (beta == 1) {
		subvolume_reactions[i].add_reaction(rate,eq);

	} else {
		subvolume_reactions[i].add_reaction(rate*pow(subvolumes.get_cell_volume(i),1-eq.lhs.get_num_reactants()),eq);
	}
}

void NextSubvolumeMethod::add_diffusion(Species &s, const double rate) {
	const int n = subvolumes.size();
		for (int i = 0; i < n; ++i) {
			const std::vector<int>& neighbrs = subvolumes.get_neighbour_indicies(i);
			const int nn = neighbrs.size();
			for (int j = 0; j < nn; ++j) {
				ReactionSide lhs;
				lhs.push_back(ReactionComponent(1.0,s,i));
				ReactionSide rhs;
				rhs.push_back(ReactionComponent(1.0,s,neighbrs[j]));
//				if (i==0) {
//					std::cout <<"adding reaction with rate = "<<rate<<" and lhs[0].compartment_index = "<<lhs[0].compartment_index<<" and rhs[0].compartment_index"<<rhs[0].compartment_index<<" to the 0th compartment"<<std::endl;
//				}
				subvolume_reactions[i].add_reaction(rate,ReactionEquation(lhs,rhs));
			}
		}
}

void NextSubvolumeMethod::add_diffusion(Species &s) {
	const int n = subvolumes.size();
		for (int i = 0; i < n; ++i) {
			const std::vector<int>& neighbrs = subvolumes.get_neighbour_indicies(i);
			const std::vector<double>& neighbrs_distances = subvolumes.get_neighbour_distances(i);
			const int nn = neighbrs.size();
			for (int j = 0; j < nn; ++j) {
				const double rate = s.D/pow(neighbrs_distances[j],2);
				ReactionSide lhs; 
				lhs.push_back(ReactionComponent(1.0,s,i));
				ReactionSide rhs; 
				rhs.push_back(ReactionComponent(1.0,s,neighbrs[j]));
//				if (i==0) {
//					std::cout <<"adding reaction with rate = "<<rate<<" and lhs[0].compartment_index = "<<lhs[0].compartment_index<<" and rhs[0].compartment_index"<<rhs[0].compartment_index<<" to the 0th compartment"<<std::endl;
//				}
				subvolume_reactions[i].add_reaction(rate,ReactionEquation(lhs,rhs));
			}
		}
}

void NextSubvolumeMethod::add_diffusion_between(Species &s, const double rate, std::vector<int>& from, std::vector<int>& to) {
	ASSERT(from.size() == to.size(), "From and To vectors must be the same length");
	const int n = from.size();
	for (int i = 0; i < n; ++i) {
		ReactionSide lhs;
		lhs.push_back(ReactionComponent(1.0,s,from[i]));
		ReactionSide rhs;
		rhs.push_back(ReactionComponent(1.0,s,to[i]));
		subvolume_reactions[from[i]].add_reaction(rate,ReactionEquation(lhs,rhs));
	}

}

void NextSubvolumeMethod::reset_all_priorities() {
	const int n = subvolumes.size();
	for (int i = 0; i < n; ++i) {
		reset_priority(i);
	}
}

void NextSubvolumeMethod::reset_priority(const int i) {
	const bool in_queue = subvolume_reactions[i].get_propensity()!=0;
	const double inv_total_propensity = subvolume_reactions[i].recalculate_propensities();

	if (inv_total_propensity != 0) {
		const double time_at_next_reaction = time - inv_total_propensity*log(uni());
		if (in_queue) {
			(*subvolume_heap_handles[i]).time_at_next_reaction = time_at_next_reaction;
			heap.update(subvolume_heap_handles[i]);
		} else {
			subvolume_heap_handles[i] = heap.push(HeapNode(time_at_next_reaction,i));
		}
	} else {
		if (in_queue) {
			heap.erase(subvolume_heap_handles[i]);
		}
	}
//	if (i==929) {
//		std::cout <<"recalcing priority for 929: time to next react = "<<h.time_at_next_reaction<<" total propensity = "<<1.0/inv_total_propensity<<std::endl;
//	}

}

void NextSubvolumeMethod::recalc_priority(const int i) {
	const double old_propensity = subvolume_reactions[i].get_propensity();
	const bool in_queue = old_propensity!=0;
	const double inv_total_propensity = subvolume_reactions[i].recalculate_propensities();
	//HeapNode& h = *(subvolume_heap_handles[i]);
	if (inv_total_propensity != 0) {
		//const double time_to_next_reaction = -inv_total_propensity*log(uni());
		if (in_queue) {
			const double old_time = (*subvolume_heap_handles[i]).time_at_next_reaction;
			const double new_time = time + old_propensity*inv_total_propensity*(old_time - time);
			(*subvolume_heap_handles[i]).time_at_next_reaction = new_time;
			heap.update(subvolume_heap_handles[i]);
		} else {
			const double new_time = time - inv_total_propensity*log(uni());
			subvolume_heap_handles[i] = heap.push(HeapNode(new_time,i));
		}

	} else {
		if (in_queue) {
			heap.erase(subvolume_heap_handles[i]);
		}
	}
//	if (i==929) {
//		std::cout <<"recalcing priority for 929: time to next react = "<<h.time_at_next_reaction<<" total propensity = "<<1.0/inv_total_propensity<<std::endl;
//	}
}

void NextSubvolumeMethod::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	const double final_time = time + dt;
	while (get_next_event_time() < final_time) {
		const int sv_i = heap.top().subvolume_index;
		time = heap.top().time_at_next_reaction;
		//std::cout << "dealing with subvolume with time = " << time << " and index = " << sv_i << std::endl;
		const double rand = uni();
		ReactionEquation r = subvolume_reactions[sv_i].pick_random_reaction(rand);
		react(r);
	}
	time = final_time;
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

void NextSubvolumeMethod::list_reactions() {
	const int n = subvolumes.size();
	for (int i = 0; i < n; ++i) {
		std::cout << "Compartment " << i << " has the following reactions:" <<std::endl;
		subvolume_reactions[i].list_reactions();
	}
}

void NextSubvolumeMethod::set_interface_reactions(Species& s,
		std::vector<int>& indicies, const double dt) {

	/*
	 * save reactions for all indicies. if no reactions assume it is already set to an interface
	 * and do not save
	 */
	BOOST_FOREACH(int i,indicies) {
		if (subvolume_reactions[i].size() > 0) {
			saved_subvolume_reactions[i] = subvolume_reactions[i];
			subvolume_reactions[i].clear();
			reset_priority(i);
		}
	}

	/*
	 * update diffusion reaction rates for neighbouring cells
	 */
	BOOST_FOREACH(int i,indicies) {
		const std::vector<int>& neighbours = subvolumes.get_neighbour_indicies(i);
		const std::vector<double>& neighbrs_distances = subvolumes.get_neighbour_distances(i);
		const int n = neighbours.size();
		for (int ii = 0; ii < n; ++ii) {
			const int j = neighbours[ii];
			const double distance = neighbrs_distances[ii];
			if (subvolume_reactions[j].size() == 0) continue;
			ReactionSide lhs;
			lhs.push_back(ReactionComponent(1.0,s,j));
			ReactionSide rhs;
			rhs.push_back(ReactionComponent(1.0,s,i));
			if (subvolume_reactions[j].delete_reaction(ReactionEquation(lhs,rhs))) {
				const double rate = (2.0*distance/sqrt(PI*s.D*dt))*(s.D/pow(distance,2));
				subvolume_reactions[j].add_reaction(rate,ReactionEquation(lhs,rhs));
				//std::cout << "new number of reactions is "<<subvolume_reactions[j].size()<< " for compartment "<<j<<std::endl;
				reset_priority(j);
			}
		}
	}
}

void NextSubvolumeMethod::unset_interface_reactions(Species& s,
		std::vector<int>& indicies) {
	/*
	 * update diffusion reaction rates for neighbouring cells
	 */
	BOOST_FOREACH(int i,indicies) {
		const std::vector<int>& neighbours = subvolumes.get_neighbour_indicies(i);
		const std::vector<double>& neighbrs_distances = subvolumes.get_neighbour_distances(i);
		const int n = neighbours.size();
		for (int ii = 0; ii < n; ++ii) {
			const int j = neighbours[ii];
			const double distance = neighbrs_distances[ii];
			if (subvolume_reactions[j].size() == 0) continue;
			ReactionSide lhs;
			lhs.push_back(ReactionComponent(1.0,s,j));
			ReactionSide rhs;
			rhs.push_back(ReactionComponent(1.0,s,i));
			if (subvolume_reactions[j].delete_reaction(ReactionEquation(lhs,rhs))) {
				const double rate = s.D/pow(distance,2);
				subvolume_reactions[j].add_reaction(rate,ReactionEquation(lhs,rhs));
				reset_priority(j);
			}
		}
	}
	/*
	 * return these cells to their original reactions list
	 */
	BOOST_FOREACH(int i,indicies) {
		//ASSERT(s.copy_numbers[i] == 0, "Can't unset interface cells that have particles in them.");
		subvolume_reactions[i] = saved_subvolume_reactions[i];
		reset_priority(i);
	}

}



void NextSubvolumeMethod::react(ReactionEquation& eq) {
	bool print = false;
	BOOST_FOREACH(ReactionComponent rc, eq.lhs) {
		rc.species->copy_numbers[rc.compartment_index] -= rc.multiplier;
//		if ((rc.compartment_index==0)||(rc.compartment_index==560)) {
//			print = true;
//			std::cout<<" compartment "<<rc.compartment_index<<" changed to have "<<rc.species->copy_numbers[rc.compartment_index]<<" molecules"<<std::endl;
//		}
	}
	if (eq.lhs.size() > 0) reset_priority(eq.lhs[0].compartment_index);
	BOOST_FOREACH(ReactionComponent rc, eq.rhs) {
		rc.species->copy_numbers[rc.compartment_index] += rc.multiplier;
//		if ((rc.compartment_index==0)||(rc.compartment_index==560)) {
//			std::cout<<" compartment "<<rc.compartment_index<<" changed to have "<<rc.species->copy_numbers[rc.compartment_index]<<" molecules"<<std::endl;
//		}
	}
	if (eq.rhs.size() > 0) recalc_priority(eq.rhs[0].compartment_index);
}

std::ostream& operator<< (std::ostream& out, NextSubvolumeMethod &b) {
	return out << "\tNext Subvolume Method";
}

}
