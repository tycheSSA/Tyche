/*
 * Reaction.cpp
 *
 *  Created on: 17 Oct 2012
 *      Author: robinsonm
 */

#include "Reaction.h"
extern "C" {
#include "rxnparam.h"
}
#include <map>

namespace Tyche {
void UniMolecularReaction::calculate_probabilities(const double dt) {
	total_probability = (1.0-exp(-dt*total_rate));
	const int n = probabilities.size();
	probabilities[0] = (rates[0]/total_rate)*total_probability;
	for (int i = 1; i < n; ++i) {
		probabilities[i] = probabilities[i-1] + (rates[i]/total_rate)*total_probability;
	}
	ASSERT(probabilities[n-1]==total_probability,
			"total_probability = "<<total_probability<<" should be equal to last entry of cummulative sum (="
			<<probabilities[n-1]<<")");
}

ReactionSide& UniMolecularReaction::get_random_reaction(const double rand) {
	const int n_minus_one = probabilities.size()-1;
	ASSERT((rand < total_probability) && (rand < probabilities[n_minus_one]),"rand is outside allowed range");
	for (int i = 0; i < n_minus_one; ++i) {
		if (rand < probabilities[i]) return product_list[i];
	}
	return product_list[n_minus_one];
}

UniMolecularReaction::UniMolecularReaction(const double rate,const ReactionEquation& eq, const double dt, const double reverse_rate):Reaction(rate) {
	CHECK((eq.lhs.size()==1) && (eq.lhs[0].multiplier == 1), "Reaction equation is not unimolecular!");
	this->add_species(*(eq.lhs[0].species));
	product_list.push_back(eq.rhs);
	probabilities.push_back(0);
	if (dt>0) {
		double difc;
		if (eq.rhs.size()==1) {
			difc = 2.0*eq.rhs[0].species->D;
		} else {
			difc = eq.rhs[0].species->D+eq.rhs[1].species->D;
		}
		const double binding_radius = bindingradius(reverse_rate,dt,difc,-1,0);
		LOG(2,"reverse binding radius for reaction "<<eq<<" with reverse rate = "<<reverse_rate<<" is calculated as r = "<<binding_radius);
//		const double unbinding_radius = unbindingradius(0.5,dt,difc,binding_radius);
//		LOG(2,"unbinding radius for reaction "<<eq<<" with reverse rate = "<<reverse_rate<<" is calculated as r = "<<unbinding_radius);

		init_radii.push_back(2.0*binding_radius);
	} else {
		init_radii.push_back(0);
	}
	total_rate = rate;
	rates.push_back(rate);
};

void UniMolecularReaction::add_reaction(const double rate, const ReactionEquation& eq, const double init_radius) {
		CHECK((eq.lhs.size()==1) && (eq.lhs[0].multiplier == 1), "Reaction equation is not unimolecular!");
		CHECK(eq.lhs[0].species == all_species[0], "Reactant is different from previous equation");
		product_list.push_back(eq.rhs);
		probabilities.push_back(0);
		rates.push_back(rate);
		init_radii.push_back(init_radius);
		total_rate += rate;
	}

void UniMolecularReaction::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
#ifdef DEBUG
	std::map<int,int> count;
	count[all_species[0]->id] = 0;
	std::map<int,int>::iterator it;
#endif
	calculate_probabilities(dt);
	boost::uniform_real<> uni_dist(0.0,1.0);
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
	Molecules &mols = all_species[0]->mols;
	const int n = mols.size();
	for (int i = 0; i < n; ++i) {
		if (mols.saved_index[i] == SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE) continue;
		const double rand = uni();
		if (rand < total_probability) {
			ReactionSide& products = get_random_reaction(rand);
			BOOST_FOREACH(ReactionComponent component, products) {
				for (int j = 0; j < component.multiplier; ++j) {
					if (j==0) {
						component.species->mols.add_molecule(mols.r[i]);
					} else {
					const double angle = 2.0*PI*uni();
					const double angle2 = 2.0*PI*uni();
					component.species->mols.add_molecule(mols.r[i] + Vect3d(init_radii[0]*sin(angle)*cos(angle2),
																			init_radii[0]*sin(angle)*sin(angle2),
																			init_radii[0]*cos(angle)));
					}
#ifdef DEBUG
					//LOG(2,"adding molecule of species "<<component.species->id);
					it = count.find(component.species->id);
					if (it==count.end()) {
						count.insert(std::pair<int,int>(component.species->id,1));
					} else {
						it->second++;
					}
#endif
				}
			}
			mols.mark_for_deletion(i);
#ifdef DEBUG
			count[all_species[0]->id]--;
#endif
		}
	}
	mols.delete_molecules();
#ifdef DEBUG
	for (it = count.begin();it!=count.end();it++) {
		LOG(2,"Created/Deleted "<<it->second<<" molecules of species ("<<it->first<<")");
	}
#endif
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

std::ostream& operator<< (std::ostream& out, UniMolecularReaction &r) {
	out << "\tUnimolecular Reaction with reactions:";
	BOOST_FOREACH(ReactionSide side, r.product_list) {
		out << "\t1("<<r.all_species[0]->id<<") >> "<<side<<" (rate = "<<r.rate<<")";
	}
	return out;
}


template<typename T>
void BiMolecularReaction<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	//std::cout << "there are " << s1.size() << " and " << s2.size() << " molecules" << ". self reaction = "<<self_reaction << std::endl;

	if (dt != binding_radius_dt) {
		LOG(2, "dt has changed, recalculating binding radius....");
		binding_radius = bindingradius(rate,dt,all_species[0]->D+all_species[1]->D,-1,0);
		binding_radius2 = binding_radius*binding_radius;
		neighbourhood_search.reset(neighbourhood_search.get_low(), neighbourhood_search.get_high(), binding_radius);
	}

	Molecules &mols1 = all_species[0]->mols;
	Molecules &mols2 = all_species[1]->mols;

	neighbourhood_search.embed_points(mols1.r);
	const int n = mols2.size();

	for (int mols2_i = 0; mols2_i < n; ++mols2_i) {
		if (!mols2.alive[mols2_i]) continue;
		if (mols2.saved_index[mols2_i] == SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE) continue;
		const Vect3d pos2 = mols2.r[mols2_i];
		const int id2 = mols2.id[mols2_i];
		std::vector<int>& neighbrs_list = neighbourhood_search.find_broadphase_neighbours(pos2);
		BOOST_FOREACH(int mols1_i, neighbrs_list) {
			if (!mols1.alive[mols1_i]) continue;
			if (mols1.saved_index[mols1_i] == SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE) continue;
			const Vect3d pos1 = mols1.r[mols1_i];
			const int id1 = mols1.id[mols1_i];
			if (self_reaction && (id1==id2)) continue;
			if ((pos2-neighbourhood_search.correct_position_for_periodicity(pos2, pos1)).squaredNorm() < binding_radius2) {
				BOOST_FOREACH(ReactionComponent component, products) {
					for (int i = 0; i < component.multiplier; ++i) {
						component.species->mols.add_molecule(0.5*(pos1+pos2));
					}
				}
				mols1.mark_for_deletion(mols1_i);
				mols2.mark_for_deletion(mols2_i);
			}
		}
	}
	mols1.delete_molecules();
	mols2.delete_molecules();
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

template<typename T>
BiMolecularReaction<T>::BiMolecularReaction(const double rate, const ReactionEquation& eq, const double dt, Vect3d low, Vect3d high, Vect3b periodic):
		Reaction(rate),
		products(eq.rhs),
		binding_radius_dt(dt),
		neighbourhood_search(low,high,periodic) {
	if (eq.lhs.size() == 1) {
		CHECK(eq.lhs[0].multiplier == 2, "Reaction equation is not bimolecular!");
		this->add_species(*(eq.lhs[0].species));
		this->add_species(*(eq.lhs[0].species));
	} else {
		CHECK((eq.lhs.size()==2) && (eq.lhs[0].multiplier == 1) && (eq.lhs[1].multiplier == 1), "Reaction equation is not bimolecular!");
		this->add_species(*(eq.lhs[0].species));
		this->add_species(*(eq.lhs[1].species));
	}

	if (all_species[0] == all_species[1]) {
		self_reaction = true;
	} else {
		self_reaction = false;
	}

	binding_radius = bindingradius(rate,dt,all_species[0]->D+all_species[1]->D,-1,0);
	binding_radius2 = binding_radius*binding_radius;

	LOG(2,"Binding radius for reaction "<<eq<<" with rate = "<<rate<<" is calculated as r = "<<binding_radius);
	LOG(2,"Ratio between rms diffusion step (D1+D2) and binding radius is rms/br = "<<sqrt(2.0*(all_species[0]->D+all_species[1]->D)*dt)/binding_radius);

	neighbourhood_search.reset(neighbourhood_search.get_low(), neighbourhood_search.get_high(), binding_radius);
};

template class BiMolecularReaction<BucketSort>;


void ZeroOrderMolecularReaction::add_species(Species& s) {
	Reaction::add_species(s);
	rates.push_back(rate);
}

void ZeroOrderMolecularReaction::add_species(const double _rate, Species& s) {
	Reaction::add_species(s);
	rates.push_back(_rate);
}

void ZeroOrderMolecularReaction::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	boost::uniform_real<> p_dist(0,1);
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > N(generator, p_dist);
	const Vect3d max_minus_min = max-min;
	const double volume = max_minus_min.squaredNorm();
	const int n_s = all_species.size();
	for (int i_s = 0; i_s < n_s; ++i_s) {
		boost::poisson_distribution<> p_dist(rates[i_s]*volume*dt);
		boost::variate_generator<base_generator_type&, boost::poisson_distribution<> > P(generator, p_dist);

		Molecules &mols = all_species[i_s]->mols;
		const int num_created = P();
		for (int i = 0; i < num_created; ++i) {
			const Vect3d new_position = min + max_minus_min.cwiseProduct(Vect3d(N(),N(),N()));
			mols.add_molecule(new_position);
		}
	}
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}


std::ostream& operator <<(std::ostream& out, ZeroOrderMolecularReaction& r) {
	out << "\tZero order molecular Reaction with reactants:";
	const int n = r.all_species.size();
	for (int i = 0; i < n; ++i) {
		out << "\t\tid = "<<r.all_species[i]->id<<" (rate = "<<r.rates[i]<<")";
	}
	return out;
}


}

