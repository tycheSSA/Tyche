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
#include <boost/math/tools/roots.hpp>

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

UniMolecularReaction::UniMolecularReaction(const double rate,const ReactionEquation& eq, const double init_radius):Reaction(rate) {
	CHECK((eq.lhs.size()==1) && (eq.lhs[0].multiplier == 1), "Reaction equation is not unimolecular!");
	this->add_species(*(eq.lhs[0].species));
	product_list.push_back(eq.rhs);
	probabilities.push_back(0);
	init_radii.push_back(init_radius);
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
			const double angle = 2.0*PI*uni();
			const double angle2 = 2.0*PI*uni();
			BOOST_FOREACH(ReactionComponent component, products) {
				for (int j = 0; j < component.multiplier; ++j) {
					if (j==0) {
						Vect3d new_pos(init_radii[0]*sin(angle)*cos(angle2),
									   init_radii[0]*sin(angle)*sin(angle2),
									   init_radii[0]*cos(angle));
						component.species->mols.add_molecule(mols.r[i]);
					} else {
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

//class reversiblef {
//public:
//	reversiblef(const double alpha, const double kappa):alpha(alpha),kappa(kappa) {}
//    double operator()(const double beta) {
//        return 4.0*PI*alpha*(beta - std::tanh(beta)) / (std::cosh(beta-beta*alpha)*std::tanh(beta) - std::sinh(beta-beta*alpha)) - kappa;
//    }
//    const double alpha;
//    const double kappa;
//};
//
//class irreversiblef {
//public:
//	irreversiblef(const double kappa):kappa(kappa) {}
//    double operator()(const double beta) {
//        return 4.0*PI*(1.0-std::tanh(beta)/beta) - kappa;
//    }
//    const double kappa;
//};

class reversiblef1 {
public:
	reversiblef1(const double unbinding,
			const double rate,
			const double lambda,
			const double difc) {
		ub = unbinding;
		lhs = rate/(4.0*PI*ub*difc);
		s = std::sqrt(lambda/difc);
	}
    double operator()(const double b) {
//        return 4.0*PI*alpha*(beta - std::tanh(beta)) / (std::cosh(beta-beta*alpha)*std::tanh(beta) - std::sinh(beta-beta*alpha)) - kappa;
    	const double bs = b*s;
    	const double bubs = (b-ub)*s;
    	return (bs - std::tanh(bs)) / (std::cosh(bubs)*std::tanh(bs) - std::sinh(bubs)) - lhs;
    }
    double lhs;
    double s;
    double ub;

};

class reversiblef {
public:
	reversiblef(const double rate,
			const double unbinding,
			const double binding,
			const double d_plus_d) {
		ub = unbinding;
		b = binding;
		alpha = ub/b;
		difc = d_plus_d;
		kappa = rate;
	}
    double operator()(const double lambda) {
    	const double beta = b * std::sqrt(lambda / difc);
    	//LOG(1,"beta = "<<beta<<" beta - tanh(beta) = "<<beta - std::tanh(beta) << "top = "<<4.0*PI*ub*difc * (beta - std::tanh(beta))<<" bottom = "<<std::cosh(beta - beta*alpha)*std::tanh(beta) - std::sinh(beta-beta*alpha));
    	//LOG(1,"arg = "<<beta-beta*alpha<<" cosh = "<<std::cosh(beta - beta*alpha)<<" tanh  = "<<std::tanh(beta)<<" sinh = "<<std::sinh(beta-beta*alpha));
    	return 4.0*PI*ub*difc * (beta - std::tanh(beta)) - (std::cosh(beta - beta*alpha)*std::tanh(beta) - std::sinh(beta-beta*alpha))*kappa;
    }
    double kappa;
    double b;
    double ub;
    double alpha;
    double difc;

};

class irreversiblef {
public:
	irreversiblef( const double kappa_constant):kappa_constant(kappa_constant) {}
    double operator()(const double beta) {
    	const double kappa = kappa_constant/beta;
        return 4.0*PI*(1.0-std::tanh(beta)/beta) - kappa;
    }
    const double kappa_constant;


};

template<typename T>
void BiMolecularReaction<T>::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);

//	if (dt != binding_radius_dt) {
//		LOG(2, "dt has changed, recalculating binding radius....");
//		const double difc = all_species[0]->D + all_species[1]->D;
//		/*
//		 * try to set lambda with unbinding radius = 0
//		 */
//		const double kapa = rate / (binding_radius*difc);
//
//		CHECK(kapa < 0.1, "failed to set bimolecular parameters with unbinding radius = 0");
//
//
//
//		neighbourhood_search.reset(neighbourhood_search.get_low(), neighbourhood_search.get_high(), binding_radius);
//	}

	Molecules &mols1 = all_species[0]->mols;
	Molecules &mols2 = all_species[1]->mols;

	boost::uniform_real<> uni_dist(0.0,1.0);
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

	const double binding_radius2 = binding_radius*binding_radius;
	const double dt_lambda = dt*lambda;
	ASSERT(dt_lambda <= 1.0, "probabilty of reaction greater than 1.0. reduce dt");
	LOG(2,"binding radius = " <<binding_radius << " dt*lambda = " << dt_lambda);

	neighbourhood_search.embed_points(mols1.r);
	const int n = mols2.size();
	for (int mols2_i = 0; mols2_i < n; ++mols2_i) {
		if (!mols2.alive[mols2_i]) continue;
		if (mols2.saved_index[mols2_i] == SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE) continue;
		const Vect3d pos2 = mols2.r[mols2_i];
		const int id2 = mols2.id[mols2_i];
		std::vector<int>& neighbrs_list = neighbourhood_search.find_broadphase_neighbours(pos2);
		for (auto mols1_i : neighbrs_list) {
			if (!mols1.alive[mols1_i]) continue;
			if (mols1.saved_index[mols1_i] == SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE) continue;
			const Vect3d pos1 = mols1.r[mols1_i];
			const int id1 = mols1.id[mols1_i];
			if (self_reaction && (id1==id2)) continue;
			if ((pos2-neighbourhood_search.correct_position_for_periodicity(pos2, pos1)).squaredNorm() < binding_radius2) {
				if (uni() < dt_lambda) {
					for (auto component : products) {
						for (int i = 0; i < component.multiplier; ++i) {
							component.species->mols.add_molecule(0.5*(pos1+pos2));
						}
					}
					mols1.mark_for_deletion(mols1_i);
					mols2.mark_for_deletion(mols2_i);
				}
			}
		}
	}
	mols1.delete_molecules();
	mols2.delete_molecules();
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

template<typename T>
void BiMolecularReaction<T>::report_dt_suitability(const double dt) {
	LOG(1,"probability of reaction (per timestep) = "<<
			     dt*lambda<<". ratio of diffusion step to binding radius = "<<
				std::sqrt(2.0*(all_species[0]->D + all_species[1]->D)*dt)/binding_radius);
}

void UniMolecularReaction::report_dt_suitability(const double dt) {
	calculate_probabilities(dt);
	LOG(1,"probability of reaction (per timestep) = "<<
			     total_probability);
}

template<typename T>
BiMolecularReaction<T>::BiMolecularReaction(const double rate, const ReactionEquation& eq,
			const double binding,
			const double unbinding,
			Vect3d low, Vect3d high, Vect3b periodic,
			const bool reversible):
		Reaction(rate),
		products(eq.rhs),
		binding_radius_dt(0),
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

	const double difc = all_species[0]->D + all_species[1]->D;

	double lambda_min =  std::sqrt(std::numeric_limits<double>::min());
	double lambda_max = 0.1*difc * pow(std::numeric_limits<double>::max_exponent,2) / pow(binding,2);
//	double binding_min = binding_max * 0.000001;
	LOG(1,"lambda_min = "<<lambda_min<<" and lambda_max = "<<lambda_max);
	unsigned int maximum_iterations = 50000;
	boost::uintmax_t max_iter = maximum_iterations;
	boost::math::tools::eps_tolerance<double> tol(60);
	std::pair<double, double> r;
	if (reversible) {
		reversiblef fop(rate,unbinding,binding,difc);

		CHECK(fop(lambda_min)*fop(lambda_max) < 0, "brackets of root not valid. f(lambda_min) = "<<fop(lambda_min)<<" and f(lambda_max) = "<<fop(lambda_max));
		while (std::isinf(fop(lambda_max))) {
			const double middle = 0.5*(lambda_max-lambda_min);
			if (fop(lambda_min)*fop(middle) < 0) {
				lambda_max = middle;
			} else {
				lambda_min = middle;
			}
		}
		LOG(2,"solving root with f(lambda_min) = "<<fop(lambda_min)<<" and f(lambda_max) = "<<fop(lambda_max));

		r = boost::math::tools::toms748_solve(
									fop,
									lambda_min, lambda_max, tol, max_iter);
		LOG(1,"f(binding) = f("<<0.5*(r.first+r.second)<<") = "<<fop(0.5*(r.first + r.second)));
		LOG(1,"f(first) = f("<<r.first<<") = "<<fop(r.first));
		LOG(1,"f(second) = f("<<r.second<<") = "<< fop(r.second));

	} else {
//		irreversiblef f(kappa_constant);
//		CHECK(f(beta_min)*f(beta_max) < 0, "brackets of root not valid. f(beta_min) = "<<f(beta_min)<<" and f(beta_max) = "<<f(beta_max));
//		while (std::isinf(f(beta_max))) {
//			const double middle = 0.5*(beta_max-beta_min);
//			if (f(beta_min)*f(middle) < 0) {
//				beta_max = middle;
//			} else {
//				beta_min = middle;
//			}
//		}
//		LOG(2,"solving root with f(beta_min) = "<<f(beta_min)<<" and f(beta_max) = "<<f(beta_max));
//		r = boost::math::tools::toms748_solve(
//									f,
//									beta_min, beta_max, tol, max_iter);
	}
	CHECK(max_iter < maximum_iterations, "could not solve for root. binding_min = "<<r.first<<" and binding_max = "<<r.second);
	lambda = 0.5*(r.first + r.second);
	binding_radius = binding;
	binding_radius_dt = binding_radius;
	unbinding_radius = unbinding;
	const double alpha = unbinding/binding_radius;
	const double beta = binding_radius*std::sqrt(lambda/difc);
	CHECK(alpha < 1.0, "unbinding radius must be less than binding radius = "<<binding_radius);

	//lambda = pow(beta / binding_radius,2)*difc;
	double germinate;
	if (reversible) {
		germinate = 1.0 - std::sinh(alpha*beta) / (alpha*beta*std::cosh(beta));
	} else {
		germinate = 0.0;
	}
	LOG(1,"created bimolecular reaction with eq: " << eq <<" binding radius = " << binding_radius <<" unbinding radius = "<<unbinding_radius<< " lambda = " << lambda <<" and germinate recombination probability = " << germinate);
	LOG(1,"kappa  = " << rate / (binding * difc) << " beta = "<< binding * std::sqrt(lambda / difc));

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

