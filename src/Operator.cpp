/*
 * Operator.cpp
 *
 *  Created on: 30 Oct 2012
 *      Author: robinsonm
 */

#include "Operator.h"
#include <sstream>
#include <string>
#include <algorithm>

namespace Tyche {
//boost::timer::cpu_timer Operator::global_timer;
double Operator::total_global_time = 0;

Operator::Operator() {
//	timer.stop();
//	global_timer.stop();
	total_time = 0;
	time = 0;
	all_species.clear();

}

void Operator::add_species(Species& s) {
	all_species.push_back(&s);
}

void Operator::resume_timer() {
//	timer.resume();
//	global_timer.resume();
	timer.restart();
}

void Operator::operator ()(const double dt) {
   time += dt;
}

void Operator::reset() {
	time = 0;
}

void Operator::stop_timer() {
	//timer.stop();
	//global_timer.stop();
	const double time = timer.elapsed();
	total_time += time;
	total_global_time += time;
}


template <typename T>
const std::string to_string(const T& data)
{
   std::ostringstream conv;
   conv << data;
   return conv.str();
}

std::string Operator::get_time() {
	//return timer.format();
	return "Time to execute: " + to_string(total_time) + " s (" + get_time_percentage() + ")";
}

std::string Operator::get_global_time() {
	//return global_timer.format();
	return "Time to execute all Operators: " + to_string(total_global_time) + " s";
}

std::string Operator::get_time_percentage() {
//	boost::timer::cpu_times percent;
//	boost::timer::cpu_times this_time = timer.elapsed();
//	boost::timer::cpu_times total = global_timer.elapsed();
//	percent.wall = total.wall / this_time.wall;
//	percent.user = total.user / this_time.user;
//	percent.system = total.system / this_time.system;
//	return boost::timer::format(percent);

	const double percent = 100*total_time/total_global_time;
	return to_string(percent) + "%";
}

void CountMolsOnGrid::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);
	const int n = all_species.size();
	for (int i = 0; i < n; ++i) {
		Species &s = *(all_species[i]);
		std::fill(s.mol_copy_numbers.begin(),s.mol_copy_numbers.end(),0);
		BOOST_FOREACH(Vect3d &r,s.mols.r) {
			if (s.grid.is_in(r)) {
				s.mol_copy_numbers[s.grid.get_cell_index(r)]++;
			}
		}
	}
}

std::ostream& operator<< (std::ostream& out, CountMolsOnGrid& b) {
	return out << "\tCount Molecules in Grid";
}


}



