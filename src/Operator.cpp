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
#include "Log.h"

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

int Operator::get_species_index(Species& s) {
	const int n = all_species.size();
	for (int i = 0; i < n; ++i) {
		if (all_species[i] == &s) {
			return i;
		}
	}
	ERROR("did not find species index");
	return -1;
}

bool Operator::add_species(Species& s) {
	for (auto i: all_species) {
		if (&s == i) {
			return false;
		}
	}
	add_species_execute(s);
	all_species.push_back(&s);
	return true;
}


void Operator::resume_timer() {
//	timer.resume();
//	global_timer.resume();
	timer.restart();
}

void Operator::operator ()(const double dt) {
	Operator::resume_timer();
	LOG(2, "Starting Operator: " << *this);

	integrate(dt);

	time += dt;
	LOG(2, "Stopping Operator: " << *this);
	Operator::stop_timer();
}

void Operator::reset() {
	reset_execute();
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

std::string Operator::get_time_string() {
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

void Operator::add_species_execute(Species &s) {

}
void Operator::reset_execute() {

}
void Operator::integrate(const double dt) {

}
void Operator::print(std::ostream& out) const {
	out << "Default Operator";
}



OperatorList operator+(Operator& arg1, Operator& arg2) {
	OperatorList result;
	result.push_back(&arg1);
	result.push_back(&arg2);
	return result;
}
OperatorList operator+(Operator& arg1, OperatorList& arg2) {
	arg2.push_back(&arg1);
	return arg2;
}
OperatorList operator+(OperatorList& arg1, Operator& arg2) {
	arg1.push_back(&arg2);
	return arg1;
}
OperatorList operator+(OperatorList& arg1, OperatorList& arg2) {
	arg1.push_back(&arg2);
	return arg1;
}

}



