/*
 * Operator.h
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
 *  Created on: 30 Oct 2012
 *      Author: robinsonm
 */

#ifndef OPERATOR_H_
#define OPERATOR_H_

#include "Species.h"
#include <boost/timer.hpp>
#include <initializer_list>

namespace Tyche {

class Operator {
public:
	Operator();
	virtual ~Operator() {};
	bool add_species(Species &s);

	void operator()(const double dt);
	std::string get_time_string() const;
	std::string get_global_time() const;
	std::string get_time_percentage() const;
	void reset();
	int get_species_index(Species& s);
	double get_time() const {return time;}
	const std::vector<Species*>& get_species() const {return all_species;}
	friend std::ostream& operator<<( std::ostream& out, const Operator& b ) {
		b.print(out);
		return out;
	}

protected:
	virtual void add_species_execute(Species &s);
	virtual void reset_execute();
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) const;

private:
	void resume_timer();
	void stop_timer();

	boost::timer timer;
	double total_time;
	static double total_global_time;
	double time;
	//boost::timer global_timer;
	std::vector<Species*> all_species;
};

class CountMolsOnGrid: public Operator {
public:
	CountMolsOnGrid() {}
	void operator()(const double dt);
};

std::ostream& operator<< (std::ostream& out, CountMolsOnGrid& b);


class OperatorList: public Operator {
public:
	OperatorList() {}
	OperatorList(Operator& o) {
		list.push_back(&o);
	}
	OperatorList(std::initializer_list<Operator*> arg) {
		for (auto i: arg) {
			list.push_back(i);
			for (auto s: i->get_species()) {
				add_species(*s);
			}
		}

	}
	virtual ~OperatorList() {}

	void push_back(Operator* const i) {
		list.push_back(i);
		for (auto s: i->get_species()) {
			add_species(*s);
		}
	}
	void push_back(const OperatorList& i) {
		for (Operator* j: i.list) {
			list.push_back(j);
			for (auto s: j->get_species()) {
				add_species(*s);
			}
		}
	}

//	OperatorList& operator+=(const OperatorList &rhs) {
//		for (Operator* j: rhs.list) {
//			list.push_back(j);
//			for (auto s: j->get_species()) {
//				add_species(*s);
//			}
//		}
//		return *this;
//	}
//	OperatorList& operator+=(Operator &rhs) {
//		list.push_back(&rhs);
//		for (auto s: rhs.get_species()) {
//			add_species(*s);
//		}
//		return *this;
//	}
	std::vector<Operator*>& get_operators() {return list;}

protected:
	virtual void print(std::ostream& out) const {
		out << "List of "<<list.size()<< " operators: ("<<get_time_string()<<")"<< std::endl;
		for (auto i : list) {
			out << "\t" << *i << " ("<<i->get_time_string()<<")"<<std::endl;
		}
	}
	virtual void integrate(const double dt) {
		for (auto i : list) {
			i->operator ()(dt);
		}
	}
	virtual void reset_execute() {
		for (auto i : list) {
			i->reset();
		}
	}

	virtual void add_species_execute(Species &s) {
		for (auto i : list) {
			i->add_species(s);
		}
	}

	std::vector<Operator*> list;
};


//OperatorList operator+(Operator& arg1, Operator& arg2);
//OperatorList operator+(Operator& arg1, OperatorList& arg2);
//OperatorList operator+(OperatorList& arg1, Operator& arg2);
//OperatorList operator+(OperatorList& arg1, OperatorList& arg2);


}
#endif /* OPERATOR_H_ */
