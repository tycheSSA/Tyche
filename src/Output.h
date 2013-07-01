/*
 * Output.h
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
 *  Created on: Nov 18, 2012
 *      Author: mrobins
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <fstream>
#include <string>
#include <iomanip>
#include <map>

#include "Operator.h"

namespace Tyche {



class Output: public Operator {
public:
	typedef std::map<std::string, std::vector<double> > DataType;
	Output(const std::string filename, const double execute_dt):filename(filename),execute_dt(execute_dt) {
		next_execute = get_time() + execute_dt;
	}
	~Output() {
		f.close();
	}

	void write();
	void write(double time);
	void write(std::string filename);
	bool is_execute_time() {
		if (get_time() > next_execute) {
			next_execute = get_time() + execute_dt;
			return true;
		}
		return false;
	}
	const std::vector<double>& get_data(const std::string name) {
		ASSERT(data.find(name)!=data.end(), "cannot find name");
		return data[name];
	}
protected:

	virtual void reset_execute() {
		next_execute = get_time() + execute_dt;
	}

	DataType data;
	double execute_dt;
	double next_execute;
	std::string filename;
	std::ofstream f;
};


class OutputConcentrations: public Output {
public:
	OutputConcentrations(const std::string filename, const double dt, const StructuredGrid& grid):
		Output(filename, dt),grid(grid) {
		//std::vector<double> tmp;
		data.insert(std::pair<std::string,std::vector<double> >("x",std::vector<double>()));
		data.insert(std::pair<std::string,std::vector<double> >("y",std::vector<double>()));
		data.insert(std::pair<std::string,std::vector<double> >("z",std::vector<double>()));
		data.insert(std::pair<std::string,std::vector<double> >("Concentration(M)",std::vector<double>()));
		data.insert(std::pair<std::string,std::vector<double> >("Concentration(C)",std::vector<double>()));
		data.insert(std::pair<std::string,std::vector<double> >("Concentration",std::vector<double>()));
	}

protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) {
		out << "\tOutput concentrations";
	}

	const StructuredGrid& grid;
};

class OutputSumConcentrations: public Output {
public:
	OutputSumConcentrations(const std::string filename, const double dt, const StructuredGrid& grid):
		Output(filename, dt),grid(grid) {
		data.insert(std::pair<std::string,std::vector<double> >("Time",std::vector<double>()));
		data.insert(std::pair<std::string,std::vector<double> >("Concentration(M)",std::vector<double>()));
		data.insert(std::pair<std::string,std::vector<double> >("Concentration(C)",std::vector<double>()));
		data.insert(std::pair<std::string,std::vector<double> >("Concentration",std::vector<double>()));
	}


protected:
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) {
		out << "\tOutput sum concentrations";
	}
	virtual void add_species_execute(Species &s);

	const StructuredGrid& grid;
};


template<typename T>
class OutputCompareWithFunction: public Output {
public:
	OutputCompareWithFunction(const std::string filename,
			const double calc_dt, const double start_time, const int average_over,
			const StructuredGrid& grid, T calc_error):
				Output(filename, calc_dt),
				grid(grid),
				start_time(start_time),
				average_over(average_over),
				calc_error(calc_error) {
		data.insert(std::pair<std::string,std::vector<double> >("Time",std::vector<double>()));
		data.insert(std::pair<std::string,std::vector<double> >("Error",std::vector<double>()));
	}


	void set_param(const std::string name, const double value);


protected:
	virtual void reset_execute();
	virtual void integrate(const double dt);
	virtual void print(std::ostream& out) {
		out << "\tCompare concentration with exact function";
	}
private:
	const double start_time;
	const int average_over;
	std::vector<double> errors;
	const StructuredGrid& grid;
	T calc_error;
	std::map<std::string, double> params;
};


template<typename T>
struct rmsError {
	rmsError(T function):function(function) {}
	double operator()(const double time, std::vector<double>& concentrations, std::vector<Vect3d>& positions, std::vector<double>& volumes) {
		const int n = concentrations.size();
		ASSERT((n == positions.size()) && (n == volumes.size()), "all vectors assumed to be equal size");
		double sum,error = 0;
		for (int i = 0; i < n; ++i) {
			sum += volumes[i];
			error += pow(concentrations[i]-function(time,positions[i]),2)*volumes[i];
		}
		return sqrt(error/sum);
	}
	T function;
};

template<typename T>
struct ScaledVarianceError {
	ScaledVarianceError(T function):function(function) {}
	double operator()(const double time, std::vector<double>& concentrations, std::vector<Vect3d>& positions, std::vector<double>& volumes) {
		const int n = concentrations.size();
		ASSERT((n == positions.size()) && (n == volumes.size()), "all vectors assumed to be equal size");
		double error = 0;
		for (int i = 0; i < n; ++i) {
			const double chat = function(time,positions[i]);
			error += pow(concentrations[i]-chat,2)*volumes[i]/chat;
		}
		return error;
	}
	T function;
};

template<typename T>
struct halfSumError {
	halfSumError(T function):function(function) {}
	double operator()(const double time, std::vector<double>& concentrations, std::vector<Vect3d>& positions, std::vector<double>& volumes) {
		const int n = concentrations.size();
		ASSERT((n == positions.size()) && (n == volumes.size()), "all vectors assumed to be equal size");
		double error = 0;
		double halfsum = 0;
		double fullsum = 0;
		for (int i = 0; i < n; ++i) {
			const double expectedN = function(time,positions[i])*volumes[i];
			fullsum += expectedN;
			if (positions[i][0] > 0.5) {
				halfsum += expectedN;
				error += concentrations[i]*volumes[i];
			}
		}
		return (error - halfsum)/fullsum;
	}
	T function;
};
}
#include "Output.impl.h"


#endif /* OUTPUT_H_ */
