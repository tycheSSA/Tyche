/* 
 * python.cpp
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of Tyche.
 *
 * Tyche is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Tyche is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PDE_BD.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Feb 9, 2013
 *      Author: mrobins
 */


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "Tyche.h"


using namespace boost::python;

namespace Tyche {

void python_init(boost::python::list& py_argv) {
	using boost::python::len;

	int argc = len(py_argv);
	char **argv = new char*[len(py_argv)];
	for (int i = 0; i < len(py_argv); ++i) {
		std::string this_argv = boost::python::extract<std::string>(py_argv[i]);
		argv[i] = new char[this_argv.size()];
		std::strcpy(argv[i],this_argv.c_str());
	}
	init(argc,argv);
	for (int i = 0; i < len(py_argv); ++i) {
		delete argv[i];
	}
	delete argv;
}



BOOST_PYTHON_MODULE(pyTyche) {
	def("init", python_init);

	/*
	 * vector
	 */
	class_<std::vector<double> >("Vectd")
	        .def(boost::python::vector_indexing_suite<std::vector<double> >())
	    ;

	/*
	 * Species
	 */
	void    (Species::*fill_uniform1)(const int) = &Species::fill_uniform;
	void    (Species::*fill_uniform2)(const Vect3d, const Vect3d, const unsigned int) = &Species::fill_uniform;
	void	(Species::*get_concentration1)(const Vect3d, const Vect3d, const Vect3i,
			std::vector<double>&) const = &Species::get_concentration;
	class_<Species,typename std::auto_ptr<Species> >("Species",boost::python::init<double>())
			.def("fill_uniform",fill_uniform1)
			.def("fill_uniform",fill_uniform2)
			.def("get_concentration",get_concentration1)
			;
	def("new_species",Species::New);

   /*
    * Operator
    */
	class_<Operator,typename std::auto_ptr<Operator> >("Operator",boost::python::no_init)
			.def("integrate",&Operator::operator())
			.def("reset",&Operator::reset)
			.def("add_species",&Operator::add_species)
			.def("get_species_index",&Operator::get_species_index)
			;

	/*
	 * Boundaries
	 */
    def("new_reflective_boundary",ReflectiveBoundary<xplane>::New);
    def("new_reflective_boundary",ReflectiveBoundary<yplane>::New);
    def("new_reflective_boundary",ReflectiveBoundary<zplane>::New);

    def("new_jump_boundary",JumpBoundary<xplane>::New);
    def("new_jump_boundary",JumpBoundary<yplane>::New);
    def("new_jump_boundary",JumpBoundary<yplane>::New);

    /*
     * Diffusion
     */
    def("new_diffusion",Diffusion::New);


    /*
     * Reactions
     */
    def("new_zero_reaction",ZeroOrderMolecularReaction::New);


    def("new_uni_reaction",UniMolecularReaction::New);

	def("new_bi_reaction",BiMolecularReaction<BucketSort>::New);

}

}
