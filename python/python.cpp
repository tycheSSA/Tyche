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
//#include <python2.7/Python.h>
#include <numpy/arrayobject.h>


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

//std::auto_ptr<Operator> new_reaction(const double rate, const boost::python::list& reactants, const boost::python::list& products) {
//	using boost::python;
//	std::auto_ptr<Operator> result;
//	const int nr = len(reactants);
//	if (nr == 0) {
//		result = ZeroOrderMolecularReaction::New(rate,min,max);
//	} else if (nr == 1) {
//
//	} else if (nr == 2) {
//
//	} else {
//		ERROR("Number of reactants must be 0, 1 or 2");
//	}
//	return result;
//}

std::auto_ptr<Operator> new_zero_reaction(const double rate, const boost::python::list& min,const boost::python::list& max) {
	Vect3d cmin,cmax;
	CHECK(len(min)==3,"length of min should be 3");
	CHECK(len(max)==3,"length of max should be 3");
	for (int i = 0; i < 3; ++i) {
		cmin[i] = extract<double>(min[i]);
		cmax[i] = extract<double>(max[i]);
	}
	return ZeroOrderMolecularReaction::New(rate,cmin,cmax);
}

std::auto_ptr<Operator> new_uni_reaction(const double rate, Species* s, const boost::python::list& products,const double init_radius=0.0) {
	ReactionSide rhs;
	const int np = len(products);
	for (int i = 0; i < np; ++i) {
		Species *s = extract<Species*>(products[i]);
		rhs = rhs + *s;
	}
	return UniMolecularReaction::New(rate,*s >> rhs,init_radius);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(new_uni_reaction_overloads, new_uni_reaction, 3, 4);

std::auto_ptr<Operator> new_bi_reaction(const double rate, Species* s1, Species* s2, const boost::python::list& products,
					const double binding,
					const double unbinding,
					const double dt,
					const boost::python::list& min, const boost::python::list& max, const boost::python::list& periodic,
					const bool reversible=false) {
	Vect3d cmin,cmax;
	Vect3b cperiodic;
	CHECK(len(min)==3,"length of min should be 3");
	CHECK(len(max)==3,"length of max should be 3");
	CHECK(len(periodic)==3,"length of periodic should be 3");
	for (int i = 0; i < 3; ++i) {
		cmin[i] = extract<double>(min[i]);
		cmax[i] = extract<double>(max[i]);
		cperiodic[i] = extract<bool>(periodic[i]);

	}
	ReactionSide rhs;
	const int np = len(products);
	for (int i = 0; i < np; ++i) {
		Species *s = extract<Species*>(products[i]);
		rhs = rhs + *s;
	}
	return BiMolecularReaction<BucketSort>::New(rate,*s1 + *s2 >> rhs,binding,unbinding,dt,cmin,cmax,cperiodic,reversible);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(new_bi_reaction_overloads, new_bi_reaction, 10, 11);


void Species_fill_uniform2(Species& self, const boost::python::list& min, const boost::python::list& max, const unsigned int n) {
	Vect3d cmin,cmax;
	CHECK(len(min)==3,"length of min should be 3");
	CHECK(len(max)==3,"length of max should be 3");
	for (int i = 0; i < 3; ++i) {
		cmin[i] = extract<double>(min[i]);
		cmax[i] = extract<double>(max[i]);
	}
	self.fill_uniform(cmin,cmax,n);
}


std::vector<double>* Species_get_concentration1(Species& self, const boost::python::list& min, const boost::python::list& max, const boost::python::list& n) {
	std::vector<double>* result = new std::vector<double>();
	Vect3d cmin,cmax;
	Vect3i cn;
	CHECK(len(min)==3,"length of min should be 3");
	CHECK(len(max)==3,"length of max should be 3");
	CHECK(len(n)==3,"length of n should be 3");

	for (int i = 0; i < 3; ++i) {
		cmin[i] = extract<double>(min[i]);
		cmax[i] = extract<double>(max[i]);
		cn[i] = extract<int>(n[i]);

	}
	self.get_concentration(cmin,cmax,cn,*result);
	return result;
//	std::cout << "finished calcing concen n = " << result.size()<<std::endl;
//	npy_intp dimens[3] = {cn[0], cn[1],cn[2]};
//
//	std::cout << "dims are " << cn[0] <<' '<< cn[1]<<' ' << Py_IsInitialized()<<std::endl;
//
//	PyObject* array = PyArray_SimpleNew(3, dimens, NPY_DOUBLE);
//	std::cout << "dims are " << cn[0] <<' '<< cn[1]<<' ' << cn[2]<<std::endl;

////	numarray.resize(cn[0],cn[1],cn[2]);
//	for (long int i = 0; i < cn[0]; ++i) {
//		for (long int j = 0; j < cn[1]; ++j) {
//			for (long int k = 0; k < cn[2]; ++k) {
//				std::cout << "i j k = "<<i<<' '<<j<<' '<<k<<std::endl;
//				const int ii = i*cn[1]*cn[2] + j*cn[2] + k;
//				double* data = (double *)PyArray_GETPTR3(array, i,  j,  k);
//				*data = result[ii];
////				numarray[make_tuple(i,j,k)] = result[ii];
//			}
//		}
//	}
//	std::cout << "finished copythig results" << std::endl;
//	return array;
////	double *arr = new double[result.size()];
////	std::copy(result.begin(), result.end(), arr);
////
////	object obj(handle<>(PyArray_SimpleNewFromData(3, dimens, PyArray_DOUBLE, arr)));
////	return extract<numeric::array>(obj);
////	object obj(handle<>(array));
////	numeric::array numarray = extract<numeric::array>(obj);
////	return numarray;
}

void Species_get_concentration2(Species& self, const boost::python::list& min, const boost::python::list& max, const boost::python::list& n, std::vector<double>& result) {
	Vect3d cmin,cmax;
	Vect3i cn;
	CHECK(len(min)==3,"length of min should be 3");
	CHECK(len(max)==3,"length of max should be 3");
	CHECK(len(n)==3,"length of n should be 3");

	for (int i = 0; i < 3; ++i) {
		cmin[i] = extract<double>(min[i]);
		cmax[i] = extract<double>(max[i]);
		cn[i] = extract<int>(n[i]);

	}
	self.get_concentration(cmin,cmax,cn,result);
}

void    (Species::*fill_uniform1)(const int) = &Species::fill_uniform;




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

	class_<Species,typename std::auto_ptr<Species> >("Species",boost::python::init<double>())
			.def("fill_uniform",fill_uniform1)
			.def("fill_uniform",Species_fill_uniform2)
			.def("get_concentration",Species_get_concentration1, return_value_policy<manage_new_object>())
			.def("get_concentration",Species_get_concentration2)
			.def(self_ns::str(self_ns::self))
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
			.def(self_ns::str(self_ns::self))
			;



	/*
	 * Geometry
	 */
	def("new_xplane",xplane::New);
	def("new_yplane",yplane::New);
	def("new_zplane",zplane::New);

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
    def("new_zero_reaction",new_zero_reaction);
    def("new_uni_reaction",new_uni_reaction, new_uni_reaction_overloads());
    def("new_bi_reaction",new_bi_reaction, new_bi_reaction_overloads());


}

}
