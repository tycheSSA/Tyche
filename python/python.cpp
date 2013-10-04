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

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>


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
		Species* s = extract<Species*>(products[i]);
		rhs = rhs + *s;
	}
	ReactionEquation eq = *s1 + *s2 >> rhs;

	return BiMolecularReaction<BucketSort>::New(rate,eq,binding,unbinding,dt,cmin,cmax,cperiodic,reversible);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(new_bi_reaction_overloads, new_bi_reaction, 10, 11);

std::auto_ptr<Operator> new_bi_reaction2(const double rate, Species* s1, Species* s2, const boost::python::list& products,
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
		Species* s = extract<Species*>(products[i]);
		rhs = rhs + *s;
	}
	ReactionEquation eq = *s1 + *s2 >> rhs;

	return BiMolecularReaction<BucketSort>::New(rate,eq,dt,cmin,cmax,cperiodic,reversible);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(new_bi_reaction_overloads2, new_bi_reaction2, 8, 9);

std::auto_ptr<Operator> new_tri_reaction(const double rate, Species* s1, Species* s2, Species* s3,
					const boost::python::list& products,
					const double dt,
					const boost::python::list& min, const boost::python::list& max, const boost::python::list& periodic) {
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
		Species* s = extract<Species*>(products[i]);
		rhs = rhs + *s;
	}
	ReactionEquation eq = *s1 + *s2 + *s3 >> rhs;

	return TriMolecularReaction::New(rate,eq,dt,cmin,cmax,cperiodic);
}

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


std::auto_ptr<Operator> group(const boost::python::list& ops) {
	OperatorList* result = new OperatorList();
	const int n = len(ops);
	for (int i = 0; i < n; ++i) {
		Operator *s = extract<Operator*>(ops[i]);
		result->push_back(s);
	}
	return std::auto_ptr<Operator>(result);
}

//template <typename T>
//std::auto_ptr<Operator> new_jump_boundary(const T& geometry, PyObject* py_jump_by){
//	Eigen::Map<Vect3d> jump_by((double *) PyArray_DATA(py_jump_by));
//	return JumpBoundary<T>::New(geometry,jump_by);
//}

template <typename T>
std::auto_ptr<Operator> new_jump_boundary(const T& geometry, const boost::python::list& py_jump_by){
	CHECK(len(py_jump_by)==3,"length of jump_by should be 3");
	Vect3d jump_by;
	for (int i = 0; i < 3; ++i) {
		jump_by[i] = extract<double>(py_jump_by[i]);
	}
	return JumpBoundary<T>::New(geometry,jump_by);
}


template<class T>
struct vtkSmartPointer_to_python {
	static PyObject *convert(const vtkSmartPointer<T> &p) {
		std::ostringstream oss;
		oss << (void*) p.GetPointer();
		std::string address_str = oss.str();

		using namespace boost::python;
		object obj = import("vtk").attr("vtkObject")(address_str);
		return incref(obj.ptr());
	}
};


BOOST_PYTHON_MODULE(pyTyche) {
	def("init", python_init);

	/*
	 * vector
	 */
	class_<std::vector<double> >("Vectd")
	        .def(boost::python::vector_indexing_suite<std::vector<double> >())
	    ;

	to_python_converter<vtkSmartPointer<vtkUnstructuredGrid>,vtkSmartPointer_to_python<vtkUnstructuredGrid> >();

	/*
	 * Species
	 */

	class_<Species,typename std::auto_ptr<Species> >("Species",boost::python::init<double>())
			.def("fill_uniform",fill_uniform1)
			.def("fill_uniform",Species_fill_uniform2)
			.def("get_concentration",Species_get_concentration1, return_value_policy<manage_new_object>())
			.def("get_concentration",Species_get_concentration2)
			.def("get_vtk",&Species::get_vtk)
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
			.def("integrate_for_time",&Operator::integrate_for_time)
			.def(self_ns::str(self_ns::self))
			;

	def("group",group);

	/*
	 * Geometry
	 */
	def("new_xplane",xplane::New);
	def("new_yplane",yplane::New);
	def("new_zplane",zplane::New);
	def("new_xrect",xrect::New);
	def("new_yrect",yrect::New);
	def("new_zrect",zrect::New);

	class_<xplane,typename std::auto_ptr<xplane> >("Xplane",boost::python::no_init);
	class_<yplane,typename std::auto_ptr<yplane> >("Yplane",boost::python::no_init);
	class_<zplane,typename std::auto_ptr<zplane> >("Zplane",boost::python::no_init);
	class_<xplane,typename std::auto_ptr<xplane> >("Xplane",boost::python::no_init);
	class_<yplane,typename std::auto_ptr<yplane> >("Yplane",boost::python::no_init);
	class_<zplane,typename std::auto_ptr<zplane> >("Zplane",boost::python::no_init);
	class_<xrect,typename std::auto_ptr<xrect> >("Xrect",boost::python::no_init);
	class_<yrect,typename std::auto_ptr<yrect> >("Yrect",boost::python::no_init);
	class_<zrect,typename std::auto_ptr<zrect> >("Zrect",boost::python::no_init);

	/*
	 * Boundaries
	 */
    def("new_reflective_boundary",ReflectiveBoundary<xplane>::New);
    def("new_reflective_boundary",ReflectiveBoundary<yplane>::New);
    def("new_reflective_boundary",ReflectiveBoundary<zplane>::New);
    def("new_reflective_boundary",ReflectiveBoundary<xrect>::New);
    def("new_reflective_boundary",ReflectiveBoundary<yrect>::New);
    def("new_reflective_boundary",ReflectiveBoundary<zrect>::New);

//	class_<ReflectiveBoundary<xplane>,typename std::auto_ptr<ReflectiveBoundary<xplane> > >("ReflectiveBoundaryXplane",boost::python::no_init);
//	class_<ReflectiveBoundary<yplane>,typename std::auto_ptr<ReflectiveBoundary<yplane> > >("ReflectiveBoundaryYplane",boost::python::no_init);
//	class_<ReflectiveBoundary<zplane>,typename std::auto_ptr<ReflectiveBoundary<zplane> > >("ReflectiveBoundaryZplane",boost::python::no_init);


    def("new_jump_boundary",new_jump_boundary<xplane>);
    def("new_jump_boundary",new_jump_boundary<yplane>);
    def("new_jump_boundary",new_jump_boundary<zplane>);
    def("new_jump_boundary",new_jump_boundary<xrect>);
    def("new_jump_boundary",new_jump_boundary<yrect>);
    def("new_jump_boundary",new_jump_boundary<zrect>);

//	class_<JumpBoundary<xplane>,typename std::auto_ptr<JumpBoundary<xplane> > >("JumpBoundaryXplane",boost::python::no_init);
//	class_<JumpBoundary<yplane>,typename std::auto_ptr<JumpBoundary<yplane> > >("JumpBoundaryYplane",boost::python::no_init);
//	class_<JumpBoundary<zplane>,typename std::auto_ptr<JumpBoundary<zplane> > >("JumpBoundaryZplane",boost::python::no_init);


    /*
     * Diffusion
     */
    def("new_diffusion",Diffusion::New);
//	class_<Diffusion,typename std::auto_ptr<Diffusion> >("Diffusion",boost::python::no_init);


    /*
     * Reactions
     */
    def("new_zero_reaction",new_zero_reaction);
//	class_<ZeroOrderMolecularReaction,typename std::auto_ptr<ZeroOrderMolecularReaction> >("ZeroOrderMolecularReaction",boost::python::no_init);

    def("new_uni_reaction",new_uni_reaction, new_uni_reaction_overloads());
//	class_<UniMolecularReaction,typename std::auto_ptr<UniMolecularReaction> >("UniMolecularReaction",boost::python::no_init);

    def("new_bi_reaction",new_bi_reaction, new_bi_reaction_overloads());
    def("new_bi_reaction",new_bi_reaction2, new_bi_reaction_overloads2());
//	class_<BiMolecularReaction<BucketSort>,typename std::auto_ptr<BiMolecularReaction<BucketSort> > >("BiMolecularReaction",boost::python::no_init);
    def("new_tri_reaction",new_tri_reaction);


}

}
