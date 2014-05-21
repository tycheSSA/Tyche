/*
 * Visualisation.cpp
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
 *  Created on: 22 Nov 2012
 *      Author: robinsonm
 */

#include "Visualisation.h"
#include <boost/foreach.hpp>
#include "Constants.h"

#include <limits>

#include <vtkVersion.h>
#include <vtkProperty.h>
#include <vtkPlaneSource.h>
#include <vtkVertexGlyphFilter.h>
#include "vtkActor2D.h"
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRegularPolygonSource.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>

#include <vtkAxis.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>


namespace Tyche {
struct Visualisation::MyVTKdata {
	MyVTKdata() {
		renderer = vtkSmartPointer<vtkRenderer>::New();
		renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

		renderWindow->SetSize(800,600);

		renderWindow->AddRenderer(renderer);
		renderWindowInteractor->SetRenderWindow(renderWindow);

		renderer->SetBackground(.1,.2,.3); // Background color dark blue

		vtkSmartPointer<vtkAxesActor> axes =
				vtkSmartPointer<vtkAxesActor>::New();

		orientationWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
		//orientationWidget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
		orientationWidget->SetOrientationMarker( axes );
		orientationWidget->SetInteractor( renderWindowInteractor );
		//orientationWidget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
		orientationWidget->SetEnabled( 1 );
		//orientationWidget->InteractiveOn();
		camera_track = 0;

		//renderWindowInteractor->Initialize();
	}
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindow> renderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer<vtkOrientationMarkerWidget> orientationWidget;
	double camera_track;
};

Visualisation::Visualisation(const double dt):vis_dt(dt),next_vis(0),low(Vect3d(0,0,0)),high(Vect3d(0,0,0)) {
	my_vtk_data = new MyVTKdata();
}

void Visualisation::integrate(const double dt) {

	if (get_time() > next_vis) {
		LOG(2, "Starting Operator: " << *this);
		my_vtk_data->renderer->RemoveAllViewProps();
		BOOST_FOREACH(Species* s, this->get_species()) {
			add_molecules_to_vis(*s);
			add_compartments_to_vis(*s);
		}
		add_planes_to_vis();
		show_vis();
		next_vis = get_time() + vis_dt;
		LOG(2, "Stopping Operator: " << *this);
	}

}


template<unsigned int DIM>
vtkSmartPointer<vtkActor> add_plane_to_vis(AxisAlignedPlane<DIM>& plane, const double xres, const double yres, const Vect3d& high, const Vect3d& low, Visualisation::MyVTKdata* my_vtk_data) {
	static const int dim_map[][2] = {{1,2}, {0,2}, {0,1}};
	double origin[3];
	origin[DIM] = plane.get_coord();
	origin[dim_map[DIM][0]] = low[dim_map[DIM][0]];
	origin[dim_map[DIM][1]] = low[dim_map[DIM][1]];
	double p0[3];
	p0[DIM] = origin[DIM];
	p0[dim_map[DIM][0]] = high[dim_map[DIM][0]];
	p0[dim_map[DIM][1]] = origin[dim_map[DIM][1]];
	double p1[3];
	p1[DIM] = origin[DIM];
	p1[dim_map[DIM][0]] = origin[dim_map[DIM][0]];
	p1[dim_map[DIM][1]] = high[dim_map[DIM][1]];

	// Create a plane
	vtkSmartPointer<vtkPlaneSource> planeSource =
			vtkSmartPointer<vtkPlaneSource>::New();

	planeSource->SetOrigin(origin);
	planeSource->SetPoint1(p0);
	planeSource->SetPoint2(p1);
	planeSource->SetXResolution(xres);
	planeSource->SetYResolution(yres);
	planeSource->Update();

	vtkPolyData* the_plane = planeSource->GetOutput();

	// Create a mapper and actor
	vtkSmartPointer<vtkPolyDataMapper> mapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInput(the_plane);
#else
	mapper->SetInputData(the_plane);
#endif

	vtkSmartPointer<vtkActor> actor =
			vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetRepresentationToWireframe();
	actor->GetProperty()->SetLineWidth(2);
	my_vtk_data->renderer->AddActor(actor);
	return actor;
}

void Visualisation::add_planes_to_vis() {
	BOOST_FOREACH(xplane p, xplanes) {
		add_plane_to_vis(p, 1, 1, high, low, my_vtk_data);
	}
	BOOST_FOREACH(yplane p, yplanes) {
		add_plane_to_vis(p, 1, 1, high, low, my_vtk_data);
	}
	BOOST_FOREACH(zplane p, zplanes) {
		add_plane_to_vis(p, 1, 1, high, low, my_vtk_data);
	}
}


void Visualisation::add_molecules_to_vis(Species& s) {
	vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();
	BOOST_FOREACH(Vect3d r, s.mols.r) {
		points->InsertNextPoint(r[0],r[1],r[2]);
	}

	vtkSmartPointer<vtkPolyData> polydata =
			vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
			vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	vertexGlyphFilter->AddInput(polydata);
#else
	vertexGlyphFilter->AddInputData(polydata);
#endif
	vertexGlyphFilter->Update();


	vtkSmartPointer<vtkPolyDataMapper> mapper =
			vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(vertexGlyphFilter->GetOutputPort());
	mapper->Update();

	vtkSmartPointer<vtkActor> actor =
			vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
	actor->GetProperty()->SetPointSize(3);
	my_vtk_data->renderer->AddActor(actor);
}

template<unsigned int DIM>
void add_compartment_slice_to_vis(Species& s, AxisAlignedPlane<DIM>& p, const Vect3d& high, const Vect3d& low, Visualisation::MyVTKdata* my_vtk_data) {
	static const int dim_map[][2] = {{1,2}, {0,2}, {0,1}};
	const StructuredGrid *grid = dynamic_cast<const StructuredGrid*>(s.grid);
	if (grid != NULL) {
		Vect3i res = grid->get_cells_along_axes();
		vtkSmartPointer<vtkActor> actor = add_plane_to_vis(p, res[dim_map[DIM][0]], res[dim_map[DIM][1]], high, low, my_vtk_data);
		vtkDataSet* data = actor->GetMapper()->GetInput();
		const int n = data->GetNumberOfCells();
		std::vector<int> compartment_indicies;
		grid->get_slice(p, compartment_indicies);
		ASSERT(n == compartment_indicies.size(),"slice of compartments has different size to visualisation plane");
		vtkSmartPointer<vtkDoubleArray> copynumber_data = vtkSmartPointer<vtkDoubleArray>::New();
		for (int i = 0; i < n; ++i) {
			copynumber_data->InsertNextValue(s.copy_numbers[compartment_indicies[i]]);
		}
		data->GetCellData()->SetScalars(copynumber_data);
	}

}

void Visualisation::add_compartments_to_vis(Species& s) {
	BOOST_FOREACH(xplane p, compartment_slices_x) {
		add_compartment_slice_to_vis(s, p, high, low, my_vtk_data);
	}
	BOOST_FOREACH(yplane p, compartment_slices_y) {
		add_compartment_slice_to_vis(s, p, high, low, my_vtk_data);
	}
	BOOST_FOREACH(zplane p, compartment_slices_z) {
		add_compartment_slice_to_vis(s, p, high, low, my_vtk_data);
	}
}


void Visualisation::show_vis() {
	// Render and interact
//	my_vtk_data->renderWindow->Render();
//	my_vtk_data->renderWindowInteractor->Initialize();
//	const Vect3d middle = 0.5*(high-low) + low;
//	my_vtk_data->renderer->GetActiveCamera()->SetFocalPoint (middle.data());
	my_vtk_data->camera_track += PI*0.0005;
	//std::cout <<"new azimuth = " << 45*std::sin(3*my_vtk_data->camera_track) <<" camera track = " << my_vtk_data->camera_track << std::endl;
	my_vtk_data->renderer->GetActiveCamera()->Azimuth(0.2*std::cos(my_vtk_data->camera_track));
	//my_vtk_data->renderer->GetActiveCamera()->Azimuth(0.01);
	my_vtk_data->renderer->GetActiveCamera()->Elevation(0.1*std::sin(2*my_vtk_data->camera_track));
	my_vtk_data->renderer->ResetCamera();
	my_vtk_data->renderWindow->Render();
//	my_vtk_data->orientationWidget->On();
	//my_vtk_data->renderWindowInteractor->Start();
}


void Visualisation::add_geometry(const xplane& geometry) {
	if (geometry.get_coord() < low[0]) {
		low[0] = geometry.get_coord();
	}
	if (geometry.get_coord() > high[0]) {
		high[0] = geometry.get_coord();
	}
	xplanes.push_back(geometry);
}
void Visualisation::add_geometry(const yplane& geometry) {
	if (geometry.get_coord() < low[1]) {
		low[1] = geometry.get_coord();
	}
	if (geometry.get_coord() > high[1]) {
		high[1] = geometry.get_coord();
	}
	yplanes.push_back(geometry);
}
void Visualisation::add_geometry(const zplane& geometry) {
	if (geometry.get_coord() < low[2]) {
		low[2] = geometry.get_coord();
	}
	if (geometry.get_coord() > high[2]) {
		high[2] = geometry.get_coord();
	}
	zplanes.push_back(geometry);
}


void Visualisation::add_compartment_slice(const xplane& geometry) {
	if (geometry.get_coord() < low[0]) {
		low[0] = geometry.get_coord();
	}
	if (geometry.get_coord() > high[0]) {
		high[0] = geometry.get_coord();
	}
	compartment_slices_x.push_back(geometry);
}
void Visualisation::add_compartment_slice(const yplane& geometry) {
	if (geometry.get_coord() < low[1]) {
		low[1] = geometry.get_coord();
	}
	if (geometry.get_coord() > high[1]) {
		high[1] = geometry.get_coord();
	}
	compartment_slices_y.push_back(geometry);
}
void Visualisation::add_compartment_slice(const zplane& geometry) {
	if (geometry.get_coord() < low[2]) {
		low[2] = geometry.get_coord();
	}
	if (geometry.get_coord() > high[2]) {
		high[2] = geometry.get_coord();
	}
	compartment_slices_z.push_back(geometry);
}

struct Plot2d::MyVTKdata {
	MyVTKdata() {
		view = vtkSmartPointer<vtkContextView>::New();
		chart = vtkSmartPointer<vtkChartXY>::New();
		line = chart->AddPlot(vtkChart::LINE);
		chart->GetAxis(1)->SetBehavior(vtkAxis::FIXED);
		chart->GetAxis(0)->SetBehavior(vtkAxis::FIXED);

		view->GetScene()->AddItem(chart);
		view->GetInteractor()->Initialize();
	}
	vtkSmartPointer<vtkContextView> view;
	vtkSmartPointer<vtkChartXY> chart;
	vtkSmartPointer<vtkPlot> line;
};

Plot2d::Plot2d(const double dt, const std::vector<double>& x, const std::vector<double>& y,
			const char* x_label, const char* y_label, const char* title):
				vis_dt(dt),x(x),y(y),x_label(x_label),y_label(y_label),title(title),
				xlow(std::numeric_limits<double>::infinity()),
				xhigh(-std::numeric_limits<double>::infinity()),
				ylow(std::numeric_limits<double>::infinity()),
				yhigh(-std::numeric_limits<double>::infinity()) {
	my_vtk_data = new MyVTKdata();
	my_vtk_data->chart->SetTitle(title);
	my_vtk_data->chart->GetAxis(1)->SetTitle(x_label);
	my_vtk_data->chart->GetAxis(0)->SetTitle(y_label);
}

void Plot2d::integrate(const double dt) {

	if (get_time() > next_vis) {
		const int numPoints = x.size();
		ASSERT(numPoints == y.size(), "x and y should have same length");
		if (numPoints == 0) return;
		//my_vtk_data->chart->ClearPlots();
		//my_vtk_data->view->GetRenderer()->RemoveAllViewProps();

		vtkSmartPointer<vtkTable> table =
				vtkSmartPointer<vtkTable>::New();

		vtkSmartPointer<vtkDoubleArray> arrX =
				vtkSmartPointer<vtkDoubleArray>::New();
		arrX->SetName(x_label.c_str());
		table->AddColumn(arrX);

		vtkSmartPointer<vtkDoubleArray> arrY =
				vtkSmartPointer<vtkDoubleArray>::New();
		arrY->SetName(y_label.c_str());
		table->AddColumn(arrY);

		//		vtkSmartPointer<vtkFloatArray> arrS =
		//				vtkSmartPointer<vtkFloatArray>::New();
		//		arrS->SetName("Sine");
		//		table->AddColumn(arrS);

		// Fill in the table with some example values

		table->SetNumberOfRows(numPoints);
		for (int i = 0; i < numPoints; ++i) {
			table->SetValue(i, 0, x[i]);
			table->SetValue(i, 1, y[i]);
			if (yhigh < y[i]) yhigh = y[i];
			if (ylow > y[i]) ylow = y[i];
			if (xhigh < x[i]) xhigh = x[i];
			if (xlow > x[i]) xlow = x[i];
		}

		// Add multiple line plots, setting the colors etc

#if VTK_MAJOR_VERSION <= 5
		my_vtk_data->line->SetInput(table, 0, 1);
#else
		my_vtk_data->line->SetInputData(table, 0, 1);
#endif
		my_vtk_data->line->SetColor(0, 255, 0, 255);
		my_vtk_data->line->SetWidth(5.0);
		my_vtk_data->chart->GetAxis(1)->SetMinimum(xlow);
		my_vtk_data->chart->GetAxis(1)->SetMaximum(xhigh);
		my_vtk_data->chart->GetAxis(0)->SetMinimum(ylow);
		my_vtk_data->chart->GetAxis(0)->SetMaximum(yhigh);

		//my_vtk_data->chart->RecalculateBounds();

//		line = chart->AddPlot(vtkChart::LINE);
//#if VTK_MAJOR_VERSION <= 5
//		line->SetInput(table, 0, 2);
//#else
//		line->SetInputData(table, 0, 2);
//#endif
//		line->SetColor(255, 0, 0, 255);
//		line->SetWidth(5.0);
//		line->GetPen()->SetLineType(2);//For dotted line, can be from 2 to 5 for different dot patterns

		my_vtk_data->chart->GetScene()->SetDirty(true);
		my_vtk_data->view->GetRenderer()->Render();
		next_vis = get_time() + vis_dt;
	}
}

void Visualisation::print(std::ostream& out) {
	out << "\tVisualisation showing:" << std::endl;
	const int nplanes = xplanes.size() + yplanes.size() + zplanes.size();
	if (nplanes > 0) {
		cout << "\t\t Planes at:" << std::endl;
		BOOST_FOREACH(xplane p, xplanes) {
			out << "\t\t\t" << p << std::endl;
		}
		BOOST_FOREACH(yplane p, yplanes) {
			out << "\t\t\t" << p << std::endl;
		}
		BOOST_FOREACH(zplane p, zplanes) {
			out << "\t\t\t" << p << std::endl;
		}
	}
	const int nslice = compartment_slices_x.size() + compartment_slices_y.size() + compartment_slices_z.size();
	if (nslice > 0) {
		out << "\t\t Compartment slices at:" << std::endl;
		BOOST_FOREACH(xplane p, compartment_slices_x) {
			out << "\t\t\t" << p << std::endl;
		}
		BOOST_FOREACH(yplane p, compartment_slices_y) {
			out << "\t\t\t" << p << std::endl;
		}
		BOOST_FOREACH(zplane p, compartment_slices_z) {
			out << "\t\t\t" << p << std::endl;
		}
	}
	BOOST_FOREACH(Species *s, get_species()) {
		if (s->mols.size() > 0) out << "\t\t "<< s->mols.size() << " molecules from Species "<< s->id << std::endl;
	}
}

}

