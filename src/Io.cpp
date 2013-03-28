/* 
 * io.cpp
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of PDE_BD.
 *
 * PDE_BD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PDE_BD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PDE_BD.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Feb 16, 2013
 *      Author: mrobins
 */

#include "Io.h"
#include "Tyche.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkSmartPointer.h"

#include <iostream>
#include <fstream>

namespace Tyche {

void write_grid(std::string filename, vtkUnstructuredGrid* grid) {
	const int my_rank = 0;
	const int num_procs = 1;
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer =
			vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
	writer->SetNumberOfPieces(num_procs);
	writer->SetStartPiece(my_rank);
	writer->SetEndPiece(my_rank);
	writer->SetInput(grid);
	writer->SetDataModeToBinary();
	writer->SetFileName(filename.c_str());
	writer->Write();
}

}
