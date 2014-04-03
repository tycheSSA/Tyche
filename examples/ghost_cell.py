#!/usr/bin/python
import pyTyche as tyche
import sys
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
#from mayavi import mlab,tools

def experiment():
    L = 1.0
    num_comps = 20
    D = 1.0
    
    num_particles = 100
    mol_dt = 0.000001
    
    A = tyche.new_species(D)
    
    xlow = tyche.new_xplane(0,1)
    xghost = tyche.new_xplane(L/2.,1);
    xhigh = tyche.new_xplane(L,-1)
    ylow = tyche.new_yplane(0,1)
    yhigh = tyche.new_yplane(L,-1)
    zlow = tyche.new_zplane(0,1)
    zhigh = tyche.new_zplane(L,-1)
    interface = tyche.new_box([0,0,0],[L/2.,L,L],True)

    xghostboundary = tyche.new_reflective_boundary(xghost)
    xminboundary = tyche.new_reflective_boundary(xlow)
    xmaxboundary = tyche.new_reflective_boundary(xhigh)
    yminboundary = tyche.new_reflective_boundary(ylow)
    ymaxboundary = tyche.new_reflective_boundary(yhigh)
    zminboundary = tyche.new_reflective_boundary(zlow)
    zmaxboundary = tyche.new_reflective_boundary(zhigh)

    boundaries = tyche.group([xghostboundary,xmaxboundary,yminboundary, ymaxboundary, zminboundary, zmaxboundary])
    boundaries.add_species(A)
    
    compartments = tyche.new_compartments([0,0,0],[L,L,L],[L*1./num_comps,L,L])
    compartments.fill_uniform(A, [0,0,0],[L/2.,L,L],num_particles/2)
    compartments.add_diffusion(A);
    A.fill_uniform([L/2.,0,0], [L,L,L], num_particles/2)

    compartments.set_ghost_cell_interface(interface)
    ghostcoupling = tyche.new_ghost_cell_boundary(interface,compartments)
    ghostcoupling.add_species(A)

    bd = tyche.new_diffusion_with_tracking(interface,compartments)
    bd.add_species(A)

    algorithm = tyche.group([bd,boundaries,compartments])
    
    compartments.list_reactions()

    output_dt = 1e-2
    time = 0;
    print algorithm
    
    ion()
    figure()
    ax = subplot(211)
    ax2 = subplot(212)

    nums = []
    while True:
        for i in range(10):
            time = algorithm.integrate_for_time(output_dt,mol_dt)
            h,b = histogram(A.get_particles()[0], range=(0,L), bins=num_comps)
            c = A.get_compartments()[:,0,0]
            c[num_comps/2] = 0
            h += c
            if h.sum()!=num_particles:
                raise ValueError("Inconsistent Ghost Cell! Number of particles is %d and not %d as it should be!" % (h.sum(), num_particles))
            nums.append(h)
        ax.clear()
        ax2.clear()
        x = (arange(0,num_comps)+0.5)*L*1./num_comps
        ax.plot(x,mean(nums,axis=0))
        ax2.plot(x,std(nums,axis=0))
        ax2.set_xlabel(str(time))

        ax.set_xlabel("x")
        ax.set_ylabel("Mean Frequency")
        ax2.set_xlabel("x")
        ax2.set_ylabel("Standard deviation")
        draw()        
        
tyche.init(sys.argv);
experiment()
