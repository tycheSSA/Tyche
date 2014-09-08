#!/usr/bin/python
import pyTyche as tyche
import sys
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

def experiment():
    # Domain size
    L = 1.0
    # Number of compartments in x direction
    num_comps = 20
    # Diffusion constant
    D = 1.0
    # Assumed average number of particles to the left of domain (fixed
    # concentration boundary conditions
    A0 = 10
    
    # Time step
    mol_dt = 0.000001
    
    # Create species
    A = tyche.new_species(D)
    
    # Create domain bounds
    xlow = tyche.new_xplane(1e-5,1)
    xghost = tyche.new_xplane(L/2.,1);
    xhigh = tyche.new_xplane(L,-1)
    ylow = tyche.new_yplane(0,1)
    yhigh = tyche.new_yplane(L,-1)
    zlow = tyche.new_zplane(0,1)
    zhigh = tyche.new_zplane(L,-1)
    # Create interface box
    interface = tyche.new_box([0,0,0],[L/2.,L,L],True)

    # Boundary conditions on the domain bounds and the interface
    # Need a reflective boundary at the ghost cell interface
    xghostboundary = tyche.new_reflective_boundary(xghost)
    xminboundary = tyche.new_reflective_boundary(xlow)
    xmaxboundary = tyche.new_reflective_boundary(xhigh)
    yminboundary = tyche.new_reflective_boundary(ylow)
    ymaxboundary = tyche.new_reflective_boundary(yhigh)
    zminboundary = tyche.new_reflective_boundary(zlow)
    zmaxboundary = tyche.new_reflective_boundary(zhigh)

    # Add all boundaries to a single operator
    boundaries = tyche.group([xghostboundary,xmaxboundary,yminboundary, ymaxboundary, zminboundary, zmaxboundary])
    boundaries.add_species(A)
    
    # Create compartments throughout the domain
    compartments = tyche.new_compartments([0,0,0],[L,L,L],[L*1./num_comps,L,L])
    compartments.add_diffusion(A);

    # Boundary reaction on leftmost boundary (particle creation)
    compartments.add_reaction_on(A0/(L/num_comps)**3, [[],[A]], xlow)
    compartments.add_reaction_on(1./(L/num_comps)**2, [[A],[]], xlow)

    # Removal reaction
    compartments.add_reaction(1., [[A],[]])
    
    # Introduce ghost cell interface at L/2
    compartments.set_ghost_cell_interface(interface)

    # Diffusion operator with tracking (tracks number of particles in ghost
    # compartment at the interface
    bd = tyche.new_diffusion_with_tracking(interface,compartments)
    bd.add_species(A)

    # Removal reaction in molecular domain
    remove = tyche.new_uni_reaction(1., [[A],[]])

    # Create algorithm
    algorithm = tyche.group([bd,boundaries,compartments, remove])
    
    # Get data every Dt
    output_dt = 1e-2
    
    # Set up plotting
    ion()
    figure()
    ax = subplot(211)
    ax2 = subplot(212)

    # List to hold number of particles in compartments
    nums = []
    # Run until steady state is reached
    time = algorithm.integrate_for_time(5.,mol_dt)
    while True:
        # Collect data
        for i in range(10):
            # Run for Dt
            time = algorithm.integrate_for_time(output_dt,mol_dt)
            # Put molecules in bins similar to compartment spacing
            h,b = histogram(A.get_particles()[0], range=(0,L), bins=num_comps)
            # Get compartment data
            c = A.get_compartments()[:,0,0]
            # Fix ghost compartment so that we do not count it twice
            c[num_comps/2] = 0
            # Add compartments to molecule histogram
            h += c
            nums.append(h)

        #P Plotting
        ax.clear()
        ax2.clear()
        x = (arange(0,num_comps)+0.5)*L*1./num_comps

        # Plot mean number of particles
        ax.plot(x,mean(nums,axis=0))
        # Solution for steady state of reaction-diffusion PDE (has an error
        # whereby the first compartment should be shifted by h/2)
        ax.plot((1-x), A0*1./cosh(1.)*cosh(x))

        # Plot standard deviation
        ax2.plot(x,std(nums,axis=0))

        ax.set_xlabel("x")
        ax.set_ylabel("Mean Frequency")
        ax2.set_xlabel("x")
        ax2.set_ylabel("Standard deviation")
        draw()        
        
tyche.init(sys.argv);
experiment()
