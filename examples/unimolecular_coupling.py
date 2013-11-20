#!/usr/bin/python
import pyTyche as tyche
import sys
import matplotlib.pyplot as plt
import numpy as np
from mayavi import mlab,tools




def experiment():
    L = 1.0
    D = 1.0
    
    num_particles = 10000.0
    max_t = 10.0;
    factor = 2.0;
    mol_dt = 0.001 
    
    A = tyche.new_species(D)
    
    bd = tyche.new_diffusion()
    bd.add_species(A)

    xlow = tyche.new_xplane(0,1)
    xhigh = tyche.new_xplane(L,-1)
    ylow = tyche.new_yplane(0,1)
    yhigh = tyche.new_yplane(L,-1)
    zlow = tyche.new_zplane(0,1)
    zhigh = tyche.new_zplane(L,-1)
    interface = tyche.new_box([0,0,0],[L/2.0,L,L],True)


    xminboundary = tyche.new_reflective_boundary(xlow)
    yminboundary = tyche.new_reflective_boundary(ylow)
    ymaxboundary = tyche.new_reflective_boundary(yhigh)
    zminboundary = tyche.new_reflective_boundary(zlow)
    zmaxboundary = tyche.new_reflective_boundary(zhigh)


    boundaries = tyche.group([xminboundary, yminboundary, ymaxboundary, zminboundary, zmaxboundary])
    boundaries.add_species(A)
    
    dr = tyche.new_uni_reaction(1.0,[[A],[]])

    compartments = tyche.new_compartments([0,0,0],[L,L,L],[L/20.0,L,L])
    compartments.add_diffusion(A);
    compartments.add_reaction(1.0,[[A],[]])
    compartments.add_reaction_on(num_particles*20.0/L,[[],[A]],xlow)
    
    corrected_interface = True 
    compartments.set_interface(interface, mol_dt, corrected_interface)
    coupling = tyche.new_coupling_boundary(interface,compartments,corrected_interface)
    coupling.add_species(A)
    
    
    algorithm = tyche.group([bd,compartments,dr,boundaries,coupling])

    
    output_dt = max_t/100.0
    time = 0;
    print algorithm
    
#    plt.figure()
#    plt.ion()
#    concentration = A.get_concentration([0,0,0],[2.0*L,L,L],[40,1,1])
#    x = np.arange(0,2.0*L,L/20.0)
#    l, = plt.plot(x,concentration[:,0,0])
#    plt.ylim([0,num_particles])

    for i in range(100):
        print A
        print 'time = ',time,' ',i,' percent done'
        time = algorithm.integrate_for_time(output_dt,mol_dt)
#        concentration = A.get_concentration([0,0,0],[2.0*L,L,L],[40,1,1])
#        l.set_ydata(concentration[:,0,0])
#        plt.draw()
#        x,y,z = A.get_particles()
#        mlab.points3d(np.array(x),np.array(y),np.array(z),scale_factor = 0.05,color=(1,1,1))
#        [0,0,0],[L,L,L],[L/20.0,L,L]
#        x, y, z = np.mgrid[0:L:L/20.0, 0:L:L, 0:L:L]
#        comps  = A.get_compartments()
#        values = tools.pipeline.scalar_field(x,y,z,comps)
#        mlab.pipeline.volume(values)
        
        
    print algorithm


        
        
tyche.init(sys.argv);
experiment()
