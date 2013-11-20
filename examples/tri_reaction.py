#!/usr/bin/python
import pyTyche as tyche
import sys
from tvtk.api import tvtk
from tvtk.tools import mlab
import numpy as np
import matplotlib.pyplot as plt

def init_vis(grid):
    g = tvtk.Glyph3D(scale_mode='data_scaling_off', input=grid)
    cs = tvtk.SphereSource(radius=0.01)
    g.source = cs.output
    m = tvtk.PolyDataMapper(input=g.output)
    a = tvtk.Actor(mapper=m)
    rw = tvtk.RenderWindow(size=(600, 600))
    ren = tvtk.Renderer(background=(0.5, 0.5, 0.5))
    rw.add_renderer(ren)
    rwi = tvtk.RenderWindowInteractor(render_window=rw)
    ren.add_actor(a)
    rwi.initialize()
    rwi.start()
    return g

def experiment(num_particles,timesteps):
    L = 1
    D = 1
    k2 = 1000.0
    k1 = k2/num_particles**3
    max_t = 1.0

    DAB = D+D;
    DABC = D + D*D/DAB;
    Rsq =  ( k1 / (4*3.14**3*(DAB*DABC)**(3.0/2.0) ) )**0.5;
    mol_dt = Rsq/10000;
    
 
    A = tyche.new_species(D)
    B = tyche.new_species(D)
    C = tyche.new_species(D)


    bd = tyche.new_diffusion()
    bd.add_species(A)
    bd.add_species(B)
    bd.add_species(C)
    
    xlow = tyche.new_xplane(0,1)
    xhigh = tyche.new_xplane(L,-1)
    ylow = tyche.new_yplane(0,1)
    yhigh = tyche.new_yplane(L,-1)
    zlow = tyche.new_zplane(0,1)
    zhigh = tyche.new_zplane(L,-1)


    xminboundary = tyche.new_jump_boundary(xlow,[L,0,0])
    xmaxboundary = tyche.new_jump_boundary(xhigh,[-L,0,0])
    yminboundary = tyche.new_jump_boundary(ylow,[0,L,0])
    ymaxboundary = tyche.new_jump_boundary(yhigh,[0,-L,0])
    zminboundary = tyche.new_jump_boundary(zlow,[0,0,L])
    zmaxboundary = tyche.new_jump_boundary(zhigh,[0,0,-L])

    boundaries = tyche.group([xminboundary, xmaxboundary, yminboundary, ymaxboundary, zminboundary, zmaxboundary])
    boundaries.add_species(A)
    boundaries.add_species(B)
    boundaries.add_species(C)

    dr3 = tyche.new_tri_reaction(k1, A, B, C, [B,C], mol_dt,
                          [0,0,0], [L,L,L], [True, True, True])
        
     
    dr = tyche.new_zero_reaction(k2,[0,0,0],[L,L,L])
    dr.add_species(A)

    algorithm = tyche.group([bd,boundaries,dr,dr3])
    
    A.fill_uniform([0,0,0],[L,L,L],int(num_particles))
    B.fill_uniform([0,0,0],[L,L,L],int(num_particles))
    C.fill_uniform([0,0,0],[L,L,L],int(num_particles))
    
    N = 1000
    output_dt = max_t/N
    time = 0;
    print algorithm
    
    grid = tvtk.to_tvtk(A.get_vtk())
    
    numA = np.zeros(N)
    for i in range(N):
        print A,B,C
        numA[i] = A.get_concentration([0,0,0],[L,L,L],[1,1,1])[0]
        print 'time = ',time,' ',i*100.0/N,' percent done'
        time = algorithm.integrate_for_time(output_dt,mol_dt)
#        grid = tvtk.to_tvtk(A.get_vtk())
#        g = init_vis(grid)
    print algorithm
    return numA

        
        
tyche.init(sys.argv);
#numA1 = experiment(100,1e5)
#numA2 = experiment(100,5e6)
numA = experiment(100,1e7)

plt.figure(figsize=(6,4.5))

        
#n, bins, patches = plt.hist(numA1, 20,range=[mmin, mmax])
n, bins, patches = plt.hist(numA, 30)

#n2, bins, patches = plt.hist(numA2, 20,histtype = 'step',range=[mmin, mmax])
#n3, bins, patches = plt.hist(numA3, 20,histtype = 'step',range=[mmin, mmax])
#plt.ylim([0,max(n.max(),n2.max(),n3.max())])
plt.xlabel("number of molecules of A")
plt.ylabel("stationary distribution")        
plt.savefig("stationary_distribution_of_A.pdf")
