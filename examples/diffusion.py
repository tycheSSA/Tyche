#!/usr/bin/python
import pyTyche as tyche
import sys
from tvtk.api import tvtk
from tvtk.tools import mlab

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

def experiment(num_particles):
    L = 1
    #L = 50e-9
    D = 1
    #D = 1e-10
    
    timesteps = 100000
    #timesteps = 10000
    max_t = 4.0/(2.0*num_particles**0.5);
    mol_dt = max_t/timesteps
    #mol_dt = 0.025*1e-8
    #max_t = mol_dt * timesteps
    
    
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


    algorithm = tyche.group([bd,boundaries])
    
    A.fill_uniform([0,0,0],[L,L,L],int(num_particles*L**3))
    B.fill_uniform([0,0,0],[L,L,L],int(num_particles*L**3))
    C.fill_uniform([0,0,0],[L,L,L],int(num_particles*L**3))

    
    output_dt = max_t/100.0
    time = 0;
    print algorithm
    
    
    for i in range(100):
        print A,B,C
        print 'time = ',time,' ',i,' percent done'
        time = algorithm.integrate_for_time(output_dt,mol_dt)
        grid = tvtk.to_tvtk(A.get_vtk())
        g = init_vis(grid)
    print algorithm


        
        
tyche.init(sys.argv);
experiment(5)