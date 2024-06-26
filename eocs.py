#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code calculates the experimental order of convergence (EOC) of computed solutions.
@author: Dr Chinedu Nwaigwe  (Associate Professor of Applied Mathematics and Scientific Computing)
Date:    2023
Email:   nwaigwe.chinedu@ust.edu.ng; nedunwaigwe@gmail.com; nwaigwe.c@ufs.ac.za
"""

import numpy as np;



#%% This function computes the EOCs if the problem has exact solution
def eocExact(solver, problem, Ni, loops, tol=1e-8): #N = initial mesh, loops no of meshes
    olderror = 0.0;
    eoc = -1;
    for N in [Ni*2**i for i in range(loops)]: 
        grid = np.linspace(0, 1, N);
        solver.initialize(problem, grid);
        uh =  solver.compute(tol);  #numerical solution
        uex = problem.exact( grid );     #exact solution
        error =  np.max( abs(uh - uex) ); #max error
        if(N>Ni):
            eoc = np.log(error/olderror)/np.log(0.5); #eoc
        print(N,  "&", error, "&", eoc, "\\\ \hline" );
        olderror = error;
#%% eoc without exact solution
 
       
#%% Grid convergence study using max-norm of solution        
def gridConvergence(solver, problem, Ni, loops, tol=1e-8): #N = initial mesh, loops no of meshes
    for N in [Ni*2**i for i in range(loops)]: 
        grid = np.linspace(0, 1, N);  
        solver.initialize(problem, grid);
        uh =  solver.compute(tol);  #numerical solution
        #h = abs( grid[1] - grid[0]);
        print(N,  np.max(abs(uh)) );
        
#%%#%% This function solves a problem1 using 100 grid-points
#The Newton solver will converge with tolerence of tol = 1*10^{-8}            
#N=100;
def run(solver, problem, a, b, N ):
    grid = np.linspace(a, b, N);  #grid;
    solver.initialize(problem, grid);
    uh = solver.compute( 1e-8 );     #compute the solution
    return uh, grid;
