#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code uses the Solver Ishikawa Iteration to solve Nonlinear Fredholm Integral Equations with Abitrary Kernels.
@author: Dr Chinedu Nwaigwe  (Associate Professor of Applied Mathematics and Scientific Computing)
Date:    2023
Email:   nwaigwe.chinedu@ust.edu.ng; nedunwaigwe@gmail.com; nwaigwe.c@ufs.ac.za
"""

from solver import *;  #fixed-point solver
from eocs import *;    #error analysis
import problems;       #the problems to solve
import pylab as plt;
import time;



#%% EOCs:
prob = problems.Problem1();
solver = FredholmFixedPointIteration();

t = time.time();
eocExact(solver, prob, 5, 10, 1e-8 );

'''
#print( time.time()-t)
#gridConvergence(prob, 2, 10, 1e-8 );        
ls = ["solid", "dashdot", "dotted", "dashed",  "solid", "dashdot", "dotted", "dashed"]
cs = ["green",    "blue",  "magenta", "black",  "red",   "blue",  "magenta"]
 

xa = 0; xb = 1;
N = [5, 10, 30, 100];
for i in range( len(N) ):
    uh, grid = run(solver, prob, xa, xb, N[i]);
    plt.plot(grid, uh, linewidth=4, \
             label="Numerical Solution with "+str(N[i])+" Points",\
             linestyle=ls[i], color=cs[i]);

plt.plot(grid, prob.exact(grid), "red", linewidth=3, linestyle="--",\
         label="Exact Solution of "+str(prob.label()) );    
plt.xlabel("$x$", fontsize=25);
plt.ylabel("$u(x)$", fontsize=25);
plt.xlim([xa,xb])
#plt.ylim([min(abs(uh)),  max(abs(uh))]);
plt.tick_params(axis='both', width=6, labelsize=25);  #change width of ticks and the size of their labels.
plt.rcParams['axes.linewidth'] = 3; 
plt.legend(fontsize=25)
plt.show();


#%% This function solves a problem1 using 100 grid-points
#The Newton solver will converge with tolerence of tol = 1*10^{-8}            
#N=100;
#def run(solver, problem, a, b, N ):
#    #N = 20;
#    #a = 0.0; b = 1.0;
#    #problem = problems.Problem2(); #problem;
#    grid = np.linspace(a, b, N);  #grid;
#    solver.initialize(problem, grid);
#    uh = solver.compute( 1e-8 );     #compute the solution
#    #plt.plot(grid, uh);              #plot the solution
#    #for i in range(N):
#    #    print(i, grid[i],   uh[i] );
#    return uh, grid;
#run();
'''
