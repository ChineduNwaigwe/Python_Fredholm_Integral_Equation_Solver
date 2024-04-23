'''
A Solver Based on Ishikawa Iteration for Nonlinear Fredholm Integral Equations with Abitrary Kernels.
Author: Dr Chinedu Nwaigwe  (Associate Professor of Applied Mathematics and Scientific Computing)
Date:   2023
Email:   nwaigwe.chinedu@ust.edu.ng; nedunwaigwe@gmail.com; nwaigwe.c@ufs.ac.za
'''
import numpy as np;
import numpy.linalg as la #for solve() and norm
import time;
#%%
class FredholmFixedPointIteration:
    def __init__(self):
        self.problem = None;  self.grid = None;
        print("\nThis is FixedPointSolver.\n")
        
    def initialize(self, prob, grid ):     
        self.problem = prob;
        self.grid = grid;
        self.h = abs( grid[1]-grid[0]);
        self.N = len( grid ); 
            
    def picard(self, uold): #Tu_n
        problem = self.problem; 
        grid = self.grid;
        N = self.N; h = self.h;
        up = np.zeros( N ); # uold.copy();        
	    
	    #compute for i = 1, 2, ...                                        
        for i in range( N ):
	        xi = grid[i];
	        #1) compute k sum = Trapezoidal Integration
	        ksum_i = 0.5*(  problem.k(xi, grid[0], uold[0]) \
                         +  problem.k(xi, grid[-1], uold[-1])    );
	        for j in range( 1, N-1): #leaving out the first and last entries
	            ksum_i += problem.k(xi, grid[j], uold[j] );            
	        
	        #2) Picard Iteration: u(x_i) = Tu(x_i) = g(x_i) + f(x_i, h*ksum_i);
	        up[i] = problem.g(xi) + h*ksum_i;
	        
        return up;        
    
    def krasnoselskij(self, uold, lam):
        return (1.0 - lam)*uold + lam*picard(uold);
	
    def ishikawa(self, uold, n):
        an = 1.0/np.sqrt( n );
        ustar = (1.0 - an)*uold  + an*self.picard(uold );
        return  (1.0 - an)*ustar + an*self.picard(ustar);
           
    #%%loop tto carry out the iterations
    def compute(self, etol ):
    ##Initialize u = u_0(x) with g(x) + \int_a^b k(x, y, g(x))
        xi = self.grid[0];     gx = self.problem.g( xi );
        k0 = 0.5*( self.problem.k(xi, self.grid[0], gx ) \
             + self.problem.k(xi, self.grid[-1], gx) ) ;
        for j in (1, self.N-1):
            k0 += self.problem.k(xi, self.grid[j], gx );   
        up = np.ones(self.N) *(gx + self.h*k0);
            	    
        error = 1e17*etol; #ensures that initial error is larger than etol.
        counter = 1;
        while( error > etol  ):
            unew  = self.ishikawa( up, counter );   #.copy();
            #unew  = self.picard( up );
            error = la.norm( unew - up );
	    #print( "Inside while loop")
            up = unew.copy();
        return unew;    
	
	
	
	
