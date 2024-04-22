#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Problem 1
This code implements the nonlinear functional mixed volterra-Fredholm Equation 
given as:
u(x) = g(x) + f(x, \int \int k(t, y, u(y))dydt  )
@author: Chinedu Nwaigwe
"""

import numpy as np;
import math;

#%% Example 1
class Problem1:
    'Example 1 '
    def __init__(self):
        print("\nThis is example 1. \n")
        
    def k(self, x, y, u):
        return x*y**3*np.cos( u )/(1.0 + x**2);
    
    def g(self, x ):
        return ( -x/6.0 + (x**2 + 1.0)*np.arccos(x*x) )/(1.0 + x**2);

    def exact(self, x):
        return np.arccos( x*x );
    def label(self):
        return "Problem 1";


#%% Example 2
class Problem2:
    'Example 2 '
    def __init__(self):
        print("\nThis is example 2. \n")
        
    def k(self, x, y, u):
        return (x**2*y**2*u)/(1.0 + u**4 + u**2);
    
    def g(self, x ):
        return x**2*(x - np.sqrt(3)*np.pi/54.0 );

    def exact(self, x):
        return x*x*x;
    def label(self):
        return "Problem 2";


class Problem3:
    'Example 3 '
    def __init__(self):
        print("\nThis is example 3. \n")
        
    def k(self, x, y, u):
        return (x**2 + y + u**3/3 - u**2)/(12);
    
    def g(self, x ):
        a = -x**2/12.0 + np.sin(x);
        a += -np.sin(2.0)/48.0 - 1.0/54.0;
        a += -np.cos(3.0)/432.0;
        a +=  np.cos(1.0)/48.0;
        return a

    def exact(self, x):
        return np.sin(x);
    def label(self):
        return "Problem 3";



















#%% These examples are from my paper on 
#Two Fixed  Methods
##Example 
class Problem4:
    'Example 1 '
    def __init__(self):
        print("\nThis is example 4. \n")
        
    def k(self, x, y, u):
        return u*np.exp(x-y-u*u);
    
    def g(self, x ):
        return np.exp(x - 1.0)/2.0 - np.exp(x - np.exp(-2.0))/2.0 + np.exp(-x)

    def exact(self, x):
        return np.exp( -x );
    def label(self):
        return "Problem 4";
    
#%%
class Problem5:
    'Example 5'
    def __init__(self):
        print("This is example 5.")
        
    def k(self, x, y, u):
        return x**3*y**5/(1.0 + u*u);
    
    def g(self, x):
        return x**3*(6.0 - np.log(2.0))/6.0;                       
    
    def exact(self, x):
        return x**3;
    def label(self):
        return "Problem 5";


class Problem6:
    'Example 6'
    def __init__(self):
        print("This is example 6.")
        
    def k(self, x, y, u):
        return x**4 + y**2 - u**3/3.0 + u**2;

    def g(self, x):
        return -x**4 + np.cos(x) - 5.0/6.0 - np.sin(2.0)/4.0 + np.sin(3.0)/36.0 + np.sin(1.0)/4.0;                
    
    def exact(self, x):
        return np.cos(x); 
    def label(self):
        return "Problem 6";

