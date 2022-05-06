#!/usr/bin/env python
# coding: utf-8

from sympy import *
from numpy import *
# Define constant values
w_B  = -16.6        #bacteria            # w_b = -16.6kbT
w_H  = -6.72        #host                # w_h = -6.72kbT
A_p = 400           #Amstrong square
v_p = 33**3         #Amstrong cube
A_B = 400           #Amstrong square (Surface area of each cells)
A_H = A_B        #Amstrong square
a_B = 71            #Amstrong square (Area of each lipids)
a_H = 74             #Angstrom square 
A_t = 12*10**8 #Ah and Ab
At = 12*10**8 #Ah and Ab

# Priminary Constant values
avo_con = 6.023*10**23 #Avogadro constant ==> 6.022 × 10^23 mol^(-1)
mu_M_unit = 10**-6*avo_con *10**-27
# ccb = np.logspace(np.log10(0.000001) , np.log10(10.0),num = 30)
unit =10**-24 #(unit conversion)
# PL_B=1/48; PL_H=1/99;
PL_H = 0.015228218882144314
PL_B = 0.17699496922894759


def MIC_competitive_membrane(CB, CH, sur_B, sur_H, multi):
    Cp, PL_h = symbols('Cp PL_h')

    MIC_0 = 1/v_p* A_B/a_B* PL_B /(1-  A_B/a_B* PL_B)*N(exp(w_B))
    MHC_0 = 1/v_p* A_H/a_H* PL_h /(1-  A_H/a_H* PL_h)*N(exp(w_H))
    #     MIC_0 = 1/v_p* A_B/a_B*sur_B * PL_B /(1-  A_B/a_B*sur_H * PL_B)*N(exp(w_B))
    #     MHC_0 = 1/v_p* A_H/a_H*sur_B * PL_h /(1-  A_H/a_H*sur_H * PL_h)*N(exp(w_H))

    eq1 = Eq( (PL_B*At/a_B *sur_B  *CB*unit + PL_h*At *multi /a_H *sur_H* CH*unit) + MIC_0, Cp*mu_M_unit)
    eq2 = Eq( (PL_B*At/a_B *sur_B  *CB*unit + PL_h*At *multi /a_H *sur_H* CH*unit) + MHC_0, Cp*mu_M_unit)
    solution = solve((eq1, eq2), (Cp, PL_h)) 

    return solution[0]


def MHC_competitive_membrane(CB, CH, sur_B, sur_H,multi):
    Cp, PL_b = symbols('Cp PL_b')
    
    MIC_0 = 1/v_p* A_B/a_B* PL_b /(1-  A_B/a_B* PL_b)*N(exp(w_B))
    MHC_0 = 1/v_p* A_H/a_H* PL_H /(1-  A_H/a_H* PL_H)*N(exp(w_H))
    
    eq1 = Eq( (PL_b*At/a_B *sur_B  *CB*unit + PL_H*At *multi /a_H *sur_H* CH*unit) + MIC_0, Cp*mu_M_unit)
    eq2 = Eq( (PL_b*At/a_B *sur_B  *CB*unit + PL_H*At *multi /a_H *sur_H* CH*unit) + MHC_0, Cp*mu_M_unit)
    solution = solve((eq1, eq2), (Cp, PL_b)) 
    return solution[0]


def MIC_competitive_cell(CB, CH, Np, multi):
    sol  = MIC_competitive_membrane(CB, CH, 1, 1,17)
    PL_h = sol[1]
    
    Cp = symbols('Cp')
    MIC_0 = 1/v_p* A_B/a_B * PL_B /(1-  A_B/a_B * PL_B)*N(exp(w_B))

    eq1 = Eq( (PL_B*At/a_B  + Np) *CB*unit + (PL_h*At *multi /a_H * CH*unit) + MIC_0, Cp*mu_M_unit)
    solution = solve((eq1), (Cp)) 
    
    return solution[0]
    

def MHC_competitive_cell(CB, CH, Np, multi):
    sol  = MHC_competitive_membrane(CB,CH, 1, 1,17)
    PL_b = sol[1]
    
    Cp = symbols('Cp')
    MHC_0= 1/v_p* A_H/a_H* PL_H /(1-  A_H/a_H* PL_H)*N(exp(w_H))

    eq1 = Eq( (PL_b*At/a_B  + Np )*CB*unit  + (PL_H*At *multi /a_H   + Np )* CH*unit  + MHC_0, Cp*mu_M_unit)

    solution = solve((eq1),(Cp)) 
    
    return solution[0]



def Multiplication_factor_MIC(CB, CH,sur_B,sur_H,Np):
    # sol  = MIC_competitive_membrane(CB,CH,sur_B,sur_H)
    # PL_h = sol[1]
    con1 = PL_B*At/a_B
    # surface = ((con1 + Np )*CB + (con2 + Np )*CH )/ (con1*CB + con2*CH)
    surface = ((con1 + Np ) )/ (con1 )

    return surface


def Multiplication_factor_MHC(CB, CH,PL_b,Np):
    # sol  = MHC_competitive_membrane(CB,CH,sur_B,sur_H)
    # PL_b = sol[1]
    
    con1 = PL_b*At/a_B
    con2 = PL_H*At/a_H
    S_B = (con1 + Np )*CB / (con1*CB)
    S_H = ((con2 + Np )*CH )/ (con2*CH)
    
    return S_B, S_H

    
