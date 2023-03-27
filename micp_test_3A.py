"""Simulation of Test 3A, based on the data from Martinez et al. 2013."""

import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

phreeqc = phreeqc_mod.IPhreeqc('/usr/local/lib/libiphreeqc.so')

initial_conditions = """
#-------------------------------------------
SOLUTION_MASTER_SPECIES
    Urea Urea  0 CO(NH2)2 60.06 
    
SOLUTION_SPECIES 
    Urea = Urea
            log_k  0 
            Vm  14.5
    Cl- + H+ = HCl
            log_k   0.7
            delta_h 0       kcal
            -gamma  0       0.06
    Cl- + Ca+2 = CaCl+
            log_k   0.3
            delta_h 0       kcal
            -gamma  0       0.06
    2Cl- + Ca+2 = CaCl2
            log_k   .66
            delta_h 0       kcal
            -gamma  0       0.06
#---------------------------------------------

RATES
Ureolysis
-START
15 initial_por = .36
16 Cb = (3.24*(cell_no*5-2.5)^3-1.47e2*(cell_no*5-2.5)^2+9.56e3*(cell_no*5-2.5)-1.79e4)/1e6
18 decay_ratio = exp(-TOTAL_TIME/288000)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))
30 moles = rate * TIME * decay_ratio # * get_por(cell_no)/initial_por
40 save moles
-END

Ureolysis_cell
-START
10 cell_D = 8 # the number of cell at the D outlet in the experiment. 
17 Cb = (3.24*(cell_D*5-2.5)^3-1.47e2*(cell_D*5-2.5)^2+9.56e3*(cell_D*5-2.5)-1.79e4)/1e6 
18 decay_ratio = exp(-5810/288000)
20 rate = Cb*soln_vol*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))
30 moles = rate * TIME *decay_ratio
40 save moles
-END


Calcite
-START
20 initial_por = .36
21 V_cell = 102 # ml
22 d_calcite = 36.9 # cm³/mol
23 pow_a = 2
25 A0 = .5
30 si_cc = SI("Calcite")
70 Area = A0 * initial_por/get_por(cell_no)^pow_a * soln_vol
80 rate = Area * (parm(1) + parm(2)*act("H+")) * (1-10^SI("Calcite"))
100 moles = rate * TIME * get_por(cell_no)/initial_por 
110 change_por(initial_por - d_calcite*kin("Calcite")/V_cell,cell_no) 
200 SAVE moles
-END
#-------------------------------------------------------------------------------

SOLUTION 0
Urea .05; Amm .0567; Ca .05; Cl .1567; Na .0038
C .0038 as Ca Calcite 0.0
-units 	mol/kgw
-pH 	7	charge_balance
water 36.9e-3 # kg
END

SOLUTION 1-10
     -units 	mol/kgw
     -pH 	7	charge_balance
water 36.9e-3 # kg
END

Kinetics 1-10 
    Ureolysis
         -m0 100  
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms  7.0e-5 0.0173
    Calcite
        -m0   0
        -parms 1.55e-6 .501
END

TRANSPORT 
   cells 10
   shifts 10 
   time_step 221 seconds 
   flow_direction forward
   boundary_conditions flux constant 
   length 0.0523
   diffusion_coefficient 1.0e-9
   porosities 10*.36
   -dispersivities 10*.005

SAVE SOLUTION 1-10
END

USE SOLUTION 8
Kinetics 
    Ureolysis_cell
         -m0 100  
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms  7.0e-5 0.0173 # 7.0e-5 0.02
    Calcite
        -m0   0# mol/L
        -parms 1.55e-6 .501
        -steps 7200 second in 100 steps

SELECTED_OUTPUT
    -step    true
     -h    false
    -sim    false
    -soln   true
     -pe     false
    -totals false
    -time   true
    -pH true
    -tot  Ca  Urea Amm
    -si Calcite
    -k Ureolysis_cell Calcite
   -dist false
   -state false
END
"""
np.set_printoptions(precision=2)
pd.set_option('display.precision', 3)

phreeqc.load_database(r"path-to-database")
phreeqc.run_string(initial_conditions)
output = phreeqc.get_selected_output_array()
table = pd.DataFrame(output[1:], columns=output[0])
# table = pd.DataFrame(data=np.transpose(output[1:]), index=output[0])
print("\n")
print(table)

n_cells = 10  # totla number of cells
V_cell = 102e-3 #l
t_step = 221  # time step
shift = 10
row_0 = 12  # cell number
t = table.loc[4:]['time']
C_u_0 = .05 # initial urea concentration molar
C_ca_0 = .05 # initial calcium concentration molar
C_amm_0 = 0.0567
por_0 = .36
l_cell = 5.23 #cm
k_step = 100


fig,ax = plt.subplots()
ax.plot(np.arange(2/k_step,2+2/k_step,2/k_step),table['Ca(mol/kgw)']/C_ca_0,'k:',label="[Calcium]: Phreeqc model")
ax.plot(np.arange(2/k_step,2+2/k_step,2/k_step),table['Urea(mol/kgw)']/C_u_0,'k-',label="[Urea]: Phreeqc model")
ax.legend()
ax.set(xlabel='time(h)',ylabel='Normalized concentration')

fig,ax = plt.subplots()
ax.plot(np.arange(2/k_step,2+2/k_step,2/k_step), table['pH'],'k-', label='pH: Phreeqc model')
ax.set(xlabel='time(h)',ylabel='pH')
ax.legend()

fig,ax = plt.subplots()
ax.plot(np.arange(2/k_step,2+2/k_step,2/k_step), table['Amm(mol/kgw)']/C_amm_0,'k-', label='Phreeqc model')
ax.set(xlabel='time(h)',ylabel='Amm')
ax.legend()

plt.show()
