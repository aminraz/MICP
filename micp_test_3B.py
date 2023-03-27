"""Simulation for the paper, based on the data from Martinez et al. 2013."""

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
# 10 if (MOL("Urea")  <= 1.0e-5) then goto 40
11 cell_l = 5.2
12 t_d = 288000
15 initial_por = .38
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))
30 moles = rate * TIME * decay_ratio
40 save moles
-END

Calcite
-START
20 initial_por = .38
21 V_cell = 102 # ml
22 d_calcite = 36.9 # cm³/mol
23 pow_a = 2
25 A0 = .5 # m2/l
30 si_cc = SI("Calcite")
70 Area = A0 * (initial_por/get_por(cell_no))^pow_a * soln_vol
80 rate = Area * (parm(1) + parm(2)*act("H+")) * (1-10^SI("Calcite"))
100 moles = rate * TIME * get_por(cell_no)/initial_por 
110 change_por(initial_por - d_calcite*kin("Calcite")/V_cell,cell_no) 
200 SAVE moles
-END
#-------------------------------------------------------------------------------

SOLUTION 0
Urea .05
Amm .0567
Ca .05
Cl .1567
C .0038 as Ca Calcite 0.0
Na .0038
-units 	mol/kgw
-pH 	7	charge_balance
water 38.7e-3 # kg
END

SOLUTION 1-10
     -units 	mol/kgw
     -pH 	7	charge_balance
water 38.7e-3 # kg
END

Kinetics 1-10 
    Ureolysis
         -m0 100  
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms 7.0e-5 0.0173
    Calcite
        -m0   0# mol/L
        -parms 1.55e-6 .501
END

TRANSPORT 
   cells 10
   shifts 178
   time_step 1055 seconds
   flow_direction forward
   boundary_conditions flux constant 
   length 0.052 
   diffusion_coefficient 1.0e-9
   porosities 10*.38
   -dispersivities 10*.005

SELECTED_OUTPUT
    -step    true
     -h    false
    -sim    false
    -soln   true
     -pe     false
    -totals false
    -time   true
    -tot  Ca  
    -k Ureolysis Calcite
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
t_step = 1055  # time steps
shift_m = 178
row_0 = 12  # cell number
t = table.loc[4:]['time']
C_u_0 = .05 # initial urea concentration molar
C_ca_0 = .05 # initial calcium concentration molar
por_0 = .38
l_cell = 5.2 #cm

fig, ax = plt.subplots()
shift = 177
ax.plot(table["k_Calcite"][row_0+shift*n_cells:row_0+n_cells*(shift+1)]/V_cell,np.arange(1,11)*l_cell,'ko', label = 'Simulation', linewidth=2)
ax.set(ylabel='Column height (cm)',xlabel='Calcite precipitation (mol/l)')
# ax.plot(C_cal,l_cell*n_cells-dis,'k',label="Experimental data", linewidth=3)
ax.grid()
ax.legend()

plt.show()

