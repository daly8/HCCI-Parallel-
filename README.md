# HCCI-Parallel-
Requires Python 3.X

Cantera-based Homogeneous Charge Compression Ignition simulator with parallel processing for large parametric studies

This code performs 0-D HCCI simulations. The user specifies the chemical mechanism, initial conditions (Temperature, pressure, fuel species, etc.) and the simulations are performed in parallel. 

If the simulation is ever interrupted, simply restart it and the cases that were previously completed will be skipped (so long as you are working in the same directory as before)

To set up initial conditions, engine and fuel parameters, heat tranfer constants, etc., open up HCCI_LTC_data.yml and input your values.  Otherwise, from the command line you can enter the argument "--interactive_input" and you can enter all the parameters in the command line.

Variables: 

{'Fuel names': ['FACE_G_KAUST'],
        'Fuel species': [['IC5H12','C7H16-2','CPT','T124MBZ','C6H12-1','IC8H18','C4H10']],
        'Species mole fractions':[[0.098,0.07,0.158,0.084,0.084,0.437,.069]],
        'Mechanism': 'KAUST_FACE_F_Reduced_Mech.cti',
        'Pressure':{'units':'Pa', 'value':100000},
        'RPM':{'Initial':800.0, 'step':100, 'final':3100.0},
        'Phi': {'Initial':0.15, 'step':0.01, 'final':0.45},
        'Temperature':{'units':'K', 'Initial':350, 'step':20, 'final':570},
        'Compression ratio':{'value':13},
        'Cylinder bore': {'units':'m', 'value':.0828},
        'Displacement volume':{'units':'m', 'value':616*10**(-6)},
        'Connecting rod to crank radius ratio':{'value':4.44},
        'Nu Const. A':{'value':.035},
        'Nu Const. B':{'value':.8},
        'Cylinder gas Const.':{'value':2.28},
        'Wall temperature':{'units':'K', 'value':430},
        'Prandtl number': {'value':0.7},
        'Nu Const. C':{'value':0.0},
        'RPM range': {'value':286},
        'Working directory': 'C:\\Users\\Daly\\Documents\\Python Scripts\\',
        'Number of cores':{'value':10},
        'Print progress to cmd':{'value':10000}}

---HCCI parametric study values---
Fuel names: The name of the fuel under study, e.g., PRF100
Fuel species: The species makup up the fuel under study, e.g., nC7H16 and iC8H18. The chemical names should correspond to the chemical nomenclature in the chemical mechanism
mechanism: mechanism file (.cti)
Pressure: Initial pressure at BDC
RPM: Engine RPM range for parametric study
Phi: Fuel / Air equivalence ratio for parametric study
Temperature: Initial temperature at BDC for parametric study  

---Engine parameters---
These are all geometric features of the ICE. Default values are a CRC research engine

---Heat transfer correlation Parameters---
"Correlations"
w_bar=C11*Sm                   #Woschni average cylinder gas velocity
Nu=aa*Re**bb*Pr**cc            #Nusselt Number
 
Nu Const. A: Nusselt number constant "aa"
Nu Const. B: Nusselt number constant "bb"
Cylinder gas Const.: Woschni average cylinder gas velocity const. C11
Nu Const. C: Nusselt number constant "cc"

---Recording Parameters---
RPM range: crank angle degree's to simulate after initiation of simulation (from BDC)
Working directory: where the chemical mechanism, *.yml and *.py file are located

---Miscellaneous Parameters---
Numer of Cores: Number of CPU cores to utilize for parallel processing
Print progress to cmd: frequency of print-to-screen values. Set to 0 to diable, higher the value means more print-to-screen values 

