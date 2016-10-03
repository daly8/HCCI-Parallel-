# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 09:51:10 2016
Writing a YAML file for HCCI LTC index data
@author: Daly
"""

import yaml

data = {'Fuel names': ['FACE_D_Reduced_Daly'],
        'Fuel species': [['CPT','NC7H16','C6H5CH3','O-XYL','IC5H12','IC6H14','IC8H18']],
        'Species mole fractions':[[0.0689388,0.0858836,0.0808812,0.300643,0.161949,0.0769644,0.22474]],
        'Mechanism': 'FACE_D_Reduced_mech_sk.cti',
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

with open('HCCI_LTC_data.yml', 'w') as outfile:
    yaml.dump(data, outfile, default_flow_style=False)
