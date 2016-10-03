# -*- coding: utf-8 -*-
"""
@author: Shane
Contact: dalys@onid.oregonstate.edu 
0-D Homogeneous Charge Compression Ignition (HCCI) Engine Simulation
Features Woschnii heat transfer (variable heat transfer coefficient)
"""
from __future__ import division       #Python version 3 features
import numpy                          #Mathematical operations
try:
    import cantera                        #Cantera library
except ImportError as e:
    print(e,'...you need to install Cantera')

from joblib import Parallel, delayed  #Parallel processsing function import 
try:
    import h5py                           #Data storage in HDF5 format
except ImportError as e:
    print(e,'...you need to install h5py')
import os                             #Operating system operations (to create working directories) 
import time
try:                           #Allows us to track processor time
    import yaml
except ImportError as e:
    print(e,'...you need to install PyYAML')
	
import sys
import ctypes    
import argparse
#--------------------->HCCI Simulation Function Code<--------------------------------
def HCCI(RPM,Ti,phi,fuel_species,X_fuel,FuelFolder,FilePrefix, P,mechanism):
    if sys.platform.startswith("win"):
            # Don't display the Windows GPF dialog if the invoked program dies.
            # See comp.os.ms-windows.programmer.win32
            SEM_NOGPFAULTERRORBOX = 0x0002 # From MSDN
            ctypes.windll.kernel32.SetErrorMode(SEM_NOGPFAULTERRORBOX)
    try:
        Pi=P[0]                             #Initial Pressure (BDC)
        "--->Engine Parameters<---"
        rc=P[1]                                #Compression Ratio [Vmax/Vmin]    
        B=P[2]                              #Cylinder Bore 
        Vd=P[3]                      #Displacement Volume of Engine [m^3] 
        c=P[4]                               #ratio of connecting rod to crank radius 
        "--->Woschni heat transfer constants<---"  
        aa=P[5]                             #Nusselt Number Correlation Constant
        bb=P[6]                               #Nusselt Number Correlation Re power constant
        C11=P[7]                            #from Chemkin
        Tw=P[8]                              #Wall (Bore) Temperature    
        Pr=P[9]
        cc=P[10]         
        
        "--->Recorded Time-step and Crank Angle Duration<---" 
        RPM_range=P[11]                        #Crank Angle Degree simulation finishes (CAD aBDC)
        #timesteps=286/.05                   #Recoreded timesteps for crank angle range    
        
        FileName = "%s_RPM=%s_Ti=%s_phi=%s.hdf5" %(FilePrefix,RPM,Ti,phi)                         #Generating a filename based on fuel name and initial parametric study values    
        FileName_failed = "%s_RPM=%s_Ti=%s_phi=%s_FAILED.hdf5" %(FilePrefix,RPM,Ti,phi)    
        if os.path.exists(FuelFolder + FileName):                                   #If the fuel working directory already exists, don't try and overwrite the folder
            #print("This simulation has already been performed...skipping to next case")
            return
        elif os.path.exists(FuelFolder + FileName_failed):                                   #If the fuel working directory already exists, don't try and overwrite the folder
            #print("This simulation has already been performed, which failed previously...skipping to next case")
            return 
        else:
            "--->Chemical Mechanism<---"  
            while True: #This simply suppresses the Exception from mechanism import (Happens with parallel processing forks trying to access file at same time), and allows a re-try without spamming alot of useless information in print window.  
                try:
                    mechanism=cantera.Solution(mechanism) #mechanism import  #Princeton_10atm Sarathy_Mech_Reduced_TRF_isopentane.cti Andrae_Xylene_mech FACE_C_mech_re MCH_mech_skeletal
                except:
                    continue
                break 
              
            #Specifying air for the ideal gas reactor (required for engine simulation) and property evaluations for viscosity, density, and thermal conductivity  
            air=cantera.Solution('air.xml')        #air (the piston model is air separated by a wall(piston) to the reactor)  
            try:    
                io2 =mechanism.species_index('o2')     #Finding species indices 
                in2 = mechanism.species_index('n2')    #Finding species indices 
                ico2 =mechanism.species_index('co2')   #Finding species indices
                ico =mechanism.species_index('co')     #Finding species indices  
                ih2o = mechanism.species_index('h2o')  #Finding species indices 
            except ValueError:
                io2 =mechanism.species_index('O2')     #Finding species indices 
                in2 = mechanism.species_index('N2')    #Finding species indices 
                ico2 =mechanism.species_index('CO2')   #Finding species indices
                ico =mechanism.species_index('CO')     #Finding species indices  
                ih2o = mechanism.species_index('H2O')  #Finding species indices
            #Engine and Heat Transfer Parameters
            s=4*Vd/(numpy.pi*B**2)                 #Stroke [m]
            Vmax=Vd*rc/(rc-1)                      #Maximum cylinder Volume [m^3]
            Ap=(numpy.pi*B**2)/4                   #Piston Area
            wdot=RPM*numpy.pi/30                   #rotation rate in rad/s 
            #Simulation Resolution     
            t_end=RPM_range/(6*RPM) #Simulation time domain of crank angle's to be resolved 
            #theta=t*wdot*180/numpy.pi                       #Simulation crank angle (Just a recorded value, not used in any calculations)
            "slider crank rule for piston velocity" #Cantera calls for piston(reactor wall) velocity, not Volume(t) or dV/dt (for the energy equation)
            def piston_speed(t):
                return (wdot*s/2*numpy.sin(t*wdot))*(1-(numpy.cos(t*wdot))/(numpy.sqrt(c**2-numpy.sin(t*wdot)**2)))
            "Piston Volume Function" #(used for heat transfer Area calculation) 
            def Volume(t):
                return Vd/(rc-1)+Vd/2*(1+c-numpy.cos(t*wdot-numpy.pi)-numpy.sqrt((c)**2-numpy.sin(t*wdot-numpy.pi)**2))
            
            "Heat Transfer Function" 
            #Chemkin Property Evaluation constants
            C1_mu= 1.1258e-05
            C2_mu= .485
            C1_k= 132.7906
            C2_k= .5
            def q(t):
                rho = r2.kinetics.density_mass    
                mu=(C1_mu*(r2.T)**(C2_mu))*0.1       #Viscosity of air
                k_g=(C1_k*(r2.T)**(C2_k))*10**-(5)   #Thermal Conductivity of air [W/m-K]         
                Sm=2.0*s*RPM/60                #Average Piston Speed
                w_bar=C11*Sm                   #Woschni average cylinder gas velocity
                Re=B*w_bar*rho/mu              #Reynolds Number
                Nu=aa*Re**bb*Pr**cc            #Nusselt Number
                h=Nu*k_g/(B)                   #convective heat transfer coefficient
                return (h*(Tw-r2.T)*(2*numpy.pi*B**3+16*Volume(t))/(4*B))/Ap 
        
                
            def net_heat_release_per_deg(t):
                return (-Volume(t)*numpy.pi*sum(r2.kinetics.partial_molar_enthalpies*r2.kinetics.net_production_rates)+q(t)*Ap)/(wdot*180)
            
                
            attempt = 0    
            while attempt <= 2:  #if the simulation integrator fails, re-try at relaxed tolerance 3 times... 
                "Fuels and composition Parameters" #(Solves for Stoich F/A micture)
                length_f=len(fuel_species)                  #Defining this value for loops below 
                ifuel=numpy.zeros(length_f)                 #pre-allocating an array for fuel indices 
                stoich_o2_fuel=numpy.zeros(length_f)        #pre-allocating an array for stoich O2 amounts for each pure component in the fuel 
                x = numpy.zeros(len(mechanism.X))           #pre-allocating an array for mole fractions to later define the relative mole fractions for the chemical mechanism
                output_variables=[]                         #pre-defining an empty list for the string list of output variables 
                for i in range(length_f):
                    "Solving for stoich O2 mixture"        
                    ifuel[i] = mechanism.species_index(fuel_species[i])  #capturing the proper indiced in the chemical mechanism for the fuel component 
                    stoich_o2_fuel[i]=X_fuel[i]*(mechanism.n_atoms(fuel_species[i],'C')+mechanism.n_atoms(fuel_species[i],'H')/4-mechanism.n_atoms(fuel_species[i],'O')/2) #solving for stoich amount of oxygen for each fuel pure component
                    "Creating fuel composition array for mechanism input" 
                    x[ifuel[i]] = X_fuel[i]*phi                              #relative amount of fuel to air 
                    "Creating first performance parameters as fuel mole fraction" 
                    output_variables.append('X'+fuel_species[i])             #the fuel mole fractions "output variable names" are being added to the list 
                "Stoich O2 Mixture"
                stoich_o2=numpy.sum(stoich_o2_fuel)                          #The above solved for the stoich mixture of the fuels components separately. The final answer is the sum of the oxygen amounts 
                x[io2] = stoich_o2                                           #stoichiometric A/F amount of Oxygen               
                x[in2] = stoich_o2*(0.79/0.21)                               #stoichiometric A/F amount of Nitrogen 
                "Setting Fuel and air mixture composition in mechanism" 
                mechanism.TPX=Ti,Pi,x[:]                                      #The chemical mechanism is set to have the user defined initial conditions 
                r2 = cantera.IdealGasReactor(mechanism,  name='HCCI Reactor') #specifying that r2 is an ideal gas reactor (varying temperature and pressure) 
                r2.volume=Vmax                                                #Specifying initial volume of ideal gas reactor (BDC volume)
                
                "Populating output variables -- columns are each performance variable and rows indicate crank angle"     
                output_variables = output_variables + ['Xo2','Xco','Xco2','Xn2' ,'Xh2o', 'pres','temp','volume','crank_rotation_angle','net_heat_release_rate_per_CA','mass','molecular_weight','mass_density','specific_heat','mixture_enthalpy','time','chem_heat_release_per_deg','net_heat_release_per_deg','heat_loss','rho']
                Performance_Parameters=numpy.zeros([1,len(output_variables)]) #The matrix of numerical values corresponding to the user selected reactor parameters    
                Performance_Parameters_add=numpy.zeros([1,len(output_variables)])
                "------>Cantera HCCI Simulation<----------"
                
                env = cantera.Reservoir(air, name='Environment (Air at the specified heat transfer correlation temperature')
                cantera.Wall(env, r2,velocity=piston_speed,A=Ap,Q=q)   #This is the "piston" in the simulation. The "wall" is moving and has heat transfer through it 
                sim = cantera.ReactorNet([r2])  #this is the completed cantera reactor that will simulate HCCI 
                time = sim.time                 #starting simulation time (default 0.0s for sim.time)   
                n=0                               
                end_time = t_end           #Simulation End-time
                if attempt == 0:
                    sim.rtol=1e-7                   #Sets maximum integration relative tolerance 
                    sim.atol=1e-11                  #Sets maximum integration absolute tolerance        
                count=0               
                heat_release = []    
                
                while time < end_time:          #Cantera Simulation integraiton loop  
                    heat_release.append(net_heat_release_per_deg(time))        
                    if count ==1:
                        Performance_Parameters = numpy.concatenate((Performance_Parameters, Performance_Parameters_add), axis=0)
                    
                    "Reactor Parameters to get"
                    for ff in range(length_f):            
                        if ff<=length_f:
                            Performance_Parameters[n,ff] = r2.kinetics.X[ifuel[ff]]
                    Performance_Parameters[n,(ff+1):] = r2.kinetics.X[io2],r2.kinetics.X[ico],r2.kinetics.X[ico2],r2.kinetics.X[in2],r2.kinetics.X[ih2o],r2.thermo.P,r2.thermo.T,r2.volume,time*wdot*180/numpy.pi,r2.kinetics.enthalpy_mass,r2.mass,r2.kinetics.mean_molecular_weight,r2.kinetics.density,r2.kinetics.cp,r2.kinetics.enthalpy_mass,time,-r2.volume*numpy.pi*sum(r2.kinetics.partial_molar_enthalpies*r2.kinetics.net_production_rates)/(wdot*180),net_heat_release_per_deg(time),q(time)*Ap,r2.kinetics.density_mass  
                 
        
                    "Advancing Simulation" 
                    try:                    
                        sim.step(end_time)
                    except Exception as e:
                        attempt +=1
                        if attempt <=2:
                            print(e, "re-trying simulation at relaxed tolerance")   
                            sim.rtol = sim.rtol*10
                            sim.atol = sim.atol*100
                            break #break "while time < end_time:", go back to top while loop 
                        else:
                            print('Simulation failed...  :(  check files with FAILED tag')
                            FileName = "%s_RPM=%s_Ti=%s_phi=%s_FAILED.hdf5" %(FilePrefix,RPM,Ti,phi)                         #Generating a filename based on fuel name and initial parametric study values 
                            f=h5py.File(FuelFolder + FileName,"w")
                            f.close()
                            return
                    time = sim.time
                    count = 1
                    n=n+1
                break #if the integrator is successful, then break the "re-try simulation" while loop
            #print('HCCI Simulations successful!')
            
            try:
                "net_heat_release analysis to determine frequency of data collection (finding appropriate d_theta)"
                        #If net heat release is greater than some threshhold, record at 0.1 CAD, otherwise collect at 2 CAD.           
                for xx in range(len(output_variables)):
                    if output_variables[xx] == 'crank_rotation_angle':
                        THETA_INDICE = xx
                dtheta=[]
                ind=[]
                p1=0
                for p in range(len(heat_release)-1):
                    dtheta.append(abs(Performance_Parameters[p1,THETA_INDICE] - Performance_Parameters[p,THETA_INDICE]))
                    if heat_release[p] >= 0.2: #Record at 0.1 CAD
                        #Get every 0.1 indice 
                        #If dtheta > 0.1 grab indice and reset beginning d_theta ind.            
                        if dtheta[p] >= 0.1:
                            p1=p #reset the initial d_theta_indice
                            ind.append(p)    #grab indice of simulation at which this happened 
                    else: #Record at 2 CAD 
                        #Get every 2 CAD Indice
                        if dtheta[p] >= 2:
                            p1=p #reset the initial d_theta_indice
                            ind.append(p)    #grab indice of simulation at which this happened            
                    
                "File Specification" 
                FileName = "%s_RPM=%s_Ti=%s_phi=%s.hdf5" %(FilePrefix,RPM,Ti,phi)                         #Generating a filename based on fuel name and initial parametric study values 
                f=h5py.File(FuelFolder + FileName,"w")                                                    #Creating file for each parametric study point #f.name is 'u/'
                "Performance Paramater Saving" 
                for v in range(len(output_variables)):
                    f.create_dataset("%s" %(output_variables[v]),data=Performance_Parameters[ind,v])        #creating hdf5 datasets for each performance parameter inside the file created above        
                f.close()                                                                                 #close the file we just created 
                return 
            except Exception as e:
                print(e)
                print('something went wrong with saving parameters...')
                FileName = "%s_RPM=%s_Ti=%s_phi=%s_FAILED.hdf5" %(FilePrefix,RPM,Ti,phi)                         #Generating a filename based on fuel name and initial parametric study values 
                f=h5py.File(FuelFolder + FileName,"w")
                f.close()
                return
    except Exception as e:
        print(e)
        print('something went wrong with initializing simulation...')
        FileName = "%s_RPM=%s_Ti=%s_phi=%s_ERROR_INIT.hdf5" %(FilePrefix,RPM,Ti,phi)                         #Generating a filename based on fuel name and initial parametric study values 
        f=h5py.File(FuelFolder + FileName,"w")
        f.close()
        return
        
def Command_Line_interface():
    print('Lets set up the fuel parameters we want to calculate LTC index for...')
    time.sleep(2)
    print(" ")
    print(" ")
    while True:
        try:
            print('Please enter a list of fuel names')
            print('e.g. "FACE A","FACE B"')
            fuel_names = [eval(input('fuel list:    '))]
            print("------------------------------------------------------------------------------------")
            print(" ")
            break
        except:
            print("------------------------------------------------------------------------------------")
            print("lets try that again...")
            time.sleep(2)
            print("------------------------------------------------------------------------------------")
            continue
    while True:
        try:
            num_fuels = len(fuel_names) #How many fuel mixtures are in parametric study
            print('Please enter a list of species for each fuel corresponding to the chemical mechanism naming convention')
            print('e.g.: for multiple fuels ["heptane","octane"],["CH4","c2h5oh"]')
            fuel_species = [eval(input('species list:    '))]
            print("------------------------------------------------------------------------------------")
            print(" ")
            break
        except:
            print("------------------------------------------------------------------------------------")
            print("lets try that again...")
            time.sleep(2)
            print("------------------------------------------------------------------------------------")
            continue
    while True:
        try:
            print('Please enter mole fraction for species list(s)')
            print('e.g.: [0.5,0.5],[0.1,0.9]')
            X_fuel=[ eval(input('mole fractions:    '))]
            print("------------------------------------------------------------------------------------")
            print(" ")
            
            break
        except:
            print("------------------------------------------------------------------------------------")
            print("lets try that again...")
            time.sleep(2)
            print("------------------------------------------------------------------------------------")
            continue
        break
    while True:
        try:
            print('Please enter the chemical mechanism file name')
            print('e.g. KAUST_FACE_G_Reduced_Mech.cti')
            mechanism=input('mechanism:    ')
            print("------------------------------------------------------------------------------------")
            print(" ")
            break
        except:
            print("------------------------------------------------------------------------------------")
            print("lets try that again...")
            time.sleep(2)
            print("------------------------------------------------------------------------------------")
            continue
        break
    

    print('################################################################################################')
    print('Would you like to use the default initial conditions ("y") or modify them (type anything)?')
    initial_conditions = input('[y/n]?:    ')
    print("------------------------------------------------------------------------------------")
    print(" ")
    print(" ")
    if initial_conditions == "y":
        "--->Parametric Study Values<---"
        Pi=100000                             #Initial Pressure (BDC)
        dRPM=100.0                            #RPM stepsize
        RPM=numpy.arange(800.0,3100.0,dRPM)   #Revolutions per minute
        dphi=0.01                             #equivalence ratio stepsize 
        phi=numpy.arange(0.15,0.45,dphi)      #Equivalence Ratio
        phi=-numpy.sort(-phi) 
        dT=20                                #Temperature stepsize
        Ti=numpy.arange(350,550+dT,dT)       #Initial Temperature (BDC)
        print('################################################################################################')
        print('------------Parametric Study Values-------------')
        print('Initial Pressure (at BDC) [Pa]:',Pi)
        print("------------------------------------------------")
        print('RPMs:',RPM)
        print("------------------------------------------------")
        print('Equivalence Ratios:',phi)
        print("------------------------------------------------")
        print('Temperatures [K]:',Ti)
        print("------------------------------------------------")
        print('Total number of cases:',len(Ti)*len(phi)*len(RPM))
        print("------------------------------------------------")
        print('################################################################################################')
    else:
        while True:
            Pi=int(input('Initial Pressure (at BDC) [Pa]:'))                             #Initial Pressure (BDC)
            RPM_input = (eval(input('RPMs (initial,final,step size):   ')))   #Revolutions per minute
            RPM=numpy.arange(RPM_input[0],RPM_input[1],RPM_input[2])
            phi_input = (eval(input('Equivalence Ratios (initial,final,step size):   ')))   #Revolutions per minute
            phi=numpy.arange(phi_input[0],phi_input[1],phi_input[2])
            Ti_input = (eval(input('Temperatures [K] (initial,final,step size):   ')))   #Revolutions per minute
            Ti=numpy.arange(phi_input[0],Ti_input[1],Ti_input[2])
            print('################################################################################################')
            print('------------Your parametric Study Values-------------')
            print('Initial Pressure (at BDC) [Pa]:',Pi)
            print("------------------------------------------------")
            print('RPMs:',RPM)
            print("------------------------------------------------")
            print('Equivalence Ratios:',phi)
            print("------------------------------------------------")
            print('Temperatures [K]:',Ti)
            print("------------------------------------------------")
            print('Total number of cases:',len(Ti)*len(phi)*len(RPM))
            print("------------------------------------------------")
            print('################################################################################################')
            print('Does this look right? [y/n]')
            answer = input()
            if answer == 'y':
                break
            else:
                continue
        
    print('')    
    print('Default engine parameters ("y") or modify them (type anything)?')
    engine_conditions = input('[y/n]?:    ')
    print("------------------------------------------------------------------------------------")
    print(" ")
    print(" ")
    if engine_conditions == "y":   
        "--->Engine Parameters<---"
        rc=13                                #Compression Ratio [Vmax/Vmin]    
        B=.0828                              #Cylinder Bore 
        Vd=616*10**(-6)                      #Displacement Volume of Engine [m^3] 
        c=4.44                               #ratio of connecting rod to crank radius 
        "--->Woschni heat transfer constants<---"  
        aa=.035                             #Nusselt Number Correlation Constant
        bb=.8                               #Nusselt Number Correlation Re power constant
        C11=2.28                            #from Chemkin
        Tw=430                              #Wall (Bore) Temperature    
        Pr=0.7
        cc=0 
        print('################################################################################################')
        print('---------Default Engine Parameters (CRC engine)----------')
        print('Compression Ratio:',rc)
        print('Cylinder Bore:',B)
        print('Displacement Volume:',Vd)
        print('Connecting rod to crank radius ratio:',c)
        print('---------Default Heat Transfer Parameters---------------')
        print('Nusselt Number Correlation "a" Constant:',aa)
        print('Nusselt Number Correlation "b" constant:',bb)
        print('Nusselt Number Correlation "c" constant:',C11)
        print('Prandtl Number:',Pr)
        print('Cylinder Wall Temperature [K]:',Tw)
        print('Prandtl Number consant:',cc)
        print('################################################################################################')
        print('')
    else:
        while True:
            rc=int(input('Compression Ratio:'))                                #Compression Ratio [Vmax/Vmin]    
            B=int(input('Cylinder Bore [m]:'))                              #Cylinder Bore 
            Vd=int(input('Displacement Volume [m^3]:'))                      #Displacement Volume of Engine [m^3] 
            c=int(input('Connecting rod to crank radius ratio:'))                               #ratio of connecting rod to crank radius 
            "--->Woschni heat transfer constants<---"  
            aa=int(input('Nusselt Number Correlation "a" Constant:'))                             #Nusselt Number Correlation Constant
            bb=int(input('Nusselt Number Correlation "b" constant:'))                               #Nusselt Number Correlation Re power constant
            C11=int(input('Cylinder gas "C11" constant:'))                            #from Chemkin
            Tw=int(input('Wall temperature:'))                              #Wall (Bore) Temperature 
            cc=int(input('Prandtl Number consant:'))
            Pr=int(input('Cylinder Wall Temperature [K]:'))
            print('################################################################################################')
            print('---------Your Engine Parameters (CRC engine)----------')
            print('Compression Ratio:',rc)
            print('Cylinder Bore:',B)
            print('Displacement Volume:',Vd)
            print('Connecting rod to crank radius ratio:',c)
            print('---------Default Heat Transfer Parameters---------------')
            print('Nusselt Number Correlation "a" Constant:',aa)
            print('Nusselt Number Correlation "b" constant:',bb)
            print('Nusselt Number Correlation "c" constant:',C11)
            print('Prandtl Number:',Pr)
            print('Cylinder Wall Temperature [K]:',Tw)
            print('Prandtl Number consant:',cc)
            print('################################################################################################')
            print('')
            print('Does this look right? [y/n]')
            answer = input()
            if answer == 'y':
                break
            else:
                continue
    
    print('Default recording parameters ("y") or modify them (type anything)?')
    record_conditions = input('[y/n]?:    ')
    print("------------------------------------------------------------------------------------")
    if record_conditions == "y":   
        "--->Recorded Time-step and Crank Angle Duration<---" 
        RPM_range=286                        #Crank Angle Degree simulation finishes (CAD aBDC)
        print('crank angle duration:',RPM_range)
    else:
        RPM_range=input('Input crank angle degrees aBDC to record for simulation:')
        #timesteps=286/.05                   #Recoreded timesteps for crank angle range    
    P = [Pi,rc,B,Vd,c,aa,bb,C11,Tw,Pr,cc,RPM_range] #collection of parameters
    
    print('################################################################################################')
    print('')
    "--->File Specifications<---" 
    while True:
        working_dir = input('Please enter your working directory    :')
        if type(working_dir) == 'str':
            break
        else:
            print('Please enter a "string"')
            continue
    print("------------------------------------------------------------------------------------")
    print(" ")
    print(" ")
    #working_dir = "C:\\Users\\Daly\\Documents\\Python Scripts\\"
    "--->Number of Processing Cores<---" 
    num_cores = int(input('Please enter the number of CPU cores you would like to allocate for this job    :'))
    
    print("------------------------------------------------------------------------------------")
    
    "--->print progress to screen<---" 
    verbosity = int(input('Please enter the frequency of print-progress-to-screen feature (0: disabled, > 1000 everything prints)   :'))
    
    print("------------------------------------------------------------------------------------")
    print(" ")
    print(" ")
    print(" ")
    
    return RPM,Ti,phi,fuel_species,X_fuel,P,mechanism,num_cores,verbosity,working_dir 
    
def arguments():
    parser = argparse.ArgumentParser(description='Performs parallel HCCI simulations') #optional help description for foo
    parser.add_argument('-ii','--interactive_input',
                        help="Use the interactive input argument to input variables in command prompt, as opposed to reading from *.yml file", 
                        action="store_true")
    args = parser.parse_args()
    return args
        
if __name__ == '__main__': 
    cmd_args = arguments()
    if cmd_args.interactive_input:  
        variables = Command_Line_interface()
        RPM=variables[0]
        Ti=variables[1]
        phi=variables[2]
        fuel_species=variables[3]
        X_fuel=variables[4]
        P=variables[5]
        mechanism=variables[6]
        num_cores=variables[7]
        Verbosity = variables[8]
        working_dir = variables[9]
    else:                                                     
        f = open('HCCI_LTC_data.yml')
        dataMap = yaml.load(f)
        f.close()
        fuel_names= dataMap['Fuel names']
        num_fuels = len(fuel_names)
        fuel_species=dataMap['Fuel species']
        X_fuel=dataMap['Species mole fractions']
        mechanism=dataMap['Mechanism']
        Pi=dataMap['Pressure']['value']
        RPM=numpy.arange(dataMap['RPM']['Initial'],dataMap['RPM']['final'],dataMap['RPM']['step'])
        phi=numpy.arange(dataMap['Phi']['Initial'],dataMap['Phi']['final'],dataMap['Phi']['step'])
        Ti=numpy.arange(dataMap['Temperature']['Initial'],dataMap['Temperature']['final'],dataMap['Temperature']['step'])
        rc=dataMap['Compression ratio']['value']
        B=dataMap['Cylinder bore']['value']
        Vd=dataMap['Displacement volume']['value']
        c= dataMap['Connecting rod to crank radius ratio']['value']
        aa=dataMap['Nu Const. A']['value']
        bb= dataMap['Nu Const. B']['value']
        C11=dataMap['Cylinder gas Const.']['value']
        Tw=  dataMap['Wall temperature']['value']
        Pr=dataMap['Prandtl number']['value']
        cc= dataMap['Nu Const. C']['value']
        RPM_range=dataMap['RPM range']['value']
        working_dir=dataMap['Working directory']
        num_cores=dataMap['Number of cores']['value']
        Verbosity = dataMap['Print progress to cmd']['value']
        P = [Pi,rc,B,Vd,c,aa,bb,C11,Tw,Pr,cc,RPM_range]
        
    print("Your HCCI cases will be saved in:   " + working_dir + '\\' + 'CR_' + str(rc) + '\\' + 'the_fuel_names')
    print("------------------------------------------------------------------------------------")
    time.sleep(3)
    print(" ")
    print(" ")
    print(" ")
    print(" ")
    print('Beginning HCCI Simulations!')

    try:
        for q in range(num_fuels):                                               #Fuel Loop 
            FuelFolder= working_dir + "CR_" + str(rc) + "\\" + fuel_names[q]                  #fuel folder to save data in 
            FilePrefix="\\"+fuel_names[q]+"_CR=%s_Pi=%.1f" %(rc,Pi*10**(-5))     #prefix to name of saved file within fuel folder                                                 #Starting processor time of simulation     
            if not os.path.exists(FuelFolder):                                   #If the fuel working directory already exists, don't try and overwrite the folder 
                os.makedirs(FuelFolder)
            #HCCI(RPM[0],Ti[4],phi[0],fuel_species[q],X_fuel[q],FuelFolder,FilePrefix,P)                                  #Creating working directory
            print("Progress will now be printed to screen:")
			#This is me entering a random change to verify the GitHub changes
            Parallel(n_jobs=num_cores, verbose=Verbosity, batch_size = 1)(delayed(HCCI)(RPM[o],Ti[z],phi[k],fuel_species[q],X_fuel[q],FuelFolder,FilePrefix,P,mechanism) for o in range(len(RPM)) for k in range(len(phi)) for z in range(len(Ti))) #parallel process executor 
    except:
        print('There was an error with the input parameters')
