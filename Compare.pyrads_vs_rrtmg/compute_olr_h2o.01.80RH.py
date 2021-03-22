import numpy as np
import sys,os

sys.path.append("..")
import pyrads
import python_rrtmg  # !

from scipy.integrate import trapz,simps
import copy

### -----------------------------------
### Helpers
class Dummy:
    pass


### -----------------------------------

# ---
## setup thermodynamic parameters
params0 = Dummy()

params0.Rv = pyrads.phys.H2O.R # moist component
params0.cpv = pyrads.phys.H2O.cp
params0.Lvap = pyrads.phys.H2O.L_vaporization_TriplePoint
params0.satvap_T0 = pyrads.phys.H2O.TriplePointT
params0.satvap_e0 = pyrads.phys.H2O.TriplePointP
params0.esat = lambda T: pyrads.Thermodynamics.get_satvps(T,params0.satvap_T0,params0.satvap_e0,params0.Rv,params0.Lvap)

params0.R = pyrads.phys.air.R  # dry component
params0.cp = pyrads.phys.air.cp
params0.ps_dry = 1e5           # surface pressure of dry component

params0.g = 9.8             # surface gravity
params0.cosThetaBar = 3./5. # average zenith angle used in 2stream eqns
params0.RH = 0.8             # relative humidity; see Kluft et al

params0.drymmw = pyrads.phys.air.MolecularWeight   # needed for RRTM!

params0.R_CO2 = pyrads.phys.CO2.R


# ---
## setup resolution (vertical,spectral)

N_press = 40       # 
wavenr_min = 0.1   # [cm^-1]
wavenr_max = 2500. #
dwavenr = 0.01     # 

Tstrat = 200.      # stratospheric temperature

## setup range of temperatures, and if/where output is saved to:

Ts_grid = np.arange(280.,320.1,5.)
filename = 'output.compute_olr_h2o.01.80RH.txt'

saveOutput = True  # Save the output/plots? [Yes/No]
if saveOutput:
    OUTDIR = "./"
    print( "Saving output to ",OUTDIR )
    if not os.path.isdir( OUTDIR ):
        os.makedirs( OUTDIR )


### -----------------------------------
## MAIN LOOP

# save resolution etc for a given loop to file:
if saveOutput:
    f = open(OUTDIR+filename,'w')
    f.write("wavenr_min,wavenr_max,dwave [cm^-1] = %.4f,%.4f,%.4f" % (wavenr_min,wavenr_max,dwavenr) )
    f.write("\n")
    f.write("N_press = %.1f" % N_press )
    f.write("\n")
    f.write("\n")
    f.write("Ts [K],\tps [bar],\tOLR pyrads [W/m2],\tdOLR/dTs pyrads [W/m2],\tCO2forcing pyrads [W/m2],\tOLR rrtmg [W/m2],\tdOLR/dTs rrtmg [W/m2],\tCO2forcing rrtmg [W/m2]")
    f.write("\n")

    f.close()
        
    ## main loop here
    for Ts in Ts_grid:
        f = open(OUTDIR+filename,'a')

        # ..
        params = copy.deepcopy(params0)
        params.ppv_CO2 = 348e-6      # value from Kluft et al
        params_forcing = copy.deepcopy(params)
        params_forcing.ppv_CO2 = 2.*348e-6   # ...
        
        # setup grid:
        dTs = 1.
        g1 = pyrads.SetupGrids.make_grid( Ts,Tstrat,N_press,wavenr_min,wavenr_max,dwavenr,params, RH=params.RH, pTOA_decade=-1.75 )
        g2 = pyrads.SetupGrids.make_grid( Ts+dTs,Tstrat,N_press,wavenr_min,wavenr_max,dwavenr,params, RH=params.RH, pTOA_decade=-1.75 )
        g1_forcing = copy.deepcopy(g1)

        # compute optical thickness:
        #   -> this is the computationally most intensive step
        g1.tau,tmp,tmp2 = pyrads.OpticalThickness.compute_tau_H2ON2_CO2dilute(g1.p,g1.T,g1.q,params.ppv_CO2,g1,params, RH=params.RH )
        g2.tau,tmp,tmp2 = pyrads.OpticalThickness.compute_tau_H2ON2_CO2dilute(g2.p,g2.T,g2.q,params.ppv_CO2,g2,params, RH=params.RH )
        g1_forcing.tau,tmp,tmp2 = pyrads.OpticalThickness.compute_tau_H2ON2_CO2dilute(g1_forcing.p,g1_forcing.T,g1_forcing.q,params_forcing.ppv_CO2,g1_forcing,params_forcing, RH=params_forcing.RH )

        # compute Planck functions etc:
        #   -> here: fully spectrally resolved!
        g1.T_2D = np.tile( g1.T, (g1.Nn,1) ).T               # [press x wave]
        g1.B_surf = np.pi* pyrads.Planck.Planck_n( g1.n,g1.Ts )     # [wave]
        g1.B = np.pi* pyrads.Planck.Planck_n( g1.wave, g1.T_2D )    # [press x wave]

        g2.T_2D = np.tile( g2.T, (g2.Nn,1) ).T               # [press x wave]
        g2.B_surf = np.pi* pyrads.Planck.Planck_n( g2.n,g2.Ts )     # [wave]
        g2.B = np.pi* pyrads.Planck.Planck_n( g2.wave, g2.T_2D )    # [press x wave]

        g1_forcing.T_2D = np.tile( g1_forcing.T, (g1_forcing.Nn,1) ).T               # [press x wave]
        g1_forcing.B_surf = np.pi* pyrads.Planck.Planck_n( g1_forcing.n,g1_forcing.Ts )     # [wave]
        g1_forcing.B = np.pi* pyrads.Planck.Planck_n( g1_forcing.wave, g1_forcing.T_2D )    # [press x wave]

        # compute OLR etc:
        olr_spec1 = pyrads.Get_Fluxes.Fplus_alternative(0,g1) # (spectrally resolved=irradiance)
        olr1 = simps(olr_spec1,g1.n)

        olr_spec2 = pyrads.Get_Fluxes.Fplus_alternative(0,g2) # (spectrally resolved=irradiance)
        olr2 = simps(olr_spec2,g2.n)

        olr_spec_forcing = pyrads.Get_Fluxes.Fplus_alternative(0,g1_forcing) # (spectrally resolved=irradiance)
        olr_forcing = simps(olr_spec_forcing,g1_forcing.n)

        lambda_pyrads = (olr2 - olr1) / dTs
        forcing_pyrads = (olr_forcing - olr1) / dTs

#         olr1 = np.nan
#         lambda_pyrads = np.nan
#         forcing_pyrads = np.nan

        print( "\n",Ts,g1.ps/1e5,olr1,forcing_pyrads, "\n" )


        # ============================
        # NOW: compare with RRTM
        python_rrtmg.make_input(g1.p,g1.T,g1.q,g1.Ts,g1.ps,params,filename="INPUT_RRTM1")
        python_rrtmg.run_rrtm("INPUT_RRTM1","OUTPUT_RRTM1")
        level,plevel,Fup,Fdown,Fnet,Q = python_rrtmg.read_rrtm_output("OUTPUT_RRTM1")
        olr_rrtm1 = Fnet[0]

        #
        python_rrtmg.make_input(g2.p,g2.T,g2.q,g2.Ts,g2.ps,params,filename="INPUT_RRTM2")
        python_rrtmg.run_rrtm("INPUT_RRTM2","OUTPUT_RRTM2")
        level,plevel,Fup,Fdown,Fnet,Q = python_rrtmg.read_rrtm_output("OUTPUT_RRTM2")
        olr_rrtm2 = Fnet[0]

        #
        python_rrtmg.make_input(g1_forcing.p,g1_forcing.T,g1_forcing.q,g1_forcing.Ts,g1_forcing.ps,params_forcing,filename="INPUT_RRTM3")
        python_rrtmg.run_rrtm("INPUT_RRTM3","OUTPUT_RRTM3")
        level,plevel,Fup,Fdown,Fnet,Q = python_rrtmg.read_rrtm_output("OUTPUT_RRTM3")
        olr_rrtm_forcing = Fnet[0]
        
        #
        lambda_rrtm = (olr_rrtm2 - olr_rrtm1)/dTs
        forcing_rrtm = (olr_rrtm_forcing - olr_rrtm1)/dTs

        # ============================
        # ...
        f.write("%.2f,\t%.4f,\t%.8f,\t%.8f,\t%.8f,\t%.8f,\t%.8f,\t%.8f" % (Ts,g1.ps/1e5,olr1,lambda_pyrads,forcing_pyrads,olr_rrtm1,lambda_rrtm,forcing_rrtm) )
        f.write("\n")

        f.close()
    
