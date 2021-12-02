import numpy
import os,copy,sys
#import phys
sys.path.append('../pyrads')
from pyrads import phys

'''
A number of utilities to run command-line RRTM with python.
dkoll (25.10.2017)
'''


### -----------------------------------
### Global definitions here
rrtm_exe = "../rrtmg_lw/build/rrtmg_lw_v4.85_gfortran"


### -----------------------------------
### Given temperature-pressure profiles, construct helper profiles (ptop,Ttop, etc),
### create input txt file for RRTM.
### By default, assume {p,T}[0] is TOA, {p,T}[-1] is near surface.
###
### INPUT: T,p,q. Ts,ps. 'params' is placeholder object.
###     Assume [p,ps] = Pa.
###
### OPTIONAL INPUT:  'Tbot,Ttop'
###     If T mid-values aren't given, construct via simple interpolation
###                  'filename'
###     If filename not given, save to 'INPUT_RRTM' text file
###
### OUTPUT: written to text file
###
### NOTE:  last value in p,T must be above surface!
###     

def make_input(p,T,q,Ts,ps,params,filename=None,flipInput=True,Tbot=None,Ttop=None):
    N_lev = len(p)

    pbot = numpy.concatenate( ( (p[0:-1]+ p[1:])/2.,[ps]) )
    ptop = numpy.concatenate( ([0.], (p[0:-1]+ p[1:])/2.) )

    if Tbot==None and flipInput:
        Tbot = numpy.concatenate( ( (T[0:-1]+ T[1:])/2.,[Ts]) )
    elif Tbot==None and ~flipInput:
        Tbot = numpy.concatenate( ([Ts],(T[0:-1]+ T[1:])/2.) )

    if Ttop==None and flipInput:
        Ttop = numpy.concatenate( ( [T[0]], (T[0:-1]+ T[1:])/2. ) )
    elif Ttop==None and ~flipInput:
        Ttop = numpy.concatenate( ( (T[0:-1]+ T[1:])/2.,[T[0]] ) )


    if flipInput:
        T = numpy.flipud(T)    ## RRTM assumes surface comes first -> flip all arrays?
        p = numpy.flipud(p)
        q = numpy.flipud(q)
        pbot = numpy.flipud(pbot)
        ptop = numpy.flipud(ptop)
        Tbot = numpy.flipud(Tbot)
        Ttop = numpy.flipud(Ttop)

    p = p/1e2        # convert all pressures to hPa
    ptop = ptop/1e2
    pbot = pbot/1e2

    ## compute molecular abundances
    ##    -> here only H2O, CO2, and N2
    ##    -> careful: this snippet makes a number of simplifications that break down at high CO2 or H2O
    R_mean = (1.-q)*params.R + q*params.Rv
    pH2O = q * params.Rv/R_mean   # !
    pCO2 = p*0. + params.ppv_CO2  # !
    pO3 = p*0.
    pN2O = p*0.
    pCO = p*0. 
    pCH4 = p*0.
    pO2 = p*0. + 0.21*(1.-q)
    N = 0.79*(1.-q) * 1e2 * (pbot - ptop)/params.g * phys.N_avogadro/params.drymmw * 0.1   # inert broadening gas
    #      above: [1e2 from hPa -> Pa, 0.1 from kg/m2 -> g/cm2; here [mmw] = g/mol!]


    ###     --- MAKE RRTM_INPUT FILE ---
    if filename==None:
        filename = "INPUT_RRTM"

    f = open(filename,'w')

    #0        1         2         3         4         5         6         7         8         9
    #123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-12345
    #$ TITLE HERE
    #                                              IATM              IXSECT        ISCAT *    &    #
    # Here:   * = ISCAT, & = IOUT, # = ICLD.      See 'rrtm_instructions' file!
    f.write("$ Custom sounding\n")
    f.write("                                                 0                   0            0 0   00    0\n")


    #       TBOUND,  IEMIS, IREFLECT, (SEMISS(IB),IB=1,16)
    #         1-10,     12,       15,        16-95
    #         E10.3 1X, I1,   2X, I1,        16E5.3
    f.write( "{:>10.3f}".format( Ts ) + " 0  0" + "\n"  )


    #          IFORM, NLAYRS, NMOL
    #            2     3-5,   6-10
    #          1X,I1    I3,    I5
    f.write( " 1" + "{:>3d}".format(N_lev) + "   7"  + "\n"  )
    
    # NOTE:
    #   IFORM      (0,1) column amount format flag
    #   = 0  read PAVE, WKL(M,L), WBROADL(L) in F10.4, E10.3, E10.3 formats (default)
    #   = 1  read PAVE, WKL(M,L), WBROADL(L) in E15.7 format


    for i in range(N_lev):
        # PAVE,  TAVE,    PZ(L-1),  TZ(L-1),   PZ(L),  TZ(L)
        # 1-10, 11-20,     44-51,    52-58,    66-73,  74-80
        # F10.4, F10.4, 23X, F8.3,    F7.2,  7X, F8.3,   F7.2
        #    --> same as 
        # p(i), T(i), pb(i), Tb(i), pt(i), Tt(i)
        f.write( "%15.7E"%p[i] + "%10.4f"%T[i] + " "*23 + "%8.3f"%pbot[i] + "%7.2f"%Tbot[i]+ " "*7 + "%8.3f"%ptop[i] + "%7.2f"%Ttop[i] + "\n" )

        #   (WKL(M,L), M=1, 7), WBROADL(L)
        #       (8E10.3)
        #    --> same as
        # pH20  pCO2  pO3   pN2O  pCO   pCH4  pO2   total
        f.write( "%15.7E"%pH2O[i]+ "%15.7E"%pCO2[i] + "%15.7E"%pO3[i] + "%15.7E"%pN2O[i] + "%15.7E"%pCO[i] + "%15.7E"%pCH4[i] + "%15.7E"%pO2[i] + "%15.7E"%N[i] +"\n" )
    
    f.write("%%\n")
    f.write("%%\n")
    f.write("%%\n")
    f.close()


### -----------------------------------
### Given input txt file, run RRTM on it.
###
### INPUT:  input.txt, output name
### OUTPUT: write to output.txt
###
def run_rrtm(file_in,file_out,clean_input=False):
    os.system("cp " + file_in + " INPUT_RRTM" )
    os.system( rrtm_exe )
    os.system("mv OUTPUT_RRTM " + file_out )

    if clean_input:
        os.remove('INPUT_RRTM')  # clean up tmp file


### -----------------------------------
### Given output txt file, extract values
###

def read_rrtm_output(file):
    data = numpy.genfromtxt(file, skip_header=3, skip_footer=18)

    level= data[:,0]
    plevel = data[:,1]
    Fup = data[:,2]
    Fdown = data[:,3]
    Fnet = data[:,4]
    Q = data[:,5]

    return level,plevel,Fup,Fdown,Fnet,Q
