import numpy as np
import matplotlib.pyplot as plt

import phys

### -----------------------------------
### Helpers



### -----------------------------------
### 

data = np.genfromtxt('output.compute_olr_h2o.01.80RH.txt',skip_header=4,delimiter=',')
Ts = data[:,0]
olr_pyrads = data[:,2]
feedback_pyrads = data[:,3]
forcing_pyrads = data[:,4]
olr_rrtmg = data[:,5]
feedback_rrtmg = data[:,6]
forcing_rrtmg = data[:,7]


### -----------------------------------
### 
saveOutput = True
if saveOutput:
    OUTDIR = "./"
    print( "Saving output to ",OUTDIR )


tmp = np.linspace( Ts.min(),Ts.max() )


# ---
plt.figure()
#
plt.plot(Ts,olr_pyrads,"-",lw=2,label="PyRADS")
plt.plot(Ts,olr_rrtmg,"-",lw=2,label="RRTMG")
#
plt.legend(loc="upper left")
plt.xlabel("Surface Temperature (K)")
plt.ylabel("OLR (W/m$^2$)")
if saveOutput: plt.savefig("plot_olr.pdf",format="pdf")


# ---
plt.figure()
#
plt.plot(Ts,-1*feedback_pyrads,"-",lw=2,label="PyRADS")
plt.plot(Ts,-1*feedback_rrtmg,"-",lw=2,label="RRTMG")
plt.plot(tmp,tmp*0.,"k--")
#
plt.legend(loc="upper left")
plt.xlabel("Surface Temperature (K)")
plt.ylabel("Feedback (W/m$^2$/K)")
plt.ylim([-2.2,0.2])
if saveOutput: plt.savefig("plot_feedback.pdf",format="pdf")


# ---
plt.figure()
#
plt.plot(Ts,-1*forcing_pyrads,"-",lw=2,label="PyRADS")
plt.plot(Ts,-1*forcing_rrtmg,"-",lw=2,label="RRTMG")
#plt.plot(tmp,tmp*0.,"k--")
#
plt.legend(loc="upper left")
plt.xlabel("Surface Temperature (K)")
plt.ylabel("Feedback (W/m$^2$/K)")
#plt.ylim([3.5,4.8])
if saveOutput: plt.savefig("plot_forcing.pdf",format="pdf")


# ---
# ECS, estimated from: ECS ~ forcing / feedback
plt.figure()
#
plt.plot(Ts,-1*forcing_pyrads/feedback_pyrads,"-",lw=2,label="PyRADS")
plt.plot(Ts,-1*forcing_rrtmg/feedback_rrtmg,"-",lw=2,label="RRTMG")
plt.plot(tmp,tmp*0.,"k--")
#
plt.legend(loc="upper left")
plt.xlabel("Surface Temperature (K)")
plt.ylabel("(K)")
plt.title("Estimated ECS = forcing/feedback")
#plt.ylim([-2.2,0.2])
if saveOutput: plt.savefig("plot_ecs.pdf",format="pdf")

#plt.show()

