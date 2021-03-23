# PyRADS versus RRTMG
Compare PyRADS and RRTMG at high surface temperatures. 

Installation, for Mac Os X with conda:
- follow the PyRADS instructions: https://github.com/ddbkoll/PyRADS
- compile RRTMG:
    - cd $PyRADS_vs_RRTMG/rrtmg_lw/build
    - make -f makefiles/make_rrtmg_lw_OS_X_gfortran_new
    - run RRTMG tests to make sure it compiled correctly; see $PyRADS_vs_RRTMG/rrtmg_lw/runs_test
- Make sure the python wrapper script for RRTMG points to the correct location.
    See: https://github.com/ddbkoll/PyRADS_vs_RRTMG/blob/main/Compare.pyrads_vs_rrtmg/python_rrtmg.py
    Set the correct executable name in 'rrtm_exe'.
    
Run comparison:
- cd $PyRADS_vs_RRTMGCompare.pyrads_vs_rrtmg
- python compute_olr_h2o.01.80RH.py    [this will take a while]
- python make_plots.py

Outputs:
- (Compare.pyrads_vs_rrtmg/plot_forcing.pdf)
- (Compare.pyrads_vs_rrtmg/plot_feedback.pdf)
- (Compare.pyrads_vs_rrtmg/plot_ecs.pdf)

