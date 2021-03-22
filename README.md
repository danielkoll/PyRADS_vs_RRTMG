Compare RRTMG and PyRADS at high temperature.

(Daniel Koll, March 2021)

Installation, on a Mac OS X with conda:
- Follow the PyRADS installation instructions here: https://github.com/ddbkoll/PyRADS
- Compile RRTMG:
  - cd $PyRADS_vs_RRTMG/rrtmg_lw/build
  - make -f makefiles/make_rrtmg_lw_OS_X_gfortran_new
  - make sure the tests in $PyRADS_vs_RRTMG/rrtmg_lw/runs_test work

To run the scripts:
- cd $PyRADS_vs_RRTMG/Compare.pyrads_vs_rrtmg
- python compute_olr_h2o.01.80RH.py
- python make_plots.py
