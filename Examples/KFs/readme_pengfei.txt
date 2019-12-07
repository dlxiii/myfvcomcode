date: 12/30/2012
3 folders
true run: generate the true field
EnKF_assimilation_run_official_source: Official EnKF assimilation test (the official code was modified based on FVCOM31_source formal version)

NOTICE: in the EnKF_assimilation_run_official_source, I modify the subroutine fvcom.F ,mod_enkf.F,mod_enkfassim.F, use cvs status to check

How to run "true run"?
go to true_run/Run: mpiexec -n 8 ./fvcom --casename=tst
this is a case of circlular region with 10m wate depth and 1m amplitude of m2 tide.
and we use its model output as true state and 20 restart files generated from  day 10 ,iint 24300 backwards with 100 step interval for assimialtion test.


How to run "assimilation run"  
go to EnKF_assimilation_run_official_source/Input:
1). run cp.sh to copy initial ensemble restart files from the true run.
2). edit the tst_assim_enkf.nml 

go to EnKF_assimilation_run_official_source/output (if not existed, mkdir output):
1) in the folder output: mkdir 2 folders: enkftemp and out_err
2) in the folder enkftemp, link file from true run as true state use ln -s ../../../true_run/output/tst_0001.nc true_field.nc

run in the "Run" folder with command: mpiexec -n 8 ./fvcom --casename=tst

How to analysis:
go to EnKF_assimilation_run_official_source/analysis/
see run.sh (it is a readme, not a shell file)
1.run for_true_analysis_contourplot.m
2.source grrider.sh
3.gri zeta.gri
4.gri err_plot.gri
5. gri grid_plot.gri
