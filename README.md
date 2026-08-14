# egusphere-2026-2142
*repo for data and scripts to generate figures for egusphere-2026-2142 manuscript*

# Contents
1. Regression_Metrics.m - MATLAB function characterizing linear and quadratic regressions
2. aps_soars_polar_plc_output_202402.xlsx - APS PLC results
  a. Size resolved efficiency and loss given as a percent, diameter is in units of µm
  b. 200x3 spreadsheet. Diameter [µm], Efficiency [%], Loss [%]
4. egusphere_2026_2142_figures.m - matlab figure script
5. ops_soars_polar_plc_output_202402.xlsx - OPS PLC results
  a. Size resolved efficiency and loss given as a percent, diameter is in units of µm
  b. 200x3 spreadsheet. Diameter [µm], Efficiency [%], Loss [%]
7. smps_soars_polar_plc_output_202402.xlsx - SMPS PLC results
  a. Size resolved efficiency and loss given as a percent, diameter is in units of µm
  b. 200x3 spreadsheet. Diameter [µm], Efficiency [%], Loss [%]
9. soars_polar_bubbles_final.mat - subsurface bubble data ("all_stats" matlab structure)
  a. SST binned results
10. soars_polar_final.mat - smps-aps data ("all_stats" matlab structure)
  a. SST binned results
  b. 1x5 cell, each cell is a 1x1 structure containing:
   1. D_mode [µm]
   2. Dg [µm]
   3. N [cm^{-3}]
   4. N_accumulation [cm^{-3}]
   5. N_aitken [cm^{-3}]
   6. N_coarse [cm^{-3}]
   7. N_submicron [cm^{-3}]
   8. N_supermicron [cm^{-3}]
   9. S [µm^2 cm^{-3}]
   10. T [°C]
   11. V [µm^2 cm^{-3}]
   12. dN [cm^{-3}] indicating dN/dlogDp
   13. dS [µm^2 cm^{-3}] indicating dS/dlogDp
   14. dV [µm^3 cm^{-3}] indicating dV/dlogDp
   15. sigma_g [unitless]
   16. D [µm] indicating dry, physical diameter
   17. num_observations
   18. notes
*All fields are 1x1 structures containing fields "mean," "std," and "SE," except D, num_observations, notes* 
11. soars_polar_sm_ops.mat - sm-ops data ("all_stats" matlab structure)
  a. SST binned results
  b. 1x25 cell, each cell is a 1x1 structure containing the same variables listed above
13. soars_polar_smops_final.mat - sm-ops data ("all_stats" matlab structure)
  a. Full sm-ops data, not SST binned
  b. 1x25 cell, each cell is a 1x1 structure containing the same variables listed above
15. soars_polar_smps_aps.mat - smps-aps data ("all_stats" matlab structure)
  a. Full smps-aps data, not SST binned
  b. 1x25 cell, each cell is a 1x1 structure containing the same variables listed above
17. spider_soars_polar_plc_output_202402.xlsx - SM PLC results
  a. Size resolved efficiency and loss given as a percent, diameter is in units of µm
  b. 200x3 spreadsheet. Diameter [µm], Efficiency [%], Loss [%]

