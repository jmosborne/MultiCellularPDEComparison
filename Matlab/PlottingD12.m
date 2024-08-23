clear; close all

data = readtable('/Users/jmosborne/testoutput/TestDiffusion/TestDiffusionSmall12.dat')

plot(data.time_minutes_,data.average_concentration_only_cells_uM_)