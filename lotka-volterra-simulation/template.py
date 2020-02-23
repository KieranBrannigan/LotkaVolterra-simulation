import pandas
import numpy as np
from pychemkin import chemkin, ChemSolver, ChemViz, simpleIO
T = 1000
x_init = np.ones(8) #initial concentration
t_max = 5.e-13 #integration end time in seconds
dt = 1.e-16 #step size in seconds
cs = ChemSolver(chemkin('rxns_reversible.xml'))
cs.solve(x_init,T,t_max,dt, algorithm='vode', method='bdf', nsteps=500)

time_array, conc_array, rxnrate_array = cs.get_results()

df = cs.to_df()
cs.save_results('simulationdata.csv')

cv = ChemViz(cs)
cv.plot_time_series('concentration', tmin=0, tmax=4.5e-13, species=['H','OH','O2','H2O'], outputfile='modeldocfig1.png')

cv.plot_time_series('reactionrate', tmin=1.e-13, tmax=4.5e-13, outputfile='modeldocfig2.png')

cv.plot_network([0, 1.5e-13, 3e-13], figsize = (8,15), outputfile= 'modeldocfig3.png')

simpleIO('cs.pkl').to_pickle(cs)
cs2 = simpleIO('cs.pkl').read_pickle()
