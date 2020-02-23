import pandas
import numpy as np
from pychemkin import chemkin, ChemSolver, ChemViz, simpleIO, InputParser
import xml.etree.ElementTree as ET
from memory_profiler import profile

@profile
def run(myFileString):
    global file_name
    T = 273.15
    x_init = np.array([100,0,1,1]) #initial concentration
    t_max = 300 #integration end time in seconds
    dt = 0.1 #step size in seconds
    #rxn_system = chemkin(filename)
    cs = ChemSolver(chemkin(file_name))
    cs.solve(x_init,T,t_max,dt, algorithm='lsoda')

    #time_array, conc_array, rxnrate_array = cs.get_results()

    cs.save_results('SimulationData/TaskDSimulationData' + myFileString + '.csv')

    cv = ChemViz(cs)
    cv.plot_time_series('concentration', tmin=0, tmax=300, species=['a','b','c','d'], outputfile='ConcentrationGraphs/TaskDConcentrationSeries' + myFileString + '.png')

    #cv.plot_time_series('reactionrate', tmin=1.e-13, tmax=300, outputfile='TaskDReactionRateSeries' + str(number) + '.png')
    #cv.html_report('HTMLReports/TaskDHTMLReport' + myFileString + '.html')
    # simpleIO('cs.pkl').to_pickle(cs)
    # cs2 = simpleIO('cs.pkl').read_pickle()
    return

    
@profile
def iterate():
    global file_name
    step_size = 1.e0
    tree = ET.parse(file_name)
    root = tree.getroot()
    constants = [i for i in root.iter(tag='k')] #Obtain elements for constants.  K values are text value of the element. 
    
    myFileString = '&#9k2='+str(constants[1].text)+',k3='+str(constants[2].text)
    
    if float(constants[1].text) < 12: # If k2 is less than 12 run, if not end.
        if float(constants[2].text) < 12: # if k3 is less than 12 iterate it, if its greater than 12 iterate k2 and set k3 to start number  (0.1)
            constants[2].text = str(float(constants[2].text) + step_size)
            tree.write(file_name)
            return myFileString
        else:
            constants[1].text = str(float(constants[1].text) + step_size)
            constants[2].text = str(float(1.e-1))
            tree.write(file_name)
            return myFileString
            
    else:
        global finished
        finished = True

if __name__ == "__main__":
    finished = False
    file_name = 'TaskD.xml'
    while not finished:
        myFileString = iterate()
        run(myFileString)
    input('SIMULATION FINISHED...PRESS ANY KEY TO EXIT')

