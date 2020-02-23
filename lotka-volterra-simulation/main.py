import pandas
import numpy as np
from pychemkin import chemkin, ChemSolver, ChemViz, simpleIO, InputParser
import xml.etree.ElementTree as ET
from memory_profiler import profile


class Simulation():
    def __init__(self, filename):
        self.filename = filename
        self.finished = False
        print("init")
    
    def runonce(self):
        global myFileString
        print("start")
        T = 273.15
        x_init = np.array([1.e10,1,1,1]) #initial concentration
        t_max = 300 #integration end time in seconds
        dt = 0.1 #step size in seconds
        cs = ChemSolver(chemkin('TaskD.xml'))
        cs.solve(x_init,T,t_max,dt, algorithm='vode', nsteps=5000)

        #time_array, conc_array, rxnrate_array = cs.get_results()

        cs.save_results('SimulationData/TaskDSimulationData' + myFileString + '.csv')

        cv = ChemViz(cs)
        cv.plot_time_series('concentration', tmin=0, tmax=300, species=['a','b','c','d'], outputfile='ConcentrationGraphs/TaskDConcentrationSeries' + myFileString + '.png')
        #cv.html_report('HTMLReports/TaskDHTMLReport' + myFileString + '.html')

        return

    def start(self):
        global myFileString
        print("start")
        T = 273.15
        x_init = np.array([100,1,1,1]) #initial concentration
        t_max = 300 #integration end time in seconds
        dt = 0.01 #step size in seconds
        cs = ChemSolver(chemkin('TaskD.xml'))
        cs.solve(x_init,T,t_max,dt, algorithm='vode', nsteps=5000)

        #time_array, conc_array, rxnrate_array = cs.get_results()

        cs.save_results('SimulationData/TaskDSimulationData' + myFileString + '.csv')

        cv = ChemViz(cs)
        cv.plot_time_series('concentration', tmin=0, tmax=300, species=['a','b','c','d'], outputfile='ConcentrationGraphs/TaskDConcentrationSeries' + myFileString + '.png')
        #cv.html_report('HTMLReports/TaskDHTMLReport' + myFileString + '.html')

        self.iterate()
        return

    def iterate(self):
        global myFileString
        print("iterate")
        step_size = float(0.5)
        tree = ET.parse('TaskD.xml')
        root = tree.getroot()
        constants = [i for i in root.iter(tag='k')] #Obtain elements for constants.  K values are text value of the element. 
        
        myFileString = '----k2='+constants[1].text+',k3='+constants[2].text
        
        if float(constants[2].text) < 5: # If k3 is less than 5 run, if not end.
            print("k3 less than 5")
            if float(constants[1].text) < 5: # if k2 is less than 5 iterate it, if its greater than 5 iterate k3 and set k2 to start number  (0.5)
                print("k2 less than 5")
                constants[1].text = str(float(constants[1].text) + step_size)
                constants[1].set('updated', 'yes')
                tree.write('TaskD.xml')
                return
            else:
                print("k2 greater than 5")
                constants[2].text = str(float(constants[2].text) + step_size)
                constants[2].set('updated', 'yes')
                constants[1].text = str(5.e-1)
                constants[1].set('updated', 'yes')
                tree.write('TaskD.xml')
                return
        
        else:
            print("k3 greater than 5")
            self.finished = True
            return


if __name__ == "__main__":
    file_name = 'TaskD.xml'
    myFileString = '----THELASTTEST'
    s = Simulation(file_name)
    s.start()
    # while not s.finished:
    #     s = Simulation(file_name, )
    #     s.start()
    input('SIMULATION FINISHED...PRESS ENTER TO EXIT')

