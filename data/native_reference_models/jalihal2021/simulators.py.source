import os
import sys
import subprocess
import pandas as pd
import PyDSTool as dst
import nutrient_signaling.modelreader as md
from nutrient_signaling.cpputils.pydstool2cpp import PyDSTool2CPP
from nutrient_signaling.utils import get_protein_abundances

class SimulatorPython:
    def __init__(self, modeldict, scale=False):
        self.modeldict = modeldict
        self.pars = {}
        self.ics = {}
        self.variables = {}
        self.createModelObject(modeldict)
        self.scale = scale
        
    def simulateModel(self):
        """
        Takes as input PyDSTool Model object
        and returns PyDSTool Pointset with
        default dt=0.01.
    
        :param model: PyDSTool Model object
        :return pts: Solution of Model
        :type pts: dict
        """
        pts = self.model.compute('test').sample()
        if self.scale:
            scaled = dict()
            scale_abundance_dict = get_protein_abundances()
            for k in self.variables: # loop over variables:
                if k in scale_abundance_dict.keys():
                    scaled[k] =[p*float(scale_abundance_dict[k])
                                for p in pts[k]]
                else:
                    scaled[k] = pts[k]
            return(scaled)
        else:
            return(pts)
    
    def simulate_and_get_points(self):
        return(self.simulateModel())
    def simulate_and_get_ss(self):
        """
        Returns steady states of all variables in the model
    
        :return SSPoints: Dictionary containing steady state values
        """
        Points = self.simulateModel()
        SSPoints={}
        for k in Points.keys():
            SSPoints[k]=Points[k][-1]
        return(SSPoints)
    
    def get_ss(self):
        """
        Returns steady states of all variables in the model
    
        :return SSPoints: Dictionary containing steady state values
        """
        Points = self.simulateModel()
        SSPoints={}
        for k in Points.keys():
            SSPoints[k]=Points[k][-1]
        return(SSPoints)

    def set_attr(self, pars={}, ics={}, tdata=[0, 90]):
        self.pars.update(pars)
        self.ics.update(ics)
        self.model.set(pars=self.pars, ics=self.ics, tdata=tdata)

    def get_attr(self):
        print("pars ", self.pars )
        print("ics ", self.ics )
        print("tend ", self.tend )
        
        
    def get_pars(self):
        return(self.pars)
    
    def get_variables(self):
        return(self.variables)
    
    def get_ics(self):
        return(self.ics)    
        
    def createModelObject(self, modeldict):
        """
        Takes dictionary object as input, and
        returns a PyDSTool Model object
    
        :param model: Model stored as dictionary. Should contain the keys `variables`, `parameters`, and `initiaconditions`
        :type model: dict
        :return ModelDS: a PyDSTool object that can be simulated to obtain the simulated trajectory
        :type ModelDS: PyDSTool Object
        """
        self.pars = modeldict['parameters']
        self.variables = modeldict['variables'].keys()
        self.ics = modeldict['initialconditions']                
        ModelArgs = dst.args(
            name='test',
            varspecs=modeldict['variables'],
            pars=self.pars,
            ics=self.ics,
            tdata=[0, 60],
            ## Define inline functions that are specific to the NutSig modelx
            ## tRNA() computes min(tRNA_total, Amino acid)
            ## pRib() computes min(Rib, eIF)
            ## shs() is the soft-heaviside function
            ## shsm() is the soft heaviside function with a maximum specified.
            fnspecs={'tRNA': (['tRNA_tot', 'AmAc'], 'min(tRNA_tot, AmAc)'),
                     'pRib': (['rib_comp', 'init_factor'],
                              'min(rib_comp,init_factor)'),
                     'shs': (['sig', 'summation'],
                             '1/(1+e^(-sig*summation))'),
                     'shsm': (['sig', 'm', 'summation'],
                              'm/(1+e^(-sig*summation))')}
        )
        self.model = dst.Vode_ODEsystem(ModelArgs)

class SimulatorCPP:
    def __init__(self, modeldict,
                 execpath='./src/',
                 executable='main.o',
                 simfilename='values.dat',
                 scale=False):
        self.execpath = execpath
        self.variables = modeldict['variables'].keys()
        self.pars = modeldict['parameters']
        self.ics = modeldict['initialconditions']
        self.tdata = [0, 90]
        self.tend = self.tdata[1]
        self.step = 0.001
        self.plot = False
        self.simfilename = simfilename
        self.executable = executable
        self.scale = scale
        
    def set_attr(self, pars={}, ics={},tdata=[0,90]):
        self.pars.update(pars)
        self.ics.update(ics)
        self.tdata = tdata
        self.tend = self.tdata[1]
        
    def get_attr(self):
        print("pars ", self.pars )
        print("ics ", self.ics)
        print("tend ", self.tend )
        
    def get_pars(self):
        return(self.pars)
    
    def get_variables(self):
        return(self.variables)
    
    def get_ics(self):
        return(self.ics)    
        
    def construct_call(self):
        """
        Returns:
        --------
        cmd : str
            Shell command to call cpp executable
        """
        cmd = self.execpath + self.executable +" --ode --tend " + str(self.tend) + " "\
              + "--step " + str(self.step) + " "
        
        if self.plot:
            cmd += "--plot "
    
        if self.simfilename is not None:
            cmd += "--out " + self.simfilename + " "
            
        if self.pars is not None:
            cmd += "--pars "
            for k,v in self.pars.items():
                cmd += k + " " + str(v)  + " "
                
        if self.ics is not None:
            cmd += "--ics "
            for k,v in self.ics.items():
                cmd += k + " " + str(v)  + " "
        return cmd
    
    def read_sim(self):
        simpoints = pd.read_csv(self.simfilename, sep="\t", index_col = False)
        return simpoints
    
    def get_ss_as_dict(self, simpoints):
        SS = simpoints.tail(1)
        SSdict = SS.to_dict(orient='record')[0]
        SSdict = {k:float(v) for k,v in SSdict.items()}
        return SSdict
    
        
    def simulate(self, cmd):
        #so = os.popen(cmd).read()
        so = subprocess.check_call(cmd.split(),
                                   stdout=subprocess.DEVNULL,
                                   stderr=subprocess.STDOUT)
        #so = subprocess.Popen(cmd.split())        
    
    #########
    ## Helpers
    def simulate_and_get_ss(self):
        """
        Returns steady states of all variables in the model
    
        :return SSPoints: Dictionary containing steady state values
        """
        Points = self.simulate_and_get_points()
        SSPoints={}
        for k in Points.keys():
            SSPoints[k]=Points[k][-1]
        return(SSPoints)
    
    # wrapper
    def get_ss(self):
        return(self.simulate_and_get_ss())
    
    def simulate_and_get_points(self):
        cmd = self.construct_call()
        self.simulate(cmd)
        D = self.read_sim()
        pts = D.to_dict(orient='list')
        if self.scale:
            scaled = {}
            scale_abundance_dict = get_protein_abundances()
            for k in pts.keys(): # loop over variables:
                if k in scale_abundance_dict.keys():
                    scaled[k] =[p*float(scale_abundance_dict[k])
                                for p in pts[k]]
                else:
                    scaled[k] = pts[k]
            return(scaled)
        else:
            return(pts)
    
    def simulateModel(self):
        return(self.simulate_and_get_points())

    
    #########
    ## Helpers
    def simulate_and_get_ss(self):
        """
        Returns steady states of all variables in the model
    
        :return SSPoints: Dictionary containing steady state values
        """
        Points = self.simulate_and_get_points()
        SSPoints={}
        for k in Points.keys():
            SSPoints[k]=Points[k][-1]
        return(SSPoints)
    
    # wrapper
    def get_ss(self):
        return(self.simulate_and_get_ss())
    
    def simulate_and_get_points(self):
        cmd = self.construct_call()
        self.simulate(cmd)
        D = self.read_sim()
        pts = D.to_dict(orient='list')
        if self.scale:
            scaled = {}
            scale_abundance_dict = get_protein_abundances()
            for k in pts.keys(): # loop over variables:
                if k in scale_abundance_dict.keys():
                    scaled[k] =[p*float(scale_abundance_dict[k])
                                for p in pts[k]]
                else:
                    scaled[k] = pts[k]
            return(scaled)
        else:
            return(pts)
    
    def simulateModel(self):
        return(self.simulate_and_get_points())        


def get_simulator(modelpath='./', simulator='py',**kwargs):
    """
    Initialize an simulator object.
    This involves reading the model defined at modelpath.
    If the simulator is 'cpp' and the corresponding executable does
    not exist, then the script compiles the C++ code and moves the 
    resulting executable to the appropriate location.
    """
    try:
        simulator in ["py","cpp","ppy"]
    except:
        print("Invalid simulator specification. Specific one of 'py' or 'cpp'.")
        sys.exit()
    model = md.readinput(modelpath) ## change
    scale = kwargs.get('scale',False)
    
    if simulator == 'py':
        simobj = SimulatorPython(model,scale=scale)
        return(simobj)
    
    elif simulator == 'cpp':
        validpath = False
        execpath = kwargs.get('execpath','./src/')
        executable = kwargs.get('executable','main.o')
        simfilename = kwargs.get('simfilename','values.dat')
        
        if not os.path.isdir(execpath):
            print('Path to executable does not exist')
            validpath = False
            
        if not os.path.isfile(execpath +  executable):
            print('Executable not found')
            validpath = False
            
        if not validpath:
            if not os.path.exists(execpath + executable):
                print(execpath + executable +' does not exist. Creating model file...')
                cwd = os.path.dirname(os.path.realpath(__file__))
                outpath = execpath
                if not os.path.exists(outpath):
                    os.mkdir(outpath)
                p2c = PyDSTool2CPP(modelpath)
                p2c.setwritepath(cwd + '/cpputils/')
                print(p2c.getwritepath())
                p2c.writecpp()
                headerpath = cwd + '/cpputils/'
                cmd = "g++ -std=c++11 {}main.cpp {}model.cpp -o {}/{}".format(headerpath,
                                                                              headerpath,
                                                                              execpath,
                                                                              executable)
                print(cmd)
                so = os.popen(cmd).read()
                #so = os.popen('cp ' + headerpath + executable+ ' ' + outpath)
        simobj = SimulatorCPP(model,
                              execpath=execpath,
                              executable=executable,
                              simfilename=simfilename,
                              scale=scale)
        return(simobj)
            
            
