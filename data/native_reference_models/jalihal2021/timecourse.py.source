import os
import sys
import copy
import yaml
import pandas as pd
import numpy as np
import matplotlib
from tqdm import tqdm
import matplotlib.pyplot as plt
import nutrient_signaling.utils as utils
from nutrient_signaling.simulators import get_simulator

class TimeCourse:
    """
    Minimal script:
    
    __plottimecourse.py__
    modelpath = "../data/2018-9-26-12-3-no-sigma"
    from nutrient_signaling import TimeCourse
    tc = TimeCourse()
    tc.setSimulator(modelpath, simulatortype="cpp")
    tc.readData()
    tc.makeSimulator()
    tc.comparison()
    tc.plotall()

    Notes:
    Some important variables:
    1. self.data_path ::  The default path is  relative to the library. Set this explicitly when
       calling readData() if a different data set has to be compared. Please ensure that
       the time course YAML file is correctly formatted!
    2. self.customFuncDef :: This is a dictionary of lambda functions to appropriately standardize
       experimental data if it is not a straighforward min-max normalization.
    """
    def __init__(self):
        self.data = []
        self.predictions = []
        self.debugflag = False
        self.parameterSetsDF = None # DataFrame        
        self.data_path = '../data/yaml/time-course-data.yaml'
        self.customFuncDef = {'mig1': (lambda mig_ts: [np.log10((m)/(1e-3+1-m)) for m in mig_ts]), # ensures denom is not 0
                 'rib': (lambda rib_ts: [rib_ts[i]/(1e-2+rib_ts[0]) for i in range(0,len(rib_ts))])} #ensures denom is not 0
        self.numToUse = 0 # Number of parameter sets to use
        self.customylabel = {'mig1' : 'log(nMig1/cMig1)',
                    'rib':'Rel Rib'}
    
    def toggledebug(self):
        self.debugflag = not self.debugflag
        
    def setSimulator(self, modelpath, simulatortype='py'):
        self.modelpath = modelpath
        self.simulatortype = simulatortype

    def makeSimulator(self):
        simobj = get_simulator(self.modelpath, self.simulatortype)
        return simobj
        
    def debug(self, prnt):
        if self.debugflag:
            print(prnt)
            
    def readData(self, data_path=None):
        cwd = os.path.dirname(os.path.realpath(__file__))
        if data_path is None:
            data_path = cwd + '/../data/yaml/time-course-data.yaml'
        self.data_path = data_path
        with open(self.data_path,'r') as infile:
            self.data = yaml.safe_load(infile)

    def loadAlternateParameterSets(self,pathToAlternateParameters):
        self.parameterSetsDF = pd.read_csv(pathToAlternateParameters, sep='\t',index_col=None)

    def setNumberOfParameterSets(self, numToUse=50):
        if self.parameterSetsDF is None:
            print("Please use readData() to first load a parameter set file")
        else:
            self.numToUse = numToUse
            self.ParameterSetsToUse = self.parameterSetsDF.sort_values(by='cost')[0:5000]
            self.ParameterSetsToUse = self.ParameterSetsToUse.sample(self.numToUse)            
            print(self.ParameterSetsToUse.shape)
        
    def compareToExperiment(self, experiment, ax=None,plot=True):
        self.debug(experiment['description'])
        if ax is None:
            f, ax = plt.subplots(1,1)
        self.plotdata(ax, experiment)
        name = experiment['description']
        if self.numToUse > 0:
            for i, psetSeries in tqdm(self.ParameterSetsToUse.iterrows()):
                traj = self.simulate(experiment,
                                     parameterSet=psetSeries.to_dict())
                if experiment['tunits'] == 's':
                    traj['t'] = [t for t in traj['t']]
                if experiment['readout'].lower() in self.customFuncDef.keys():
                    traj['v'] = self.customFuncDef[experiment['readout'].lower()](traj['v'])
                self.plotpred(ax, traj, c='k',alpha=0.2,lw=1)
        traj = self.simulate(experiment)
        if experiment['tunits'] == 's':
            traj['t'] = [t for t in traj['t']]
        if experiment['readout'].lower() in self.customFuncDef.keys():
            traj['v'] = self.customFuncDef[experiment['readout'].lower()](traj['v'])        
        self.plotpred(ax, traj)        
        ax.set_ylabel(experiment['readout'])
        ax.set_xlabel('time')

            
    def comparison(self): 
        """
        Call simulate on each item in yaml file
        TODO: This is getting messy again. 
        3. cAMP should always have 0-1 range.
        """
        print("Starting Comparison")
        #store_attributes = ['value','time','readout']
        for simid, experiment in enumerate(self.data):
            self.debug(simid)
            traj = self.simulate(experiment)
            if experiment['tunits'] == 's':
                traj['t'] = [t*60. for t in traj['t']]
            if experiment['readout'].lower() in self.customFuncDef.keys():
                traj['v'] = self.customFuncDef[experiment['readout'].lower()](traj['v'])
            self.predictions.append(traj)

    def plotall(self):
        for simid, (experiment, pred) in enumerate(zip(self.data, self.predictions)):
            f, ax = plt.subplots(1,1)
            self.plotdata(ax, experiment)
            self.plotpred(ax, pred, experiment)
            name = experiment['description']
            plt.savefig(f'img/{name}.png')
            plt.close()
        
    def plotdata(self, ax, experiment):
        self.debug(experiment['title'])

        if experiment['tunits'] == 's':
            time = [t/60. for t in experiment['time']]
        else:
            time = experiment['time']            
        if experiment['readout'] !='Mig1':
            lo = experiment['normalize']['lower']
            hi = experiment['normalize']['upper']
            if np.isnan(lo) and len(experiment['value']) > 0 :
                lo = min(experiment['value'])
            if np.isnan(hi) and len(experiment['value']) > 0 :
                hi = max(experiment['value'])
            ax.set_ylim([0,1.1])
            plotdata = utils.minmaxnorm_special(experiment['value'], lo, hi)
        else:
            plotdata = experiment['value']
        ax.plot(time, plotdata, 'ko',markerfacecolor="none",ms=15)
        ax.set_title(experiment['title'])
        
    def plotpred(self, ax, traj, experiment, c='r',alpha=1.0,linestyle='-',lw=3):
        if experiment['tunits'] == 's':
            time = [t/60. for t in traj['t']]
        else:
            time = traj['t']        
        ax.plot(time, traj['v'],
                c=c,
                alpha=alpha,
                ls=linestyle,
                lw=lw)
            
    def simulate(self,experiment, parameterSet=None):
        """
        Simulate a time course using the preshift and post shift specifications defined in
        the experiment
        """
        self.debug(experiment['readout'])
        model = self.makeSimulator()
        if parameterSet is not None:
            model.set_attr(pars=parameterSet)
        preshift = experiment['spec']['preshift']
        postshift = experiment['spec']['postshift']
        self.debug(preshift)
        self.debug(postshift)
        tmax = 90
        if len(experiment['time']) >0:
            tmax = max(experiment['time'])
        
        if experiment['tunits'] == 's':
            tmax = tmax/60.

        # Set up preshift
        preics = model.get_ss()
        #preics.update(preshift['ics'])
        model.set_attr(ics=preics, pars=preshift['pars'])
        # set up post shift
        postics = model.get_ss()
        #postics.update(postshift['ics'])

        model.set_attr(ics=postics, pars=postshift['pars'], tdata=[0,tmax])
        y = model.simulate_and_get_points()
        traj = {'t':y['t'],
                'v':y[experiment['readout']]}

        return traj
            
    def plot_predicted(self, ax, traj, plotargs=None):
        """
        traj should be dictionary with keys 't' and 'v' with values being lists
        containing the time points and the values respectively
        """
        ax.plot(traj['t'], traj['v'])
