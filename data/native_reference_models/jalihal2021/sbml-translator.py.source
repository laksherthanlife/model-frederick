"""
Author: Amogh Jalihal
Usage: At the commandline, pass the path to the directory where the model 
       is defined in three tab separated files.
       1. variables.txt contains the variable and the corresponding ODE RHS
       2. parameters.txt contains kinetic parameters and their values
       3. initialconditions.txt contains the variable names and initial values
"""
import sys
from libsbml import *
import re
import sys
from libsbml import *
import re

def checkstring(inputfile,line):
    if line[0]=='#':
        return 1
    elif line[0]=='\n':
        return 1
    else:
        K,V=line.split('\n')[0].split('\t')
        if '#' in V:
            V=V.split('#')[0]
        if inputfile!='variables':
            return((K,float(V)))
        else:
            return((K,V))

def readinput(PATH):

    modeldefinition={
    'variables':{},
    'parameters':{},
    'initialconditions':{}
    }
    
    for inputfile in modeldefinition.keys():
        print("Reading " + inputfile)
        with open(PATH+'/'+inputfile+".txt",'r') as infile:
            for line in infile.readlines():
                if checkstring(inputfile,line)==1:
                    continue
                else:
                    key,value=checkstring(inputfile,line)
                    modeldefinition[inputfile][key]=value
    return(modeldefinition)

def create_sbml_model(Model):
    try:
        document = SBMLDocument(3,2)
    except ValueError:
        raise SystemExit('Could not create SBML Documentation object')

    model = document.createModel()
    model.setTimeUnits("second")
    model.setExtentUnits("mole")
    model.setSubstanceUnits('mole')
    # There is a single compartment
    c = model.createCompartment()
    c.setId('c')
    c.setConstant(True)
    c.setSize(1)
    c.setSpatialDimensions(3)
    c.setUnits('litre')
    list_of_vars = list(Model['variables'].keys())
    ## Parameters
    for par in Model['parameters']:
        k = model.createParameter()
        k.setId(par)
        k.setConstant(True)
        k.setValue(Model['parameters'][par])
    ## Species aka variables
    for var in list_of_vars:
        s = model.createSpecies()
        s.setId(var)
        s.setConstant(False)
        s.setCompartment('c')
        s.setInitialAmount(Model['initialconditions'][var]) # ics
        s.setBoundaryCondition(False)
        s.setHasOnlySubstanceUnits(False)

    # Second pass. define reactions

    trna_func = re.compile("tRNA\(([\w\s\_\(\)\+\-]*),([\w\s\_\+\*\-\(\)]*)\)")
    prib_func = re.compile("pRib\(([\w\s\_\(\)]*),([\w\s\_\+\*\-\(\)]*)\)")    
    shs_func = re.compile("shs\(([\w\s\_]*),([\w\s\_\+\*\-\(\)\,\.]*)[\s]*-([\s\w\_\)]*)")

    for var in list_of_vars:
        r1 = model.createReaction()
        r1.setId(var + "_r")
        r1.setReversible(False)
        r1.setFast(False)

        equation = Model['variables'][var]
        # react = r1.createReactant()
        # react.setSpecies(var)
        # react.setConstant(True)

        for v_modifier in list_of_vars:
            if v_modifier in equation:
                modifier = r1.createModifier()
                modifier.setSpecies(v_modifier)
                #modifier.setConstant(True)
        if "tRNA(" in equation:
            reout = trna_func.search(equation)
            terms = reout.groups()
            equation = equation.replace("tRNA(" + terms[0] + "," + terms[1] +")",
                                        "min(" + terms[0] + "," + terms[1] +")")
        if "pRib(" in equation:
            reout = prib_func.search(equation)
            terms = reout.groups()
            equation = equation.replace("pRib(" + terms[0] + "," + terms[1] +")",
                                        "min(" + terms[0] + "," + terms[1] +")")

        if 'shs' in equation:
            reout = shs_func.search(equation)
            if reout is None:
                print(equation)
                break
            terms = reout.groups()
            terms = list(terms)
            terms[1]  = terms[1][:-1]
            terms = tuple(terms)
            exp = str(equation)
            exp = exp.replace("shs("+terms[0] + "," + terms[1],
                        "(1/(1+exp(-("+ terms[0] + ")*("+terms[1]+")))")
            product = r1.createProduct()
            product.setSpecies(var)
            product.setConstant(True)
            kinetic_law = r1.createKineticLaw()
            mathAST = parseL3Formula(exp)
            kinetic_law.setMath(mathAST)

        else:
            kinetic_law = r1.createKineticLaw()
            kinetic_law.setMath(parseL3Formula(equation))
            product = r1.createProduct()
            product.setSpecies(var)
            product.setConstant(True)
    return writeSBMLToString(document)

def main(inputpath):
    return(readinput(inputpath))

if __name__=='__main__':

    from optparse import OptionParser
    
    parser = OptionParser()
    parser.add_option('-p','--path',dest='path',default="./",
                  help='Path to model definition')
    (options,args)=parser.parse_args()
    
    Model = main(options.path)
    SBMLdoc = create_sbml_model(Model)
    with open("model.xml",'w') as outfile:
        outfile.write(SBMLdoc)
