
from importlib import reload
import BuildLP
reload(BuildLP)

def CheckUnitFlux(m,reac,neg=False,lp=None):
    """ pre: reac in m.sm.cnames
       post: the lp solution for a unit flux in reac
             or {} if none exists"""

    if lp ==None:
        lp = BuildLP.BuildLP(m)

    j = -1 if neg else 1
    lp.SetFixedFlux({reac:j})
    lp.Solve(False)
    lp.ClearFluxConstraint(reac)
    return lp.GetPrimSol()

def EnergyCheck(m,lp=None,check='ATPASE-RXN'):
    """Pre: model
	Post: LP solution if reaction can carry flux in absence of any media uptake"""
    
    if lp ==None:
        lp = BuildLP.BlockUptakeLP(m)
    sol = CheckUnitFlux(m,check,lp=lp)
    return sol

def NoFeasible(sols):
	rv = []
	for r in list(sols.keys()):
		if sols[r] == {}:
			rv.append(r)
	return rv

def CheckBMProduction(m,neg=True,lp=None):
    """ Pre: model
            Post: dictionary of biomass comp that can be produced with LP solution as values"""

    rv = {}
    not_produced=[]
    produced=[]

    if lp == None:
        lp = BuildLP.BuildLP(m)
    for bm in list(filter(lambda s: s.endswith('_bm_tx'), m.smx.cnames)):
        sol=CheckUnitFlux(m,bm,neg,lp=lp)
        rv[bm] = sol

    for k in rv.keys():
        if rv[k] == {}:
            not_produced.append(k)
    for k in rv.keys():
        if rv[k] != {}:
            produced.append(k)

    print('Not produced:', not_produced)
    print('Produced:', produced)
    
    return rv

def CheckBMProduction_NoAAs(m,neg=True,lp=None):
	"""Pre: model
	Post: dictionary of biomass that can be produced in absence of amino acids with values as LP solution"""
	rv = {}
	if lp ==None:
        	lp = BuildLP.BlockAAUptakeLP(m)
	rv = CheckBMProduction(m,neg,lp)
	return rv
    
def CheckEproduction_NoAAs(m,lp=None,check='ATPASE-RXN'): # Maite
    """ Pre: model
            Post: LP solution for E production in absence of media AAs"""
    
    if lp==None:
        lp = BuildLP.BlockAAUptakeLP(m) # Generate LP block AA influx
    sol=CheckUnitFlux(m,check,neg=False,lp=lp)
    return sol

def CheckEprod_NoAAs_blocks(m, Obj=[], check='ATPASE-RXN', blocks=[]): # Maite
    """Pre: model
        Post: LP solution for E production in absence of AA and extra constraints"""

    if Obj==[]:
        lp = BuildLP.BlockAAUptakeLP(m) 
    else:
        lp = BuildLP.Block_AAUptake_LPobjective(m, Obj) # objective: minimize flux through reacs in Obj
    
    for r in blocks:
        lp.SetFixedFlux({r:0})
    sol=CheckUnitFlux(m, check, neg=False, lp=lp)
    return sol

####### Fxs for curation = identify reacs to be modified to allow BMP

def CheckProdRev(m, targ, neg=True, lp=None):
    """Pre: model
        Post: LP solution for production of BM components with all model reacs made reversible"""

    orig_rp = dict(m.sm.RevProps) # store a copy of dic of reversibility values of stoich matrix of m
    for r in filter(lambda s: not s.endswith("_tx"), m.sm.cnames): # for reac in m excluding transporters
        m.sm.RevProps[r] = m.smx.RevProps[r] = "<>" # make reversible (in sm and smx to avoid internal inconsistency)

    rv = CheckUnitFlux(m,reac=targ,neg=True,lp=None)
##        lp = m.GetLP()
##        lp.SetObjective(m.sm.cnames)
##        lp.SetFixedFlux({targ:-1})  # check out if our 'targeted' compound can be made
##        lp.Solve()

    m.sm.RevProps = m.smx.RevProps = orig_rp # revert changes in directionallity of reacs in stoich matrix of m

    return rv
##        if lp.IsStatusOptimal(): # if solution is optimal (target can be made)
##                return lp.GetPrimSol() # print solution
        

def RevCands(m, sol):
    """Pre: sol = CheckProdRev(m, targ)
        Post: list of reactions to be checked for reversibility"""
    rv = []
    for reac in sol:
        if sol[reac] <0 and m.sm.RevProps[reac] == "->": # If flux of reac is <0 and in model direc originally '->'
            rv.append(reac)

    return rv

def DictRevCands(m):
    """ Pre: model
            Post: dictionary with BM compound as keys and reactions candidates for reversibility
            for their production as values"""
    rv={}
    no_sol=[]
    for bmtx in filter(lambda s: s.endswith('bm_tx'), m.smx.cnames):
        rev_sol= CheckProdRev(m, bmtx, neg=True, lp=None)
        rev_cands= RevCands(m, rev_sol)
        rv[bmtx]=rev_cands
        
        if rv[bmtx]==[]:
            no_sol.append(bmtx)

    print('BM components with no sol', no_sol)

    return  rv

#############################

#####   DO NOT TOUCH BELOW CODE ORIGINAL ##########################################################
###################################################################################################

##from ScrumPy.Util import Set
##from importlib import reload
##import BuildLP
##reload(BuildLP)
##import re
##import FluxDic
##reload(FluxDic)
##
##def CheckUnitFlux(m,reac,neg=False,lp=None):
##    """ pre: reac in m.sm.cnames
##       post: the lp solution for a unit flux in reac
##             or {} if none exists"""
##
##    if lp ==None:
##        lp = BuildLP.BuildLP(m)
##
##    j = -1 if neg else 1
##    lp.SetFixedFlux({reac:j})
##    lp.Solve()
##    lp.ClearFluxConstraint(reac)
##    return lp.GetPrimSol()
## 
##
##def CheckMultiUnitFluxes(m, reacs, neg=False, lp=None):
##    """ pre: reac in m.sm.cnames
##       post: the lp solution for a unit flux in reactions
##             or {} if none exists"""
##    if lp ==None:
##        lp =  BuildLP.BuildLP(m)
##
##    j = -1 if neg else 1
##    for reac in reacs:
##        lp.SetFixedFlux({reac:j})
##
##    lp.Solve()
##    lp.ClearFluxConstraints()
##    return lp.GetPrimSol()
##    
##
##def EnergyCheck(m,lp=None,check='ATPASE-RXN'):
##    """Pre: model
##	Post: LP solution if reaction can carry flux in absence of any media uptake"""
##    if lp ==None:
##        lp = BuildLP.BlockUptakeLP(m)
##    sol = CheckUnitFlux(m,check,lp=lp)
##    return sol
##
##def NoFeasible(sols):
##	rv = []
##	for r in list(sols.keys()):
##		if sols[r] == {}:
##			rv.append(r)
##	return rv
##			
##def CheckBMProduction(m,neg=True,lp=None,fd=None):
##	"""Pre: model
##	Post: dictionary of biomass that can be produced with values as LP solution"""
##	rv = {}
##	if fd == None:
##		bm = list(filter(lambda s:'_bm_tx' in s, m.sm.cnames))
##	else:
##		bm = list(fd.keys())
##	if lp ==None:
##        	lp = BuildLP.BuildLP(m)
##	for reac in bm:
##		sol = CheckUnitFlux(m,reac,neg,lp)
##		rv[reac] = sol
##	print ('Not feasible', NoFeasible(rv))
##	return rv
##	
##def CheckBMProduction_NoAA(m,neg=True,lp=None,fd=None):
##	"""Pre: model
##	Post: dictionary of biomass that can be produced in absence of amino acids with values as LP solution"""
##	rv = {}
##	if lp ==None:
##        	lp = BuildLP.BlockAAUptakeLP(m)
##	rv = CheckBMProduction(m,neg,lp,fd)
##	return rv
##
##
