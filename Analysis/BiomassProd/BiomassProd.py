######
########
##########   This module is for LPA of production of single BM components in isolation
##########    and identification of reac candidates for reversibility corrections
########    NOTE: remember, there is no NITRATE (NO3) in MMc media !!! and there is no 
######            GLN (glutamine), so block their importers accordingly during analysis

#import sys
#sys.path.append('../Analysis/BiomassProd/')
import BuildLP

## Production of single biomass components

def Check(m, n="ALA_AA_bm_tx"): # lp object checks for production of single biomass component (n)

    """pre: m = the model ; n = biomass transporter of compound of interest
       post: if biomass compound can be produced, prints flux through tx in
       solution and number of reacs in solution"""

    sol = {}
    lp = BuildLP.BuildLP(m)
    lp.SetFixedFlux({n:-1}) # Set flux of n to -1
    lp.Solve()
    sol = lp.GetPrimSol()
    for r in sol:
        if '_tx' in r:   # for transporters in solution
            print (r, sol[r])  # print reac name and flux

    print ("This is the number of reacs in the solution: " , "\n", len(sol))
    return sol

def CheckProdBlocks(m, n="ALA_AA_bm_tx", blocks=["GLN_AA_mm_tx", "NO3_mm_tx"]): # lp object checks for production of single biomass component (n)

    """pre: m = the model ;
            n = biomass transporter of compound of interest
            blocks = reactions to be blocked
       post: if biomass compound can be produced, prints flux through tx in
       solution and number of reacs in solution"""

    sol = {}
    lp = BuildLP.BuildLP(m)
    lp.SetFixedFlux({n:-1}) # Set flux of n to -1

    for r in blocks:
        if r in m.smx.cnames:
            lp.SetFixedFlux({r:0.0})

    lp.Solve()
    sol = lp.GetPrimSol()
    for r in sol:
        if '_tx' in r:   # for transporters in solution
            print (r, sol[r])  # print reac name and flux

    print ("This is the number of reacs in the solution: " , "\n", len(sol))
    return sol


def BMtx(m):

    """pre: m = the model
       post: returns a list of biomass transporters"""
    
    rv = []

    for r in filter(lambda s: "_bm_tx" in s, m.smx.cnames):
        rv.append(r)

    return rv

def CheckProds(m, AA=None, blocks=["GLN_AA_mm_tx", "NO3_mm_tx"]):

    # "GLN_AA_mm_tx", "NO3_mm_tx" are blocked cause are not present in MMc media

    """pre: m = the model
            AA = flux constraint through AA tx
            blocks = reactions to be blocked
       post: returns a dictionary with transporters of biomass products as keys
       and solution leading to production and export of 1 unit of each product as values"""

    rv = {}

    #lp = m.GetLP() # LP analysis of m
    #lp.SetObjective(m.sm.cnames) # minimize flux in all functions

    lp = BuildLP.BuildLP(m)

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    for r in blocks:
        if r in m.smx.cnames:
            lp.SetFixedFlux({r:0.0})

    for prod in BMtx(m):
        lp.SetFixedFlux({prod:-1})
        lp.Solve(False)                 # Dont print solution
        if lp.IsStatusOptimal():
            sol = lp.GetPrimSol()
            rv[prod] = sol
        else:
            rv[prod] = 0
        lp.ClearFluxConstraint(prod)    # Clear constraint so m back to normal

    Produced = []
    Targets = []

    for prod in list(rv.keys()):
        if rv[prod] != 0:
            Produced.append(prod)
        else:
            Targets.append(prod)

    print ("These are the exporters of BM compounds that can be produced and excreted:", "\n", Produced, "\n")
    print ("These are the exporters of BM compounds that cannot be produced or excreted:", "\n", Targets)

    return rv

def CheckProdRev(m, targ, AA=None, blocks=["GLN_AA_mm_tx", "NO3_mm_tx"]):

    """pre: m = the model ;
            targ = tx for product that currently cannot be made and excreted
            AA = flux constraint through AA tx
            blocks = reactions to be blocked
       post: if production is possible after making all reacs reversible, returns lp solution """

    orig_rp = dict(m.sm.RevProps)                                       # stores copy of dic of rev values of reacs in m
    for reac in filter(lambda s: not s.endswith("_tx"), m.sm.cnames):
            m.sm.RevProps[reac] = m.smx.RevProps[reac] = "<>"           # make rev in sm and smx to avoid internal inconsistency

    lp = BuildLP.BuildLP(m)

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    for r in blocks:
        if r in m.smx.cnames:
            lp.SetFixedFlux({r:0.0})

    lp.SetFixedFlux({targ:-1})  # check out if 'targeted' prod can be made
    lp.Solve()

    m.sm.RevProps = m.smx.RevProps = orig_rp    # revert changes in directionallity of reacs

    if lp.IsStatusOptimal():                    # if solution is optimal (target can be made)
        sol = lp.GetPrimSol()                   # print alternative solution
        return sol
    else:
        print ('No solution')

##    return sol

def RevCandsRL(m,altsol): # RIGHT TO LEFT

    """pre: altsol = CheckProdRev(m, targ...) ; alternative solution
        post: returns a list of canditate reacs to re-evaluate directionality R-L"""

    rv = []
    for r in list(altsol.keys()):
        if altsol[r] <0 and m.sm.RevProps[r] == "->": # If flux of reac is <- and in m was originally '->'
            rv.append(r)

    print ("These reactions have switched flux to negative in alternative solution:" , "\n", rv)
    return rv

def RevCandsLR(m,altsol): # LEFT TO RIGHT

    """pre: altsol = CheckProdRev(m, targ...)  ; alternative solution
        post: returns a list of canditate reacs to re-evaluate directionality L-R"""

    rv = []
    for r in list(altsol.keys()):
        if altsol[r] >0 and m.sm.RevProps[r] == "<-": # If flux of reac is -> and in m was originally '<-'
            rv.append(r)

    print ("These reactions have switched flux to positive in alternative solution:" , "\n", rv)
    return rv

def CheckCandReacsREV(m, targ, AA=None, blocks=["GLN_AA_mm_tx", "NO3_mm_tx"]):

    """pre: m = the model
            targ = tx for product that currently cannot be made and excreted
            AA = flux constraint through AA tx
            blocks = reactions to be blocked
        post: returns a list of canditate reacs to re-evaluate directionality"""
    
    altsol = CheckProdRev(m, targ, AA, blocks)
    
    CandReacsRevRL = RevCandsRL(m,altsol)
    
    CandReacsRevLR = RevCandsLR(m,altsol)

##def RevCandsRL(m, altsol): # RIGHT TO LEFT
##
##    """pre: altsol = CheckProdRev(m, targ)
##        post: returns a list of canditate reacs to re-evaluate directionality"""
##
##    rv = []
##    for r in altsol:
##        if altsol[r] <0 and m.sm.RevProps[r] == "->": # If flux of reac is <- and in m was originally '->'
##            rv.append(r)
##
##    print "These reactions have switched flux to negative in alternative solution:" , "\n", rv
##    return rv
##
##def RevCandsLR(m, altsol): # LEFT TO RIGHT
##
##    """pre: altsol = CheckProdRev(m, targ)
##        post: returns a list of canditate reacs to re-evaluate directionality"""
##
##    rv = []
##    for r in altsol:
##        if altsol[r] >0 and m.sm.RevProps[r] == "<-": # If flux of reac is -> and in m was originally '<-'
##            rv.append(r)
##
##    print "These reactions have switched flux to positive in alternative solution:" , "\n", rv
##    return rv
