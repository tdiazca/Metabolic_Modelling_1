##
###
#### Finds LP solutions for ATP synthesis in different conditions
###     Note: no NITRATE (no3) or GLN in MHHW media
##          but both present in synovial fluid
#

import ScrumPy
from ScrumPy.Util import Set
import BuildLP

# Import DataForAnalysis for dicts of no. C and N atoms in media comps

def MediaTx(m):
    '''Pre: model
    Post: freturn list of media transporters'''
    return list(filter(lambda s: '_mm_tx' in s, m.smx.cnames))

def AAMediaTx(m):
    '''Pre: model
    Post: return list of amino acids media transporters'''
    return list(filter(lambda s: '_AA_mm_tx' in s, m.smx.cnames))

def Tx(m):
    '''Pre: model
    Post: returns a list of transporters'''
    return list(filter(lambda s: '_tx' in s, m.smx.cnames))

def BMTx(m):
    '''Pre: model
    Post: returns a list of biomass transporters'''
    return list(filter(lambda s: '_bm_tx' in s, m.smx.cnames))

def BPsTx(m):
    '''Pre: m = model
    Post: returns a list of byproducts transporters'''
    return(filter(lambda s: '_bp_tx' in s, m.smx.cnames))

### Analyse system for production of ATP in presence/absence O2 and NO3

def ProdATP(m, ATPdemand=1.0, o2=None, no3=None, blocks=[]):
    '''Pre: m, ATP demand, O2, NO3, reacs to block
    Post: lp solution for ATP production'''

    lp =BuildLP.BuildLP(m) # lp=m.GetLP(); lp.SetObjective(m.sm.cnames)
    lp.SetFixedFlux({"ATPASE-RXN": ATPdemand})

    if o2 != None:
        lp.SetFixedFlux({"OXYGEN-MOLECULE_tx":o2})

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    lp.Solve(False)
    sol = lp.GetPrimSol()
    ObjVal = lp.GetObjVal()

    Txs = Tx(m)
    BM = BMTx(m)

    NoTxsInSol = Set.Complement(sol.keys(), Txs)

    print ("no. of reacs",len(sol.keys()))
    print ("no. of reacs excluding transporters", len(NoTxsInSol))
    print ("Objective value: ", ObjVal)
    print ("\n","Net import of media components: ", list(filter(lambda s: "_mm_tx" in s, sol.keys())))
##    print ("Byproducts excreted:")
##    for bp in BPs:
##        if bp in sol.keys():
##            print(bp, sol[bp])
    print("\n","Compounds imported/exported","\n")
    for tx in Set.Complement(Txs, BM):
        if tx in sol.keys():
            print(tx, sol[tx])

    return sol

def MakeSubModelFromSolution(m, sol, File='SubModel.spy'):
    '''Pre: sol
    Post: generates submodel containing reactions in lp solution sol'''

    reacs = list(sol.keys())
    sm = m.smx.SubSys(reacs)
    sm.ToScrumPy(File) # creates file describing model containing reacs

def ElementaryModesSubmodel(m, sol, File='SubModel.spy'):
    '''Pre: sol, submodel file
    Post: net stoichiometries for the elementary modes ion submodel'''

    submodel = MakeSubModelFromSolution(m,sol,File)
    mm = ScrumPy.Model(File)
    #k = m.sm.NullSpace()
    #k.IntegiseC()
    #print(k)

    elmo = mm.ElModes()
    elmo.Integise()

    print('These are the elementary modes in the subnetwork: ',len(elmo),"\n")
    print ('These are the net stoichiometries of the elementary modes:', "\n",elmo.Stos())

    return elmo          


    










