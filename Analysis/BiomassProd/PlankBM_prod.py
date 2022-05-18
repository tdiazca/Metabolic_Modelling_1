##
###
#### Finds LP solutions for synthesis of planktonic biomass under different conditions
###    while meeting the cells energy demand (NGAM + GAM ATP demand)  
##        Note: no NITRATE (no3) or GLN in MHHW media
#         but both present in synovial fluid

import ScrumPy
from ScrumPy.Util import Set
import BuildLP
from ScrumPy.Bioinf import Biomass
import SepiBiomass

#DefaultBiomass=SepiBiomass.All

# Import DataForAnalysis for dicts of no. C and N atoms in media comps

def Tx(m):
    '''Pre: model
    Post: returns a list of transporters'''
    return list(filter(lambda s: '_tx' in s, m.smx.cnames))

def BMTx(m):
    '''Pre: model
    Post: returns a list of biomass transporters'''
    return list(filter(lambda s: '_bm_tx' in s, m.smx.cnames))

## Analyse system for production of 1 unit of cell biomass (BMP)
    # while meeting ATP cell maintenance requirements

def PlanktonicBMprod(m, AMaint=0.0, o2=None, no3=None,
                     fd=None, blocks=[]):
    '''Pre: model
    Post: LP solution for production of planktonic biomass while meeting ATP cell maintenace demand'''

    lp = BuildLP.BiomassLP(m,AMaint,fd)

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
    print ("\n","Compounds imported/exported","\n")
    for tx in Set.Complement(Txs, BM):
        if tx in sol.keys():
            print(tx, sol[tx])

    return sol

def PlanktonicBMprodFreeBM(m, AMaint=0.0, o2=None, no3=None,
                     fd=None, blocks=[]):
    '''Pre: model
    Post: LP solution for production of planktonic biomass
    while meeting ATP cell maintenace demand with free lower flux bound
    for bm transporters'''

    lp = BuildLP.BiomassLPfreeBM(m,AMaint,fd) 

    if o2 != None:
        lp.SetFixedFlux({"OXYGEN-MOLECULE_tx":o2})

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    lp.Solve(False)
    sol = lp.GetPrimSol()
    ObjVal = lp.GetObjVal()

    fd=BuildLP.GetFluxDic(m) # Fixed BM tx fluxes (1 gram BM)
    BMcompsExcreted = {}
     #for r in filter(lambda s: "bm_tx" in s, list(sol.keys())):
    for r in Set.Intersect(list(sol.keys()),list(fd.keys())):
        # (PIA1_bf_bm_tx, PIA2_bf_bm_tx, PIA3_bf_bm_tx with flux 0 if PropBF = 0.0)
        if sol[r] != fd[r]:
            comp = r.split("_")[0]
            NetExp = sol[r] - fd[r] 
            BMcompsExcreted[comp+"_NetExport"] = NetExp

    Txs = Tx(m)
    BM = BMTx(m)
    NoTxsInSol = Set.Complement(sol.keys(), Txs)

    print ("no. of reacs",len(sol.keys()))
    print ("no. of reacs excluding transporters", len(NoTxsInSol))
    print ("Objective value: ", ObjVal)
    print ("\n","Net import of media components: ", list(filter(lambda s: "_mm_tx" in s, sol.keys())))
    print ("\n","Compounds imported/exported","\n")
    for tx in Set.Complement(Txs, BM):
        if tx in sol.keys():
            print(tx, sol[tx])

    print ("These are the net amounts of BM compounds excreted (mmol): ", "\n",)
    for comp in sorted(list(BMcompsExcreted.keys())):
        print (comp, ":  ", round(BMcompsExcreted[comp],6))

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
