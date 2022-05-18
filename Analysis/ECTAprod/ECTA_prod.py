##
###
#### Finds LP solutions for ECTA synthesis under different conditions
###     Note: no NITRATE (no3) or GLN in MHHW media
##          but both present in synovial fluid
#

import ScrumPy
from ScrumPy.Util import Set
import BuildLP

# Import DataForAnalysis for dicts of no. C and N atoms in media comps

##### n of C atoms per comp in medium (MMc or synovial fluid)
##CatomsAll = {"ASP_AA_mm_tx":4,"SER_AA_mm_tx":3,"THR_AA_mm_tx":4,"ASN_AA_mm_tx":4,"ILE_AA_mm_tx":6, 
##           "PHE_AA_mm_tx":9,"LYS_AA_mm_tx":6,"TYR_AA_mm_tx":9,"CYS_AA_mm_tx":3,"GLY_AA_mm_tx":2,
##           "HIS_AA_mm_tx":6,"GLT_AA_mm_tx":5,"MET_AA_mm_tx":5,"TRP_AA_mm_tx":11,"GLN_AA_mm_tx":5,
##           "LEU_AA_mm_tx":6,"ARG_AA_mm_tx":6,"VAL_AA_mm_tx":5,"ALA_AA_mm_tx":3,"PRO_AA_mm_tx":5,
##           "GLC_mm_tx":6,"CITRULLINE_mm_tx":6,"ORNIT_mm_tx":5,"4-AMINO-BUTYRATE_mm_tx":4,
##           "UREA_mm_tx":1}
##
##### n of N atoms per comp in medium (MMc or synovial fluid)                                        
##NatomsAll = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
##           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
##           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
##           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
##           "NO3_mm_tx":1,"NH4_mm_tx":1,"CITRULLINE_mm_tx":3,"ORNIT_mm_tx":2,"4-AMINO-BUTYRATE_mm_tx":1,
##           "UREA_mm_tx":2}

def Tx(m):
    '''Pre: model
    Post: returns a list of transporters'''
    return list(filter(lambda s: '_tx' in s, m.smx.cnames))

def BMTx(m):
    '''Pre: model
    Post: returns a list of biomass transporters'''
    return list(filter(lambda s: '_bm_tx' in s, m.smx.cnames))

## Analyse system for production of ECTA (extracellular teichoic acids = WTA)
    ## in presence / absence O2 and NO3

def ProdECTA(m, ECTAdemand=-1.0, ATPdemand=None, o2=None, no3=None, blocks=[]):
    '''Pre: model
    Post: lp solution for synthesis of ECTA (extracellular teichoic acid)'''

    lp = BuildLP.BuildLP(m) # lp = m.GetLP(); # lp.SetObjective(m.sm.cnames)

    lp.SetFixedFlux({"ECTA_bf_bm_tx": ECTAdemand})

    if ATPdemand!= None:
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
