######
########
##########   This module is for LPA of a response for biomass production while
##########   meeting the cells energy demand (NGAM + GAM ATP demand)   
########    NOTE: remember, there is no NITRATE (NO3) or GLN in MMc media !!!
######              but there is GLN and NO3 in synovial fluid

import ScrumPy
import sys
from ScrumPy.Util import Set
import BuildLP
from ScrumPy.Bioinf import Biomass
import SepiBiomass
import DataForAnalysis # so it can recognise CatomsAll and NatomsAll

#DefaultBiomass=SepiBiomass.All

def Tx(m):
    """pre: m = the model
       post: returns a list of transporters"""
    rv=[]
    for r in filter(lambda s: "_tx" in s, m.smx.cnames):
        rv.append(r)
        
    return rv

def BMtx(m):
    """pre: m = the model
       post: returns a list of biomass transporters"""
    rv = []
    for r in filter(lambda s: "_bm_tx" in s, m.smx.cnames):
        rv.append(r)
        
    return rv

## Generate LP for production of whole biomass while meeting
    ## the ATP cell maintenace reqirements

def GetFluxDic(m, PropBF=0.0, AMaint=45.0): # set flux bm_tx for BMP and ATP demand
    """ lp object constrained to biomass fluxes defined by BM composition """
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
       post: generates a flux dictionary for biomass exporters according to description of
       biomass composition and sets flux through ATPase to the ATP maintenace cost"""
    ##bm = SepiBiomass.All.Copy() # Copy() doesnt work!
    bm =  Biomass.Composition(
        1.0,{
            "WholeCell": SepiBiomass.WholeCell,
            "Biofilm" : SepiBiomass.Biofilm
        }
    )
    
    bm.SetAmmount(PropBF, "Biofilm")
    bm.SetAmmount(1.0-PropBF, "WholeCell")
    
    fd = SepiBiomass.FluxDic(m, bm)
    fd["ATPase"] = AMaint

    return fd

def GetLP(m, PropBF=0.0, AMaint=45.0):
    """pre: GetFluxDic(m, PropBF=0.0, AMaint=45.0)     
       post: lp object constrained to biomass fluxes defined by BM composition
       and ATP maintenace demand"""

    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)
    fd = GetFluxDic(m, PropBF, AMaint)
    for r in list(fd.keys()):
        if 'bm_tx' in r:
            lp.SetFixedFlux(fd)
    lp.SetFluxBounds({"ATPase": (AMaint, None)}) # If we dont want ATPase upper flux to be fixed

    return lp

def GetLPfreeBM(m, PropBF=0.0, AMaint=45.0):
    """pre: GetFluxDic(m, PropBF=0.0, AMaint=45.0)     
       post: lp object constrained to biomass fluxes defined by BM composition
       and ATP maintenace demand"""
    
    orig_rp = dict(m.sm.RevProps)
    for r in filter(lambda s: "bm_tx" in s, m.sm.cnames):
        m.sm.RevProps[r] = m.smx.RevProps[r] = "<>"

    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)

    fd = GetFluxDic(m, PropBF, AMaint)
        # if PropBF=1, still all other bm_tx are keys in dic, but flux in fd = 0
    
    for r in list(fd.keys()):
        if 'bm_tx' in r:    # otherwise allows ATPase to go to zero flux! (since also in fd)
            lp.SetFluxBounds({r:(None,fd[r])}) # Bug = takes val of '0' if r not declared REV
#	lp.SetFluxBounds({r:(fd[r],None)})  # Bug = takes positive and negative values! if r not declared REV
    m.sm.RevProps = m.smx.RevProps = orig_rp
    
    lp.SetFluxBounds({"ATPase": (AMaint, None)}) # If we dont want ATPase upper flux to be fixed
                    
    return lp

## Analyse system for production of 1 unit of cell biomass (BMP)
    # while meeting ATP cell maintenance requirements

def ProdUnitBM(m, PropBF=0.0, AMaint=45.0, AA=None, o2=None, no3=None,
               targets=["GLN_AA_mm_tx"], FreeBM=None):
    # "GLN_AA_mm_tx", "NO3_mm_tx" are not present in MMc media
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            targets = targeted reactions to be blocked
            FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unbound
       post: returns a dictionary with solution for production of one unit of cell biomass"""

    sol = {}

    if FreeBM != None:
        lp = GetLPfreeBM(m, PropBF, AMaint)
    else:
        lp = GetLP(m, PropBF, AMaint)

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})

    for r in targets:
        lp.SetFixedFlux({r:0.0})

    lp.Solve(False)
    sol = lp.GetPrimSol()
    ObjVal = lp.GetObjVal()

    ByProds = {}
    for r in sol:
        if r in filter(lambda s: "_bp_tx" in s, m.smx.cnames):
            ByProds[r] = round(sol[r],3)

    CompUsed = {}
    for r in sol:
        if r in Set.Complement(Tx(m), BMtx(m)):
            CompUsed[r] = round(sol[r],3)

    TotalSolReac = []
    for r in sol:
        TotalSolReac.append(r)

    SolReac = []
    for r in sol:
        if r not in Tx(m):
            SolReac.append(r)

    if FreeBM != None:

        fd=GetFluxDic(m, PropBF, AMaint) # Fixed BM tx fluxes (1 gram BM)

        BMcompsExcreted = {}
        for r in filter(lambda s: "bm_tx" in s, list(sol.keys())):
            for r in Set.Intersect(list(sol.keys()),list(fd.keys())):
                # ATPase and (PIA1_bf_bm_tx, PIA2_bf_bm_tx, PIA3_bf_bm_tx with flux 0 if PropBF = 0.0)
                if sol[r] != fd[r]:
                    comp = r.split("_")[0]
                    NetExp = sol[r] - fd[r] 
                    BMcompsExcreted[comp+"_NetExport"] = NetExp

    print ("Total no of reactions in solution: ", len(TotalSolReac))
    print ("Total no of reactions in solution without transporters: ", len(SolReac))
    print ("Objective value of the solution in mmol/gDW/h: ", round(ObjVal,3))
    if "GLC_mm_tx" in list(sol.keys()): # OJO, 'GLC_PTS_mm_tx' could also be in sol!
        if "GLC_PTS_mm_tx" in list(sol.keys()): 
            print ("This is the ammount of GLC utilised (mmol): ", round(sol["GLC_mm_tx"] + sol["GLC_PTS_mm_tx"],3), "\n",)
        else:
            print ("This is the ammount of GLC utilised (mmol): ", round(sol["GLC_mm_tx"],3), "\n",)
    elif "GLC_PTS_mm_tx" in list(sol.keys()):
            print ("This is the ammount of GLC utilised (mmol): ", round(sol["GLC_PTS_mm_tx"],3), "\n",)
    else:
        print ("No GLC is consumed", "\n",)
##    if "GLC_mm_tx" in sol.keys():
##        print "This is the ammount of GLC utilised (mmol): ", round(sol["GLC_mm_tx"],3)
    if "NH4_mm_tx" in list(sol.keys()):
        print ("This is the ammount of NH4 utilised/excreted (mmol): ", round(sol["NH4_mm_tx"],3))
    print ("These are the ammounts of by-products exported (mmol): ", "\n",)
    for bp in sorted(list(ByProds.keys())):
        print (bp, ":  ", round(ByProds[bp],3))
    print ("These are the amounts of compounds uptaken/excreted (mmol): ", "\n",)
    for comp in sorted(list(CompUsed.keys())):
        print (comp, ":  ", round(CompUsed[comp],3))

    if FreeBM != None:
        print ("These are the net amounts of BM compounds excreted (mmol): ", "\n",)
        for comp in sorted(list(BMcompsExcreted.keys())):
            print (comp, ":  ", round(BMcompsExcreted[comp],6))

    if FreeBM != None:
        return sol, ObjVal, ByProds, CompUsed, TotalSolReac, SolReac, BMcompsExcreted

    else:
        return sol, ObjVal, ByProds, CompUsed, TotalSolReac, SolReac

####### NEW FX 9 sept 2020

def MakeSubModelFromSolution(m,res,File='SubModel.spy'):
    """pre: res = analysis solution,
        post: sub-model containing all reacs involved in solution"""
    reacs = list(res.keys()) #res.keys() python 2
    # e.g. if res = running an LP analysis such as ProdUnitBM[0] above
    sm = m.smx.SubSys(reacs)
    sm.ToScrumPy(File,Externs=m.Externals())
    #sm.ToScrumPy(m.md.QuoteMap,File,Externs=m.Externals()) #ScrumPy2
        # writes a File describing a model containing reacs

def ElementaryModesSubmodel(m,res,File='SubModel.spy'):
    """pre: m = model,
        res = analysis solution,
        MakeSubModelWholeDataSet(...) = generates a file containing all reactions in a solution (res)
     post: returns the net stoichiometries of the elementary modes of the submodel from reacs in res"""
    submodel = MakeSubModelFromSolution(m,res,File='SubModel.spy')
    mm=ScrumPy.Model(File)
    #k=m.sm.NullSpace()
    #k.IntegiseC()          # show entries in k matrix as integers
    #print k

    elmo = mm.ElModes()
    elmo.Integise()         # show elmo coeficients as integers

    print ('There are n elementary modes in the subnetwork: ',len(elmo),"\n")
    print ('These are the net stoichiometries of the elementary modes of the submodel: ',"\n","\n",elmo.Stos()) 
    return elmo

        #print elmo.ReacsOf('ElMo_0')

####### NEW FX 9 sept 2020 END - above

def ProdBMcharact(m, PropBF=0.0, AMaint=45.0, AA=None, o2=None, no3=None,
                  targets=["GLN_AA_mm_tx"], FreeBM=None, 
                  CatomsAll=DataForAnalysis.Catoms_synovial, NatomsAll=DataForAnalysis.Natoms_synovial):

    """pre: m = the model
            Objective = reaction which flux must be minimised as LP objective
            AMaint = ATP maintenance demand (NGAM plus GAM),growth rate 1gDW/h = 45 mmol/gDW = 45.0
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            targets = targeted reactions to be blocked
            CatomsAll = dic of n of C atoms in media components
            NatomsAll = dic of n of N atoms in media components
       post: returns a dictionary with solution for production of ATP"""

    res = ProdUnitBM(m, PropBF, AMaint, AA, o2, no3, targets, FreeBM)[0] 

    CompContribCuptake = ContribToTotalCuptake(res, CatomsAll)

    CompContribNuptake = ContribToTotalNuptake(res, NatomsAll)

    return res

def TotalCconsumed(res, CatomsAll=DataForAnalysis.Catoms_synovial):
    """pre: res = LP sol (e.g. ProdUnitBM(m))
            CatomsAll = dic of n of C atoms in media components
       post: returns a dic of total C atoms imported per medium components
            and total C atoms imported"""

    CperComp = {}     # C uptaken
    for comp in list(CatomsAll.keys()):
        if comp in list(res.keys()):
            if res[comp] > 0.0: 
                CperComp[comp] = CatomsAll[comp]*res[comp]
        else:
            CperComp[comp] = 0.0

    TotalC = sum(CperComp.values())

    print ("This is the total C consumed (mmol):", round(TotalC,3), "\n",)

    return CperComp, TotalC

def ContribToTotalCuptake(res, CatomsAll=DataForAnalysis.Catoms_synovial):
    """pre: res = LP sol (e.g. ProdUnitBM(m))
            CatomsAll = dic of n of C atoms in media components
       post: returns a dic of percentage of contribution of each medium
       component to the total C uptaken"""

    CperComp = TotalCconsumed(res, CatomsAll)[0]

    TotalC = TotalCconsumed(res, CatomsAll)[1]

    propTdicC = {}
    for comp in CatomsAll:
        if comp not in list(res.keys()): # if comp in medium is not uptaken   
            propTdicC[comp] = 0.0
        elif res[comp] < 0.0:
            propTdicC[comp] = 0.0
        else:
            propTdicC[comp] = CperComp[comp]/TotalC
            
    TotalpercentC = round((sum(propTdicC.values())*100),3)

    print ("\n", "These are the contributions of medium components to total C uptake (percentage): ", "\n")
    for comp in sorted(list(CatomsAll.keys())):
        print (comp, "","","","", round(propTdicC[comp]*100,3))

    print ("\n", "Sum: ", TotalpercentC, "\n")

    return propTdicC 

def TotalNconsumed(res, NatomsAll=DataForAnalysis.Natoms_synovial): # NatomsAll=NatomsAll
    """pre: res = LP sol (e.g. ProdUnitBM(m))
            NatomsAll = dic of n of N atoms in media components
       post: returns a dic of total N atoms imported per medium components
            and total N atoms imported"""

    NperComp = {}     # N uptaken
    for comp in list(NatomsAll.keys()):
        if comp in list(res.keys()):
            if res[comp] > 0.0: 
                NperComp[comp] = NatomsAll[comp]*res[comp]
        else:
            NperComp[comp] = 0.0

    TotalN = sum(NperComp.values())

    print ("These is the total N consumed (mmol):", round(TotalN,3), "\n",)

    return NperComp, TotalN

def ContribToTotalNuptake(res, NatomsAll=DataForAnalysis.Natoms_synovial): #NatomsAll=NatomsAll
    """pre: res = LP sol (e.g. ProdUnitBM(m))
            NatomsAll = dic of n of N atoms in media components
       post: returns a dic of percentage of contribution of each medium
       component to the total N uptaken"""

    NperComp = TotalNconsumed(res, NatomsAll)[0]

    TotalN = TotalNconsumed(res, NatomsAll)[1]

    propTdicN = {}

    for comp in NatomsAll:
        if comp not in list(res.keys()): # if comp in medium is not uptaken   
            propTdicN[comp] = 0.0
        else:
            if res[comp] > 0.0: 
                propTdicN[comp] = NperComp[comp]/TotalN
            
    TotalpercentN = round((sum(propTdicN.values())*100),3)

    print ("\n", "These are the contributions of medium components to total N uptake (percentage): ", "\n")
    for comp in sorted(list(NatomsAll.keys())):
        if comp in list(propTdicN.keys()):
            print (comp, "","","","", round(propTdicN[comp]*100,3))
    print ("\n", "Sum: ", TotalpercentN, "\n")

    return propTdicN

#################################################################

####
####def PrintSol(ProdUnitBM):
####
####    unit = ProdUnitBM.ProdUnitBM(m)
####
####    for r in unit[2]:
####        print r, ": ", unit[2][r] 

##    for r in unit[2]:
##	print r, ": ", unit[2][r]
##
##ACET_bp_tx :  -1.10255623366
##NH4_mm_bp_tx :  -7.84871627868
##CO2_bp_tx :  -13.8542883583
##SUC_bp_tx :  -7.3565056005
##
##    for r in unit[3]:
##	print r, ": ", unit[3][r]
##
##Pi_mm_tx :  1.30183469267
##PRO_AA_mm_tx :  0.116176691276
##TYR_AA_mm_tx :  0.118867158965
##LEU_AA_mm_tx :  0.281918600027
##HS_tx :  -0.00744487634041
##CYS_AA_mm_tx :  0.0192293230272
##Autoinducer_2_tx :  -0.00744487634041
##ARG_AA_mm_tx :  0.261341408021
##ACET_bp_tx :  -1.10255623366
##HIS_AA_mm_tx :  0.0733605812859
##ASP_AA_mm_tx :  0.627151220528
##VAL_AA_mm_tx :  0.207370510749
##GLT_AA_mm_tx :  7.716035946
##ASN_AA_mm_tx :  0.172716375919
##PHE_AA_mm_tx :  0.137417216616
##NH4_mm_bp_tx :  -7.84871627868
##GLC_mm_tx :  5.77727066344
##SER_AA_mm_tx :  0.197630611268
##CO2_bp_tx :  -13.8542883583
##THR_AA_mm_tx :  1.03555913191
##MET_AA_mm_tx :  0.0911185817296
##ILE_AA_mm_tx :  0.269180268803
##TRP_AA_mm_tx :  0.0229736703097
##SUC_bp_tx :  -7.3565056005
##O2_mm_tx :  12.6980778213
##ALA_AA_mm_tx :  6.67696826638
##NIACINE_mm_tx :  0.00201261420041
##LYS_AA_mm_tx :  0.338245197349
##H2O_mm_tx :  -14.2398582548
