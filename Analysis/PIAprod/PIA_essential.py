######
########
##########   This module is for defining essential reactions for 
##########   PIA synthesis 
########
######            

##    NOTE: remember, there is no NITRATE (NO3) or GLN in MMc media !!!
##              but there is GLN and NO3 in synovial fluid

import sys
from ScrumPy.Util import Set
import BuildLP
from ScrumPy.Bioinf import Biomass
import SepiBiomass

#DefaultBiomass=SepiBiomass.All

###############################################################################################
## Analyse system to define reactions in the same subset as
## reactions of interest

    ## PIA1 producing reaction : 'PIA1_synth', key for PIA production

def TargetEnzSubsets(m, target=['PIA1_synth']):
    """pre: m = the model
            target = list of reactions of interest
       post: returns the enzyme subsets containing the targeted reactions"""
    
    enzsubs=m.EnzSubsets()
    targetEnzSubs={}
    for r in target:
        for k in list(enzsubs.keys()):
            if r in list(enzsubs[k].keys()):
                targetEnzSubs[r]=k

    for r in target:
        print ('These are the reactions in the same enz subset as: ', r, "\n",)
        print (targetEnzSubs[r], enzsubs[targetEnzSubs[r]], "\n",)
        
    return enzsubs, targetEnzSubs

###############################################################################################
## Analyse system for production of 1 unit of PIA1
    ## while generating single and double reaction flux blocks

###############################################################################################

def CheckProdPIA1singleR(m, AMaint=0.0, AA=None, o2=None, no3=None, blocks=[]):
    # "GLN_AA_mm_tx", "NO3_mm_tx" are not present in MMc media
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
       post: returns a dictionary with solution for production of one unit of PIA1
       upon blocking flux of a single reaction at a time and list of essential
       reactions for PIA1 synthesis"""

    ProdPIA = {}

    constraints=[]

    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)

    lp.SetFixedFlux({"PIA1_bf_bm_tx": -1.0})
    constraints.append("PIA1_bf_bm_tx")

    if AMaint != None:
        lp.SetFluxBounds({"ATPase": (AMaint, None)})
        constraints.append("ATPase")

    if AA != None:
        AAtx=filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames)
        for aa in AAtx:
            lp.SetFixedFlux({aa:AA})
            constraints.append(aa)

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})
        constraints.append("O2_mm_tx")

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})
        constraints.append("NO3_mm_tx")

    for r in blocks:
        lp.SetFixedFlux({r:0.0})
        constraints=constraints+blocks

    print ('constraints: ', len(constraints), constraints)

    for r in Set.Complement(m.smx.cnames, constraints):
        # not to overwrite constraints!
        lp.SetFixedFlux({r:0.0})                              
        lp.Solve(False)
        print (".",) 
        if lp.IsStatusOptimal():
            ProdPIA[r] = lp.GetPrimSol()
        else:
            ProdPIA[r] = "EssentialForPIA"
        lp.ClearFluxConstraint(r) # now, this doesnt remove constraints either

    EssentialForPIA = []
    for r in ProdPIA.keys():
        if ProdPIA[r] == "EssentialForPIA":
            EssentialForPIA.append(r)

    NonEssentialForPIA = []
    for r in list(ProdPIA.keys()):
        if ProdPIA[r] != "EssentialForPIA":
            NonEssentialForPIA.append(r)

    print ("These reactions are essential for production of PIA1: ", "\n", len(EssentialForPIA), sorted(EssentialForPIA), "\n")
    print ("This is the number of reactions not essential for PIA1: ", "\n", len(NonEssentialForPIA))
    
    return ProdPIA, EssentialForPIA, NonEssentialForPIA

#### I THINK THIS FX BELOW IS NOT RIGHT! (YET)

def CheckProdPIA1doubleR(m, AMaint=0.0, AA=None, o2=None, no3=None, blocks=[]):
    # "GLN_AA_mm_tx", "NO3_mm_tx" are not present in MMc media
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
       post: returns a dictionary with solution for production of one unit of PIA1
       upon blocking flux of pairs of reactions at a time and list of essential
       reactions pairs for PIA1 synthesis"""

    ReacsEssentialForPIA = CheckProdPIA1singleR(m,AMaint,AA,o2,no3,blocks)[1]

    PIAprod1 = {}

    constraints=[]
    
    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)

    lp.SetFixedFlux({"PIA1_bf_bm_tx": -1.0})
    
    constraints.append("PIA1_bf_bm_tx")
    
    if AMaint != None:
        lp.SetFluxBounds({"ATPase": (AMaint, None)})
        constraints.append("ATPase")

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})
            constraints.append(aa)

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})
        constraints.append("O2_mm_tx")

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})
        constraints.append("NO3_mm_tx")

    for r in blocks:
        lp.SetFixedFlux({r:0.0})
        constraints=constraints+blocks

    print ('constraints: ', len(constraints), constraints)

    ReacsToKnockOut=[]
    for r in Set.Complement(m.smx.cnames, ReacsEssentialForPIA):
        # no need to scan pair already including 1 essential reac!
        if r in Set.Complement(m.smx.cnames, constraints): # not to overwrite constraints!
            ReacsToKnockOut.append(r)

    for r in ReacsToKnockOut:
        PIAprod2={}
        target1=[]
        lp.SetFixedFlux({r:0.0})
        target1.append(r)
        for r in Set.Complement(ReacsToKnockOut, target1):
            target2=[]
            lp.SetFixedFlux({r:0.0})
            target2.append(r)
            lp.Solve(False)
            print (".",)
            if lp.IsStatusOptimal():
                sol = lp.GetPrimSol()
##              PIAprod2[r] = lp.GetObjVal() # store ObjVal 
            else:
                PIAprod2[r] = "StopsPIA"
            lp.ClearFluxConstraint(r)      
        PIAprod1[r]=PIAprod2
        lp.ClearFluxConstraint(r)

    StopsPIAproduction = {}
    for r1 in list(PIAprod1.keys()):
        stops=[]
        for r2 in list(PIAprod1[r1].keys()):
            if PIAprod1[r1][r2] == "StopsPIA":
                stops.append(r2)
        StopsPIAproduction[r1] = stops

##    NoEffect = {}
##    for r1 in PIAprod1.keys():
##        noneStops=[]
##        for r2 in PIAprod1[r1].keys():
##            if PIAprod1[r1][r2] != "StopsPIA":
##                noneStops.append(r2)
##        NoEffect[r1] = noneStops

    print ("These pairs of reactions are essential for production of PIA1: ", "\n",)
    for r in sorted(list(StopsPIAproduction.keys())):
        print (r, ":  ", StopsPIAproduction[r])

    return StopsPIAproduction

#############################################################################################
## Analyse system for production of 1 unit of PIA biomass
        # and allowing BM prods to be exported at a higher proportion
            # than needed to produce BM!
            
######      USE WITH CAUTION!!!!!!!!! STILL IN DEVELOPMENT !!!!!

## Generate LP for production of biomass while meeting
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
            lp.SetFluxBounds({r:(None,fd[r])}) # Bug = takes val of '0' if r not declared REV first
#	lp.SetFluxBounds({r:(fd[r],None)})  # Bug = takes positive and negative values! if r not declared REV
    m.sm.RevProps = m.smx.RevProps = orig_rp
    
    lp.SetFluxBounds({"ATPase": (AMaint, None)}) # If we dont want ATPase upper flux to be fixed
                    
    return lp

def CheckProdPIA_BM_singleR(m, PropBF=1.0, AMaint=0.0, AA=None,
        o2=None, no3=None, blocks=[], FreeBM=None):
    # "GLN_AA_mm_tx", "NO3_mm_tx" are not present in MMc media
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
            FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unbound
       post: returns a dictionary with solution for production of one unit of PIA
       (PIA1 plus PIA2) upon blocking flux of a single reaction at a time and dic
       of essential reactions pairs for PIA1 synthesis"""

    ProdPIA = {}

    if FreeBM != None:
        lp = GetLPfreeBM(m, PropBF, AMaint)
    else:
        lp = GetLP(m, PropBF, AMaint)

    fd = GetFluxDic(m, PropBF, AMaint)

    constraints=[]

    constraints=constraints+list(fd.keys()) ## All bm_tx are here...
    # constraints=constraints+list(fd.keys()) ## All bm_tx are here...

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})
            constraints.append(aa)

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})
        constraints.append("O2_mm_tx")

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})
        constraints.append("NO3_mm_tx")

    for r in blocks:
        lp.SetFixedFlux({r:0.0})
        constraints=constraints+blocks

    print ('constraints: ', len(constraints), constraints)

    for r in Set.Complement(m.smx.cnames, constraints):
        # careful not to overwrite constraints!
        lp.SetFixedFlux({r:0.0})                              
        lp.Solve(False)
        print ("."), 
        if lp.IsStatusOptimal():
            ProdPIA[r] = lp.GetPrimSol()
        else:
            ProdPIA[r] = "EssentialForPIA"
        lp.ClearFluxConstraint(r) # now, this doesnt remove constraints either

    EssentialForPIA = []
    for r in list(ProdPIA.keys()):
        if ProdPIA[r] == "EssentialForPIA":
            EssentialForPIA.append(r)

    NonEssentialForPIA = []
    for r in list(ProdPIA.keys()):
        if ProdPIA[r] != "EssentialForPIA":
            NonEssentialForPIA.append(r)

    print ("These reactions are essential for production of PIA: ", "\n", len(EssentialForPIA), sorted(EssentialForPIA), "\n")
    print ("This is the number of reactions not essential for PIA: ", "\n", len(NonEssentialForPIA))
    
    return ProdPIA, EssentialForPIA, NonEssentialForPIA

def CheckProdPlank_BM_singleR(m, PropBF=0.0, AMaint=0.0, AA=None,
        o2=None, no3=None, blocks=[], FreeBM=None):
    # "GLN_AA_mm_tx", "NO3_mm_tx" are not present in MMc media
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
            FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unbound
       post: returns a dictionary with solution for production of one unit of PIA
       (PIA1 plus PIA2) upon blocking flux of a single reaction at a time and dic
       of essential reactions pairs for PIA1 synthesis"""

    ProdPlank = {}

    if FreeBM != None:
        lp = GetLPfreeBM(m, PropBF, AMaint)
    else:
        lp = GetLP(m, PropBF, AMaint)

    fd = GetFluxDic(m, PropBF, AMaint)

    constraints=[]

    constraints=constraints+fd.keys() ## All bm_tx are here...

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})
            constraints.append(aa)

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})
        constraints.append("O2_mm_tx")

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})
        constraints.append("NO3_mm_tx")

    for r in blocks:
        lp.SetFixedFlux({r:0.0})
        constraints=constraints+blocks

    print ('constraints: ', len(constraints), constraints)

    for r in Set.Complement(m.smx.cnames, constraints):
        # careful not to overwrite constraints!
        lp.SetFixedFlux({r:0.0})                              
        lp.Solve(False)
        print (".",)
        if lp.IsStatusOptimal():
            ProdPlank[r] = lp.GetPrimSol()
        else:
            ProdPlank[r] = "EssentialForPlank"
        lp.ClearFluxConstraint(r) # now, this doesnt remove constraints either

    EssentialForPlank = []
    for r in list(ProdPlank.keys()):
        if ProdPlank[r] == "EssentialForPlank":
            EssentialForPlank.append(r)

    NonEssentialForPlank = []
    for r in ProdPlank.keys():
        if ProdPlank[r] != "EssentialForPlank":
            NonEssentialForPlank.append(r)

    print ("These reactions are essential for production of planktonic BM: ", "\n", len(EssentialForPlank), sorted(EssentialForPlank), "\n")
    print ("This is the number of reactions not essential for planktonic BM: ", "\n", len(NonEssentialForPlank))
    
    return ProdPlank, EssentialForPlank, NonEssentialForPlank


#### I THINK THIS FX BELOW IS NOT RIGHT! (YET)

def CheckProdPIA_BM_doubleR(m, PropBF=1.0, AMaint=0.0, AA=None,
        o2=None, no3=None, blocks=[], FreeBM=None):
    # "GLN_AA_mm_tx", "NO3_mm_tx" are not present in MMc media
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
            FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unbound
       post: returns a dictionary with solution for production of one unit of PIA
       (PIA1 plus PIA2) upon blocking flux of pairs of reactions at a time and dic
       of essential reactions pairs for PIA1 synthesis"""

    ReacsEssentialForPIA = CheckProdPIA_BM_singleR(m,PropBF,AMaint,AA,o2,no3,blocks,FreeBM)[1]

    PIAprod1 = {}

    if FreeBM != None:
        lp = GetLPfreeBM(m, PropBF, AMaint)
    else:
        lp = GetLP(m, PropBF, AMaint)

    fd = GetFluxDic(m, PropBF, AMaint)

    constraints=[]

    constraints=constraints+list(fd.keys()) ## All bm_tx are here...

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})
            constraints.append(aa)

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})
        constraints.append("O2_mm_tx")

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})
        constraints.append("NO3_mm_tx")

    for r in blocks:
        lp.SetFixedFlux({r:0.0})
        constraints=constraints+blocks

    print ('constraints: ', len(constraints), constraints)

    ReacsToKnockOut=[]
    for r in Set.Complement(m.smx.cnames, ReacsEssentialForPIA):
        # no need to scan pair already including 1 essential reac!
        if r in Set.Complement(m.smx.cnames, constraints): # not to overwrite constraints!
            ReacsToKnockOut.append(r)

    for r in ReacsToKnockOut:
        PIAprod2={}
        target1=[]
        lp.SetFixedFlux({r:0.0})
        target1.append(r)
        for r in Set.Complement(ReacsToKnockOut, target1):
            target2=[]
            lp.SetFixedFlux({r:0.0})
            target2.append(r)
            lp.Solve(False)
            print (".",)
            if lp.IsStatusOptimal():
                sol = lp.GetPrimSol()
##              PIAprod2[r] = lp.GetObjVal() # store ObjVal 
            else:
                PIAprod2[r] = "StopsPIA"
            lp.ClearFluxConstraint(r)      
        PIAprod1[r]=PIAprod2
        lp.ClearFluxConstraint(r)

    StopsPIAproduction = {}
    for r1 in list(PIAprod1.keys()):
        stops=[]
        for r2 in list(PIAprod1[r1].keys()):
            if PIAprod1[r1][r2] == "StopsPIA":
                stops.append(r2)
        StopsPIAproduction[r1] = stops

##    NoEffect = {}
##    for r1 in PIAprod1.keys():
##        noneStops=[]
##        for r2 in PIAprod1[r1].keys():
##            if PIAprod1[r1][r2] != "StopsPIA":
##                noneStops.append(r2)
##        NoEffect[r1] = noneStops

    print ("These pairs of reactions are essential for production of PIA1: ", "\n",)
    for r in sorted(list(StopsPIAproduction.keys())):
        print (r, ":  ", StopsPIAproduction[r])

    return StopsPIAproduction
