import ScrumPy
from ScrumPy.Bioinf import Biomass
import SepiBiomass

def GetFluxDic(m):
    '''Pre: modellp object constrained to biomass fluxes defined by BM composition
    Post: dictionary of fluxes for biomass transporters according to biomass composition
    in SepiBiomass.py and fixing here amount of biofilm plus ATPASE-RXN for cell maintenance
    '''
    fd = SepiBiomass.FluxDic(m)

    return fd

##def GetFluxDic(m, PropBF=0.0):
##    '''Pre: modellp object constrained to biomass fluxes defined by BM composition
##    Post: dictionary of fluxes for biomass transporters according to biomass composition
##    in SepiBiomass.py and fixing here amount of biofilm plus ATPASE-RXN for cell maintenance
##    '''
##    ##bm = SepiBiomass.All.Copy() # Copy() doesnt work!
##    bm =  Biomass.Composition(
##        1.0,{
##            "WholeCell": SepiBiomass.WholeCell,
##            "Biofilm" : SepiBiomass.Biofilm
##        }
##    )
##    
##    bm.SetAmmount(PropBF, "Biofilm")
##    bm.SetAmmount(1.0-PropBF, "WholeCell")
##    
##    fd = SepiBiomass.FluxDic(m, bm)
##    #fd["ATPASE-RXN"] = AMaint
##
##    return fd

def BuildLP(m):
    """Pre: model
	Post: lp object, with minimisation of total flux as default objective function"""
    lp = m.GetLP()
    lp.SetObjective(m.sm.cnames)
    return lp

def ClosedLP(m): # Maite
    """Pre:
        Post: lp object constraining flux through tx reacs to zero"""

    lp = BuildLP(m)
    
    for tx in list(filter(lambda s: "_tx" in s, m.sm.cnames)): # all tx in model
        lp.SetFixedFlux({tx:0.0}) # set flux to 0.0

    return lp

def BiomassLP(m, AMaint, fd=None):
    
    '''Pre: model
    Post: lp object constrained to biomass fluxes defined by BM composition
    and cell maintenance cost'''

    if fd==None:
        fd = GetFluxDic(m)

    lp = BuildLP(m)
    lp.SetFixedFlux(fd)
    lp.SetFixedFlux({'ATPASE-RXN':AMaint})
    return lp

def BiomassLPfreeBM(m, AMaint, fd=None):

    '''Pre:
    Post: lp object constrained to biomass fluxes defined by BM composition
    and cell maintenance cost with unconstrained lower flux for bm tx'''

    orig_rp = dict(m.sm.RevProps)
    for r in filter(lambda s: "bm_tx" in s, m.sm.cnames):
        m.sm.RevProps[r] = m.smx.RevProps[r] = "<>"

    if fd==None:
        fd = GetFluxDic(m)

    lp = BuildLP(m)

    for r in list(fd.keys()):
        #if 'bm_tx' in r:    # otherwise allows ATPase to go to zero flux! (since also in fd)
        lp.SetFluxBounds({r:(None,fd[r])}) # Bug = takes val of '0' if r not declared REV
    m.sm.RevProps = m.smx.RevProps = orig_rp

    lp.SetFixedFlux({'ATPASE-RXN':AMaint})
    return lp

def BlockUptakeLP(m):
    """ Pre: model
	Post: lp object, all media transporters constrained to 0 fluxes """

    lp = BuildLP(m)
    
    for tx in list(filter(lambda x: x.endswith("_mm_tx"), m.sm.cnames)):
        lp.SetFixedFlux({tx:0})

    return lp

def BlockAAUptakeLP(m):
    """ Pre: model
	Post: lp object, all amino acid media transporters constrained to 0 fluxes """
    lp = BuildLP(m)

    for tx in list(filter(lambda x: x.endswith("_AA_mm_tx"), m.sm.cnames)):
        lp.SetFixedFlux({tx:0})

    return lp


def BuildLPobjective(m, Obj=[]):
    """Pre: model
        Post: lp object with minimisation of flux through reacs in Obj as objective value"""

    if Obj == []:
        lp = BuildLP(m)
    else:
        lp = m.GetLP()
        for r in Obj:
            lp.SetObjective(r)
        
    return lp

def BlockUptakeLPobjective(m, Obj=[]):
    """ Pre: model
	Post: lp object, all media transporters constrained to 0 fluxes """

    if Obj == []:
        lp = BlockAAUptakeLP(m)
    else:
        lp = BuildLPobjective(m, Obj)
    for tx in list(filter(lambda x: x.endswith("_mm_tx"), m.sm.cnames)):
        lp.SetFixedFlux({tx:0})

    return lp

def Block_AAUptake_LPobjective(m, Obj=[]):
    """ Pre: model
	Post: lp object, all amino acid media transporters constrained to 0 fluxes
            plus extra flux constraints"""

    if Obj == []:
        lp =  BlockAAUptakeLP(m)
    else:
        lp = BuildLPobjective(m,Obj)

    for tx in list(filter(lambda x: x.endswith("_AA_mm_tx"), m.sm.cnames)):
        lp.SetFixedFlux({tx:0})

    return lp


#####   DO NOT TOUCH BELOW CODE ORIGINAL ##########################################################
###################################################################################################


##from ScrumPy.Util import Set
##from importlib import reload
##import FluxDic
##reload(FluxDic)
##
##
##def BuildLP(m):
##    """Pre: model
##	Post: lp object, with minimisation of total flux as default objective function"""
##    lp =m.GetLP()
##    lp.SetObjective(m.sm.cnames)
##    return lp

##def BiomassLP(m,fd=None):
##    """ Pre: model
##	Post: lp object, constrained to biomass fluxes defined in FluxDic"""
##    if fd == None:
##    	fd = FluxDic.fd
##
##    lp = BuildLP(m)        
##    lp.SetFixedFlux(fd)
##
##    unwanted = Set.Complement(filter(lambda x: x.endswith("_bm_tx"), m.sm.cnames), fd.keys())
##    print ('biomass not in fd', unwanted)
##    for tx in unwanted: 
##        lp.SetFixedFlux({tx:0.0})
##    
##    return lp
##
##
##def NullBiomassLP(m):
##    """ Pre: model
##	Post: lp object, all biomass transporters constrained to 0 fluxes"""
##    lp = BuildLP(m)
##
##    for tx in list(filter(lambda x: x.endswith("_bm_tx"), m.sm.cnames)):
##        lp.SetFixedFlux({tx:0})
##
##    return lp
##
##def BlockUptakeLP(m):
##    """ Pre: model
##	Post: lp object, all media transporters constrained to 0 fluxes """
##    lp = BuildLP(m)
##
##    for tx in list(filter(lambda x: x.endswith("_mm_tx"), m.sm.cnames)):
##        lp.SetFixedFlux({tx:0})
##
##    return lp
##
##def BlockAAUptakeLP(m):
##    """ Pre: model
##	Post: lp object, all amino acid media transporters constrained to 0 fluxes """
##    lp = BuildLP(m)
##
##    for tx in list(filter(lambda x: x.endswith("_AA_mm_tx"), m.sm.cnames)):
##        lp.SetFixedFlux({tx:0})
##
##    return lp
