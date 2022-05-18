##Python source code for ScrumPy3

from importlib import reload
import BuildLP
reload(BuildLP)
import ModelResults
reload(ModelResults)


def GetSol(m,lp,mediabound=True,substrate=None):
	'''Pre: m, lp and substrate to be blocked
	Post: return lp solution, obj value, no of reac with blocked substrate transporters '''
	
	if substrate != None:
		transporters = ModelResults.SubTx(m,substrate)
		for tx in transporters:
			lp.SetFixedFlux({tx:0.0})
		lp.Solve(False)
	
	else:
		lp.Solve(False)
	if lp.IsStatusOptimal():
		sol = lp.GetPrimSol()
		return (sol,lp.GetObjVal(),len(sol.keys()))
	else:
		return None

def SolMM(m,ion='AMMONIUM',Biomass=True):
	'''Pre: True
	Post:  LP Sol for biomass production on inorganic N source'''
	rv = {}
	lp = BuildLP.BlockAAUptakeLP(m,Biomass)
	lp.Solve()
	if lp.IsStatusOptimal():
		sol = lp.GetPrimSol()
		rv[ion] = (sol,lp.GetObjVal(),len(sol.keys()))
	print ("Inputs on minimal media",list(filter(lambda s:'_mm_tx' in s, sol.keys())))
	print ("Obj val and no. of reacs",rv[ion][1],rv[ion][2])
	return rv

def GetSingleNLP(m,ion='AMMONIUM',Biomass=True):
	'''Pre: True
	Post: retuns lp with block all N source and LP Sol for biomass production on inorganic N source'''
	lp = BuildLP.BlockAAUptakeLP(m,Biomass)
	for tx in ModelResults.SubTx(m, ion): # block NH4 media tx, does not include NH4 exporter
		lp.SetFixedFlux({tx:0.0})
	return lp

def SingleNSource(m,ion='AMMONIUM',Biomass=True,lp=None): #To generate results in section Single nitrogen sources
	'''Pre: True
	Post: returns dict with keys as single N source that can lead to Biomass production'''

	substrates = ModelResults.GetAllAASubstrates(m)
	rv = {}
	lp = GetSingleNLP(m,ion,Biomass)
	nf = []
	
	for subs in substrates:
		transporters = ModelResults.SubTx(m, subs)
		lp.ClearFluxConstraints(transporters)
		lp.Solve()
		if lp.IsStatusOptimal():
			sol = lp.GetPrimSol()
			rv[subs] = (sol,lp.GetObjVal(),len(sol.keys()))
		else:
			nf.append(subs)
		for tx in transporters:
			lp.SetFixedFlux({tx:0.0})
	print ('Substrates not utilised as single N source', nf)
	return rv

def SingleSubsRemoval(m):#To generate results in section Identification of substrate auxotrophies
	'''Pre: True
	Post: return dict with keys as substrate and values as tuple of LP sol, obj value and no. of reacs on individual substrate removal'''
	rv = {}
	aux = []
	substrates = ModelResults.GetAllSubstrates(m)
	lp = ModelResults.MediaBound(m)
	lp.Solve(False)
	sol = lp.GetPrimSol()
	rv['MHHW'] = (sol,lp.GetObjVal(),len(sol.keys()))
	
	for subs in substrates:
		#print (subs)
		lp = ModelResults.MediaBound(m)
		for tx in ModelResults.SubTx(m,subs):
			lp.SetFixedFlux({tx:0.0})
		lp.Solve(False)
		if lp.IsStatusOptimal():
			sol = lp.GetPrimSol()
			rv[subs] = (sol,lp.GetObjVal(),len(sol.keys()))
		else:
			aux.append(subs)
	print ('Auxotrophy substrates', aux)
	print ("Change in obj val and no. of reacs w.r.t MHHW")
	for subs in rv.keys():
		print (subs,rv[subs][1]-rv['MHHW'][1],rv[subs][2]-rv['MHHW'][2])
	
	return rv
		



	




