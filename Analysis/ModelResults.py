##Python source code for ScrumPy3
from importlib import reload
from ScrumPy.Util import Set
import BuildLP
reload(BuildLP)
import Media
reload(Media)
media = Media.media

empform = {'ACETOIN':{'C': 4, 'H': 8, 'O': 2},
	   'BUTANEDIOL':{'C': 4, 'H': 10, 'O': 2}
	  }

def MediaTx(m):
	'''Pre: model 
	Post: return list of media transporters'''
	return list(filter(lambda s:'_mm_tx' in s, m.sm.cnames))

def AAMediaTx(m):
	'''Pre: model 
	Post: return list of amion acids media transporters'''
	return list(filter(lambda s:'_AA_mm_tx' in s, m.sm.cnames))

def GetAllSubstrates(m):
	'''Pre: model 
	Post: return list of substrates in media'''
	rv = []
	for tx in MediaTx(m):
		subs = tx.split('_')[0]
		if subs not in rv:
			rv.append(subs)
	return rv

def GetAllAASubstrates(m):
	'''Pre: model 
	Post: return list of amino acids in media'''
	rv = []
	for tx in AAMediaTx(m):
		subs = tx.split('_')[0]
		if subs not in rv:
			rv.append(subs)
	return rv

def GetAllByProd(m):
	'''Pre: model 
	Post: return list of by-products'''
	rv = []
	for tx in list(filter(lambda s:"_bp_tx" in s,m.sm.cnames)):
		subs = tx.split('_bp_tx')[0]
		if subs not in rv:
			rv.append(subs)
	return rv

def SubTx(m, substrate):
	'''Pre: model and substrate
	Post: return all transporters for a given substrate. Note that there can be more than one transporter for a substrate'''
	rv = []
	for tx in MediaTx(m):
		if substrate == tx.split('_')[0]:
			rv.append(tx)
	for tx in list(filter(lambda s:"_bp_tx" in s,m.sm.cnames)):
		if substrate == tx.split('_bp_tx')[0]:
			rv.append(tx)
	return rv 

def SubBound(m,lp,substrate):
	'''Pre: model, lp and substrate
	Post: return lp with upper bound on substrate uptake'''
	tx = SubTx(m,substrate)
	if len(tx)==1.0:
		lp.SetFluxBounds({tx[0]:(0.0,media[substrate])})
	else:
		lp.SetSumFluxConstraint(tx,media[substrate],"Total"+substrate)
		lp.SetRowBounds("Total"+substrate,0,media[substrate])
	return lp

def MediaBound(m, AMaint, lp=None, fd=None):
	'''Pre: True
	Post: return lp with upper bound in media uptake'''
	substrates = GetAllSubstrates(m)
	if lp == None:
		lp = BuildLP.BiomassLP(m,AMaint,fd)
	for substrate in substrates:
		lp = SubBound(m,lp,substrate)
	return lp


def EssRxn(m,lp=None, AMaint=45.0,fd=None): # To generate essential reaction result in section Model simulation on defined rich media
	'''pre:True
	post: list of reactions essential for biomass production'''
	rv = []
	if lp == None:
		lp = BuildLP.BiomassLP(m,AMaint,fd)
	for r in list(filter(lambda s:'_tx' not in s,m.sm.cnames))+['ATPASE-RXN']:
		lp.SetFixedFlux({r:0})
		lp.Solve(False)
		if lp.GetStatusMsg() != 'optimal':
			rv.append(r)
		lp.ClearFluxConstraint(r)
	print ("No. of essential reactions",len(rv))	
	return rv

def ATPSol(m, AMaint=45.0, aerobic=True, mediabound=True, fd=None): #To generate result in Energy metabolism section 
	"""Pre:True
	Post: lp sol for ATP synthesis"""
	lp = BuildLP.BuildLP(m)
	lp.SetFixedFlux({'ATPASE-RXN':AMaint})
	if aerobic == False:
		lp.SetFixedFlux({'OXYGEN-MOLECULE_tx':0.0})
	if mediabound == True:
		lp = MediaBound(m,AMaint,lp,fd)
	lp.Solve()
	if lp.IsStatusOptimal():
		sol = lp.GetPrimSol()
	print ("no. of reacs",len(sol.keys()),"\n","Net import of aa",list(filter(lambda s:"_AA_mm_tx" in s,sol.keys())))
	return sol

def LPSol(m,mediabound=True,lp=None,AMaint=45.0,fd=None): # To generate results in section Model simulation on defined rich media
	"""Pre:True
	Post: lp sol for biomass synthesis"""
	if mediabound == True:
		lp = MediaBound(m,AMaint,lp,fd)
	else:
		lp = BuildLP.BiomassLP(m,AMaint,fd)
	lp.Solve()
	sol = lp.GetPrimSol()
	print ("Obj val",lp.GetObjVal(),"\n","no. of reactions excluding transporters",len(list(filter(lambda s:'_tx' not in s, sol.keys()))))
	return sol
	
def TotalSubFlux(m,substrate,sol):
	"""Pre:True
	Post: Total substrate uptake flux"""
	transporters = SubTx(m,substrate)
	influx = 0.0
	for tx in transporters:
		if tx in sol.keys():
			influx += sol[tx]
	return influx 

def NetAAtxFlux(m,sol=None):# To generate net substrate import results in section Model simulation on defined rich media
	'''Pre: True
	Post: net AA substrate influx (after substration flux in respective biomass transporter)'''
	rv = {}
	aa = []
	substrates = GetAllAASubstrates(m)
	if sol == None:
		sol = LPSol(m)
	for subs in substrates:
		influx = TotalSubFlux(m,subs,sol)
		bm = subs+'_AA_bm_tx'
		if bm in sol.keys():
			print (subs, influx+sol[bm])
			rv[subs] = influx+sol[bm]
	for subs in list(rv.keys()):
		if round(rv[subs],8) == 0:
			aa.append(subs)
	print ("Amino acids net import to biomass",aa)
	return rv

def ElementalCont(m,db,sol=None):
	'''Pre: model, database, lp solution
	Post: return carbon and nitrogen contribution by each substrate in the lp solution'''
	rv = {}
	byprod = GetAllByProd(m)
	substrates = GetAllSubstrates(m)
	if sol == None:
		sol = LPSol(m)
	for subs in substrates+byprod:
		emp = {}
		if subs in db.dbs['Compound'].keys(): 
			if 'C' in db[subs].EmpForm.keys():
				emp['C'] = TotalSubFlux(m,subs,sol)*db[subs].EmpForm['C']
			if 'N' in db[subs].EmpForm.keys():
				emp['N'] = TotalSubFlux(m,subs,sol)*db[subs].EmpForm['N']

		else:
			if subs in empform.keys(): 
				if 'C' in empform[subs].keys():
					emp['C'] = TotalSubFlux(m,subs,sol)*empform[subs]['C']
				if 'N' in empform[subs].keys():
					emp['N'] = TotalSubFlux(m,subs,sol)*empform[subs]['N']
		rv[subs]=emp		
	return rv

def PercentCont(m,db,Input='Total',sol=None): #To generate results, C and N contribution by each substrate, in section Model simulation on defined rich media and table 2
	'''Pre: model, database, lp solution
	Post: return percentage carbon and nitrogen contribution by each substrate in the lp solution'''
	rv = {}
	elecont = ElementalCont(m,db,sol)
	totalCin = 0.0
	totalbpC = 0.0
	totalNin = 0.0
	totalbpN = 0.0
	if "AMMONIUM_tx" in sol.keys():
		totalbpN = abs(sol["AMMONIUM_tx"])
	for subs in elecont.keys():
		if 'C' in elecont[subs].keys(): 
			if elecont[subs]['C']< 0.0:
				totalbpC+= abs(elecont[subs]['C'])
			else:
				totalCin += elecont[subs]['C']
		if 'N' in elecont[subs].keys():
			totalNin += elecont[subs]['N']

	totalbmN = totalNin-totalbpN
	totalbmC = totalCin-totalbpC
	print ('Total C and N input', totalCin, totalNin)		
	print ('Total C and N in biomass', totalbmC, totalbmN)
	#print ("****percentage C and N contribution******")
	for subs in GetAllSubstrates(m):
		pemp = {}
		if 'C' in elecont[subs]:
			if Input == 'Total':
				pemp['C'] = (elecont[subs]['C']/totalCin)*100
			if Input == 'Biomass':
				pemp['C'] = (elecont[subs]['C']/totalbmC)*100
		if 'N' in elecont[subs]:
			if Input == 'Total':
				pemp['N'] = (elecont[subs]['N']/totalNin)*100
			if Input == 'Biomass':
				pemp['N'] = (elecont[subs]['N']/totalbmN)*100
		rv[subs] = pemp
		print (subs,pemp)
	print ('****total C and N excretion',totalbpC,totalbpN)
	return rv
