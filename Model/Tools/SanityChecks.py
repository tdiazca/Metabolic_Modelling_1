#
##
### This module checks for reacs with atomic imbalance and finds unknown compounds
##
#


from ScrumPy.Bioinf import PyoCyc
from ScrumPy.Util import Set

path = '.'
def RemoveSuffix(reac): # rm suffix NADP or NAD added by ScrumPy
	"""Pre: a reaction;
	Post: Returns reaction name after removing suffix and -(NADP)/-(NAD)"""
	r = reac.rsplit('_',1)[0].split('-(NADP)')[0].split('-(NAD)')[0]
	return r

def Reactions(m): # rm suffix for reacs excluding tx
	rv = []
	tx = list(filter(lambda s:"_tx" in s,m.sm.cnames))
	for reac in Set.Complement(m.sm.cnames,tx):
		rv.append(RemoveSuffix(reac))
	return rv
    
def Updatedb(db,path=path,ExtraCompounds='ExtraCompounds.dat'):
	'''pre: True
	Post: Update db with mol formula of compounds in ExtraCompounds.dat'''

	if ExtraCompounds !=None:
		updatedb = PyoCyc.Compound.DB(path,ExtraCompounds)
		for met in list(updatedb.keys()):
			db.dbs['Compound'][met] = updatedb[met]
	
	return db

def CheckImbals(m,db,path=path,ExtraCompounds='ExtraCompounds.dat',reacs=[]):
	'''Pre:True
	Post: return dictionary of imbalanced reactions'''

	rv = {}
	db = Updatedb(db,path,ExtraCompounds)

	if reacs == []:
		reacs = list(filter(lambda x: not "_tx" in x, m.smx.cnames))
	for reac in reacs:
		stod = m.smx.InvolvedWith(reac) #mets involved in reac
		imbal =  db.dbs['Compound'].AtomImbal(stod)
		if len(imbal) > 0:
			rv[reac] = imbal
	return rv
                
def FindUnknownCompounds(imbaldic):
	""" pre: imbaldic = CheckImbals(..)
	post: list of Compound names in imbaldic that were prefixed as Unknown compound"""

	rvd = {}
	for ib in imbaldic.values():
		names = list(ib.keys())
		for name in names:
			if "Unknown" in name:
				name = name.rsplit(" ",1)[1]
				rvd[name] =1
	return rvd.keys()

def  SeperateUnknowns(imbaldic):
	'''Pre:imbaldic = CheckImbals(..)
	Post: remove reactions involedwith Unknown compound from imbalanced dictionary'''
	rv = {}
	for reac in list(imbaldic.keys()):
		imbal = imbaldic[reac]
		names = "".join(imbal.keys())
		if "Unknown" in names:
			rv[reac] = imbal
			del imbaldic[reac]
	return rv
   

