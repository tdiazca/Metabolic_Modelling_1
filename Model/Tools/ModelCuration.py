from ScrumPy.Util import Set
import re

def ETCReacs(db,Type='Electron-Transfer'):
	'''pre:True
	post:list of electron transfer reactions in database'''
	rv = []
	for p in db.dbs['Pathway'].keys():
		if 'TYPES' in list(db[p].keys()) and Type in db[p]['TYPES']:
			for r in db[p]['REACTION-LIST']:
				if r not in rv:
					rv.append(r)
	return rv


def FindtRNA(m,pattern='RNA'):
	'''pre:True
	post:list of met associated with RNA in database'''
	rv = []
	for met in m.sm.rnames:
		z = re.search(pattern, met)
		if z != None:
			rv.append(met)
	return rv



