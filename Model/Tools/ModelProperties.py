
######
########
##########   This module is used to define the general properties of the model
##########    
########        
######

from ScrumPy.Util import Set


## Gene - reaction association

def RemoveSplit(reac):
    
    """pre: r = a reaction
       post: returns reaction name after removing -(NADP)/-(NAD)"""
    
    r = reac.rsplit('-(NADP)')[0].rsplit('-(NAD)')[0] #reverse split
    
    return r

def Reactions(m):

    """pre: True
       post: returns list of reactions excluding transporters"""
    
    rv = []
    
    tx = list(filter(lambda s: "_tx" in s, m.sm.cnames))
    
    for reac in Set.Complement(m.sm.cnames, tx):
            r = RemoveSplit(reac) # removes -NAD/NADP
            rv.append(r)

    return rv

def SponRxn(db,reacs=[]):
	'''pre:True
	post:list of spontaneous reactions'''
	rv = []	
	if reacs == []:
		reacs = list(db.dbs['REACTION'].keys())
	for r in reacs:
		if r in list(db.dbs['REACTION'].keys()) and 'SPONTANEOUS?' in list(db[r].keys()) and 'T' in db[r]['SPONTANEOUS?']:
			rv.append(r)
	return rv

def ReacAsso(m,db,reactions=[]):
    
    """pre: True
       post: returns dictionary with key as all reacs in model excluding transporters
              and values are lists of associated genes
              "No Gene" if reaction is not associated to genes
              "Reaction UID not found" if reac not present in database"""
        
    rv = {}
        
    if reactions == []:
        reactions = Reactions(m)

    for r in reactions:
        if r in db.dbs['REACTION'].keys():
            genes=[]
            genes = db[r].GetGenes()
            if len(genes) != 0:
                rv[r] = genes
            else:
                rv[r] = "No genes"
        else:
            rv[r] = "Reaction UID not found"

    return rv

def Genes(m,db, reactions=[]):

    """pre: True
       post: returns list of genes associated with reactions in the model"""

    rv = []
    ra = ReacAsso(m,db, reactions)
    for r in ra.keys():
        if type(ra[r]) is list:
            for g in ra[r]:
                if g not in rv:
                    rv.append(g)
    return rv


def ReacToGene(m,db, reactions=[]):
    
    """pre: True
       post: returns dictionary with reacs with associated genes as keys
             and list of associated genes as values"""
        
    rv = {}
    r2g = ReacAsso(m,db,reactions)

    NoGene = []
    NoUID = []

    tx = list(filter(lambda s: "_tx" in s, m.sm.cnames))
    for reac in Set.Complement(m.sm.cnames, tx):
        r = RemoveSplit(reac)
        if r in list(r2g.keys()): # remember reacs names here have -NADP/NAD split!
            if r2g[r] == 'No genes': 
                NoGene.append(reac) # so add full reac name to list (before split occured)
            elif r2g[r] == "Reaction UID not found":
                NoUID.append(reac)
            else:
                rv[reac] = r2g[r]
                    
    #print ("Total no of reactions with gene association excluding tx:",len(rv.keys()))
    #print ("Total no of reactions with no gene association excluding tx:",len(NoGene))
    #print ("Total no of reactions with no UID in db excluding tx:",len(NoUID))

    return rv, NoGene, NoUID # make the dict and list available

## General model properties

def ModelProperty(m,db):

    """pre: m = model
       post: returns lists of total reactions, reactions excluding tx,
       transporters, metabolites, orphan metabolites, dead reactions,
       and reactions subsets"""
    sponrxn = SponRxn(db,m.sm.cnames)
    Transporters = list(filter(lambda s: "_tx" in s, m.smx.cnames))   # list
    gene,nogene,nouid = ReacToGene(m,db)
    print ('Total no of reactions including transporters:', len(m.sm.cnames))
    print ('Total no of transporters:', len(Transporters))
    print ('Total no of reactions excluding transporters:', len(Set.Complement(m.sm.cnames, Transporters)))
    print ('Total no of internal metabolites:', len(m.sm.rnames))
    print ('Total no of orphan metabolites:', len(m.sm.OrphanMets()))
    print ('Total no of dead reactions:', len(m.DeadReactions()))
    print ('Total no of spontaneous reaction:', len(sponrxn)) 
    print ('Total no of reaction with gene asso, no gene asso and no UID resp (update for new db):', len(gene.keys()),len(nogene)-len(sponrxn), len(nouid))




