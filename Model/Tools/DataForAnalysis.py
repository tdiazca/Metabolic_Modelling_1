
## ProdUnitsBM.py

## ATPprod.py

### n of N and C atoms per AA in biomass

Natoms_BM_AA = {"ASP_AA_bm_tx":1,"SER_AA_bm_tx":1,"THR_AA_bm_tx":1,"ASN_AA_bm_tx":2,"ILE_AA_bm_tx":1, 
        "PHE_AA_bm_tx":1,"LYS_AA_bm_tx":2,"TYR_AA_bm_tx":1,"CYS_AA_bm_tx":1,"GLY_AA_bm_tx":1,
        "HIS_AA_bm_tx":3,"GLT_AA_bm_tx":1,"MET_AA_bm_tx":1,"TRP_AA_bm_tx":2,"GLN_AA_bm_tx":2,
        "LEU_AA_bm_tx":1,"ARG_AA_bm_tx":4,"VAL_AA_bm_tx":1,"ALA_AA_bm_tx":1,"PRO_AA_bm_tx":1}

Catoms_BM_AA = {"ASP_AA_bm_tx":4,"SER_AA_bm_tx":3,"THR_AA_bm_tx":4,"ASN_AA_bm_tx":4,"ILE_AA_bm_tx":6, 
        "PHE_AA_bm_tx":9,"LYS_AA_bm_tx":6,"TYR_AA_bm_tx":9,"CYS_AA_bm_tx":3,"GLY_AA_bm_tx":2,
        "HIS_AA_bm_tx":6,"GLT_AA_bm_tx":5,"MET_AA_bm_tx":5,"TRP_AA_bm_tx":11,"GLN_AA_bm_tx":5,
        "LEU_AA_bm_tx":6,"ARG_AA_bm_tx":6,"VAL_AA_bm_tx":5,"ALA_AA_bm_tx":3,"PRO_AA_bm_tx":5}

### n of C atoms per comp in MinimalMedium
Catoms_MinMed = {
           "GLC_mm_tx":6,
           "GLC_PTS_mm_tx":6, # PTS system
           }

### n of N atoms per comp in MinimalMedium                                       
Natoms_MinMed = {
           "NO3_mm_tx":1,"NH4_mm_tx":1}

### n of C atoms per AA in MM medium (AAs only)
Catoms_MM_AA = {"ASP_AA_mm_tx":4,"SER_AA_mm_tx":3,"THR_AA_mm_tx":4,"ASN_AA_mm_tx":4,"ILE_AA_mm_tx":6, 
           "PHE_AA_mm_tx":9,"LYS_AA_mm_tx":6,"TYR_AA_mm_tx":9,"CYS_AA_mm_tx":3,"GLY_AA_mm_tx":2,
           "HIS_AA_mm_tx":6,"GLT_AA_mm_tx":5,"MET_AA_mm_tx":5,"TRP_AA_mm_tx":11,"GLN_AA_mm_tx":5,
           "LEU_AA_mm_tx":6,"ARG_AA_mm_tx":6,"VAL_AA_mm_tx":5,"ALA_AA_mm_tx":3,"PRO_AA_mm_tx":5,
           "Ala_gly_AA_mm_tx":5,"ARG_ORN_antip_AA_mm_tx":1,"GLT_NA_AA_mm_tx":5
             }  # Includes GLN

### n of N atoms per AA in MM medium (AAs only)                                    
Natoms_MM_AA = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
           "Ala_gly_AA_mm_tx":2,"ARG_ORN_antip_AA_mm_tx":2,"GLT_NA_AA_mm_tx":1
           }  # Includes GLN

### n of C atoms per comp in MM medium
Catoms_MM = {"ASP_AA_mm_tx":4,"SER_AA_mm_tx":3,"THR_AA_mm_tx":4,"ASN_AA_mm_tx":4,"ILE_AA_mm_tx":6, 
           "PHE_AA_mm_tx":9,"LYS_AA_mm_tx":6,"TYR_AA_mm_tx":9,"CYS_AA_mm_tx":3,"GLY_AA_mm_tx":2,
           "HIS_AA_mm_tx":6,"GLT_AA_mm_tx":5,"MET_AA_mm_tx":5,"TRP_AA_mm_tx":11,"GLN_AA_mm_tx":5,
           "LEU_AA_mm_tx":6,"ARG_AA_mm_tx":6,"VAL_AA_mm_tx":5,"ALA_AA_mm_tx":3,"PRO_AA_mm_tx":5,
           "GLC_mm_tx":6,"GLC_PTS_mm_tx":6,"Ala_gly_AA_mm_tx":5,"ARG_ORN_antip_AA_mm_tx":1,
           "GLT_NA_AA_mm_tx":5,
             }  # Includes GLC ; includes GLN

### n of N atoms per comp in MM medium                                      
Natoms_MM = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
           "NO3_mm_tx":1,"NH4_mm_tx":1,"Ala_gly_AA_mm_tx":2,"ARG_ORN_antip_AA_mm_tx":2,
           "GLT_NA_AA_mm_tx":1,}  # Includes NO3; includes NH4 ; includes GLN


### n of C atoms per AA in synovial fluid (AAs only) 
Catoms_synovial_AA = {"ASP_AA_mm_tx":4,"SER_AA_mm_tx":3,"THR_AA_mm_tx":4,"ASN_AA_mm_tx":4,"ILE_AA_mm_tx":6, 
           "PHE_AA_mm_tx":9,"LYS_AA_mm_tx":6,"TYR_AA_mm_tx":9,"CYS_AA_mm_tx":3,"GLY_AA_mm_tx":2,
           "HIS_AA_mm_tx":6,"GLT_AA_mm_tx":5,"MET_AA_mm_tx":5,"TRP_AA_mm_tx":11,"GLN_AA_mm_tx":5,
           "LEU_AA_mm_tx":6,"ARG_AA_mm_tx":6,"VAL_AA_mm_tx":5,"ALA_AA_mm_tx":3,"PRO_AA_mm_tx":5,
           "CITRULLINE_mm_tx":6,"ORNIT_mm_tx":5,"4-AMINO-BUTYRATE_mm_tx":4,"UREA_mm_tx":1,
           "Ala_gly_AA_mm_tx":5,"ARG_ORN_antip_AA_mm_tx":1,"GLT_NA_AA_mm_tx":5,"Urate_mm_tx":5}
            # Includes N sources as ORNITHINE etc

### n of N atoms per AA in synovial fluid (AAs only)                                   
Natoms_synovial_AA = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
           "CITRULLINE_mm_tx":3,"ORNIT_mm_tx":2,"4-AMINO-BUTYRATE_mm_tx":1,"UREA_mm_tx":2,
           "Ala_gly_AA_mm_tx":2,"ARG_ORN_antip_AA_mm_tx":2,"GLT_NA_AA_mm_tx":1,"Urate_mm_tx":4
            }  # Includes other N sources like ORNITHINE etc

### n of C atoms per comp in synovial fluid
Catoms_synovial = {"ASP_AA_mm_tx":4,"SER_AA_mm_tx":3,"THR_AA_mm_tx":4,"ASN_AA_mm_tx":4,"ILE_AA_mm_tx":6, 
           "PHE_AA_mm_tx":9,"LYS_AA_mm_tx":6,"TYR_AA_mm_tx":9,"CYS_AA_mm_tx":3,"GLY_AA_mm_tx":2,
           "HIS_AA_mm_tx":6,"GLT_AA_mm_tx":5,"MET_AA_mm_tx":5,"TRP_AA_mm_tx":11,"GLN_AA_mm_tx":5,
           "LEU_AA_mm_tx":6,"ARG_AA_mm_tx":6,"VAL_AA_mm_tx":5,"ALA_AA_mm_tx":3,"Ala_gly_AA_mm_tx":5,
           "PRO_AA_mm_tx":5,"GLC_mm_tx":6,"GLC_PTS_mm_tx":6,"CITRULLINE_mm_tx":6,"ORNIT_mm_tx":5,
           "4-AMINO-BUTYRATE_mm_tx":4,"UREA_mm_tx":1,"Ala_gly_AA_mm_tx":5,"ARG_ORN_antip_AA_mm_tx":1,
           "GLT_NA_AA_mm_tx":5,"Urate_mm_tx":5
            }  # Includes GLC; includes other C sources like ORNITHINE etc

### n of N atoms per comp in synovial fluid                                    
Natoms_synovial = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
           "NO3_mm_tx":1,"NH4_mm_tx":1,"CITRULLINE_mm_tx":3,"ORNIT_mm_tx":2,"4-AMINO-BUTYRATE_mm_tx":1,
           "UREA_mm_tx":2,"Ala_gly_AA_mm_tx":2,"ARG_ORN_antip_AA_mm_tx":2,"GLT_NA_AA_mm_tx":1,
           "Urate_mm_tx":4
            }  # Includes NO3; includes NH4; includes other N sources like ORNITHINE etc

Natoms_synovial_noNH4 = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
           "NO3_mm_tx":1,"CITRULLINE_mm_tx":3,"ORNIT_mm_tx":2,"4-AMINO-BUTYRATE_mm_tx":1,
           "UREA_mm_tx":2,"Ala_gly_AA_mm_tx":2,"ARG_ORN_antip_AA_mm_tx":2,"GLT_NA_AA_mm_tx":1,
           "Urate_mm_tx":4
            }  # Includes NO3; EXCLUDES NH4; includes other N sources like ORNITHINE etc

CatomByProds = {"CO2_bp_tx":1,"D-LACTATE_bp_tx":3,"L-LACTATE_bp_tx":3,"SUC_bp_tx":4,"ACET_bp_tx":2,
                "FORMATE_bp_tx":1,"ETOH_bp_tx":2,"2-KETOGLUTARATE_bp_tx":5,"ACETOIN_bp_tx":4,
                "BUTYRIC_ACID_bp_tx":4,"BUTANOL_bp_tx":4,"BUTANEDIOL_bp_tx":4}

## N transporters potentially available to our model

NTxs_MinMed = [
    "NH4_mm_tx",
    "NO3_mm_tx"    ## important that we turn this off when we don't need it !!!
]

NTxs_MM = [
    #"NO3_mm_tx",    ## NOT in medium. Important that we turn this off when we don't need it !!!
    "NH4_mm_tx",
    "THR_AA_mm_tx",
    "LEU_AA_mm_tx",
    "PHE_AA_mm_tx",
    "TYR_AA_mm_tx",
    "HIS_AA_mm_tx",
    #"GLN_AA_mm_tx", ## NOT in medium. Reinstate for completeness if needed.
    "GLT_AA_mm_tx",
    "SER_AA_mm_tx",
    "MET_AA_mm_tx",
    "ASN_AA_mm_tx",
    "PRO_AA_mm_tx",
    "ALA_AA_mm_tx",
    "ILE_AA_mm_tx",
    "ASP_AA_mm_tx",
    "GLY_AA_mm_tx",
    "CYS_AA_mm_tx",
    "LYS_AA_mm_tx",
    "TRP_AA_mm_tx",
    "VAL_AA_mm_tx",
    "ARG_AA_mm_tx",
    "Ala_gly_AA_mm_tx",
    "ARG_ORN_antip_AA_mm_tx",
    "GLT_NA_AA_mm_tx",
]

NTxs_synovial = [
    "NO3_mm_tx",    ## important that we turn this off when we don't need it !!!
    "NH4_mm_tx",
    "THR_AA_mm_tx",
    "LEU_AA_mm_tx",
    "PHE_AA_mm_tx",
    "TYR_AA_mm_tx",
    "HIS_AA_mm_tx",
    "GLN_AA_mm_tx",
    "GLT_AA_mm_tx",
    "SER_AA_mm_tx",
    "MET_AA_mm_tx",
    "ASN_AA_mm_tx",
    "PRO_AA_mm_tx",
    "ALA_AA_mm_tx",
    "ILE_AA_mm_tx",
    "ASP_AA_mm_tx",
    "GLY_AA_mm_tx",
    "CYS_AA_mm_tx",
    "LYS_AA_mm_tx",
    "TRP_AA_mm_tx",
    "VAL_AA_mm_tx",
    "ARG_AA_mm_tx",
    "CITRULLINE_mm_tx", # "L-CITRULLINE"
    "ORNIT_mm_tx", # "L-ORNITHINE"
    "4-AMINO-BUTYRATE_mm_tx", # "4-AMINO-BUTYRATE"
    "UREA_mm_tx" # "UREA"
    "Ala_gly_AA_mm_tx",
    "ARG_ORN_antip_AA_mm_tx",
    "GLT_NA_AA_mm_tx",
    "Urate_mm_tx"
]

# Extra n sources in synovial fluid vs MM medium =
## ExtraNTxs = [
##    'NO3_mm_tx', # "NITRATE"
##    'CITRULLINE_mm_tx', # "L-CITRULLINE"
##    'ORNIT_mm_tx', # "L-ORNITHINE"
##    '4-AMINO-BUTYRATE_mm_tx', # "4-AMINO-BUTYRATE"
##    'UREA_mm_tx' # "UREA"
##]

# C transporters potentially available to our model

CTxs_MinMed = [
    'GLC_mm_tx'
    'GLC_PTS_mm_tx'  # PTS system
]

CTxs_MM = [
    'THR_AA_mm_tx',
    'LEU_AA_mm_tx',
    'PHE_AA_mm_tx',
    'TYR_AA_mm_tx',
    'HIS_AA_mm_tx',
    #'GLN_AA_mm_tx', ## NOT in medium. Reinstate for completeness if needed.
    'GLT_AA_mm_tx',
    'SER_AA_mm_tx',
    'MET_AA_mm_tx',
    'ASN_AA_mm_tx',
    'PRO_AA_mm_tx',
    'ALA_AA_mm_tx',
    'ILE_AA_mm_tx',
    'ASP_AA_mm_tx',
    'GLY_AA_mm_tx',
    'CYS_AA_mm_tx',
    'LYS_AA_mm_tx',
    'TRP_AA_mm_tx',
    'VAL_AA_mm_tx',
    'ARG_AA_mm_tx',
    'GLC_PTS_mm_tx',  # PTS system
    'GLC_mm_tx',
    'Ala_gly_AA_mm_tx',
    'ARG_ORN_antip_AA_mm_tx',
    'GLT_NA_AA_mm_tx',
]

CTxs_synovial = [
    'THR_AA_mm_tx',
    'LEU_AA_mm_tx',
    'PHE_AA_mm_tx',
    'TYR_AA_mm_tx',
    'HIS_AA_mm_tx',
    'GLN_AA_mm_tx',
    'GLT_AA_mm_tx',
    'SER_AA_mm_tx',
    'MET_AA_mm_tx',
    'ASN_AA_mm_tx',
    'PRO_AA_mm_tx',
    'ALA_AA_mm_tx',
    'ILE_AA_mm_tx',
    'ASP_AA_mm_tx',
    'GLY_AA_mm_tx',
    'CYS_AA_mm_tx',
    'LYS_AA_mm_tx',
    'TRP_AA_mm_tx',
    'VAL_AA_mm_tx',
    'ARG_AA_mm_tx',
    'GLC_mm_tx',
    'GLC_PTS_mm_tx', # PTS system
    'CITRULLINE_mm_tx', # "L-CITRULLINE" # C 6
    'ORNIT_mm_tx', # "L-ORNITHINE" # C 5
    '4-AMINO-BUTYRATE_mm_tx', # "4-AMINO-BUTYRATE" # C 4
    'UREA_mm_tx', # "UREA" # C 1
    'Ala_gly_AA_mm_tx',
    'ARG_ORN_antip_AA_mm_tx',
    'GLT_NA_AA_mm_tx',
    'Urate_mm_tx'
]

## Extra n sources in synovial fluid vs MM medium =
## ExtraNTxs = [
##    'NO3_mm_tx', # "NITRATE" # C 0
##    'CITRULLINE_mm_tx', # "L-CITRULLINE" # C 6
##    'ORNIT_mm_tx', # "L-ORNITHINE" # C 5
##    '4-AMINO-BUTYRATE_mm_tx', # "4-AMINO-BUTYRATE" # C 4
##    'UREA_mm_tx' # "UREA" # C 1
##]

# AA defined with AAutilN.UtilSingleAA(m) as non-viable single N sources.
NonViabAA = [
    #'THR_AA_mm_tx',    # THR is viable after curation
    'LEU_AA_mm_tx',
    'PHE_AA_mm_tx',
    'TYR_AA_mm_tx',
    'MET_AA_mm_tx',
    'ILE_AA_mm_tx',
    'LYS_AA_mm_tx',
    'TRP_AA_mm_tx'
]

############# Dictionaries below are for ols Txs (before July 2020)
####
###### ProdUnitsBM.py
####
###### ATPprod.py
####
####### n of N and C atoms per AA in biomass
####
####Natoms_BM_AA = {"ASP_AA_bm_tx":1,"SER_AA_bm_tx":1,"THR_AA_bm_tx":1,"ASN_AA_bm_tx":2,"ILE_AA_bm_tx":1, 
####        "PHE_AA_bm_tx":1,"LYS_AA_bm_tx":2,"TYR_AA_bm_tx":1,"CYS_AA_bm_tx":1,"GLY_AA_bm_tx":1,
####        "HIS_AA_bm_tx":3,"GLT_AA_bm_tx":1,"MET_AA_bm_tx":1,"TRP_AA_bm_tx":2,"GLN_AA_bm_tx":2,
####        "LEU_AA_bm_tx":1,"ARG_AA_bm_tx":4,"VAL_AA_bm_tx":1,"ALA_AA_bm_tx":1,"PRO_AA_bm_tx":1}
####
####Catoms_BM_AA = {"ASP_AA_bm_tx":4,"SER_AA_bm_tx":3,"THR_AA_bm_tx":4,"ASN_AA_bm_tx":4,"ILE_AA_bm_tx":6, 
####        "PHE_AA_bm_tx":9,"LYS_AA_bm_tx":6,"TYR_AA_bm_tx":9,"CYS_AA_bm_tx":3,"GLY_AA_bm_tx":2,
####        "HIS_AA_bm_tx":6,"GLT_AA_bm_tx":5,"MET_AA_bm_tx":5,"TRP_AA_bm_tx":11,"GLN_AA_bm_tx":5,
####        "LEU_AA_bm_tx":6,"ARG_AA_bm_tx":6,"VAL_AA_bm_tx":5,"ALA_AA_bm_tx":3,"PRO_AA_bm_tx":5}
####
####### n of C atoms per comp in MinimalMedium
####Catoms_MinMed = {
####           "GLC_mm_tx":6}
####
####### n of N atoms per comp in MinimalMedium                                       
####Natoms_MinMed = {
####           "NO3_mm_tx":1,"NH4_mm_tx":1}
####
####### n of C atoms per AA in MM medium (AAs only)
####Catoms_MM_AA = {"ASP_AA_mm_tx":4,"SER_AA_mm_tx":3,"THR_AA_mm_tx":4,"ASN_AA_mm_tx":4,"ILE_AA_mm_tx":6, 
####           "PHE_AA_mm_tx":9,"LYS_AA_mm_tx":6,"TYR_AA_mm_tx":9,"CYS_AA_mm_tx":3,"GLY_AA_mm_tx":2,
####           "HIS_AA_mm_tx":6,"GLT_AA_mm_tx":5,"MET_AA_mm_tx":5,"TRP_AA_mm_tx":11,"GLN_AA_mm_tx":5,
####           "LEU_AA_mm_tx":6,"ARG_AA_mm_tx":6,"VAL_AA_mm_tx":5,"ALA_AA_mm_tx":3,"PRO_AA_mm_tx":5
####             }  # Includes GLN
####
####### n of N atoms per AA in MM medium (AAs only)                                    
####Natoms_MM_AA = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
####           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
####           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
####           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
####           }  # Includes GLN
####
####### n of C atoms per comp in MM medium
####Catoms_MM = {"ASP_AA_mm_tx":4,"SER_AA_mm_tx":3,"THR_AA_mm_tx":4,"ASN_AA_mm_tx":4,"ILE_AA_mm_tx":6, 
####           "PHE_AA_mm_tx":9,"LYS_AA_mm_tx":6,"TYR_AA_mm_tx":9,"CYS_AA_mm_tx":3,"GLY_AA_mm_tx":2,
####           "HIS_AA_mm_tx":6,"GLT_AA_mm_tx":5,"MET_AA_mm_tx":5,"TRP_AA_mm_tx":11,"GLN_AA_mm_tx":5,
####           "LEU_AA_mm_tx":6,"ARG_AA_mm_tx":6,"VAL_AA_mm_tx":5,"ALA_AA_mm_tx":3,"PRO_AA_mm_tx":5,
####           "GLC_mm_tx":6}  # Includes GLC ; includes GLN
####
####### n of N atoms per comp in MM medium                                      
####Natoms_MM = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
####           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
####           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
####           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
####           "NO3_mm_tx":1,"NH4_mm_tx":1}  # Includes NO3; includes NH4 ; includes GLN
####
####
####### n of C atoms per AA in synovial fluid (AAs only) 
####Catoms_synovial_AA = {"ASP_AA_mm_tx":4,"SER_AA_mm_tx":3,"THR_AA_mm_tx":4,"ASN_AA_mm_tx":4,"ILE_AA_mm_tx":6, 
####           "PHE_AA_mm_tx":9,"LYS_AA_mm_tx":6,"TYR_AA_mm_tx":9,"CYS_AA_mm_tx":3,"GLY_AA_mm_tx":2,
####           "HIS_AA_mm_tx":6,"GLT_AA_mm_tx":5,"MET_AA_mm_tx":5,"TRP_AA_mm_tx":11,"GLN_AA_mm_tx":5,
####           "LEU_AA_mm_tx":6,"ARG_AA_mm_tx":6,"VAL_AA_mm_tx":5,"ALA_AA_mm_tx":3,"PRO_AA_mm_tx":5,
####           "CITRULLINE_mm_tx":6,"ORNIT_mm_tx":5,"4-AMINO-BUTYRATE_mm_tx":4,
####           "UREA_mm_tx":1}  # Includes N soureces as ORNITHINE etc
####
####### n of N atoms per AA in synovial fluid (AAs only)                                   
####Natoms_synovial_AA = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
####           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
####           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
####           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
####           "CITRULLINE_mm_tx":3,"ORNIT_mm_tx":2,"4-AMINO-BUTYRATE_mm_tx":1,
####           "UREA_mm_tx":2}  # Includes other N sources like ORNITHINE etc
####
####### n of C atoms per comp in synovial fluid
####Catoms_synovial = {"ASP_AA_mm_tx":4,"SER_AA_mm_tx":3,"THR_AA_mm_tx":4,"ASN_AA_mm_tx":4,"ILE_AA_mm_tx":6, 
####           "PHE_AA_mm_tx":9,"LYS_AA_mm_tx":6,"TYR_AA_mm_tx":9,"CYS_AA_mm_tx":3,"GLY_AA_mm_tx":2,
####           "HIS_AA_mm_tx":6,"GLT_AA_mm_tx":5,"MET_AA_mm_tx":5,"TRP_AA_mm_tx":11,"GLN_AA_mm_tx":5,
####           "LEU_AA_mm_tx":6,"ARG_AA_mm_tx":6,"VAL_AA_mm_tx":5,"ALA_AA_mm_tx":3,"PRO_AA_mm_tx":5,
####           "GLC_mm_tx":6,"CITRULLINE_mm_tx":6,"ORNIT_mm_tx":5,"4-AMINO-BUTYRATE_mm_tx":4,
####           "UREA_mm_tx":1}  # Includes GLC; includes other C sources like ORNITHINE etc
####
####### n of N atoms per comp in synovial fluid                                    
####Natoms_synovial = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
####           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
####           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
####           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
####           "NO3_mm_tx":1,"NH4_mm_tx":1,"CITRULLINE_mm_tx":3,"ORNIT_mm_tx":2,"4-AMINO-BUTYRATE_mm_tx":1,
####           "UREA_mm_tx":2}  # Includes NO3; includes NH4; includes other N sources like ORNITHINE etc
####
####Natoms_synovial_noNH4 = {"ASP_AA_mm_tx":1,"SER_AA_mm_tx":1,"THR_AA_mm_tx":1,"ASN_AA_mm_tx":2,"ILE_AA_mm_tx":1, 
####           "PHE_AA_mm_tx":1,"LYS_AA_mm_tx":2,"TYR_AA_mm_tx":1,"CYS_AA_mm_tx":1,"GLY_AA_mm_tx":1,
####           "HIS_AA_mm_tx":3,"GLT_AA_mm_tx":1,"MET_AA_mm_tx":1,"TRP_AA_mm_tx":2,"GLN_AA_mm_tx":2,
####           "LEU_AA_mm_tx":1,"ARG_AA_mm_tx":4,"VAL_AA_mm_tx":1,"ALA_AA_mm_tx":1,"PRO_AA_mm_tx":1,
####           "NO3_mm_tx":1,"CITRULLINE_mm_tx":3,"ORNIT_mm_tx":2,"4-AMINO-BUTYRATE_mm_tx":1,
####           "UREA_mm_tx":2}  # Includes NO3; EXCLUDES NH4; includes other N sources like ORNITHINE etc
####
####CatomByProds = {"CO2_bp_tx":1,"D-LACTATE_bp_tx":3,"L-LACTATE_bp_tx":3,"SUC_bp_tx":4,"ACET_bp_tx":2,
####                "FORMATE_bp_tx":1,"ETOH_bp_tx":2,"2-KETOGLUTARATE_bp_tx":5,"ACETOIN_bp_tx":4,
####                "BUTYRIC_ACID_bp_tx":4,"BUTANOL_bp_tx":4,"BUTANEDIOL_bp_tx":4}
####
###### N transporters potentially available to our model
####
####NTxs_MinMed = [
####    "NH4_mm_tx",
####    "NO3_mm_tx"    ## important that we turn this off when we don't need it !!!
####]
####
####NTxs_MM = [
####    #"NO3_mm_tx",    ## NOT in medium. Important that we turn this off when we don't need it !!!
####    "NH4_mm_tx",
####    "THR_AA_mm_tx",
####    "LEU_AA_mm_tx",
####    "PHE_AA_mm_tx",
####    "TYR_AA_mm_tx",
####    "HIS_AA_mm_tx",
####    #"GLN_AA_mm_tx", ## NOT in medium. Reinstate for completeness if needed.
####    "GLT_AA_mm_tx",
####    "SER_AA_mm_tx",
####    "MET_AA_mm_tx",
####    "ASN_AA_mm_tx",
####    "PRO_AA_mm_tx",
####    "ALA_AA_mm_tx",
####    "ILE_AA_mm_tx",
####    "ASP_AA_mm_tx",
####    "GLY_AA_mm_tx",
####    "CYS_AA_mm_tx",
####    "LYS_AA_mm_tx",
####    "TRP_AA_mm_tx",
####    "VAL_AA_mm_tx",
####    "ARG_AA_mm_tx"
####]
####
####NTxs_synovial = [
####    "NO3_mm_tx",    ## important that we turn this off when we don't need it !!!
####    "NH4_mm_tx",
####    "THR_AA_mm_tx",
####    "LEU_AA_mm_tx",
####    "PHE_AA_mm_tx",
####    "TYR_AA_mm_tx",
####    "HIS_AA_mm_tx",
####    "GLN_AA_mm_tx",
####    "GLT_AA_mm_tx",
####    "SER_AA_mm_tx",
####    "MET_AA_mm_tx",
####    "ASN_AA_mm_tx",
####    "PRO_AA_mm_tx",
####    "ALA_AA_mm_tx",
####    "ILE_AA_mm_tx",
####    "ASP_AA_mm_tx",
####    "GLY_AA_mm_tx",
####    "CYS_AA_mm_tx",
####    "LYS_AA_mm_tx",
####    "TRP_AA_mm_tx",
####    "VAL_AA_mm_tx",
####    "ARG_AA_mm_tx",
####    "CITRULLINE_mm_tx", # "L-CITRULLINE"
####    "ORNIT_mm_tx", # "L-ORNITHINE"
####    "4-AMINO-BUTYRATE_mm_tx", # "4-AMINO-BUTYRATE"
####    "UREA_mm_tx" # "UREA"
####]
####
##### Extra n sources in synovial fluid vs MM medium =
###### ExtraNTxs = [
######    'NO3_mm_tx', # "NITRATE"
######    'CITRULLINE_mm_tx', # "L-CITRULLINE"
######    'ORNIT_mm_tx', # "L-ORNITHINE"
######    '4-AMINO-BUTYRATE_mm_tx', # "4-AMINO-BUTYRATE"
######    'UREA_mm_tx' # "UREA"
######]
####
##### C transporters potentially available to our model
####
####CTxs_MinMed = [
####    'GLC_mm_tx'
####]
####
####CTxs_MM = [
####    'THR_AA_mm_tx',
####    'LEU_AA_mm_tx',
####    'PHE_AA_mm_tx',
####    'TYR_AA_mm_tx',
####    'HIS_AA_mm_tx',
####    #'GLN_AA_mm_tx', ## NOT in medium. Reinstate for completeness if needed.
####    'GLT_AA_mm_tx',
####    'SER_AA_mm_tx',
####    'MET_AA_mm_tx',
####    'ASN_AA_mm_tx',
####    'PRO_AA_mm_tx',
####    'ALA_AA_mm_tx',
####    'ILE_AA_mm_tx',
####    'ASP_AA_mm_tx',
####    'GLY_AA_mm_tx',
####    'CYS_AA_mm_tx',
####    'LYS_AA_mm_tx',
####    'TRP_AA_mm_tx',
####    'VAL_AA_mm_tx',
####    'ARG_AA_mm_tx',
####    'GLC_mm_tx'
####]
####
####CTxs_synovial = [
####    'THR_AA_mm_tx',
####    'LEU_AA_mm_tx',
####    'PHE_AA_mm_tx',
####    'TYR_AA_mm_tx',
####    'HIS_AA_mm_tx',
####    'GLN_AA_mm_tx',
####    'GLT_AA_mm_tx',
####    'SER_AA_mm_tx',
####    'MET_AA_mm_tx',
####    'ASN_AA_mm_tx',
####    'PRO_AA_mm_tx',
####    'ALA_AA_mm_tx',
####    'ILE_AA_mm_tx',
####    'ASP_AA_mm_tx',
####    'GLY_AA_mm_tx',
####    'CYS_AA_mm_tx',
####    'LYS_AA_mm_tx',
####    'TRP_AA_mm_tx',
####    'VAL_AA_mm_tx',
####    'ARG_AA_mm_tx',
####    'GLC_mm_tx',
####    'CITRULLINE_mm_tx', # "L-CITRULLINE" # C 6
####    'ORNIT_mm_tx', # "L-ORNITHINE" # C 5
####    '4-AMINO-BUTYRATE_mm_tx', # "4-AMINO-BUTYRATE" # C 4
####    'UREA_mm_tx' # "UREA" # C 1
####]
####
###### Extra n sources in synovial fluid vs MM medium =
###### ExtraNTxs = [
######    'NO3_mm_tx', # "NITRATE" # C 0
######    'CITRULLINE_mm_tx', # "L-CITRULLINE" # C 6
######    'ORNIT_mm_tx', # "L-ORNITHINE" # C 5
######    '4-AMINO-BUTYRATE_mm_tx', # "4-AMINO-BUTYRATE" # C 4
######    'UREA_mm_tx' # "UREA" # C 1
######]
####
##### AA defined with AAutilN.UtilSingleAA(m) as non-viable single N sources.
####NonViabAA = [
####    #'THR_AA_mm_tx',    # THR is viable after curation
####    'LEU_AA_mm_tx',
####    'PHE_AA_mm_tx',
####    'TYR_AA_mm_tx',
####    'MET_AA_mm_tx',
####    'ILE_AA_mm_tx',
####    'LYS_AA_mm_tx',
####    'TRP_AA_mm_tx'
####]
####
