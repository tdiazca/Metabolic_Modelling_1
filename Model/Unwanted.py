# several of these reacs are not in current vs of PGDB

Reacs = [
"3.6.3.12-RXN", # Transporter K+,  EC Number: 3.6.3.12 (K+-transporting ATPase)
"3.6.3.16-RXN", # Transporter Arsenite (detoxification) , import
"ABC-27-RXN", # ABC phosphate transporter (ATPase), EC Number: 3.6.3.27 , import
"ATPSYN-RXN", # ATP-synthase, EC-3.6.3.14
"ALDOSE-1-EPIMERASE-RXN", # 'ALPHA-GLUCOSE', 'GLC'
"RXN-1685", # replaced by "GLUCOKIN-RXN" , since "ALPHA-GLUCOSE" : "GLC"
"GLUCOSE-6-PHOSPHATE-1-EPIMERASE-RXN", # 'ALPHA-GLC-6-P', 'GLC-6-P'
"RXN-9615", #  'AMMONIUM' #
"RXN-13519", # 'PROTON' # part of ETC: NADH-quinox oxreductase
"RXN-14904", # RIBOSE to RIBOSE
"RXN0-4082", # 'ec 3.2.1.14 - chitinase'
"RXN0-5304", # RIBOSE to RIBOSE
"RXNJJB-4", # Predicted ABC transporter of phosphonate
"RXNJJB-5", # Transporter MO+2
"RXNJJB-6", # Transporter CO+2
"RXNJJB-7", # Transporter PUTRESCINE
"RXNJJB-8", # Transporter SPERMIDINE 
"TRANS-RXN-1", # 'FORMATE' , rev
"TRANS-RXN-2", # Transporter Potassium (K(+)), EC-3.6.3.12
"TRANS-RXN-29", # 'PROTON', 'PRO'
"TRANS-RXN-104", # 'PROTON', 'L-LACTATE', rev PGDB 24.0
"TRANS-RXN-108", # 'PROTON', '|Nucleosides|', rev
"TRANS-RXN-132", # 'PROTON', 'URACIL'
"TRANS-RXN-137", # 'NITRITE', import PGDB 24.0
"TRANS-RXN-141", # 'MG+2', import
"TRANS-RXN-141A", # 'CO+2', import
"TRANS-RXN-249", # ATP-synthase, EC-3.6.3.14 
"TRANS-RXN0-209", # 'PROTON', 'GLUCONATE' import PGDB 24.0
"TRANS-RXN0-470", # '|Pi|'
"TRANS-RXN0-490", # 'BETAINE'
"TRANS-RXNJJB-1", # '|Anions|'
"TRANS-RXNJJB-2", # 'CARNITINE'
"TRANS-RXNJJB-3", # 'CHOLINE'
"TRANS-RXNJJB-4", # ''25-DIAMINOPENTANOATE', 'ARG'
"TRANS-RXNJJB-5", # '|Purines|', 'CYTOSINE''
"TRANS-RXNJJB-6", # '|C4-dicarboxylates|'
"TRANS-RXNJJB-7", # 'PROTON', '|Amino-Acids|'
"TRANS-RXNJJB-8", # 'GLYCEROL-3P'
"TRANS-RXNJJB-9", # '|Cations|'
"TRANS-RXNJJB-10", # 'PROTON', 'NA+'
"TRANS-RXNJJB-11", # 'NA+', '|Glutamates|'
"TRANS-RXNJJB-12", # '|Sugar|'
"TRANS-RXNJJB-13", # 'URACIL', 'XANTHINE'
"TRANS-RXNJJB-14", # 'PROTON', 'GLYCEROL'
"TRANS-RXNJJB-15", # 'PROTON', '|Cations|'
"TRANS-RXNJJB-16", # 'NA+', '|dicarboxylate|'
"TRANS-RXNJJB-17", # '|Amino-Acids|'
"TRANS-RXNJJB-18", # 'PROTON', 'SULFATE'
"RXN-14963", # EC-1.3.5.4 ,rm from latest BioCyc version
"RXN-13482", # replaced by "ORNCARBAMTRANSFER-RXN"
"RXN0-5265", # COE ;  4 E- + 1 OXYGEN-MOLECULE + 4 PROTON <> 2 WATER 
"RXN0-5248", # COE ;  2 E- + 1 PROTON + 1 NAD <> 1 NADH
"RXN0-5245", # COE;  2 E- + 1 FUM + 2 PROTON <> 1 SUC
"TRANS-RXN0-237",# COE; 2 E- + 1 Ubiquinones + 2 PROTON <> 1 Ubiquinols
"RXN0-6372", # COE; 2 E- + 1 GLC-D-LACTONE + 2 PROTON <> 1 GLC
"RXN0-6370", # COE; 1 NITRATE + 2 E- + 2 PROTON <> 1 WATER + 1 NITRITE
"RXN0-5244", # COE; 2 E- + 1 Ubiquinones + 2 PROTON <> 1 Ubiquinols
"RXN0-5242",# COE; 1 NITRATE + 2 E- + 2 PROTON <> 1 WATER + 1 NITRITE
"RXN0-5259", # COE; 2 E- + 1 ETR-Quinones + 2 PROTON <> 1 ETR-Quinols
"RXN-14451", # COE; 1 CPD-9612 + 2 E- + 2 PROTON <> 1 CPD-15301
"RXN0-6490", # COE;  2 E- + 1 OROTATE + 2 PROTON <> 1 DI-H-OROTATE
"NADH-DEHYDROG-A-RXN", # ETC/COE; replaced in ETC_core.spy by "NADH_DH_mena"
"RXN0-5330", # ETC/COE; replaced in ETC.spy
"RXN0-5266",  # ETC/COE; replaced in ETC.spy
"CYTOCHROME-C-OXIDASE-RXN", # ETC, replaced in ETC.spy by "Cytochrome_B_oxidase"
"1.10.2.2-RXN", # ETC, replaced in ETC.spy by 'Mena_oxid'
#"NITRATE-REDUCTASE-CYTOCHROME-RXN", # ETC, replaced in ETC.spy by "Nitrate_reductase"
#"1.7.2.2-RXN", # ETC; nitrite to NH3, replaced in ETC.spy by "Nitrite_reductase"
"RXN0-6369", # ETC, nitrate reduction with Ubiquinones (should be menaquinones), Nitrate_reductase: #(RXN0-3501) in ETC.spy
"RXN0-3501", # ETC, nitrate reduction with menaquinones, replaced in ETC.spy by "Nitrate_mena"
#"RXN0-6377", # ETC, nitrite reduction assim, replaced in ETC.spy by "Nitrite_red_nadh/nadph"
"RXN-12352", # ETC, this looks like a condensed reac (ubiquinones reac with O2)
#"L-LACTATE-DEHYDROGENASE-RXN", #Replaced by "DLACTDEHYDROGNAD-RXN" in model, since L-LACTATE: D-LACTATE
"SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN", #Replaced in Extras.spy "SUCCINATE-DEHYDROGENASE-MENAQUINONE-RXN" , since "Ubiquinones": "Menaquinones"
"R601-RXN", # replaced by "SUCCINATE-DEHYDROGENASE-MENAQUINONE-RXN" in Extras.spy
#"RXN-14971", # replaced by "SUCCINATE-DEHYDROGENASE-MENAQUINONE-RXN" in Extras.spy
"RXN-13720", # replaced by "PGLUCISOM-RXN" in Extras.spy. RXN-13720 no longer in MetaCyc; PGLUCISOM-RXN is the official ID now.
"RXN-11039", # COE, Only described in acetic acid bact (G-s) and in opposite direction (ETOH to ACETALD)
"RXN66-521", # same as "GLUCOKIN-RXN" after  "Glucopyranose" : "GLC"
"2.8.1.6-RXN", # replaced by in Extras.spy (|Sulfurated-Sulfur-Acceptors| removed)
"RXN-8769", # inconsist subset: lump reac (2 subreacs are already in model)
"RXN-13403", # inconsist subset: lump reac (2 subreacs are already in model)
"RXN-7614", # inconsist subset: not in any pathway in dbs, similar reac same EC in model 
"RXN-14047", # lump r of subreacs: "ACONITATEDEHYDR-RXN" and "ACONITATEHYDR-RXN" in model
"RXN-12994", # inco subset: 2 sub-reacs only present in lower eukaryotes and mammals.
"RXN-6182", #Its the same as PGLUCISOM-RXN
"RXN-12611", # r from Ecoli; we already have RXN-12610, which has been reported in S.aureus
"ARGINASE-RXN", # r with no gene associated in dbs (PRO real auxotrophy)
"ORNITHINE-RACEMASE-RXN", # no gene associated in dbs (aa auxotrophies)
"RXN-14903", # we have RXN0-7008 instead
"NMNNUCLEOSID-RXN", # EC-3.2.2.14, niacine met (no genes in RP62A)
"NMNAMIDOHYDRO-RXN", # EC-3.5.1.42, niacine met (no genes in RP62A)
"QUINOPRIBOTRANS-RXN", # EC-2.4.2.19, niacine met (no genes in RP62A)
"QUINOLINATE-SYNTHA-RXN",# EC-2.5.1.72, niacine met (no genes in RP62A),
"BENZALDEHYDE-DEHYDROGENASE-NADP+-RXN", # inco subset, only reacs in model involved with the two compounds, IRREV L-R EC-1.2.1.7, OK
"BENZALDEHYDE-DEHYDROGENASE-NAD+-RXN", # inco subset (as above) EC-1.2.1.28, REV, BRENDA IRREV, L-R, OK
"COA-DISULFIDE-REDUCTASE-NADH-RXN-(NAD)", # inconsist subset:(enz strong preference NADPH in S. aureus), OK
"24-DIAMINOPENTANOATE-DEHYDROGENASE-RXN", # inco subset (NADH more imp in pathway)
"GLUTAMATE-SYNTHASE-FERREDOXIN-RXN", # EC-1.4.7.1 not found in KEGG or BRENDA for RP62A
"MALTODEXGLUCOSID-RXN", # glycogen degradation 1 glycogen -> 1 GLC, model doesnt uptake -> in Corr_Rev
"PHOSPHOKETOLASE-RXN", # No present on Staph.No act in Sepi.Gap-filled by PathwayTools.REF:Sivakanesan1980.StrstersAndWinkler1963.
"ISOCITRATE-DEHYDROGENASE-NAD+-RXN", # ISOCITDEH-RXN in instead (NADPH) in TCA prokariotes. Removed from vs Db 23.0, was present before.
"RXN-13697", # same reac, EC and gene as RXN-13697
"ACETATE--COA-LIGASE-ADP-FORMING-RXN", # no genes, "ACETATE--COA-LIGASE-RXN" in model with genes 
"GLUTAMATE-SYNTHASE-NADH-RXN" # no genes, NADPH reac in model with genes "GLUTAMATESYN-RXN"
##"MALYL-COA-LYASE-RXN", # no gene assoc; leave in for consistency since added to Sepi_846 to allow GLT and AAs BS!
##"MALATE--COA-LIGASE-RXN" # wrong gene association in db, gene for Suc-CoA ligase (6.2.1.5);  leave in for consistency since added to Sepi_846 to allow GLT and AAs BS!
#"RXN0-6377", # EC-1.7.1.4, most orgs in BRENDA pref for NADH. REFs Scarnosus with NADH
#"PYRROLINECARBREDUCT-RXN", # EC preference for NADPH
#"BUTANAL-DEHYDROGENASE-RXN", # inco subset: stronger pref for NADPH (Stp8... notes)
#"GLUCOSE-1-DEHYDROGENASE-RXN", # EC similar preference NAD/NADP but this reac in PPP its useful to generate NADPH
#"GLUTDEHYD-RXN", # EC-1.4.1.4 , same gene as 1.4.1.2, whihc is for a NAD-dependent enz,so enz is bifunctionsl OK
#"GLUTAMATE-SYNTHASE-NADH-RXN", # EC-1.4.1.14, same gene as 1.4.1.13, whihc is for a NADP-dependent enz, so enz is bifunctionsl OK
#"GLUTAMATE-DEHYDROGENASE-NADP+-RXN", # EC-1.4.1.3 , no genes dbs (same r as GLUTDEHYD-RXN),so enz is bifunctionsl OK
##"RXN-14189", # dead, involved with "HYDROGEN-PEROXIDE" (unconserved)
##"RXN-8635", # dead, involved with "HYDROGEN-PEROXIDE" (unconserved)
##"L-ASPARTATE-OXID-RXN", # dead, involved with "HYDROGEN-PEROXIDE" (unconserved)
######"RXN-8673", # dead, involved with "HYDROGEN-PEROXIDE" (unconserved), only reac involved with CPD-8890 in MetaCyc.
######"RXN-8674", # dead, involved with "HYDROGEN-PEROXIDE" (unconserved), only reac involved with CPD-10490 in MetaCyc.
######"RXN-14240", # dead, involved with "HYDROGEN-PEROXIDE"(unconserved), only reac with CPD-15172 in MetaCyc
######"FHLMULTI-RXN", # Unconserved met = 'HYDROGEN-MOLECULE' now in Unwanted ; no needed (system can export formate)
]

Mets = [
#'CPD-8890', # betanidin quinone : only involved in 1 reac in MetaCyc and its dead (RXN-8635)
#'CPD-10490', # N-ethylglycine : only involved in 1 reac in MetaCyc and its dead (RXN-8674)
#'CPD-12377', # All reacs are dead, involved with O and H (hydroxyl radical OH-), formula OH
#'HYDROGEN-MOLECULE', # Unconserved met (StoiCons), all reacs are dead
#'OH', # Unconserved met (StoiCons), all reacs are dead, formula OH
#'SUPER-OXIDE', # Unconserved met (StoiCons), all reacs are dead, OO-, formula O2
'1-2-Diglycerides',
'2-3-4-Saturated-L-Phosphatidates',
'2-Acylglycero-Phosphocholines',
'2-Hexadecenoyl-ACPs',
'2-Me-Branched-234-Sat-FA',
'2-Me-Branched-234-Sat-FALD',
'2-Me-Branched-234-Sat-Fatty-Acyl-CoA',
'2-Octenoyl-ACPs',
'2-Oxo-carboxylates',
'23S-rRNA-pseudouridine1911-1915-1917',
'23S-rRNA-pseudouridine2457',
'23S-rRNA-pseudouridine2604',
'23S-rRNA-pseudouridine2605',
'23S-rRNA-pseudouridine955-2504-2580',
'23S-rRNA-uridine1911-1915-1917',
'23S-rRNA-uridine2457',
'23S-rRNA-uridine2604',
'23S-rRNA-uridine2605',
'23S-rRNA-uridine955-2504-2580',
'2Fe-2S-proteins',
'3-HYDROXY-DOCOSAPENTAENOYL-ACP',
'3-Hydroxy-octanoyl-ACPs',
'3-Hydroxyglutaryl-ACP-methyl-ester',
'3-KETOACYL-COA',
'3-Ketoglutaryl-ACP-methyl-ester',
'3-Ketopimeloyl-ACP-methyl-esters',
'3-Methyl-Saturated-Fatty-Acids',
'3-Methyl-Saturated-Fatty-Acyl-CoA',
'3-OXO-EICOSAPENTAENOYL-ACP',
'3-Oxo-octanoyl-ACPs',
'3-Prime-Ribonucleoside-Monophosphates',
'3-hydroxy-cis-D7-tetraecenoyl-ACPs',
'3-hydroxy-cis-D9-hexaecenoyl-ACPs',
'3-hydroxypimeloyl-ACP-methyl-esters',
'3-oxo-arachidoyl-ACPs',
'3-oxo-behenoyl-ACPs',
'3-oxo-cerotoyl-ACPs',
'3-oxo-cis-D7-tetradecenoyl-ACPs',
'3-oxo-cis-D9-hexadecenoyl-ACPs',
'3-oxo-cis-vaccenoyl-ACPs',
'3-oxo-decanoyl-ACPs',
'3-oxo-dodecanoyl-ACPs',
'3-oxo-hexanoyl-ACPs',
'3-oxo-lignoceroyl-ACPs',
'3-oxo-myristoyl-ACPs',
'3-oxo-palmitoyl-ACPs',
'3-oxo-petroselinoyl-ACPs',
'3-oxo-stearoyl-ACPs',
'3-phosphooligonucleotides',
'3-terminal-unsaturated-sugars',
'5-L-GLUTAMYL-AMINO-ACID',
'5-L-GLUTAMYL-L-AMINO-ACID',
'5-L-GLUTAMYL-PEPTIDE',
'5-Methylcytosine-DNA',
'5-OXOPROLYL-PEPTIDE',
'5-Phosphopolynucleotides',
'5-phosphooligonucleotides',
'6-Dimethylallyladenosine37-tRNAs',
'6-Phospho-b-D-galactosides',
'A-LIPID-HYDROPEROXIDE',
'ACETYL-ACP',
'ACP',
'ACYL-ACP',
'ACYL-SN-GLYCEROL-3P',
'ADP-D-ribosyl-nitrogen-reductases',
'ALA-tRNAs',
'AMINOMETHYLDIHYDROLIPOYL-GCVH',
'ARSENATE',
'ARG-tRNAs',
'ASN-tRNAs',
'ASP-tRNAs',
'Acetoacetyl-ACPs',
'Acids',
'Actinorhodin-Intermediate-1',
'Actinorhodin-Intermediate-2',
'Actinorhodin-Intermediate-4',
'Acyl-Phosphates',
'Adenylated-ThiS-Proteins',
'Alcohols',
'Aldehydes',
'Aliphatic-Alpha-Omega-Diamines',
'Aliphatic-N-Acetyl-Diamines',
'Alkyl-Hydro-Peroxides',
'All-tRNAs',
'All-trans-Retinyl-Esters',
'Amino-Acids',
'Amino-Acids-20',
'Apo-FeS-cluster-proteins',
'Apo-GcvH',
'Arachidoyl-ACPs',
'Aryl-Alcohol',
'Aryl-Sulfate',
'B-KETOACYL-ACP',
'BCAA-dehydrogenase-2MB-DH-lipoyl', # LEU / ILE degrad
'BCAA-dehydrogenase-2MP-DH-lipoyl',  # LEU / ILE degrad
'BCAA-dehydrogenase-DH-lipoyl',   # LEU / ILE degrad
'BCAA-dehydrogenase-lipoyl',  # LEU / ILE degrad
'BCCP-biotin-monomers',
'BCCP-dimers',
'BCCP-monomers',
'Behenoyl-ACPs',
'Beta-3-hydroxybutyryl-ACPs',
'Beta-D-Fructofuranosides',
'Beta-Lactams',
'Beta-hydroxydecanoyl-ACPs',
'Branched-Amino-Acids',
'Butanoyl-ACPs',
'CD-2S-SP-Complex',
'CD-S-SP-Complex',
'CD-SP-2Fe2S-Complex',
'CDP-2-3-4-Saturated-Diacylglycerols',
'CDPDIACYLGLYCEROL',
'CHITIN',
'CHITOBIOSE',
'CIS-DELTA3-ENOYL-COA',
'CPD-70',
'CPD-1084',
'CPD-12443',
'CPD-15317',
'CPD-15709',
'CPD-15972',
'CPD-15973',
'CPD-16758',
'CPD-409',
'CPD-504',
'CPD-8180', 
'CPD-8534' ,
'CPD-8550',
'CPD-8624',
'CPD-8625',
'CPD-9569',
'CPD-9998',
'CPD-11592',
'CPD-15360', 
'CPD0-2350',
'CPD0-2351',
'CPD0-2352',
'CPD0-2353',
'CPD0-2354',
'CPD0-971',
'CPD66-39',
'CYS-tRNAs',
'Carboxyadenylated-MPT-synthases',
'Carboxybiotin-BCCP',
'Carboxylates',
'Carboxylic-esters',
'Cations',
'Cerotoyl-ACPs',
'Chap-ADP-apo-SP-Complex',
'Chap-ATP-Co-chaperone-SP-2Fe2S-Complex',
'Charged-ALA-tRNAs',
'Charged-ARG-tRNAs',
'Charged-ASN-tRNAs',
'Charged-ASP-tRNAs',
'Charged-CYS-tRNAs',
'Charged-GLN-tRNAs',
'Charged-GLT-tRNAs',
'Charged-GLY-tRNAs', # Pep_synthesis Cell wall
'Charged-HIS-tRNAs',
'Charged-ILE-tRNAs',
'Charged-LEU-tRNAs',
'Charged-LYS-tRNAs',
'Charged-MET-tRNAs',
'Charged-PHE-tRNAs',
'Charged-PRO-tRNAs',
'Charged-SER-tRNAs',
'Charged-THR-tRNAs',
'Charged-TRP-tRNAs',
'Charged-TYR-tRNAs',
'Charged-VAL-tRNAs',
'Chitodextrins',
'Cis-Delta5-dodecenoyl-ACPs',
'Cis-Delta7-tetradecenoyl-ACPs',
'Cis-delta-3-decenoyl-ACPs',
'Cleaved-Repressor-LexA',
'Co-chaperone-SP-2Fe2S-Complex',
'Corrinoid-Fe-S-proteins',
'Crotonyl-ACPs',
'D-3-HYDROXYACYL-COA', 
'D-form-FeS-Cluster-Scaffold-Proteins',
'DEOXYNUCLEOTIDESM',
'DIACYLGLYCEROL',
'DIACYLGLYCEROL-PYROPHOSPHATE',
'DIHYDROLIPOYL-GCVH', # LEU / ILE degradation
'DIPEPTIDES',
'DNA-6-O-Methyl-Guanines',
'DNA-Adjacent-Pyrimidines',
'DNA-Combined-With-Exogenous-DNA',
'DNA-Cytosines',
'DNA-Guanines',
'DNA-Holder',
'DNA-N',
'DNA-Segment-Placeholder',
'DNA-Segment-in-Reverse-Orientations',
'DNA-With-Pyrimidine-Dimers',
'DNA-containing-a-Apyrimidinic-Sites',
'DNA-containing-aPurinic-Sites',
'DNA-containing-abasic-Sites',
'DNA-containing-diamino-hydro-formamidops',
'DNA-with-Alkylated-Adenines',
'DNA-with-Uracils',
'Damaged-DNA-Pyrimidine',
'Decanoyl-ACPs',
'Demethylated-Ubiquinols',
'Demethylated-methyl-acceptors',
'Demethylmenaquinols',
'Deoxy-Ribonucleoside-Diphosphates',
'Deoxy-Ribonucleoside-Monophosphates',
'Deoxy-Ribonucleoside-Triphosphates',
'Deoxy-Ribonucleosides',
'Deoxynucleotides',
'DIACYLGLYCEROL',
'Diacylglycerol-Prolipoproteins',
'Dietary-retinyl-esters',
'Dihydro-Lipoyl-Proteins', # LEU / ILE degradation
'Dipeptides-With-Asp-At-N-Terminal',
'Dipeptides-With-Proline-Carboxy',
'Dodec-2-enoyl-ACPs',
'Dodecanoyl-ACPs',
'Donor-H2',
'Double-Stranded-DNAs',
'E-', # COE
'FORMYL-L-METHIONYL-PEPTIDE',
'FORMYL-THF-GLU-N',
'FRU',
'Fatty-Acids',
'Fatty-Aldehydes',
'FeS-Cluster-Chaperones-ATP',
'FeS-Cluster-Co-Chaperones',
'Flavodoxins-Semiquinones',
'Formate-acetyltrans-glycine-radical',
'Formate-acetyltransferase-glycine',
'GLN-tRNAs',
'GLT-tRNAs',
'GLY-tRNAs', # Pep_synthesis Cell wall
'Gcv-H',
'General-Protein-Substrates',
'Glutaryl-ACP-methyl-esters',
'Glycerophosphodiesters',
'Guanine34-in-tRNAs',
'Guanine37-in-tRNA',
'Guanine46-in-tRNA',
'Guanine9-in-tRNA',
'HEME_O',
'HIS-tRNAs',
'Hex-2-enoyl-ACPs',
'Hexanoyl-ACPs',
'ILE-tRNAs',
'IS30-Insertion-Sequences',
'IS30-with-Integrated-Transposon',
'Jasmonic-Acids',
'Kanamycin-3-phosphates',
'Kanamycins',
'L-1-PHOSPHATIDYL-ETHANOLAMINE',
'L-1-PHOSPHATIDYL-GLYCEROL',
'L-1-PHOSPHATIDYL-GLYCEROL-P',
'L-3-HYDROXYACYL-COA', 
'L-Cysteine-Desulfurases',
'L-Glutamyl-Peptides',
'L-PHOSPHATIDATE',
'L-glutamyl-tRNAGln',
'L-methionyl-tRNAfmet',
'L-seryl-SEC-tRNAs',
'LEU-tRNAs',
'LONG-CHAIN-KETONE',
'LYS-tRNAs',
'Leader-Sequences',
'Light',
'Lignoceroyl-ACPs',
'Lipids',
'Lipoproteins',
'Lipoyl-Protein',
'Lipoyl-Protein-N6-lipoyllysine', # LEU / ILE degradadation
'Long-Chain-3S-Hydroxyacyl-CoAs',
'Long-Chain-Acyl-CoAs',
'Long-Chain-Fatty-Acids',
'Long-Chain-Polyphosphate',
'Long-Chain-Trans-23-Dehydroacyl-CoA',
'Long-Chain-oxoacyl-CoAs',
'LysW-C-Terminal-L-Glutamate',
'LysW-L-glutamate',
'LysW-L-glutamate-5-phosphate',
'LysW-L-glutamate-5-semialdehyde',
'LysW-L-ornithine',
'MALONYL-ACP',
'MET-tRNAs',
'METHIONYL-PEPTIDE',
'METHYLENE-THF-GLU-N',
'MPT-Synthases',
'Malonyl-acp-methyl-ester',
'Methyl-Jasmonates',
'Methylated-Ribosomal-Protein-L11s',
'Methylated-corrinoid-Fe-S-Proteins',
'Methylated-methyl-acceptors',
'Myristoyl-ACPs',
'N-Substituted-Amino-Acids',
'N-Substituted-Aminoacyl-tRNA',
'N-formyl-L-methionyl-tRNAfmet',
'NAcMur-4Peptide-NAcGlc-Undecaprenols',
'NAcMur-Peptide-NAcGlc-Undecaprenols',
'NAcMur-Peptide-Undecaprenols',
'Negatively-super-coiled-DNAs',
'Nitrogen-reductases',
'Nonmethylated-Ribosomal-Protein-L11s',
'Nucleoside-Diphosphates',
'Nucleoside-Monophosphates',
'Nucleoside-Triphosphates',
'OCTANOYL-ACP',
'OH-ACYL-ACP',
'OLIGOPEPTIDES',
'Octadec-2-enoyl-ACPs',
'Octanoyl-ACPs',
'Octanoylated-Gcv-H',
'Octanoylated-domains',
'Odd-Saturated-Fatty-Acyl-CoA',
'Odd-Straight-Chain-234-Sat-FA',
'Odd-Straight-Chain-234-Sat-FALD',
'Oleoyl-ACPs',
'Orthophosphoric-Monoesters',
'Oxo-glutarate-dehydro-suc-DH-lipoyl',
'Oxo-glutarate-dehydrogenase-DH-lipoyl', # LEU / ILE degrad
'Oxo-glutarate-dehydrogenase-lipoyl', # LEU / ILE degrad 
'PHE-tRNAs',
'PHOSPHATIDYLCHOLINE',
'POLYPEPTIDE',
'POLYRIBITOL-PHOSPHATE',
'PRO-tRNAs',
'PROT-CYS',
'PROTEIN-LIPOYLLYSINE', # LEU / ILE degrad
'PTS-I-Histidines',
'PTS-I-pi-phospho-L-histidines',
'Palmitoleoyl-ACPs',
'Palmitoyl-ACPs',
'Peptides-holder',
'Peptides-with-Leader-Sequence',
'Peptidoglycan-with-DAP-pentapeptide',
'Peptidoglycan-with-L-lysine-pentapeptide', # Pep_synthesis Cell wall
'Peptidoglycans',
'Persulfurated-L-cysteine-desulfurases',
'Petrosel-2-enoyl-ACPs',
'Petroselinoyl-ACPs',
'Phenolic-Donors',
'Phenoxyl-rad-of-phenolic-donors',
'Polyketide-ACP-Proteins',
'Polynucleotide-Holder',
'Primary-Alcohols',
'Prolipoproteins',
'Protein-Histidines',
'Protein-L-Methionine-S-Oxides',
'Protein-L-methionine',
'Protein-L-serine-or-L-threonine',
'Protein-L-serines',
'Protein-Phosphoserines',
'Protein-S-methyl-L-cysteine',
'Protein-Ser-or-Thr-phosphate',
'Protein-Tyrosines',
'Protein-With-N-Terminal-Met',
'Protein-pi-phospho-L-histidines',
'Protein-tyrosine-phosphates',
'Purine-Bases',
'Purine-Deoxyribonucleosides',
'Purine-Ribonucleosides',
'Pyruvate-dehydrogenase-acetylDHlipoyl',
'Pyruvate-dehydrogenase-dihydrolipoate', # LEU / ILE degrad
'Pyruvate-dehydrogenase-lipoate',  # LEU / ILE degrad
'R-3-Hydroxypalmitoyl-ACPs',
'R-3-hydroxy-cis-vaccenoyl-ACPs',
'R-3-hydroxyarachidoyl-ACPs',
'R-3-hydroxybehenoyl-ACPs',
'R-3-hydroxycerotoyl-ACPs',
'R-3-hydroxydodecanoyl-ACPs',
'R-3-hydroxyhexanoyl-ACPs',
'R-3-hydroxylignoceroyl-ACPs',
'R-3-hydroxymyristoyl-ACPs',
'R-3-hydroxypetroselinoyl-ACPs',
'R-3-hydroxystearoyl-ACPs',
'R2-2OH-Straight-Chain-234-Sat-FA-CoA',
'R2OH-Straight-Chain-234-Sat-FA',
'RNA-DNA-hybrids',
'RNA-Holder',
'RNASE-III-MRNA-PROCESSING-SUBSTRATE',
'RNASE-III-PROCESSING-PRODUCT-MRNA',
'Relaxed-DNAs',
'Resolution-of-Recombinational-Junction',
'Retinols',
'Ribonuc-tri-P-reductases-active',
'Ribonuc-tri-P-reductases-inactive',
'Ribonucleoside-Diphosphates',
'Ribonucleoside-Monophosphates',
'Ribonucleoside-Triphosphates',
'Ribonucleosides',
'Ribosomal-protein-L3-L-glutamine',
'Ribosomal-protein-L3-N5M-L-glutamine',
'Ribulose-phosphates',
'Ribuloses',
'S-CD-Apo-SP-Complex',
'S-CD-S-SP-Complex',
'SEC-tRNAs',
'SER-tRNAs',
'SS-Oligodeoxyribonucleotides',
'Saturated-2-Lysophosphatidates',
'Saturated-Fatty-Acyl-ACPs',
'Saturated-Fatty-Acyl-CoA',
'Secondary-Alcohols',
'Short-RNA-Fragments',
'Single-Stranded-DNAs',
'Stearoyl-ACPs',
'Sulfur-Carrier-Proteins-ThiI',
'Sulfurated-Sulfur-Acceptors',
'Sulfurylated-ThiI',
'Supercoiled-Duplex-DNAs',
'THF-GLU-N',
'THR-tRNAs',
'TRANS-D2-ENOYL-ACP',
'TRANS-D2-ENOYL-COA',
'TRIPEPTIDES',
'TRP-tRNAs',
'TYR-tRNAs',
'Teichoic-P-Gro',
'Teichoic-P-Gro-Glc',
'Teichoic-peptidoglycan',
'Tetradec-2-enoyl-ACPs',
'Thi-S',
'Thiocarboxyadenylated-ThiS-Proteins',
'Thiocarboxylated-MPT-synthases',
'Trans-D2-decenoyl-ACPs',
'Trans-D2-hexacos-2-enoyl-ACPs',
'Trans-D3-cis-D5-dodecenoyl-ACPs',
'Trans-D3-cis-D7-tetradecenoyl-ACPs',
'Trans-D3-cis-D9-hexadecenoyl-ACPs',
'Triacylglycerides',
'Unsulfurated-Sulfur-Acceptors',
'VAL-tRNAs',
'a-2-oxoglutarate-dehydrogenase-E2-protei',
'a-pyruvate-dehydrogenase-E2-protein-Nsup',
'a-radical-of-luteolin-7-iOi-diglucuronid',
'apo-ACP',
'b-Hydroxy-cis-D5-dodecenoyl-ACPs',
'b-Keto-cis-D5-dodecenoyl-ACPs',
'cis-cis-D11-23-3-hydroxyC42-2-ACPs',
'cis-cis-D11-23-3-oxo-C42-2-ACPs',
'cis-cis-D11-29-3-hydroxyC48-2-ACPs',
'cis-cis-D11-29-3-oxo-C48-2-ACPs',
'cis-cis-D13-25-3-hydroxyC44-2-ACPs',
'cis-cis-D13-25-3-oxo-C44-2-ACPs',
'cis-cis-D13-31-3-hydroxyC50-2-ACPs',
'cis-cis-D13-31-3-oxo-C50-2-ACPs',
'cis-cis-D15-27-3-hydroxyC46-2-ACPs',
'cis-cis-D15-27-3-oxo-C46-2-ACPs',
'cis-cis-D15-33-3-hydroxyC52-2-ACPs',
'cis-cis-D15-33-3-oxo-C52-2-ACPs',
'cis-cis-D17-29-3-hydroxyC48-2-ACPs',
'cis-cis-D17-29-3-oxo-C48-2-ACPs',
'cis-cis-D17-35-3-hydroxyC54-2-ACPs',
'cis-cis-D17-35-3-oxo-C54-2-ACPs',
'cis-cis-D19-31-3-hydroxyC50-2-ACPs',
'cis-cis-D19-31-3-oxo-C50-2-ACPs',
'cis-cis-D19-37-3-hydroxyC56-2-ACPs',
'cis-cis-D19-37-3-oxo-C56-2-ACPs',
'cis-cis-D21-39-3-hydroxyC58-2-ACPs',
'cis-cis-D21-39-3-oxo-C58-2-ACPs',
'cis-cis-D5-17-3-hydroxyC36-2-ACPs',
'cis-cis-D5-17-3-oxo-C36-2-ACPs',
'cis-cis-D5-23-3-hydroxyC42-2-ACPs',
'cis-cis-D5-23-3-oxo-C42-2-ACPs',
'cis-cis-D7-19-3-hydroxyC38-2-ACPs',
'cis-cis-D7-19-3-oxo-C38-2-ACPs',
'cis-cis-D7-25-3-hydroxyC44-2-ACPs',
'cis-cis-D7-25-3-oxo-C44-2-ACPs',
'cis-cis-D9-21-3-hydroxyC40-2-ACPs',
'cis-cis-D9-21-3-oxo-C40-2-ACPs',
'cis-cis-D9-27-3-hydroxyC46-2-ACPs',
'cis-cis-D9-27-3-oxo-C46-2-ACPs',
'cis-delta11-3-hydroxymelissoyl-ACPs',
'cis-delta11-3-oxo-melissoyl-ACPs',
'cis-delta11-melissoyl-ACPs',
'cis-delta13-3-hydroxylacceroyl-ACPs',
'cis-delta13-3-oxo-lacceroyl-ACPs',
'cis-delta13-lacceroyl-ACPs',
'cis-delta15-3-hydroxygheddoyl-ACPs',
'cis-delta15-3-oxo-gheddoyl-ACPs',
'cis-delta15-gheddoyl-ACPs',
'cis-delta17-3-hydroxyC36-ACPs',
'cis-delta17-3-oxo-C36-ACPs',
'cis-delta17-C36-ACPs',
'cis-delta19-3-hydroxyC38-ACPs',
'cis-delta19-3-oxo-C38-ACPs',
'cis-delta19-C38-ACPs',
'cis-delta21-3-hydroxyC40-ACPs',
'cis-delta21-3-oxo-C40-ACPs',
'cis-delta3-behenoyl-ACPs',
'cis-delta5-3-hydroxylignoceroyl-ACPs',
'cis-delta5-3-oxo-lignoceroyl-ACPs',
'cis-delta5-lignoceroyl-ACPs',
'cis-delta7-3-hydroxycerotoyl-ACPs',
'cis-delta7-3-oxo-cerotoyl-ACPs',
'cis-delta7-cerotoyl-ACPs',
'cis-delta9-3-hydroxymontanoyl-ACPs',
'cis-delta9-3-oxo-montanoyl-ACPs',
'cis-delta9-montanoyl-ACPs',
'cis-vaccen-2-enoyl-ACPs',
'cleaved-lipoprotein-signal-peptide',
'repressor-LexA',
'ssRNAs',
'tRNA-Adenosines-37',
'tRNA-Containing-5AminoMe-2-ThioUrdines',
'tRNA-Containing-5MeAminoMe-2-ThioU',
'tRNA-Containing-N1-Methylguanine-37',
'tRNA-Containing-N1-Methylguanine-9',
'tRNA-Containing-N7-Methylguanine-46',
'tRNA-precursors',
'tRNA-pseudouridine-38-40',
'tRNA-pseudouridine32',
'tRNA-uridine-38-40',
'tRNA-uridine32',
'tRNA-with-7-aminomethyl-7-deazaguanine',
'tRNAs',
'tRNAs-containing-epoxy-quenosine',
'tRNAs-with-CCA',
'tRNAs-with-queuine',
'trans-delta2-arachidoyl-ACPs',
'trans-delta2-behenoyl-ACPs',
'trans-delta2-lignoceroyl-ACPs',
]

def BadMetMatches(met):
    return False

def BadReacMatches(met):
    return False
