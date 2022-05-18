



MetNames = {
    "L-arabinopyranose"        : "ARABINOSE",
    "CPD-10353"                : "ACETOIN",
    "CPD-13469"                : "D-GLUCOSAMINE-6-P",
    "N-acetyl-D-glucosamine"   : "N-ACETYL-D-GLUCOSAMINE",
    "Glucose"                  : "GLC",
    "Glucopyranose"            : "GLC",
    "ALPHA-GLUCOSE"            : "GLC",
    "ALPHA-GLC-6-P"            : "GLC-6-P",
    "D-glucopyranose-6-phosphate": "GLC-6-P",
    "MALTOSE"                  : "ALPHA-MALTOSE",                  
    "NAD-P-OR-NOP"             : "NAD(P)",
    "NADH-P-OR-NOP"            : "NAD(P)H",
    "CPD0-1108"                : "RIBOSE",
    "CPD-10330"                : "RIBOSE",
    "CPD0-1110"                : "RIBOSE",
    "D-Ribofuranose"           : "RIBOSE", 
    "D-Ribopyranose"	       : "RIBOSE",
    "D-Xylopyranose"           : "XYLOSE",
    "25-DIAMINOPENTANOATE"     : "L-ORNITHINE",
    "Menaquinones"             : "MENAQUINONE",
    "Menaquinols"              : "REDUCED-MENAQUINONE",
    "Ubiquinones"              : "MENAQUINONE",
    "Ubiquinols"               : "REDUCED-MENAQUINONE",
    "ETR-Quinones"             : "MENAQUINONE", 
    "ETR-Quinols"              : "REDUCED-MENAQUINONE",
    "Reduced-Quinones"         : "REDUCED-MENAQUINONE",
    "Quinones"                 : "MENAQUINONE",
    "Rhodoquinols"             : "REDUCED-MENAQUINONE",
    "Rhodoquinones"            : "MENAQUINONE",
    "Reduced-cytochromes"      : "Cytochromes-C-Reduced", 
    "Oxidized-cytochromes"     : "Cytochromes-C-Oxidized",
    }
##    "Mannose-6-phosphate"      : "MANNOSE-6P"
##    "ETF-Reduced"              : "REDUCED-MENAQUINONE",
##    "ETF-Oxidized"             : "MENAQUINONE", 

##  Not substituted here, but in Corr_stoich becasue also needed to adjust protons =
##  "ETF-Reduced"  = "ETF-Oxidized" + 3H
##  while "REDUCED-MENAQUINONE" = "MENAQUINONE" + 2H

##Either all "ETF-Oxidized" to FAD for continuity of the network, or for MENAQUINONES.
##
##Example of problem:
##
##    MEPROPCOA-FAD-RXN: # EC-1.3.8- ; valine degradation
##            1 ISOBUTYRYL-COA + 1 MENAQUINONE <> 1 REDUCED-MENAQUINONE + 1 METHACRYLYL-COA
##            ~
##            
##MetaCyc 22.6 = FAD had to be susbtituted by MENAQUINONE etc. Since not any other reactions
##in model were involved with FAD and so, it could not be regenerated. In Corr_stoich, cause
## also needed to adjust protons. (FADH2 = FAD + 2H)
##
##The electron transferring-flavoproteins (ETF) serve as specific electron carriers for
##other dehydrogenases. They transfer the electrons to the main respiratory chain via
##ETF-ubiquinone oxidoreductase (ETF dehydrogenase).
##
##Most of the characterized electron-transferring flavoproteins belong to a single group.
##These proteins, including those from mammals [Roberts96], plants [Ishizaki06], and some
##bacteria, such as the methylotrophic bacterium Methylophilus methylotrophus W3A1 [Chen94]
##and the soil bacterium Paracoccus denitrificans [Roberts99], are heterodimeric, contain
##one non-covalently bound FAD and one AMP per molecule, and function as electron carriers
##between other flavoproteins (such as acyl-CoA dehydrogenases, sarcosine dehydrogenase,
##and trimethylamine dehydrogenase) and a quinone, depending on the source species.

##UNIQUE-ID - ETF-Oxidized
##CHEMICAL-FORMULA - (ETF 1)
##//
##UNIQUE-ID - ETF-Reduced
##CHEMICAL-FORMULA - (ETF 1)
##CHEMICAL-FORMULA - (H 3)
##//
##UNIQUE-ID - ETR-Quinols
##CHEMICAL-FORMULA - (ETR-Q )
##CHEMICAL-FORMULA - (H 2)
##//
##UNIQUE-ID - ETR-Quinones
##CHEMICAL-FORMULA - (ETR-Q)

       
