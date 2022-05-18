#
##
### Build the S. epi RP62A model
##
#



OrgName = 'Staph_epi_RP62A_20.1'   # name of the database passed to PyoCyc.Organism

from ScrumPy.BuildModel.UsrBuildOrg import Init, BuildModel, RebuildModel, ShowCorrec, HideCorrec

Init(OrgName,NoCompID="Cytosol")

from ScrumPy.BuildModel.UsrBuildOrg import  orgdb as sepidb

