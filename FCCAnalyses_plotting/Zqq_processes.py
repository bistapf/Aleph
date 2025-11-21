from collections import namedtuple
import ROOT

# Apranik's colors:
# colors = [
#     "#c7e9f1", "#9ed9de", "#72c7c7",
#     "#4bb0ae", "#2a8b9b"
# ]


processSpecs = namedtuple('processSpecs', ['sample_list', 'colour_key', 'title'])

########################################### MC PROCESSES ########################################### 

# merging uds jets into one process called light jets
MC_group_light_jets ={
	"Zbb":processSpecs(sample_list=["Zbb"], colour_key=ROOT.TColor.GetColor("#3fa9a9"), title="b-jets" ), 
	"Zcc":processSpecs(sample_list=["Zcc"], colour_key=ROOT.TColor.GetColor("#80d1d1"), title="c-jets" ), 
	"Zudsuds":processSpecs(sample_list=["Zss", "Zuu", "Zdd"], colour_key=ROOT.TColor.GetColor("#1f77b4"), title="light jets" )
}

# TODO: add spilt Zuds

########################################### DATA ########################################### 

zqq_data = {
	"data":processSpecs(sample_list=["Znn"], colour_key=ROOT.kBlack, title="data" ),
}

