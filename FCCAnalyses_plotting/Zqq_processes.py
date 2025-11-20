from collections import namedtuple
import ROOT


processSpecs = namedtuple('processSpecs', ['sample_list', 'colour_key', 'title'])

########################################### MC PROCESSES ########################################### 

# merging uds jets into one process called light jets
MC_group_light_jets ={
	"Zbb":processSpecs(sample_list=["Zbb"], colour_key=ROOT.TColor.GetColor("#3f7eb3"), title="b-jets" ), 
	"Zcc":processSpecs(sample_list=["Zcc"], colour_key=ROOT.TColor.GetColor("#2E9094"), title="c-jets" ), 
	"Zudsuds":processSpecs(sample_list=["Zss", "Zuu", "Zdd"], colour_key=ROOT.TColor.GetColor("#975b95"), title="light jets" )
}

# TODO: add spilt Zuds

########################################### DATA ########################################### 

zqq_data = {
	"data":processSpecs(sample_list=["Znn"], colour_key=ROOT.kBlack, title="data" ),
}

