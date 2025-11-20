from collections import namedtuple
import ROOT

#variables to plot:
default_binning = 25 
PlotSpecs = namedtuple('PlotSpecs', ['name', 'xmin', 'xmax', 'label', 'nbins'])

Zqq_data_MC_vars = {

	# event kinematics
	"m_jj":PlotSpecs(name="event_invariant_mass", xmin=50, xmax=120, label="m_{jj} in GeV", nbins=70),

	# all input features of the tagger

	
}
