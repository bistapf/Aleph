from collections import namedtuple
import ROOT

#variables to plot:
default_binning = 25 
PlotSpecs = namedtuple('PlotSpecs', ['name', 'xmin', 'xmax', 'label', 'nbins'])

Zqq_data_MC_vars = {

	# event multiplicities & kinematics
	"n_jets":PlotSpecs(name="event_njet", xmin=0, xmax=5, label="Number of jets", nbins=5),
	"m_jj":PlotSpecs(name="event_invariant_mass", xmin=50, xmax=120, label="m_{jj} in GeV", nbins=70),

	# jet properties
	"jet_mass":PlotSpecs(name="jet_mass", xmin=0, xmax=60, label="Jet mass in GeV", nbins=60),
	"jet_pT":PlotSpecs(name="jet_pT", xmin=0, xmax=60, label="p_{T} jet in GeV", nbins=60),
	"jet_eta":PlotSpecs(name="jet_eta", xmin=-2, xmax=2, label="#eta jet", nbins=30),
	"jet_phi":PlotSpecs(name="jet_phi", xmin=0, xmax=7, label="#phi jet", nbins=30),

	# jet constituents multiplicities
	"jet_nconst":PlotSpecs(name="jet_nconst", xmin=0, xmax=50, label="Number of constituents in jet", nbins=50),
	"jet_nmu":PlotSpecs(name="jet_nmu", xmin=0, xmax=5, label="Number of muons in jet", nbins=5),
	"jet_nel":PlotSpecs(name="jet_nel", xmin=0, xmax=5, label="Number of electrons in jet", nbins=5),
	"jet_ngamma":PlotSpecs(name="jet_ngamma", xmin=0, xmax=20, label="Number of photons in jet", nbins=20),
	"jet_nchad":PlotSpecs(name="jet_nchad", xmin=0, xmax=30, label="Number of charged hadrons in jet", nbins=30),
	"jet_nnhad":PlotSpecs(name="jet_nnhad", xmin=0, xmax=20, label="Number of neutral hadrons in jet", nbins=20),

	# all input features of the tagger

	
}
