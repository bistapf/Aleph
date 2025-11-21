# Default config for plotting from the stage 1 ntuples, grouping uds = light jets

import Zqq_plots 
import Zqq_processes 

class PlottingConfig:

    # Input and output paths 
    inputs_path = "/eos/user/h/hfatehi/D0fliped-good/"
    outputs_path = "./plots_data_mc_zqq_test/"

    # Dictionary of all variables to plot: info about binning, range, labels etc - edit or add in Zqq_plots.py
    plots_dict = Zqq_plots.Zqq_data_MC_vars

    # Dictionary of processes: info about samples, colours, labels - edit or add in Zqq_processes.py
    data = Zqq_processes.zqq_data["data"]
    mc_processes = Zqq_processes.MC_group_light_jets

    # Meta-data about the data plotted:
    year = "1994"
    sel_tag = "Selected events"
    lumi = 57.89 # in pb-1
    ecm = 91. # in GeV

    # Plotting style settings:
    do_log_y = True
    add_overflow = False
    ratio_range = (0., 2.)
    weighted = False # whether to apply generator level weights or not, currently not needed keep on false

    # Output file settings:
    out_format = ".png"
    store_root_file=False #if true a ROOT file of the plotted histograms ins written that can be used by e.g. combine for fitting