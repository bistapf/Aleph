# Script to plot data and MC from ALEPH

import ROOT
import os 
from collections import namedtuple
from array import array
import copy
import numpy
from argparse import ArgumentParser
import json

from plotting_config_stage1 import PlottingConfig

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

def addOverflowToLastBin(hist):
	# print("Adding overflow to last bin for histogram:", hist.GetName())
	overflow_index = hist.GetNbinsX()+1
	if hist.IsBinOverflow(overflow_index):
		# print ("Found overflow bin with index", overflow_index, "merging content with previous bin!")
		overflow = hist.GetBinContent(overflow_index)
		lastbin = hist.GetBinContent(overflow_index-1)
		hist.SetBinContent(overflow_index-1, overflow+lastbin)
		# print ("Added together overflow =", overflow, "and last bin content", lastbin, "in histogram of name", hist.GetName())



def is_valid_rdf(input_filepath, tree_name="events", do_chunks=False):

    #check if all chunks are ok
    if do_chunks:
        files = glob.glob(os.path.join(input_filepath, "chunk*"))
        if not files:
            print(f"No files found in {input_filepath}. Skipping.")
            return False

        for f in files:
            tf = ROOT.TFile.Open(f)
            if not tf or tf.IsZombie():
                print(f"File {f} is not a valid ROOT file. Skipping.")
                continue
            if not tf.Get(tree_name):
                print(f"Tree '{tree_name}' not found in file {f}. Skipping.")
                continue
    
    else:
        if not os.path.exists(input_filepath):
            print("Filepath ", input_filepath, "not found")
            return False 

        tf = ROOT.TFile.Open(input_filepath)

        if not tf or tf.IsZombie():
            print(f"File {input_filepath} is not a valid ROOT file. Skipping.")
            return False 

        if not tf.Get(tree_name):
            print(f"Tree '{tree_name}' not found in file {input_filepath}. Skipping.")
            return False 
    
    return True

def get_rdf(input_filepath):

    # print("Getting rdf from:", input_filepath)
    rdf = None
    if input_filepath.endswith(".root"):
        
        if not is_valid_rdf(input_filepath):
            print("Process ", input_filepath, " is empty. Skipping.")
            return None
        
        rdf = ROOT.RDataFrame("events", input_filepath)

        #catching exception DOES NOT WORK
        # try:
        #     rdf = ROOT.RDataFrame("events", input_filepath)
        # except:
        #     print("File {} appears to be empty.".format(input_filepath))
        #     return
    else:
        if not is_valid_rdf(input_filepath):
            print("Process ", input_filepath, " is empty. Skipping.")
            return None

        rdf = ROOT.RDataFrame("events", input_filepath+"/chunk*")

    if not rdf:
        print("Empty file for:", input_filepath, " Skipping.")
        return None

    # print(rdf.GetColumnNames())

    return rdf

def getNormFactor(sample, proc_dict_path = "/afs/cern.ch/user/b/bistapf/ALEPH_data_project/Aleph/FCCAnalyses_plotting/normalisation.json", 
                  weighted=False, lumi=57.89, print_debug=False):
    norm_factor = -1.

    #read x-section etc from the json file:
    if '.json' in proc_dict_path:
        with open(proc_dict_path, 'r') as dict_file:
            proct_dict = json.load(dict_file)
    else:
        raise Exception("Error in getNormFactor: Only support .json files for the cross-section values currently!")

    if weighted:
        norm_factor = lumi*proct_dict[sample]["crossSection"]*proct_dict[sample]["kfactor"]*proct_dict[sample]["matchingEfficiency"]/proct_dict[sample]["sumOfWeights"]
        print("used sum of weights", proct_dict[sample]["sumOfWeights"])
    else:
        norm_factor = lumi*proct_dict[sample]["crossSection"]*proct_dict[sample]["kfactor"]*proct_dict[sample]["matchingEfficiency"]/proct_dict[sample]["numberOfEvents"]
    
    if print_debug:
        print(f"{sample} - built normfactor from:")
        print("#evts = ", proct_dict[sample]["numberOfEvents"])
        print("x-sec = ", proct_dict[sample]["crossSection"])
        print("k-factor = ", proct_dict[sample]["kfactor"])
        print("match-eff = ", proct_dict[sample]["matchingEfficiency"])
        print("lumi = ", lumi)
        print("norm_factor = ", norm_factor)

    if norm_factor < 0:
        raise Exception("Error in getNormFactor: Probably missing a factor in the normalisation ...")

    return norm_factor


def get_hist_from_tree(proc_name, input_filepath, hist_var, hist_nbins, hist_xmin, hist_xmax, is_data= False):

    #going to need a TH1Model to fill:
    has_variable_binning = False
    if not isinstance(hist_nbins, int):
        has_variable_binning = True
        hist_binEdges = array("d", hist_nbins)
        hist_nBins = len(hist_nbins)-1
        #init the histogram with variable bin widths:
        hist_model = ROOT.RDF.TH1DModel(f"hist_{hist_var}", f"hist_{hist_var}", hist_nBins, hist_binEdges)

    else:
        hist_model = ROOT.RDF.TH1DModel(f"hist_{hist_var}", f"hist_{hist_var}", hist_nbins, hist_xmin, hist_xmax)

    # print(f"Reading from file {input_filepath}")

    rdf = get_rdf(input_filepath)
    if not rdf:
        print("File ", input_filepath, " is empty, error ...")
        exit()

    tmp_hist = rdf.Histo1D(hist_model, hist_var).GetValue() # DO NOT USE WEIGHTS

    #scale MC with xsec and lumi!
    if not is_data:
        norm_factor = getNormFactor(proc_name)
        tmp_hist.Scale(norm_factor)

    print("{} \t {:.2f}".format(proc_name, tmp_hist.Integral()))

    hist_to_add = copy.deepcopy(tmp_hist)
    hist_to_add.SetTitle(proc_name)
    return hist_to_add

# function to draw all the things #todo: add back support for applying extra cut?
def make_plot(plot, input_dir, data_proc, mc_processes, out_dir_base, 
              year="", sel_tag ="", lumi=1., ecm=1.,
              do_log_y=True, add_overflow=False, fix_ratio_range=(), weighted=False,
              out_format = ".png", store_root_file=False,   
              ):

    print("Plotting", plot.name)

    if weighted:
        print("Applying MC events weights")
    
    #get data histogram
    data_hist = None
    for data_file in data_proc.sample_list:
        data_filepath = os.path.join(input_dir, f"{data_file}.root")
        data_tmp_hist = get_hist_from_tree(data_file, data_filepath, plot.name, plot.nbins, plot.xmin, plot.xmax, is_data= True)

        #skip if nothing passes selection:
        if not data_tmp_hist.GetEntries():
            print("Empty histogram for:", data_file, " Skipping.")
            continue

        #move into get_hsito function!
        if add_overflow:
            addOverflowToLastBin(data_tmp_hist)
        
        #set plotting properties:
        data_tmp_hist.SetTitle(data_proc.title)
        data_tmp_hist.SetLineColor(data_proc.colour_key)
        data_tmp_hist.SetMarkerStyle(ROOT.kFullCircle)

        # make summed up histogram of all samples in the process
        if not data_hist: 
            data_hist = copy.deepcopy(data_tmp_hist)
        else:
            data_hist.Add(data_tmp_hist)


    # get MC prediction
    mc_hists_list = []

    for proc in mc_processes:
        #loop over files in our mc process and get each histogram, then add them together
        mc_proc_hist = None
        for sample in mc_processes[proc].sample_list:
            print(f"Looking for MC sample: {sample}")
            filepath = os.path.join(input_dir, f"{sample}.root")
            print(filepath)
            
            tmp_hist = get_hist_from_tree(sample, filepath, plot.name, plot.nbins, plot.xmin, plot.xmax , is_data= False)

            #skip if nothing passes selection:
            if not tmp_hist.GetEntries():
                print("Empty histogram for:", sample, " Skipping.")
                continue

            #move into get_hsito function!
            if add_overflow:
                addOverflowToLastBin(tmp_hist)

            #set plotting properties:
            tmp_hist.SetTitle(mc_processes[proc].title)
            tmp_hist.SetLineColor(mc_processes[proc].colour_key)
            tmp_hist.SetFillColor(mc_processes[proc].colour_key)

            # make summed up histogram of all samples in the process
            if not mc_proc_hist: 
                mc_proc_hist = copy.deepcopy(tmp_hist)
            else:
                mc_proc_hist.Add(tmp_hist)
            
        # add the total histogram for this process to the MC stack 
        mc_hists_list.append(mc_proc_hist)
    

    #sort the MC by size:
    mc_hists_list = sorted(mc_hists_list, key=lambda hist: hist.Integral())
    mc_stack = ROOT.THStack("mc_stack","")
    for mc_hist in mc_hists_list:
        mc_stack.Add(mc_hist)

    #draw stuff:
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800) 
    #add ratio:
    pad_up = ROOT.TPad("pad_up", "pad_up", 0., 0., 1., 1.)
    pad_up.SetFillStyle(0)
    pad_up.SetBottomMargin(0.32)
    pad_up.SetTopMargin(0.03)
    pad_up.SetLeftMargin(0.13)
    pad_up.SetRightMargin(0.05)
    pad_up.SetLogy(do_log_y)
    pad_up.Draw()

    pad_low = ROOT.TPad("pad_low", "pad_low", 0., 0., 1., 1.);
    pad_low.SetFillStyle(0)
    pad_low.SetBottomMargin(0.12)
    pad_low.SetTopMargin(0.72)
    pad_low.SetLeftMargin(0.13)
    pad_low.SetRightMargin(0.05)
    pad_low.SetGrid()
    pad_low.Draw()

    #
    pad_up.cd()

    mc_stack.Draw("hist")
    data_hist.Draw("E0P same")

    #axes:
    mc_stack.SetMinimum(1.e0)
    mc_stack.SetMaximum(1.e6)
    mc_stack.GetYaxis().SetTitle("Events")
    mc_stack.GetXaxis().SetLabelSize(0)
    data_hist.GetXaxis().SetTitle(plot.label)

    #labels:
    Text = ROOT.TLatex()

    Text.SetNDC() 

    Text.SetTextAlign(12);
    Text.SetNDC(ROOT.kTRUE) 
    Text.SetTextSize(0.025) 
    Text.DrawLatex(0.17, 0.95, "#it{ALEPH data and simulation}") 
    lumi_tag = "#bf{{ #sqrt{{s}} = {:.0f} GeV, L = {:.2f} pb^{{-1}} }}".format(ecm, lumi)
    Text.DrawLatex(0.17, 0.9, lumi_tag)
    ana_tag = "#bf{{ {} , {} }}".format(year, sel_tag)
    Text.DrawLatex(0.17, 0.85, ana_tag)

    #legend
    legsize = 0.05*((len(mc_hists_list)+1)/2)
    # leg = ROOT.TLegend(0.58,0.86 - legsize,0.86,0.88)
    # leg.Draw()
    leg = ROOT.TLegend(0.6, 0.95 - legsize, 0.95, 0.95)
    for mc_hist in reversed(mc_hists_list):
        leg.AddEntry(mc_hist, mc_hist.GetTitle(), "f")
    leg.AddEntry(data_hist, data_hist.GetTitle(), "ep")

    leg.SetFillStyle( 0 )
    leg.SetBorderSize( 0 )
    leg.SetTextFont( 43 )
    leg.SetTextSize( 22 )
    leg.SetNColumns( 2 )
    leg.SetColumnSeparation(-0.05)
    leg.Draw()

    #ratio panel
    hist_ratio = data_hist.Clone()
    hist_ratio.Divide( mc_stack.GetStack().Last() )
    hist_ratio.GetYaxis().SetTitle("data/MC")
    hist_ratio.GetYaxis().SetTitleOffset(1.95)
    hist_ratio.GetYaxis().SetNdivisions(6)

    #set some fixed ratio: 
    if fix_ratio_range:
        hist_ratio.GetYaxis().SetRangeUser(fix_ratio_range[0], fix_ratio_range[1])


    pad_low.cd()
    pad_low.Update()
    hist_ratio.Draw("E0P")
    pad_low.RedrawAxis()

    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()	

    if not os.path.exists(out_dir_base):
        os.makedirs(out_dir_base)

    if weighted:
        filename = "hists_"+plot.name+"_weighted"+out_format
    else:
        filename = "hists_"+plot.name+out_format
        
    fileout = os.path.join(out_dir_base, filename)


    canvas.SaveAs(fileout)


    #write out a root file for use in WS building if requested:
    if store_root_file:
        print("Do we need this")
        # filename = "histograms_"+variable+"_"+cut_name+".root"
        # fileoutpath = os.path.join(out_dir_base, "WS_Input")
        # if not os.path.exists(fileoutpath):
        # 	os.makedirs(fileoutpath)
        # fileoutpath = os.path.join(fileoutpath, filename)
        # fileout = ROOT.TFile(fileoutpath, "RECREATE")
        # fileout.cd()

        # #subdirectories for categories
        # # if "_dphi" in extra_cut:
        # if extra_cut:
        # 	fileout.mkdir(extra_cut_name.strip("_"))
        # 	fileout.cd(extra_cut_name.strip("_"))



        # signal_hist.Write()

        # print("total signal evts = ", signal_hist.Integral())

        # for bkg_hist in bkg_hists_list:
        # 	proc_name = bkg_hist.GetName().split(variable)[0].rstrip("_")
        # 	# print(proc_name)
        # 	bkg_hist.SetName(proc_name)
        # 	bkg_hist.SetTitle(proc_name)
        # 	bkg_hist.SetLineColor(ROOT.kBlack)
        # 	bkg_hist.Write()
        # 	print("total {} evts = {}".format(proc_name, bkg_hist.Integral()) )

        # #add an observation

        # #this somehow stopped working??
        # # obs_hist = bkg_stack.GetStack().Last().Clone("data_obs")
        # # obs_hist.SetName("data_obs")
        # # obs_hist.SetTitle("data_obs")
        # # obs_hist.Write()

        # data_hist.Write()

        # print("Wrote root file:", fileoutpath)

        # fileout.Close()
     

if __name__ == "__main__":


    for plot_name, plot_specs in PlottingConfig.plots_dict.items():
        print(plot_name, plot_specs)

        make_plot(plot_specs, PlottingConfig.inputs_path, PlottingConfig.data, PlottingConfig.mc_processes, PlottingConfig.outputs_path, 
              year=PlottingConfig.year, sel_tag =PlottingConfig.sel_tag, lumi=PlottingConfig.lumi, ecm=PlottingConfig.ecm,
              do_log_y=PlottingConfig.do_log_y, add_overflow=PlottingConfig.add_overflow, fix_ratio_range=PlottingConfig.ratio_range, 
              weighted=PlottingConfig.weighted, out_format=PlottingConfig.out_format, store_root_file=PlottingConfig.store_root_file
           )
    
    exit()


    #use the makeplot
    
    
    #get the data histogram 
    data_proc = 'Znn'

    histo_name = 'event_invariant_mass'
    histo_nbins = 70
    histo_xmin = 50.
    histo_xmax = 120.
    histo_xlabel = "m_{jj} in GeV"
    
    data_filepath = os.path.join(input_filebase, f"{data_proc}_{sel_level}.root")
    data_hist = get_hist_from_tree("data", data_filepath, histo_name, histo_nbins, histo_xmin, histo_xmax, is_data= True)

    mc_dict={
        #name:(sample_list, label)
        'Zbb':(['Zbb'], "b-jets"),
        'Zcc':(['Zcc'], "c-jets"),
        'Zudsuds':(['Zss', 'Zuu', 'Zdd'], "light jets"),
    }
    mc_stack = ROOT.THStack("mc_stack", "mc_stack")

    for mc_proc, (mc_sample_list, mc_label) in mc_dict.items():
        print(mc_label, mc_sample_list)

        hist_mc_proc = None
        tot_hist_template = ROOT.TH1D(f"tot_hist_{mc_proc}", f"tot_hist_{mc_proc}", histo_nbins, histo_xmin, histo_xmax)

        for mc_sample in mc_sample_list:
            mc_sample_filepath = os.path.join(input_filebase, f"{mc_sample}_{sel_level}.root")
            mc_hist_proc = get_hist_from_tree(mc_sample, mc_sample_filepath, histo_name, histo_nbins, histo_xmin, histo_xmax, is_data= False)
            tot_hist_template.Add(mc_hist_proc)

        hist_mc_proc = copy.deepcopy(tot_hist_template)
        hist_mc_proc.SetTitle(mc_label)

        mc_stack.Add(hist_mc_proc)
    
    #test draw:
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800) 
    canvas.SetLogy(True)
    canvas.cd()

    mc_stack.Draw("HIST")
    data_hist.Draw("HIST SAME")

    canvas.SaveAs("test.png")




    #print total data and MC
    data_evts = data_hist.Integral()
    mc_evts = mc_stack.GetStack().Last().Integral()

    print(f"Data events = {data_evts}")
    print(f"MC events = {mc_evts}")
    print(f"Average data/MC ratio = {data_evts/mc_evts}")
    print(f"Average MC/data ratio = {mc_evts/data_evts}")
    



        #     # print("Adding hists to stack:", hist_ttyy_bkg.GetTitle(), hist_vyy_bkg.GetTitle() )
        
        
        # hs_bkgs.Add(hist_vyy_bkg)
        # hs_bkgs.Add(hist_ttyy_bkg)

