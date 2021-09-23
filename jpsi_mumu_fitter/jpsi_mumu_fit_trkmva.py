#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
from ROOT import RooRealVar, RooPolynomial, RooCategory, RooFit,\
    RooDataSet, RooArgSet, RooCBShape, RooGaussian,\
    RooExponential, RooAbsRealLValue, RooAddPdf, RooArgList 
import numpy as np
import uproot
import matplotlib.pyplot as plt

from tdrstyle import *

from CMS_lumi import *

setTDRStyle()

# In[2]:


from ROOT import RooChebychev


# In[7]:


plt.style.use('../notebooks/style')


# In[3]:


c = ROOT.TCanvas("c","Tau Mass Fit",660,660)
#c.Divide(1,2)
pad1 = ROOT.TPad("pad1", "Pad for histogram", 0,0.3,1,1.0)
pad2 = ROOT.TPad("pad2", "Pad for the ratio plot", 0,0.0,1,0.3)




bujpsi_filename = '../Ntuples/2018/Muon_Validation/JPsiSelection_combined_mc_BuJPsi_private_trackerMuonId_v3.root'
bparking_filename = '../Ntuples/2018/Muon_Validation/JPsiSelection_combined_data_BScut_trackerMuonId_v3.root'

# In[11]:


def t3mselection_reader(filename, br_list=['var_Muon1_Pt', 'var_Muon2_Pt', 'var_Muon3_Pt'], cuts=''):
    tree_list = ['TreeB', 'TreeS_Ds', 'TreeS_Bu', 'TreeS_Bd']
    tree_dict = {}
    for tree_ in tree_list:
        _ = root2array(filename, tree_, br_list, cuts)
        tree_dict[tree_] = rec2array(_)
    return tree_dict


# In[12]:


def select_tracker_muons(df):
    return df[(df['Muon_IsGlobalMuon']==0) &
              (df['Muon_IsPFMuon']==1) &
              (df['Muon_IsTrackerMuon']==1)]


# In[13]:


file_bujpsi = uproot.open(bujpsi_filename)
file_bparking = uproot.open(bparking_filename)


# In[14]:


tree_bujpsi = file_bujpsi['tree']
tree_bparking = file_bparking['tree']


# In[15]:


validation_data = tree_bparking.pandas.df(['var_trackerMuonId_without_calo','Muon_IsTrackerMuon','Muon_IsGlobalMuon','Muon_IsPFMuon','InvMass'])
validation_mc = tree_bujpsi.pandas.df(['var_trackerMuonId_without_calo','Muon_IsTrackerMuon','Muon_IsGlobalMuon','Muon_IsPFMuon','InvMass'])


# In[ ]:

trkmva_bins = [-0.36, -0.10, -0.08, -0.05, -0.02, 0.0, 0.03, 0.06, 0.09, 0.12, 0.36]


# In[17]:


jpsi_mass_categories = {'data': {}, 'mc': {}}
for ibin in xrange(len(trkmva_bins)-1):
    jpsi_mass_categories['data']['trkmva_bin_'+str(ibin+1)] = select_tracker_muons(validation_data[(validation_data['var_trackerMuonId_without_calo']>trkmva_bins[ibin])&
                                                                       (validation_data['var_trackerMuonId_without_calo']<trkmva_bins[ibin+1])])['InvMass'].to_numpy()
    jpsi_mass_categories['mc']['trkmva_bin_'+str(ibin+1)] = select_tracker_muons(validation_mc[(validation_mc['var_trackerMuonId_without_calo']>trkmva_bins[ibin])&
                                                                       (validation_mc['var_trackerMuonId_without_calo']<trkmva_bins[ibin+1])])['InvMass'].to_numpy()


# In[18]:


rds_data = {}


# In[19]:


x = RooRealVar("x","x",3.0,2.95,3.3)
x_args = RooArgSet(x)
x.setBins(70)


# In[20]:


for cat_ in jpsi_mass_categories['data']:
    rds_data[cat_] = RooDataSet("rds_data_"+cat_, "", x_args)
    for d in jpsi_mass_categories['data'][cat_]:
        if d<2.95: continue
        RooAbsRealLValue.__assign__(x, d)
        rds_data[cat_].add(x_args)


# ---------------------------------------
# set ranges
# ---------------------------------------
x.setRange("R1",3.05,3.15) #first peak JPsi(3.097GeV)
x.setRange("R2",2.95,3.05) #background
x.setRange("R3",3.15,3.3) #background
x.setRange("R6",2.95,3.3) #full range


# In[32]:


class jpsi_fitter():
    
    def __init__(self):       
        self.x = RooRealVar("x","M(#mu#mu)",3.0,2.95,3.3)
        self.x_args = RooArgSet(self.x)
        self.x.setBins(70)
        
        # parameters for the jpsi crystalball fit
        self.jpsi_crystal_mean = RooRealVar("jpsi_crystal_mean","mean of crystal ball (jpsi)",3.097,3.09,3.1)
        self.jpsi_crystal_sigma = RooRealVar("jpsi_crystal_sigma","width of crystal ball (jpsi)",0.015,0.00,0.030)
        self.jpsi_n = RooRealVar("jpsi_n","power parameter in the crystal ball (jpsi)", 1.0, 0.0, 10.0)
        self.jpsi_alpha = RooRealVar("jpsi_alpha","boundry in the crystal ball (jpsi)",1.0, -10.0, 10.0)
        self.jpsi_crystal_constant = RooRealVar("jpsi_crystal_constant", "constant (cb, jpsi)",500, 0, 5.0e5)

        #parameters for the jpsi gaussian fit
        self.jpsi_gaus_constant = RooRealVar("jspi_gaus_constant","constant (gaus, jpsi)",500, 0, 5.0e5)
        self.jpsi_gaus_mean = RooRealVar("jpsi_gaus_mean","mean of the gaussian",3.097,3.09,3.1)
        self.jpsi_gaus_sigma = RooRealVar("jpsi_gaus_sigma","sigma of the gaussian",0.15,0.00,0.030)

        #parameters for the exponential fit
        self.bkg_rate = RooRealVar("bkg_rate", "falling rate of the background",1.0, -10.0, 10.0)
        self.bkg_constant = RooRealVar("bkg_constant","Normalization constant",250, 0, 5.0e5)

        #parameters for chebyshev function fit
        self.bkg_cheb_c1 = RooRealVar("bkg_cheby_c1","coeff 1 of chebychev",-0.1,-2,1)
        self.bkg_cheb_c2 = RooRealVar("bkg_cheby_c2","coeff 2 of chebychev",0.04,-2.0,1.0)
        self.bkg_cheb_constant = RooRealVar("bkg_cheb_constant","chebychev constant",2500, 0, 5.0e5)

        self.cb_jpsi = RooCBShape("jpsi_shape", "Cystal Ball Function (jpsi)",
                                  self.x,
                                  self.jpsi_crystal_mean,
                                  self.jpsi_crystal_sigma,
                                  self.jpsi_alpha,
                                  self.jpsi_n)

        self.gaus_jpsi = RooGaussian("gaus_chape", "Gaussian Function (jpsi)",
                                     self.x,
                                     self.jpsi_gaus_mean,
                                     self.jpsi_gaus_sigma)
        
        self.exp_bkg = RooExponential("bkg_shape", "Exponential Function (Bkg)",
                                      self.x,
                                      self.bkg_rate)
        self.cheb_bkg = RooChebychev("bkg_cheb_shape", "2nd order Chebychev Function (Bkg)",
                                     self.x,
                                     RooArgList(self.bkg_cheb_c1, self.bkg_cheb_c2))

        self.model = RooAddPdf("model", "model", RooArgList(self.gaus_jpsi, self.cb_jpsi, self.cheb_bkg, self.exp_bkg),
                          RooArgList(self.jpsi_gaus_constant, self.jpsi_crystal_constant, self.bkg_cheb_constant, self.bkg_constant))

        self.fitResults = 0
    
    def reset_params(self):
        
        self.jpsi_crystal_mean.setVal(3.097)
        self.jpsi_crystal_sigma.setVal(0.015)
        self.jpsi_n.setVal(1.0)
        self.jpsi_alpha.setVal(1.0)
        self.jspi_crystal_constant.setVal(500.0)
        
        self.bkg_rate.setVal(1.0)
        self.bkg_constant.setVal(250.0)
        self.bkg_cheb_c1.setVal(-0.1)
        self.bkg_cheb_c2.setVal(0.04)
        self.bkg_cheb_constant(250.0)

        self.fitResults = 0
    
    
    def fit_data(self, rds):
        # ------------------------------
        # Fit model on MC (signal)
        # ------------------------------
        self.cb_jpsi.fitTo(rds, RooFit.Range("R2"))
        self.gaus_jpsi.fitTo(rds, RooFit.Range("R2"))
        self.exp_bkg.fitTo(rds, RooFit.Range("R1,R3"))
        self.cheb_bkg.fitTo(rds, RooFit.Range("R1,R3"))
        self.fitResult = self.model.fitTo(rds, RooFit.Save(True))
        self.fitResult.Print()
    
    def plot(self, rds, c, pad1, pad2, outputFile='test_fit.png'):

        #------
        # Draw
        #------
        # Construct plot frame in 'x'
        xframe = self.x.frame()
        ras_bkg = RooArgSet(self.exp_bkg, self.cheb_bkg)
        ras_signal = RooArgSet(self.gaus_jpsi, self.cb_jpsi)

        rds.plotOn(xframe)
        self.model.plotOn(xframe, RooFit.Components(ras_signal),
                     RooFit.LineColor(ROOT.kRed),
                     RooFit.LineStyle(ROOT.kDashed))
        self.model.plotOn(xframe, RooFit.Components(ras_bkg),
                     RooFit.LineColor(ROOT.kBlue),
                     RooFit.LineStyle(ROOT.kDashed))
        self.model.plotOn(xframe, RooFit.LineColor(ROOT.kBlue))

        #print(model.getMaxVal(RooArgSet(x)))
        hresid = xframe.residHist()
        
        frame2 = self.x.frame()
        frame2.addPlotable(hresid,"Pe")

        c.cd()
        xframe.Draw()

        fitChi2 = xframe.chiSquare()
        latex = ROOT.TLatex()
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.SetTextSize(0.03)
        latex.DrawLatex(3.2,self.jpsi_crystal_constant.getVal()/3,"#chi^{2}/ndf = %0.2f" % fitChi2)

        hresid.GetYaxis().SetTitle("Events/(0.01)")
        hresid.GetXaxis().SetTitle("Residual")
        ROOT.TGaxis().SetMaxDigits(3)
        ROOT.TGaxis().SetExponentOffset(-0.075, -0.01, "y")

        #pad1.cd()
        #xframe.SetTitle("")
        #xframe.Draw()
        #pad2.cd()
        #frame2.SetTitle("")
        #frame2.Draw()

        #c.cd()
        #pad1.Draw()
        #pad2.Draw()

        #CMS_lumi(c, 4, 0)
        c.SaveAs(outputFile)
    
    def get_integral(self):
        
        #Compute integrals
        x.setRange("signal",2.95,3.3)
        signal_range = RooFit.Range("signal")
        sideband_range = RooFit.Range("sideband")
        
        l = RooArgSet(self.x)
        nx_total = RooFit.NormSet(l)

        jpsi_crystal_constant_err = self.jpsi_crystal_constant.getPropagatedError(self.fitResult)
        jpsi_gaus_constant_err = self.jpsi_crystal_constant.getPropagatedError(self.fitResult)

        #fraction of total events in 1.93,2.01 (n_signal_region_events/n_total_events)
        fsignal_model_cb = self.cb_jpsi.createIntegral(RooArgSet(self.x), nx_total, signal_range)
        fs_cb = fsignal_model_cb.getVal()
        fs_err_cb = fsignal_model_cb.getPropagatedError(self.fitResult)

        fsignal_model_gaus = self.gaus_jpsi.createIntegral(RooArgSet(self.x), nx_total, signal_range)
        fs_gaus = fsignal_model_gaus.getVal()
        fs_err_gaus = fsignal_model_gaus.getPropagatedError(self.fitResult)

        ns = fs_cb*(self.jpsi_crystal_constant.getVal())+fs_gaus*(self.jpsi_gaus_constant.getVal())
        ns_err = np.sqrt((fs_err_cb*self.jpsi_crystal_constant.getVal())**2+(jpsi_crystal_constant_err*fs_cb)**2+
                        (fs_err_gaus*self.jpsi_gaus_constant.getVal())**2+(jpsi_gaus_constant_err*fs_gaus)**2)

        return ns, ns_err

yields_data = {}

# In[34]:

for i in range(1, len(trkmva_bins)):
    binstr = 'trkmva_bin_'+str(i)
    fitter_ = jpsi_fitter()
    fitter_.fit_data(rds_data[binstr])
    fitter_.plot(rds_data[binstr], c, pad1, pad2, 'fits/jpsi_mumu_tracker_muon_data_trkmva_'+binstr+'.png')
    yields_data[binstr] = fitter_.get_integral()
for i in range(1, len(trkmva_bins)):
    binstr = 'trkmva_bin_'+str(i)
    print("%.2f\t%.2f\t%.2f" % (yields_data[binstr][0], yields_data[binstr][1], len(jpsi_mass_categories['mc'][binstr])))

