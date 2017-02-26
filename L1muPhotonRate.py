#!/usr/bin/env python

"""
Example L1TNtuple analysis program
"""

import ROOT, os
from copy import deepcopy as dc

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptStat(0)



def acceptL1_Dimuon6_OS(muons):
    #print muons
    checker = False
    if (len(muons) >= 2):
        if(muons[0].lorentzVector.Pt() >= 6.0 and muons[1].lorentzVector.Pt() >= 6.0):
            if (muons[0].chg*muons[1].chg < 1.0):
                checker = True
    return checker  
    #return True





class L1Tau(object):
    def __init__(self, et = -99., eta = -99., phi = -99. , iso = -1, bx = -99):
        self.et  = et 
        self.eta = eta
        self.phi = phi
        self.iso = iso
        self.bx  = bx
        self.lorentzVector = ROOT.TLorentzVector()
        self.lorentzVector.SetPtEtaPhiM(et, eta, phi, 1.77682)

    
    def __str__(self):
        return 'Et = %.2f \t eta %.2f \t phi %.2f \t iso %d \t bx %d' %(self.et , self.eta, self.phi, self.iso, self.bx)

class L1Muon(object):
    def __init__(self, et = -99., eta = -99., phi = -99. , chg = -99, qual = -99, iso = -1, bx = -99):
        self.et  = et 
        self.eta = eta
        self.phi = phi
        self.chg = chg
        self.qual = qual
        self.iso = iso
        self.bx  = bx
        self.lorentzVector = ROOT.TLorentzVector()
        self.lorentzVector.SetPtEtaPhiM(et, eta, phi, 0.1056583745)
    
    def __str__(self):
        return 'Et = %.2f \t eta %.2f \t phi %.2f \t Chg = %.2f \t Qual = %.2f \t iso %d \t bx %d' %(self.et , self.eta, self.phi, self.chg, self.qual, self.iso, self.bx)

class L1Eg(object):
    def __init__(self, et = -99., eta = -99., phi = -99., iso = -1, bx = -99):
        self.et  = et 
        self.eta = eta
        self.phi = phi
        self.iso = iso
        self.bx  = bx
        self.lorentzVector = ROOT.TLorentzVector()
        self.lorentzVector.SetPtEtaPhiM(et, eta, phi, 0.0)

    def __str__(self):
        return 'Et = %.2f \t eta %.2f \t phi %.2f \t iso %d \t bx %d' %(self.et , self.eta, self.phi, self.iso, self.bx)



def eventLoop(files, total_rate = 1., nevents = -1, verbose = False):
    
    countEvts = 0
    nPassedTrigBitL1 = 0
    lumiMin = -1
    lumiMax = -1
    nThresholds = 80;
    EtMin = 0.;
    EtMax = 20.;
    lumiBinMax = 700

    # load file
    outfile = ROOT.TFile.Open('rate_histograms.root', 'recreate')
    
    differential_rate = ROOT.TH1F('ditau_differential_rate', 'ditau_differential_rate', 80, 0., 80.)
    cumulative_rate = dc(differential_rate)
    cumulative_rate.SetTitle('ditau_cumulative_rate')
    cumulative_rate.SetName ('ditau_cumulative_rate')
     
    h_muon1Pt = ROOT.TH1F('h_muon1Pt', 'Leading Muon Pt; Pt (GeV);Rate (Hz)', 80, 0., 80.)
    h_muon2Pt = dc(h_muon1Pt)
    h_muon2Pt.SetTitle('Trailing Muon Pt; Pt (GeV);Rate (Hz)')
    h_muon2Pt.SetName ('h_muon2Pt')
    h_egPt = dc(h_muon1Pt)
    h_egPt.SetTitle('EG Pt; Pt (GeV);Rate (Hz)')
    h_egPt.SetName ('h_egPt')
    h_dimuonPtSpread = ROOT.TH1F('h_dimuonPtSpread', 'Dimuon Pt Spread; Pt (GeV);Rate (Hz)', 160, -80., 80.)
    h_dimuonMass = ROOT.TH1F('h_dimuonMass', 'Dimuon Mass; Pt (GeV);Rate (Hz)', 80, 0., 80.)

    h_RateVsegPt = ROOT.TH1F ("h_RateVsegPt", "h_RateVsegPt", 100, 0, 100)
    cumulative_rateVsegPt   = dc(h_muon1Pt)
    cumulative_rateVsegPt.SetTitle('egPt_cumulative_rate')
    cumulative_rateVsegPt.SetName ('egPt_cumulative_rate')

    h_NevtsL1VsLumi = ROOT.TH1F ("h_NevtsL1BitVsLumi","h_NevtsL1BitVsLumi",lumiBinMax,0,lumiBinMax)
    h_NevtsVsThreshold = ROOT.TH1F("h_NevtsVsThreshold","h_NevtsVsThreshold",nThresholds,EtMin,EtMax)
    
    
    print 'start loading %d files' %len(files)
    
    treeEvent = ROOT.TChain("l1EventTree/L1EventTree")
    treeL1up  = ROOT.TChain("l1UpgradeEmuTree/L1UpgradeTree")

    for file in files:
        treeEvent.Add(file)
        treeL1up.Add(file)

    print 'loaded %d files containing %d events' %(len(files), treeEvent.GetEntries())

    treeEvent.AddFriend(treeL1up)

    for jentry, event in enumerate(treeEvent):
        if jentry%10000 == 0:
            print '=============================> event %d / %d' %(jentry, treeEvent.GetEntries())
        if nevents > 0 and jentry >= nevents:
            break
        
        ev  = treeEvent.Event
        sim = treeL1up.L1Upgrade
        lumiBlock = ev.lumi
        nBunches = ev.bx

        ### Taus
        if verbose: print '\nevent: %d, nTaus: %d' %(ev.event, sim.nTaus)
        
        taus = []
        
        for i in range(sim.nTaus):
            taus.append(L1Tau(sim.tauEt[i], sim.tauEta[i], sim.tauPhi[i], sim.tauIso[i], sim.tauBx[i]))
        
        taus.sort(key=lambda x: x.et, reverse=True)
                
        for i, tau in enumerate(taus):
            if verbose: print '\ttau %d \t' %i, tau
        
        isotaus = [tau for tau in taus if tau.iso == 1 and abs(tau.eta) < 2.17 and tau.bx == 0]
        #############

        ### Egs
        if verbose: print '\nevent: %d, nEGs: %d' %(ev.event, sim.nEGs)
        
        egs = []
        
        for i in range(sim.nEGs):
            egs.append(L1Eg(sim.egEt[i], sim.egEta[i], sim.egPhi[i], sim.egIso[i], sim.egBx[i]))
        
        egs.sort(key=lambda x: x.et, reverse=True)
                
        for i, eg in enumerate(egs):
            if verbose: print '\teg %d \t' %i, eg
        
        isoegs = [eg for eg in egs if eg.iso == 1]
        #############

        ### Muons
        if verbose: print '\nevent: %d, nMuons: %d' %(ev.event, sim.nMuons)
        
        muons = []
        
        for i in range(sim.nMuons):
            muons.append(L1Muon(sim.muonEt[i], sim.muonEta[i], sim.muonPhi[i], sim.muonChg[i], sim.muonQual[i], sim.muonIso[i], sim.muonBx[i]))
        
        muons.sort(key=lambda x: x.et, reverse=True)

                
        for i, muon in enumerate(muons):
            if verbose: print '\tmuon %d \t' %i, muon
        
        isomuons = [muon for muon in muons if muon.iso == 1]
        #############

        countEvts+=1
        if(lumiMin < 0):
          lumiMin = lumiBlock 
        if(lumiMin >= 0 and lumiBlock < lumiMin) :
          lumiMin = lumiBlock

        if(lumiMax < 0):
          lumiMax = lumiBlock 
        if(lumiMax >= 0 and lumiBlock > lumiMax):
          lumiMax = lumiBlock
          
        if (acceptL1_Dimuon6_OS(muons)):
            #print "Passou: L1_Dimuon6_OS"
            nPassedTrigBitL1=+1
            h_NevtsL1VsLumi.Fill(lumiBlock)

            h_muon1Pt.Fill(muons[0].lorentzVector.Pt())
            h_muon2Pt.Fill(muons[1].lorentzVector.Pt())
            if(len(egs) > 0):
                h_egPt.Fill(egs[0].lorentzVector.Pt())
                h_RateVsegPt.Fill(egs[0].lorentzVector.Pt())
                threshold=[]
                bin_center=[]
                for i_bin in range(h_NevtsVsThreshold.GetNbinsX()):
                    threshold.append(float(h_NevtsVsThreshold.GetXaxis().GetBinLowEdge(i_bin+1)))
                    bin_center.append(float(h_NevtsVsThreshold.GetXaxis().GetBinCenter(i_bin+1)))
                    if(egs[0].lorentzVector.Pt()<threshold):
                        h_NevtsVsThreshold.Fill(bin_center[i_bin],279.160)
            h_dimuonPtSpread.Fill(muons[0].lorentzVector.Pt() - muons[1].lorentzVector.Pt())
            h_dimuonMass.Fill((muons[0].lorentzVector + muons[1].lorentzVector).M())
            #h_RateVsegPt.Fill(egs[0].lorentzVector.Pt())


    os.system("rm -rf l1Plots ; mkdir l1Plots")


    outfile.cd()
    h_muon1Pt.Write()     
    h_muon1Pt.Draw('HIST')
    ROOT.gPad.SetLogy()
    ROOT.gPad.SaveAs('l1Plots/h_muon1Pt_L1_Dimuon6_OS.png')

    outfile.cd()
    h_muon2Pt.Write()     
    h_muon2Pt.Draw('HIST')
    ROOT.gPad.SetLogy()
    ROOT.gPad.SaveAs('l1Plots/h_muon2Pt_L1_Dimuon6_OS.png')

    outfile.cd()
    h_egPt.Write()     
    h_egPt.Draw('HIST')
    ROOT.gPad.SetLogy()
    ROOT.gPad.SaveAs('l1Plots/h_egPt_L1_Dimuon6_OS.png')

    outfile.cd()
    h_dimuonPtSpread.Write()     
    h_dimuonPtSpread.Draw('HIST')
    ROOT.gPad.SetLogy()
    ROOT.gPad.SaveAs('l1Plots/h_dimuonPtSpread_L1_Dimuon6_OS.png')

    outfile.cd()
    h_dimuonMass.Write()     
    h_dimuonMass.Draw('HIST')
    ROOT.gPad.SetLogy(0)
    ROOT.gPad.SaveAs('l1Plots/h_dimuonMass_L1_Dimuon6_OS.png')


    allevents = float(treeEvent.GetEntries() if nevents<0 else nevents)
    print "allevents: ", allevents
    bins_eg = []
    for bin in range(h_RateVsegPt.GetNbinsX()):
        bins_eg.append(float(h_RateVsegPt.GetBinContent(bin+1)))
    rates_eg = []
    for i, bin in enumerate(bins_eg):
        rates_eg.append(total_rate * sum(bins_eg[i:])/ allevents)
        print "total_rate:",total_rate,"; ",bins_eg[i:]
    for bin in range(cumulative_rateVsegPt.GetNbinsX()):
        cumulative_rateVsegPt.SetBinContent(bin+1, rates_eg[bin])
    ROOT.gPad.SetLogy()
    cumulative_rateVsegPt.Write()
    cumulative_rateVsegPt.Draw('HIST')
    ROOT.gPad.SaveAs('rate_egPt.pdf')
    # bins = []

    # for bin in range(differential_rate.GetNbinsX()):
    #     bins.append( float(differential_rate.GetBinContent(bin+1)) )

    # rates = []
    
    # for i, bin in enumerate(bins):
    #     rates.append(total_rate * sum(bins[i:]) / allevents)
    
    # for bin in range(cumulative_rate.GetNbinsX()):
    #     cumulative_rate.SetBinContent(bin+1, rates[bin])

    # ROOT.gPad.SetLogy()
    
    # outfile.cd()
    # cumulative_rate.Write()     
    # cumulative_rate.Draw('HIST')
    # ROOT.gPad.SaveAs('rate.pdf')
    
    print "Analyzed ", countEvts, " events" 
    print "Min and max lumi section: (", lumiMin , "," , lumiMax , ")" 
    print "Nr. events passing L1 bit acceptL1_Dimuon6_OS ", ": " , nPassedTrigBitL1 
    
    
    effTrigBitL1 = nPassedTrigBitL1/countEvts
    float(effTrigBitL1)
    print "Efficiency of TrigBit L1 :", float(effTrigBitL1)
    #func_effTrigBitL1.SetParameter(0,float(effTrigBitL1))
    
    lumiLength = 23.3 #s
    deltaLumi = lumiMax - lumiMin
    Lum = 20. #10^33cms-2s-1
    lumiRef = 8.5 #10^33cms-2s-1 
    #rate do HLT da sample = ???
    #prescale of the reference trigger =??
    if nBunches < 0 :
      scale = (-1.*nBunches)/(deltalumi*23.31)  
    else:
       scale = 11246. #ZB per bunch in Hz
       scale /= countEvts # in Hz
       scale *= nBunches
       
    print "Scale by ", " 11246 / nZeroBiasevents * NumberOfBunches, with nZeroBiasevents = " 
    print countEvts , " NumberOfBunches = " , nBunches 
    
    h_RateVsThreshold = h_NevtsVsThreshold.Clone("h_RateVsThreshold")
    print "deltaLumi ", deltaLumi, "lumiBlock:", lumiBlock
    print "nBunches:", nBunches
    print "lumiLength ", lumiLength
    print "Lum " , Lum
        
    h_RateVsThreshold.Scale((1./(deltaLumi*lumiLength)*(Lum/lumiRef)))
    #h_RateVsThreshold.Scale(scale)
    
    outfile.cd()
    h_RateVsThreshold.Write()     
    h_RateVsThreshold.Draw('HIST')
    #ROOT.gPad.SetLogy()
    ROOT.gPad.SaveAs('l1Plots/h_RateVsThreshold.png')

    h_RateVsL1VsLumi = h_NevtsL1VsLumi.Clone('h_RateVsL1VsLumi')
    h_RateVsL1VsLumi.Scale(1./lumiLength)

    outfile.cd()
    h_RateVsL1VsLumi.Write() 
    h_RateVsL1VsLumi.Draw('HIST')
    ROOT.gPad.SaveAs('l1Plots/h_RateVsL1VsLumi.png')
   
    outfile.cd()
    h_NevtsVsThreshold.Write() 
    h_NevtsVsThreshold.Draw('HIST')
    ROOT.gPad.SaveAs('l1Plots/h_NevtsVsThreshold.png')


if __name__ == "__main__":
    #eventLoop(['root://xrootd.unl.edu//store/group/dpg_trigger/comm_trigger/L1Trigger/L1Menu2016/Stage2/l1-tsg-v4/ZeroBias2/crab_l1-tsg-v4__259721_ZeroBias2/160315_144909/0000/L1Ntuple_%d.root' % i for i in range(1,35)])
    eventLoop(['file:/uscmst1b_scratch/lpc1/3DayLifetime/ftorresd/L1Ntuples/L1Ntuple_%d.root' % i for i in range(1,500)])
    #eventLoop(['file:L1Ntuple_99.root' % i for i in range(1,2)])
    
