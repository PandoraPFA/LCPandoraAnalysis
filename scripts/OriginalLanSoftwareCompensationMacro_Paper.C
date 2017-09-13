#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TH2F.h>
#include <TString.h>
#include <TLeaf.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <vector>
#include <TROOT.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TStyle.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <iostream>

using namespace std;

bool useEini = false;
bool useSumEhit = false;

//const int nE = 14;
const int nE = 10;
//double Ebeam[nE] = {1, 3, 5, 7};
double Ebeam[nE] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 95};
//double Ebeam[nE] = {1, 3, 5, 7, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95};

const int NBIN = 10;
const int NPAR = 9;

//float lowMIP[NBIN]  = {0.3,  2, 5.5,  8, 10, 14, 17, 21, 25, 30};//First bin has to start at 0.3 now as the cut in MarlinPandora changed
//float highMIP[NBIN] = {  2, 5.5,   8, 10, 14, 17, 21, 25, 30, 1e6};

//const float bins[nBins] = {0.f, 2.f, 5.f, 7.5f, 9.5f, 13.f, 16.f, 20.f, 23.5f, 28.f, 1e6f};

double lowDensity[NBIN]  = {0,  2.,  5., 7.5, 9.5,  13., 16.,  20., 23.5,  28.};
double highDensity[NBIN] = {2., 5., 7.5, 9.5,  13., 16., 20., 23.5,  28.,  1e6};
//float lowDensity[NBIN]  = {0,   5, 9.5, 13, 16,  20};
//float highDensity[NBIN] = {5, 9.5,  13, 16, 20, 1e6};
//float lowDensity[NBIN]  = {0,   5, 7.5, 9.5, 13, 16, 20, 24, 28,  32};
//float highDensity[NBIN] = {5, 7.5, 9.5,  13, 16, 20, 24, 28, 32,  1e6};
float binDensity[NBIN];

bool m_Debug = false;

typedef std::vector<float> FloatVector;
typedef std::vector<int> IntVector;

double m_p[NBIN*3];


void MakeBinDensity(){
  for (int ibin = 0; ibin < NBIN; ibin++){
    binDensity[ibin] = (lowDensity[ibin]+highDensity[ibin])/2;
    if (ibin==(NBIN-1))
      binDensity[ibin] = 30;
  }
}

/*-------------------------------------------------------------------*/

int FindHitBin(float hitEnergy, float cellvolume){  
  bool m_debugDensity = false;
  
  float hitEnergyDensity = hitEnergy/cellvolume;
  
  //std::cout << "FindHitBin: hitEnergy " << hitEnergy << std::endl;

  int hitBin = -1;
  for (int ibin = 0; ibin < NBIN; ibin++){
    if (hitEnergyDensity>=lowDensity[ibin] && hitEnergyDensity<highDensity[ibin]){
      //std::cout << "hitEnergyDensity " << hitEnergyDensity << std::endl;
      hitBin = ibin;
    }
  }
  
  if (hitBin<0)
    cout << "Can't find hit bin" << endl;
  else
    return hitBin;   
}

/*-------------------------------------------------------------------*/

float FindHitBinDensity(float hitEnergy, float cellvolume){  
  bool m_debugDensity = false;
  
  float hitEnergyDensity = hitEnergy/cellvolume;
  
  float hitBinDensity = 0;
  for (int ibin = 0; ibin < NBIN; ibin++){
    if (hitEnergyDensity>=lowDensity[ibin] && hitEnergyDensity<highDensity[ibin])
      hitBinDensity = binDensity[ibin];
  }
  
  if (hitBinDensity<=0)
    cout << "Can't fin hit bin density" << endl;
  else
    return hitBinDensity;   
}



/*----------------------------------------------------------------------------------------------------------- Define class EventC -----------------------------------------------------------------------------------------------------------*/

class EventC {
  
public:
  EventC() { 
    m_HCalBinEnergy->clear();
    m_HCalHitEnergy->clear();
    m_ECalEnergy = 0;
    m_Ereco = 0;
    m_Ebeam = 0;
  }
  ~EventC(){}
  
  EventC(double Ebeam, double Ereco, FloatVector *HitEnergies, IntVector *HitType, FloatVector *CellSize0, FloatVector *CellSize1, FloatVector *CellThickness){ 
    m_HCalBinEnergy = new FloatVector;
    m_HCalBinEnergy->reserve(NBIN);
    m_HCalHitEnergy = new FloatVector;
    m_HCalHitEnergy->reserve(1000);
    m_Ebeam = Ebeam;
    m_Ereco = (double)Ereco;    
    m_ECalEnergy = 0;
    
    double EnergyPerBin[NBIN];
    double nHitsPerBin[NBIN];
    for (int i = 0; i < NBIN; i++){
      EnergyPerBin[i] = 0;
      nHitsPerBin[i] = 0;
    }
    
    for (unsigned int ihit = 0; ihit < HitEnergies->size(); ihit++){
        //if (m_Debug) cout << "ihit " << ihit << " energy " << HitEnergies->at(ihit) << endl;
        if (HitEnergies->at(ihit)<=0 || HitEnergies->at(ihit)>=10000)
          continue;

        if (HitType->at(ihit)==2)//HCAL hits      
      	{
        	  m_HCalHitEnergy->push_back(HitEnergies->at(ihit));
        	  float cellvolume = (CellSize0->at(ihit))*(CellSize1->at(ihit))*(CellThickness->at(ihit))/1000000.;
        	  int hitBin = FindHitBin(HitEnergies->at(ihit),cellvolume);
        	  EnergyPerBin[hitBin] += HitEnergies->at(ihit);
        	  nHitsPerBin[hitBin]++;
      	}
        else //ECAL hits 
      	{
        	  m_ECalEnergy += HitEnergies->at(ihit);
      	}
    }

    for (int i = 0; i < NBIN; i++){
      if (m_Debug) cout << "Bin " << i << " number of Hits " << nHitsPerBin[i] << endl;
      if (m_Debug) cout << "Energy of bin " << EnergyPerBin[i] 
			<< " inMIP (roughly) " << EnergyPerBin[i]/0.025 << endl;      
      m_HCalBinEnergy->push_back(EnergyPerBin[i]);
    }
    if (m_Debug) cout << "HCalBinEnergy size " << m_HCalBinEnergy->size() << endl;

  }

  double Ebeam() { return m_Ebeam; }
  double Ereco() { return m_Ereco; }
  double ECalEnergy() { return m_ECalEnergy; }
  FloatVector *HCalBinEnergy() { return m_HCalBinEnergy; }
  FloatVector *HCalHitEnergy() { return m_HCalHitEnergy; }

private:
  double m_Ebeam;
  double m_Ereco;
  double m_ECalEnergy;
  FloatVector *m_HCalBinEnergy;
  FloatVector *m_HCalHitEnergy;
};

/*----------------------------------------------------------------------------------------------------------- Get list of events -----------------------------------------------------------------------------------------------------------*/
std::vector<EventC> m_EventList;

void GetListOfEvents(int cellsize = 30, bool FitWithAllFiles = true, int FileNumber = -1){
  m_EventList.clear();

  //Get tree
  TString path = "/nfs/dust/ilc/user/huonglan/HCAL_Optimisation/SoftwareCompensation/RECO/Output/";

  for (int ifi = 0; ifi < nE-2; ifi++){    
    //if (ifi==4||ifi==5||ifi==6||ifi==7) continue;//For 10mm cell size
    //if (ifi==4||ifi==6) continue;//For 30mm cell size
    //if (!FitWithAllFiles && ifi!=FileNumber) continue;
    //if (ifi==6) continue;//For 40mm cell size
    if (ifi!=FileNumber) continue;

    TString fname = path + "HitsForSC_"; fname += Ebeam[ifi]; 
    fname += "GeV_Cell"; fname += cellsize; fname += ".root";

    cout << "input file " << fname << endl;

    TFile *f = TFile::Open(fname);

    TTree *m_tree = (TTree*)f->Get("HitEnergyTree");
    
    int          nPfos = 0;
    FloatVector *pfoE = 0;
    int          nClusters = 0;
    FloatVector *HitEnergies = 0;
    FloatVector *CellSize0 = 0;
    FloatVector *CellSize1 = 0;
    FloatVector *CellThickness = 0;
    FloatVector *RawEnergyOfCluster = 0;
    IntVector   *HitType = 0;
    
    TBranch *b_nPfos;
    TBranch *b_pfoE;
    TBranch *b_nClusters;
    TBranch *b_HitEnergies;
    TBranch *b_CellSize0;
    TBranch *b_CellSize1;
    TBranch *b_CellThickness;
    TBranch *b_RawEnergyOfCluster;
    TBranch *b_HitType;

    m_tree->SetBranchAddress("numberOfPfos",       &nPfos,         &b_nPfos);
    m_tree->SetBranchAddress("EnergyOfPfos",       &pfoE,          &b_pfoE);
    m_tree->SetBranchAddress("numberOfClusters",   &nClusters,     &b_nClusters);
    m_tree->SetBranchAddress("HitEnergies",        &HitEnergies,   &b_HitEnergies);
    m_tree->SetBranchAddress("CellSize0",          &CellSize0,     &b_CellSize0); 
    m_tree->SetBranchAddress("CellSize1",          &CellSize1,     &b_CellSize1);
    m_tree->SetBranchAddress("CellThickness",      &CellThickness, &b_CellThickness);
    m_tree->SetBranchAddress("RawEnergyOfCluster", &RawEnergyOfCluster, &b_RawEnergyOfCluster);
    m_tree->SetBranchAddress("HitType",            &HitType,       &b_HitType);

    for (long int iev = 0; iev < m_tree->GetEntries(); iev++){
      m_tree->GetEntry(iev);
      
      //cout << "nPfos " << nPfos << " nClusters " << nClusters << endl;
      //cout << "pfoE->size() " << pfoE->size() << endl;

      /*
      if (nPfos==1&&nClusters>1) {
	cout << "One PFO with energy " << pfoE->at(0) << endl;
	cout << "RawEnergyOfCluster->size " << RawEnergyOfCluster->size() << endl;
	for (unsigned int icl = 0; icl < RawEnergyOfCluster->size(); icl++){
	  cout << "cluster " << icl << " cluster energy " << RawEnergyOfCluster->at(icl) << endl;
	}
	}*/

      if (nPfos!=1) continue;
      if (nClusters!=1) continue;
      if (pfoE->at(0)<=0) continue;

      //std::cout << "pfoE " << pfoE << std::endl;
      
      bool hasMuonHits = false;
      for (unsigned int ic = 0; ic < HitType->size(); ++ic){
	
	if (HitType->at(ic)==3)
	  hasMuonHits = true;
      }
      if (hasMuonHits) {
	if (m_Debug) std::cout << "Cluster has muon hits" << std::endl;
	continue;
      }
      if (HitEnergies->size()==0) continue;
      
      EventC eventC(Ebeam[ifi],pfoE->at(0), HitEnergies, HitType, CellSize0, CellSize1, CellThickness);
      
      m_EventList.push_back(eventC);
    }
  }

  cout << "m_EventList size " << m_EventList.size() << endl;
  
}

/*--------------- MakeHitDensityHisto ---------------------------*/

void MakeBinDensityHisto_OneEnergy(double cellsize = 30, int fileNumber = 2){
  gStyle->SetOptStat(0);
  
  GetListOfEvents(cellsize,false,fileNumber);

  TH1D *hAll = new TH1D("All", "", 7000, 0, 35);
  TH1D *hAlter = new TH1D("Alter", "", 7000, 0, 35);
  TString hname = ""; 

  for (int ibin = 0; ibin < NBIN; ibin++){
    lowDensity[ibin] = lowDensity[ibin]*0.025/0.3/0.3/0.265;//Thickness was 0.05 before
    highDensity[ibin] = highDensity[ibin]*0.025/0.3/0.3/0.265;
  }

  for (unsigned int iev = 0; iev < m_EventList.size(); iev++){
    EventC event = m_EventList[iev];
    double ebeam = event.Ebeam();
    double ECalEnergy = event.ECalEnergy();
    FloatVector *HCalHitEnergy = event.HCalHitEnergy();

    if (ebeam<=0) continue;
    if (ECalEnergy<0) continue;
    if (HCalHitEnergy->size()<=0) continue;

    for (unsigned int i = 0; i < HCalHitEnergy->size(); i++){
      if (HCalHitEnergy->at(i)<=0) continue;

      int binIndex = -1;
      
      double hitEnergy = HCalHitEnergy->at(i)/0.3/0.3/0.265;

      if (hitEnergy>=lowDensity[9])
	binIndex=9;
      else if (hitEnergy>=lowDensity[8] && hitEnergy<highDensity[8])
	binIndex=8;
      else if (hitEnergy>=lowDensity[7] && hitEnergy<highDensity[7])
	binIndex=7;
      else if (hitEnergy>=lowDensity[6] && hitEnergy<highDensity[6])
	binIndex=6;
      else if (hitEnergy>=lowDensity[5] && hitEnergy<highDensity[5])
	binIndex=5;
      else if (hitEnergy>=lowDensity[4] && hitEnergy<highDensity[4])
	binIndex=4;
      else if (hitEnergy>=lowDensity[3] && hitEnergy<highDensity[3])
	binIndex=3;
      else if (hitEnergy>=lowDensity[2] && hitEnergy<highDensity[2])
	binIndex=2;
      else if (hitEnergy>=lowDensity[1] && hitEnergy<highDensity[1])
	binIndex=1;
      else 
	binIndex=0;

	
      if (binIndex<0) {
	std::cout << "binIndex is wrong!!!" << hitEnergy << std::endl;
	continue;
      }

      hAll->Fill(HCalHitEnergy->at(i)/0.3/0.3/0.265);
      if (binIndex%2==0) hAlter->Fill(HCalHitEnergy->at(i)/0.3/0.3/0.265);

      if ( (HCalHitEnergy->at(i)/0.3/0.3/0.265 < 38) && binIndex==3)
	cout << "HCalHitEnergy->at(i) " << HCalHitEnergy->at(i)
	     << " HCalHitEnergy->at(i)/0.3/0.3/0.265 " << HCalHitEnergy->at(i)/0.3/0.3/0.265 
	     << endl;
    }
  }

  hAll->GetXaxis()->SetTitle("Hit energy density #rho [GeV/1000 cm^{3}]");
  hAll->GetXaxis()->SetTitleSize(0.065);
  hAll->GetXaxis()->SetTitleOffset(1.1);
  hAll->GetXaxis()->SetLabelSize(0.06);
  hAll->GetYaxis()->SetTitle("N_{hits}");
  hAll->GetYaxis()->SetTitleSize(0.065);
  hAll->GetYaxis()->SetTitleOffset(0.9);
  hAll->GetYaxis()->SetLabelSize(0.055);  

  hAll->SetLineColor(kAzure+5);
  hAll->SetFillColor(kAzure+5);
  hAlter->SetLineColor(kAzure-6);
  hAlter->SetFillColor(kAzure-6);  

  TString legname = "K^{0}_{L} and neutrons at "; 
  legname += Ebeam[fileNumber]; legname += " GeV";
  TLegend *leg = new TLegend(0.35, 0.7, 0.88, 0.87,legname);
  leg->SetTextColor(2);
  
  //gPad->SetLogy();
  TCanvas *c = new TCanvas("c" , "", 900, 650);
  gPad->SetBottomMargin(0.18);
  gPad->SetLeftMargin(0.14);
  c->cd();
  hAll->Draw();
  hAlter->Draw("same");
  leg->Draw();
  //for (int ibin = 0; ibin < NBIN; ibin++){
  //hEnergyDensity[ibin]->Draw("same");
  //}

}


void MakeBinDensityHisto(double cellsize = 30){
  gStyle->SetOptStat(0);
  
  GetListOfEvents(cellsize);

  //gStyle->SetOptStat(0);

  TH1F *hECal[nE];
  TProfile *hHCal[nE];
  TH1F *hHCalHitEnergy[nE];
  TH1F *hsumEnergyHCal[nE];
  TString hname = ""; 
  
  for (int index = 0; index < nE; index++){
    hname = "EnergyEcal"; hname += Ebeam[index]; hname += "GeV";
    hECal[index] = new TH1F(hname, hname, 100, 0, 100);
    
    hname = "EnergyHcal"; hname += Ebeam[index]; hname += "GeV";
    hHCal[index] = new TProfile(hname, hname, NBIN, 0, NBIN);  

    hname = "HcalHitEnergy"; hname += Ebeam[index]; hname += "GeV";
    hHCalHitEnergy[index] = new TH1F(hname, hname, 120, 0, 120);  

    hname = "sumEnergyHCal"; hname += Ebeam[index]; hname += "GeV";
    hsumEnergyHCal[index] = new TH1F(hname, hname, 100, 0, 100);
  }

  TLegend *leg = new TLegend(0.6, 0.65, 0.9, 0.85);

  for (unsigned int iev = 0; iev < m_EventList.size(); iev++){
    EventC event = m_EventList[iev];
    double ebeam = event.Ebeam();
    int index = -1;
    for (int i = 0; i < nE; i++){
      if (ebeam == Ebeam[i])
	index = i;
    }

    double ECalEnergy = event.ECalEnergy();
    FloatVector *HCalBinEnergy = event.HCalBinEnergy();
    FloatVector *HCalHitEnergy = event.HCalHitEnergy();

    if (ECalEnergy>0)
      hECal[index]->Fill(ECalEnergy);

    for (unsigned int i = 0; i < HCalHitEnergy->size(); i++){
      //std::cout << "HCalHitEnergy " << HCalHitEnergy->at(i) << std::endl;

      if (HCalHitEnergy->at(i)<=0) continue;
      hHCalHitEnergy[index]->Fill(HCalHitEnergy->at(i)/0.025);
    }
    
    double sumAllEHCal = 0;
    for (int i = 0; i < NBIN; i++){
      hHCal[index]->Fill(i,HCalBinEnergy->at(i));      
      
      //std::cout << "Energy bin " << i << " HCalBinEnergy " << HCalBinEnergy->at(i) << std::endl;
      sumAllEHCal += HCalBinEnergy->at(i);
    }
    sumAllEHCal += ECalEnergy;

    hsumEnergyHCal[index]->Fill(sumAllEHCal);

    if (sumAllEHCal<50) {    
    //if (m_Debug){
    cout << "ECalEnergy " << ECalEnergy << endl;
    cout << "ereco " << event.Ereco() << " sumAllEHCal " << sumAllEHCal << endl;
      //}
    }
  }  
  
  TCanvas *c1 = new TCanvas("c1", "", 1100, 600);
  c1->Divide(2);
  c1->cd(1);
  hECal[4]->SetLineColor(kBlue+2);
  //hECal[4]->SetMaximum(2); 
  hECal[4]->DrawNormalized(); 
  for (int i = 0; i < nE; i++){
    hECal[i]->SetLineWidth(2);
    hECal[i]->SetLineColor(i+1);
    if (i==9)
      hECal[i]->SetLineColor(kOrange+1);
    TString legname = ""; legname += Ebeam[i]; legname += "GeV";
    leg->AddEntry(hECal[i],legname,"l");
    hECal[i]->DrawNormalized("same");      
  }
  c1->cd(2);
  hsumEnergyHCal[4]->SetLineColor(kBlue+2);
  //hsumEnergyHCal[4]->SetMaximum(2); 
  hsumEnergyHCal[4]->DrawNormalized(); 
  for (int i = 0; i < nE; i++){
    hsumEnergyHCal[i]->SetLineWidth(2);
    hECal[i]->SetLineColor(i+1);
    if (i==9)
      hsumEnergyHCal[i]->SetLineColor(kOrange+1);
    TString legname = ""; legname += Ebeam[i]; legname += "GeV";
    hsumEnergyHCal[i]->DrawNormalized("same");      
  }
  leg->Draw();

  TCanvas *c2 = new TCanvas("c2", "", 1100, 600);
  c2->Divide(2);
  c2->cd(1);
  hHCalHitEnergy[0]->SetMinimum(0);
  //hHCalHitEnergy[0]->SetMaximum(500000);
  hHCalHitEnergy[4]->SetLineColor(kBlue+2);
  hHCalHitEnergy[4]->SetFillColor(kAzure-5);
  hHCalHitEnergy[4]->GetXaxis()->SetTitle("Hit energy [MIP]");
  hHCalHitEnergy[4]->GetXaxis()->SetTitleSize(0.045);
  //hHCalHitEnergy[4]->SetMaximum(2); 
  hHCalHitEnergy[4]->DrawNormalized();
  for (int i = 0; i < nE; i++){
    hHCalHitEnergy[i]->SetLineWidth(2);
    hHCalHitEnergy[i]->SetLineColor(i+1);
    if (i==9)
      hHCalHitEnergy[i]->SetLineColor(kOrange+1);
    hHCalHitEnergy[i]->DrawNormalized("same");      
  }
  leg->Draw();
  c2->cd(2);
  hHCal[0]->SetMinimum(0);
  //hHCal[0]->SetMaximum(2);
  //hHCal[0]->DrawNormalized();
  hHCal[4]->SetLineColor(kBlue+2);
  hHCal[4]->SetFillColor(kAzure-5);
  hHCal[4]->GetXaxis()->SetTitle("Bin energy");
  hHCal[4]->GetXaxis()->SetTitleSize(0.045);
  //hHCal[4]->SetMaximum(2); 
  hHCal[4]->DrawNormalized();
  for (int i = 0; i < nE; i++){
    hHCal[i]->SetLineWidth(2);
    hHCal[i]->SetLineColor(i+1);
    if (i==9)
      hHCal[i]->SetLineColor(kOrange+1);
    hHCal[i]->DrawNormalized("same");      
  }
  leg->Draw();

}

/*----------------------------------------------------------------------------------------------------------- Fitting  -----------------------------------------------------------------------------------------------------------*/

double chi2_Minuit2(const double *par){
  bool m_DebugChi2 = false;

  double Chi2 = 0;

  for (unsigned int iev = 0; iev < m_EventList.size(); iev++){
    EventC event = m_EventList[iev];

    double ebeam = event.Ebeam();
    double ereco = event.Ereco();

    double ECalEnergy = event.ECalEnergy();
    FloatVector *HCalBinEnergy = event.HCalBinEnergy();

    if (m_DebugChi2) cout << "ebeam " << ebeam << " ereco " << ereco << endl;
    if (m_DebugChi2) cout << "HCalBinEnergy size " << HCalBinEnergy->size() << endl;

    double E_SC = 0;

    //Add energy from ECal
    E_SC += ECalEnergy;

    if (m_DebugChi2) cout << "ECalEnergy " << ECalEnergy << endl;

    //Sum bins in HCal
    for (unsigned int ibin = 0; ibin < HCalBinEnergy->size(); ibin++){
      double BinEnergy = HCalBinEnergy->at(ibin);
      
      if (BinEnergy<=0) continue;

      if (m_DebugChi2) cout << "ibin " << ibin << endl;
      if (m_DebugChi2) cout << "HCalBinEnergy " << HCalBinEnergy->at(ibin) << endl;
      
      double Weight = 0;
      if (NPAR==3)
	Weight = par[0]*exp(par[1]*binDensity[ibin])+par[2];
      else if (NPAR==9){
	float p1 = par[0] + par[1]*ebeam + par[2]*ebeam*ebeam;
	float p2 = par[3] + par[4]*ebeam + par[5]*ebeam*ebeam;
	float p3 = par[6]/(par[7]+exp(par[8]*ebeam));

	Weight = p1*exp(p2*binDensity[ibin])+p3;
      }

      if (m_DebugChi2) cout << "binDensity " << binDensity[ibin] << " Weight = " << Weight << endl;      
      
      double BinEnergyWeighted = BinEnergy*Weight;

      if (m_DebugChi2) cout << "bin " << ibin << " BinEnergy " << BinEnergy << " Weight " << Weight << " BinEnergyWeighted " << BinEnergyWeighted << endl;

      E_SC += BinEnergyWeighted;
    }

    if (m_DebugChi2) cout << "E_SC " << E_SC << endl;

    double diff = (E_SC - ebeam)*(E_SC - ebeam)/(0.5*ebeam);
    Chi2 += diff;
  }

  return Chi2;
}

/*----------------------------------------------------------------------------------------------------------- Fitting  -----------------------------------------------------------------------------------------------------------*/

void Fit_Minuit2(int cellsize = 30, bool FitWithAllFiles = true, int FileNumber = -1){
  //Make the list of density values
  MakeBinDensity();

  //Get list of events
  GetListOfEvents(cellsize,FitWithAllFiles,FileNumber);

  int npar = NPAR;
  double step = 0.01;

  double fitPar0 = 2.5;
  double fitPar1 = -0.02;
  double fitPar2 = 0.5;
  double fitPar3, fitPar4, fitPar5, fitPar6, fitPar7, fitPar8;

  const char *minName = "Minuit2";
  const char *algoName = "Migrad";
  
  //ROOT::Minuit2::Minuit2Minimizer
  ROOT::Math::Minimizer* myMinuit = ROOT::Math::Factory::CreateMinimizer(minName,algoName);

  // set tolerance , etc...
  myMinuit->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
  myMinuit->SetMaxIterations(1000000);  // for GSL
  myMinuit->SetTolerance(0.001);
  myMinuit->SetPrintLevel(1);
  
  // create function wrapper for minmizer
  ROOT::Math::Functor f(&chi2_Minuit2,npar);

  myMinuit->SetFunction(f);

  if (NPAR==3){
    myMinuit->SetVariable(0, "p0", fitPar0, step);
    myMinuit->SetVariable(1, "p1", fitPar1, step);
    myMinuit->SetVariable(2, "p2", fitPar2, step);
  } else if (NPAR==9){
    fitPar0 = 2.4;
    myMinuit->SetVariable(0, "p0", fitPar0, step);
    fitPar1 = -0.06;
    myMinuit->SetVariable(1, "p1", fitPar1, step);
    fitPar2 = 0.0008;
    myMinuit->SetVariable(2, "p2", fitPar2, step);    
    fitPar3 = -0.09;
    myMinuit->SetVariable(3, "p3", fitPar3, step);    
    fitPar4 = -0.004;
    myMinuit->SetVariable(4, "p4", fitPar4, step);    
    fitPar5 = -0.00008;
    myMinuit->SetVariable(5, "p5", fitPar5, step);    
    fitPar6 = 0.05;
    myMinuit->SetVariable(6, "p6", fitPar6, step);    
    fitPar7 = 0.07;
    myMinuit->SetVariable(7, "p7", fitPar7, step);    
    fitPar8 = -0.1;
    myMinuit->SetVariable(8, "p8", fitPar8, step);    
  }

  // do the minimization
  myMinuit->Minimize();

  const double *xs = myMinuit->X();
  std::cout << "Minimum: chi2 = " << myMinuit->MinValue()  << std::endl;

  // expected minimum is 0
  if ( myMinuit->MinValue()  < 1.E-4  && f(xs) < 1.E-4)
    std::cout << "Minimizer " << minName << " - " << algoName
	      << "   converged to the right minimum" << std::endl;
  else {
    std::cout << "Minimizer " << minName << " - " << algoName
	      << "   failed to converge !!!" << std::endl;
    std::cout << "MinValue " << myMinuit->MinValue() << std::endl;
  }
  
  for (int i = 0; i < npar; i++){
    cout << "Parameter " << i << " is " << xs[i] << endl;
    m_p[i] = xs[i];
  }

  for (int i = 0; i < npar; i++){
    cout << xs[i] << " ";
    //if (i%3==2) cout << endl;
  }
  cout << endl;
}

/*----------------------------------------------------------------------------------------------------------- Fitting  -----------------------------------------------------------------------------------------------------------*/

void chi2_Minuit1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  bool m_DebugChi2 = false;

  double Chi2 = 0;

  for (unsigned int iev = 0; iev < m_EventList.size(); iev++){
    EventC event = m_EventList[iev];

    double ebeam = event.Ebeam();
    double ereco = event.Ereco();

    double ECalEnergy = event.ECalEnergy();
    FloatVector *HCalBinEnergy = event.HCalBinEnergy();

    if (m_DebugChi2) cout << "ebeam " << ebeam << " ereco " << ereco << endl;
    if (m_DebugChi2) cout << "HCalBinEnergy size " << HCalBinEnergy->size() << endl;

    double E_SC = 0;

    //Add energy from ECal
    E_SC += ECalEnergy;

    if (m_DebugChi2) cout << "ECalEnergy " << ECalEnergy << endl;

    //Sum bins in HCal
    for (unsigned int ibin = 0; ibin < HCalBinEnergy->size(); ibin++){
      double BinEnergy = HCalBinEnergy->at(ibin);
      
      //if (BinEnergy<=0) continue;

      if (m_DebugChi2) cout << "ibin " << ibin << endl;
      if (m_DebugChi2) cout << "HCalBinEnergy " << HCalBinEnergy->at(ibin) << endl;
      
      double Weight = 0;
      if (NPAR==3)
	Weight = par[0]*exp(par[1]*binDensity[ibin])+par[2];
      else if (NPAR==9){
	float p1 = par[0] + par[1]*ebeam + par[2]*ebeam*ebeam;
	float p2 = par[3] + par[4]*ebeam + par[5]*ebeam*ebeam;
	float p3 = par[6]/(par[7]+exp(par[8]*ebeam));

	Weight = p1*exp(p2*binDensity[ibin])+p3;
      }

      if (m_DebugChi2) cout << "binDensity " << binDensity[ibin] << " Weight = " << Weight << endl;      
      
      double BinEnergyWeighted = BinEnergy*Weight;

      //cout << "bin " << ibin << " BinEnergy " << BinEnergy << " Weight " << Weight << " BinEnergyWeighted " << BinEnergyWeighted << endl;

      E_SC += BinEnergyWeighted;
    }

    if (m_DebugChi2) cout << "E_SC " << E_SC << endl;
    //cout << "E_SC " << E_SC << endl;

    double diff = (E_SC - ebeam)*(E_SC - ebeam)/(0.5*ebeam);
    Chi2 += diff;
  }

  f = Chi2;
}

/*----------------------------------------------------------------------------------------------------------- Fitting  -----------------------------------------------------------------------------------------------------------*/

void Fit_Minuit1(int cellsize = 30){
  //Make the list of density values
  MakeBinDensity();

  //Get list of events
  GetListOfEvents(cellsize);

  const int npar = NPAR;
  int ifl = 0;
  double step = 0.001;

  double fitPar0 = 0, fitPar1 = 0, fitPar2 = 0, 
    fitPar3 = 0, fitPar4 = 0, fitPar5 = 0, 
    fitPar6 = 0, fitPar7 = 0, fitPar8 = 0;
  double min = -5, max = 5;

  TMinuit *myMinuit = new TMinuit(npar);
  myMinuit->Command("SET NOW");
  myMinuit->Command("SET STR 2");//try to improve minimum search
  myMinuit->SetFCN(chi2_Minuit1);
  //myMinuit->mnparm(0,"p1",fitPar0,step,min0,max0,ifl);
  //myMinuit->mnparm(1,"p2",fitPar1,step,min1,max1,ifl);
  //myMinuit->mnparm(2,"p3",fitPar2,step,min2,max2,ifl);
  fitPar0 = 3;//2.4;
  myMinuit->mnparm(0, "p0", fitPar0, step, min, max, ifl);
  fitPar1 = -0.1;//-0.06;
  myMinuit->mnparm(1, "p1", fitPar1, step, min, max, ifl);
  fitPar2 = 0.08;//0.0008;
  myMinuit->mnparm(2, "p2", fitPar2, step, min, max, ifl);
  fitPar3 = -0.2;//-0.09;
  myMinuit->mnparm(3, "p3", fitPar3, step, min, max, ifl);
  fitPar4 = -0.03;//-0.004;
  myMinuit->mnparm(4, "p4", fitPar4, step, min, max, ifl);
  fitPar5 = -0.00034;//-0.00008;
  myMinuit->mnparm(5, "p5", fitPar5, step, min, max, ifl);
  fitPar6 = 0.01;//0.05;
  myMinuit->mnparm(6, "p6", fitPar6, step, min, max, ifl);
  fitPar7 = 0.02;//0.07;
  myMinuit->mnparm(7, "p7", fitPar7, step, min, max, ifl);
  fitPar8 = -0.001;//-0.1;
  myMinuit->mnparm(8, "p8", fitPar8, step, min, max, ifl);

  myMinuit->Command("MIGRAD 10000 1");

  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  myMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  myMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  double Par0Best = 0, Par1Best = 0, Par2Best = 0;
  double err0Best = 0, err1Best = 0, err2Best = 0;

  double BestFit[npar], BestFitErr[npar];
  for (int i = 0; i < npar; i++){
    myMinuit->GetParameter(i,BestFit[i],BestFitErr[i]);
  }

  for (int i = 0; i < npar; i++){
    cout << BestFit[i] << " ";
  }
  cout << endl;
}


