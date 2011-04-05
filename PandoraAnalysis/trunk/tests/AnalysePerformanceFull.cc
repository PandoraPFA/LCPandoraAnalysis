/**
 *  @file   PandoraAnalysis/tests/AnalysePerformanceFull.cc
 * 
 *  @brief  Implementation of pandora analyse performance full binary.
 * 
 *  $Log: $
 */

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include "AnalysisHelper.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

using namespace pandora_analysis;

void AnalysePerformance(TFile *pTFile, const std::string &outputRootFileName);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    const int nArgs(argc - 1);

    if ((nArgs < 1) || (nArgs > 2))
    {
        std::cout << std::endl
                  << "Usage: ./analysePerformance inputFileName [outputFileName]" << std::endl << std::endl
                  << "  inputFileName  : file containing pandora pfo analysis tree" << std::endl
                  << "  outputFileName : optional output root file, for histogram output" << std::endl << std::endl;
        return 1;
    }

    const std::string inputFileName(argv[1]);
    const std::string outputRootFileName((nArgs == 2) ? argv[2] : "");

    TFile *pTFile = new TFile(inputFileName.c_str(), "READ");
    AnalysePerformance(pTFile, outputRootFileName);
    pTFile->Close();

    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysePerformance(TFile *pTFile, const std::string &outputRootFileName)
{
    // Pfo analysis tree details
    TTree *pTTree = (TTree*)(pTFile->Get("PfoAnalysisTree"));
    const unsigned int nTreeEntries(pTTree->GetEntries());

    int nPfosTotal(0), nPfosNeutralHadrons(0), nPfosPhotons(0), nPfosCharged(0), qPdg(0);
    float pfoEnergyTotal(0.f), pfoEnergyNeutralHadrons(0.f), pfoEnergyPhotons(0.f), pfoEnergyCharged(0.f), pfoMassTotal(0.f),
        mcEnergyTotal(0.f), mcEnergyENu(0.f), mcEnergyFwd(0.f), eQQ(0.f), eQ1(0.f), eQ2(0.f), costQQ(0.f), costQ1(0.f), costQ2(0.f),
        mQQ(0.f), thrust(0.f), netEnergyChange(0.f), sumModulusEnergyChanges(0.f), sumSquaredEnergyChanges(0.f);

    pTTree->SetBranchAddress("nPfosTotal", &nPfosTotal);
    pTTree->SetBranchAddress("nPfosNeutralHadrons", &nPfosNeutralHadrons);
    pTTree->SetBranchAddress("nPfosPhotons", &nPfosPhotons);
    pTTree->SetBranchAddress("nPfosCharged", &nPfosCharged);
    pTTree->SetBranchAddress("pfoEnergyTotal", &pfoEnergyTotal);
    pTTree->SetBranchAddress("pfoEnergyNeutralHadrons", &pfoEnergyNeutralHadrons);
    pTTree->SetBranchAddress("pfoEnergyPhotons", &pfoEnergyPhotons);
    pTTree->SetBranchAddress("pfoEnergyCharged", &pfoEnergyCharged);
    pTTree->SetBranchAddress("pfoMassTotal", &pfoMassTotal);
    pTTree->SetBranchAddress("mcEnergyTotal", &mcEnergyTotal);
    pTTree->SetBranchAddress("mcEnergyENu", &mcEnergyENu);
    pTTree->SetBranchAddress("mcEnergyFwd", &mcEnergyFwd);
    pTTree->SetBranchAddress("eQQ", &eQQ);
    pTTree->SetBranchAddress("eQ1", &eQ1);
    pTTree->SetBranchAddress("eQ2", &eQ2);
    pTTree->SetBranchAddress("costQQ", &costQQ);
    pTTree->SetBranchAddress("costQ1", &costQ1);
    pTTree->SetBranchAddress("costQ2", &costQ2);
    pTTree->SetBranchAddress("mQQ", &mQQ);
    pTTree->SetBranchAddress("thrust", &thrust);
    pTTree->SetBranchAddress("qPdg", &qPdg);
    pTTree->SetBranchAddress("netEnergyChange", &netEnergyChange);
    pTTree->SetBranchAddress("sumModulusEnergyChanges", &sumModulusEnergyChanges);
    pTTree->SetBranchAddress("sumSquaredEnergyChanges", &sumSquaredEnergyChanges);

    // Book histograms
    TH1F *pNPFO = new TH1F("fNPFO", "number of pfos ", 200, 0., 200.);
    TH1F *pPFAMZ = new TH1F("fPFAMZ", "vector boson mass", 200, 50., 150.);
    TH1F *pPFAMW = new TH1F("fPFAMW", "vector boson mass", 200, 50., 150.);
    TH1F *pPFAMZa = new TH1F("fPFAMZa", "vector boson mass", 200, 50., 150.);
    TH1F *pPFAMWa = new TH1F("fPFAMWa", "vector boson mass", 200, 50., 150.);
    TH1F *pPFA = new TH1F("fPFA", "total PFA energy", 10000, 0., 5000.);
    TH1F *pPFAnu = new TH1F("fPFAnu", "total energy + nu", 10000, 0., 5000.);
    TH1F *pPFAnufwd = new TH1F("fPFAnufwd", "total energy + nu + fwd", 10000, 0., 5000.);
    TH1F *pPFAudscb = new TH1F("fPFAudscb", "total energy", 10000, 0., 5000.);
    TH1F *pPFAuds = new TH1F("fPFAuds", "total energy", 10000, 0., 5000.);
    TH1F *pPFAudsHM20 = new TH1F("fPFAudsHM20", "total energy", 10000, 0., 5000.);
    TH1F *pPFAudsHM10 = new TH1F("fPFAudsHM10", "total energy", 10000, 0., 5000.);
    TH1F *pPFAudsHP10 = new TH1F("fPFAudsHP10", "total energy", 10000, 0., 5000.);
    TH1F *pPFAudsHP20 = new TH1F("fPFAudsHP20", "total energy", 10000, 0., 5000.);
    TH1F *pPFAFudsHM20 = new TH1F("fPFAFudsHM20", "total energy",5000, 0., 250.);
    TH1F *pPFAFudsHM10 = new TH1F("fPFAFudsHM10", "total energy",5000, 0., 250.);
    TH1F *pPFAFudsHP10 = new TH1F("fPFAFudsHP10", "total energy",5000, 0., 250.);
    TH1F *pPFAFudsHP20 = new TH1F("fPFAFudsHP20", "total energy",5000, 0., 250.);
    TH1F *pPFAcb = new TH1F("fPFAcb", "total energy", 10000, 0., 5000.);
    TH1F *pPFA1 = new TH1F("fPFA1", "total energy 0.0-0.1", 10000, 0., 5000.);
    TH1F *pPFA2 = new TH1F("fPFA2", "total energy 0.1-0.2", 10000, 0., 5000.);
    TH1F *pPFA3 = new TH1F("fPFA3", "total energy 0.2-0.3", 10000, 0., 5000.);
    TH1F *pPFA4 = new TH1F("fPFA4", "total energy 0.3-0.4", 10000, 0., 5000.);
    TH1F *pPFA5 = new TH1F("fPFA5", "total energy 0.4-0.5", 10000, 0., 5000.);
    TH1F *pPFA6 = new TH1F("fPFA6", "total energy 0.5-0.6", 10000, 0., 5000.);
    TH1F *pPFA7 = new TH1F("fPFA7", "total energy 0.6-0.7", 10000, 0., 5000.);
    TH1F *pPFAL7A = new TH1F("fPFAL7A", "total energy <0.7", 10000, 0., 5000.);
    TH1F *pPFAL7Aud= new TH1F("fPFAL7Aud", "total energy <0.7 ud", 10000, 0., 5000.);
    TH1F *pPFAL7As = new TH1F("fPFAL7As", "total energy <0.7 s", 10000, 0., 5000.);
    TH1F *pPFAL7Ac = new TH1F("fPFAL7Ac", "total energy <0.7 c", 10000, 0., 5000.);
    TH1F *pPFAL7Ab = new TH1F("fPFAL7Ab", "total energy <0.7 b", 10000, 0., 5000.);
    TH1F *pPFAL7B = new TH1F("fPFAL7B", "total energy <0.7 - bad", 10000, 0., 5000.);
    TH1F *pPFAFL7A = new TH1F("fPFAFL7A", "total energy <0.7",5000, 0., 250.);
    TH1F *pPFAFL7Aud = new TH1F("fPFAFL7Aud", "total energy <0.7 ud",5000, 0., 250.);
    TH1F *pPFAFL7As = new TH1F("fPFAFL7As", "total energy <0.7 s",5000, 0., 250.);
    TH1F *pPFAFL7Ac = new TH1F("fPFAFL7Ac", "total energy <0.7 c",5000, 0., 250.);
    TH1F *pPFAFL7Ab = new TH1F("fPFAFL7Ab", "total energy <0.7 b",5000, 0., 250.);
    TH1F *pPFAQQ   = new TH1F("fPFAQQ", "total energy - true E", 2000, -500., 500.);
    TH1F *pPFAQQ8  = new TH1F("fPFAQQ8", "total energy - true E 8", 2000, -500., 500.);
    TH1F *pPFA8  = new TH1F("fPFA8", "total energy 0.7-0.8", 10000, 0., 5000.);
    TH1F *pPFA9  = new TH1F("fPFA9", "total energy 0.8-0.9", 10000, 0., 5000.);
    TH1F *pPFA10 = new TH1F("fPFA10","total energy 0.9-1.0", 10000, 0., 5000.);
    TH1F *pPFA11 = new TH1F("fPFA11","total energy 0.9-0.925", 10000, 0., 5000.);
    TH1F *pPFA12 = new TH1F("fPFA12","total energy 0.925-0.95", 10000, 0., 5000.);
    TH1F *pPFA13 = new TH1F("fPFA13","total energy 0.95-0.975", 10000, 0., 5000.);
    TH1F *pPFA14 = new TH1F("fPFA14","total energy 0.975-1.0", 10000, 0., 5000.);
    TH1F *pPFADMZ = new TH1F("fPFADMZ", "Delta Mz", 200, -50., 50.);
    TH1F *pPFADMZ8 = new TH1F("fPFADMZ8", "Delta Mz", 200, -50., 50.);
    TH1F *pPFADMZQQ8 = new TH1F("fPFADMZQQ8", "Delta Mz", 200, -50., 50.);
    TH1F *pPFADMZP8 = new TH1F("fPFADMZP8", "Delta Mz", 200, -50., 50.);
    TH1F *pPFADMZOMZ = new TH1F("fPFADMZOMZ" , "Delta Mz / Mz", 200, -25., 25.);
    TH1F *pPFADMZOMZQQ8 = new TH1F("fPFADMZOMZQQ8" , "Delta Mz / Mz", 200, -25., 25.);

    // Fill histograms with event information
    for (unsigned int iTree = 0; iTree < nTreeEntries; ++iTree)
    {
        pTTree->GetEntry(iTree);

        pNPFO->Fill(nPfosTotal);

        if (std::fabs(costQQ) < 0.9f)
        {
            pPFAMWa->Fill(pfoMassTotal);
            pPFAMZa->Fill(pfoMassTotal);

            if (std::fabs(mQQ - 80.3) < 4.f)
                pPFAMW->Fill(pfoMassTotal);

            if (std::fabs(mQQ - 91.2f) < 4.f)
                pPFAMZ->Fill(pfoMassTotal);

            if (mQQ > 0.1f)
            {
                pPFAQQ->Fill(pfoEnergyTotal - eQQ);
                pPFADMZ->Fill(pfoMassTotal - mQQ);

                if (mQQ > 75.f)
                    pPFADMZOMZ->Fill(100.f * (pfoMassTotal - mQQ) / pfoMassTotal);

                if (std::fabs(costQQ) < 0.8f)
                    pPFADMZ8->Fill(pfoMassTotal - mQQ);

                if (std::fabs(costQQ) > 0.8f)
                    pPFADMZP8->Fill(pfoMassTotal - mQQ);

                if (std::fabs(costQ1) < 0.8f && std::fabs(costQ2) < 0.8f)
                {
                    pPFADMZQQ8->Fill(pfoMassTotal - mQQ);
                    pPFADMZOMZQQ8->Fill(100.f * (pfoMassTotal - mQQ) / pfoMassTotal);
                    pPFAQQ8->Fill(pfoEnergyTotal - eQQ);
                }
            }
        }

        pPFA->Fill(pfoEnergyTotal);
        pPFAnu->Fill(pfoEnergyTotal + mcEnergyENu);
        pPFAnufwd->Fill(pfoEnergyTotal + mcEnergyENu + mcEnergyFwd);

        if (qPdg >= 1 && qPdg <= 3)
            pPFAuds->Fill(pfoEnergyTotal);

        if (qPdg >= 4 && qPdg <= 5)
            pPFAcb->Fill(pfoEnergyTotal);

        if (qPdg >= 1 && qPdg <= 5)
            pPFAudscb->Fill(pfoEnergyTotal);

        if (qPdg >= 1 && qPdg <= 3)
        {
            if (thrust <= 0.1f)
                pPFA1->Fill(pfoEnergyTotal);

            if (thrust > 0.1f && thrust <= 0.2f)
                pPFA2->Fill(pfoEnergyTotal);

            if (thrust > 0.2f && thrust <= 0.3f)
                pPFA3->Fill(pfoEnergyTotal);

            if (thrust > 0.3f && thrust <= 0.4f)
                pPFA4->Fill(pfoEnergyTotal);

            if (thrust > 0.4f && thrust <= 0.5f)
                pPFA5->Fill(pfoEnergyTotal);

            if (thrust > 0.5f && thrust <= 0.6f)
                pPFA6->Fill(pfoEnergyTotal);

            if (thrust <= 0.7f)
            {
                pPFAL7A->Fill(pfoEnergyTotal + mcEnergyENu);
                pPFAFL7A->Fill(pfoEnergyTotal + mcEnergyENu);

                if (qPdg <= 2)
                    pPFAL7Aud->Fill(pfoEnergyTotal + mcEnergyENu);

                if (qPdg == 3)
                    pPFAL7As->Fill(pfoEnergyTotal + mcEnergyENu);

                if (qPdg == 4)
                    pPFAL7Ac->Fill(pfoEnergyTotal + mcEnergyENu);

                if (qPdg == 5)
                    pPFAL7Ab->Fill(pfoEnergyTotal + mcEnergyENu);

                if (qPdg <= 2)
                    pPFAFL7Aud->Fill(pfoEnergyTotal + mcEnergyENu);

                if (qPdg == 3)
                    pPFAFL7As->Fill(pfoEnergyTotal + mcEnergyENu);

                if (qPdg == 4)
                    pPFAFL7Ac->Fill(pfoEnergyTotal + mcEnergyENu);

                if (qPdg == 5)
                    pPFAFL7Ab->Fill(pfoEnergyTotal + mcEnergyENu);

                pPFAudsHM20->Fill(pfoEnergyTotal + mcEnergyENu - pfoEnergyNeutralHadrons * 0.1f);
                pPFAudsHM10->Fill(pfoEnergyTotal + mcEnergyENu - pfoEnergyNeutralHadrons * 0.05f);
                pPFAudsHP10->Fill(pfoEnergyTotal + mcEnergyENu + pfoEnergyNeutralHadrons * 0.05f);
                pPFAudsHP20->Fill(pfoEnergyTotal + mcEnergyENu + pfoEnergyNeutralHadrons * 0.1f);

                pPFAFudsHM20->Fill(pfoEnergyTotal + mcEnergyENu  -pfoEnergyNeutralHadrons * 0.1f);
                pPFAFudsHM10->Fill(pfoEnergyTotal + mcEnergyENu - pfoEnergyNeutralHadrons * 0.05f);
                pPFAFudsHP10->Fill(pfoEnergyTotal + mcEnergyENu + pfoEnergyNeutralHadrons * 0.05f);
                pPFAFudsHP20->Fill(pfoEnergyTotal + mcEnergyENu + pfoEnergyNeutralHadrons * 0.1f);
            }

            if (thrust <= 0.7f && mcEnergyFwd / pfoEnergyTotal < 0.01f)
                pPFAL7B->Fill(pfoEnergyTotal + mcEnergyENu);

            if (thrust > 0.6f && thrust <= 0.7f)
                pPFA7->Fill(pfoEnergyTotal);

            if (thrust > 0.7f && thrust <= 0.8f)
                pPFA8->Fill(pfoEnergyTotal);

            if (thrust > 0.8f && thrust <= 0.9f)
                pPFA9->Fill(pfoEnergyTotal);

            if (thrust > 0.9f)
                pPFA10->Fill(pfoEnergyTotal);

            if (thrust > 0.9f  && thrust <= 0.925f)
                pPFA11->Fill(pfoEnergyTotal);

            if (thrust > 0.925f && thrust <= 0.95f)
                pPFA12->Fill(pfoEnergyTotal);

            if (thrust > 0.95f && thrust <= 0.975f)
                pPFA13->Fill(pfoEnergyTotal);

            if (thrust > 0.975f)
                pPFA14->Fill(pfoEnergyTotal);
        }
    }

    // Examine histograms and investigate variation of resolution with cos(theta)
    float cosThetaBins[14] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 1.0};
    TH1F *pEvsCHist = new TH1F("SigmaEvsCosTheta", "#sigma(E) vs cos(#theta)", 13, cosThetaBins);
    pEvsCHist->SetYTitle("#sigma(E)");
    pEvsCHist->SetXTitle("cos(#theta)");

    float sigma(0.f), sigmasigma(0.f);
    AnalysisHelper::CalculatePerformance(pPFA,  sigma, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA1, sigma, sigmasigma);  pEvsCHist->SetBinContent( 1, sigma); pEvsCHist->SetBinError( 1, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA2, sigma, sigmasigma);  pEvsCHist->SetBinContent( 2, sigma); pEvsCHist->SetBinError( 2, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA3, sigma, sigmasigma);  pEvsCHist->SetBinContent( 3, sigma); pEvsCHist->SetBinError( 3, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA4, sigma, sigmasigma);  pEvsCHist->SetBinContent( 4, sigma); pEvsCHist->SetBinError( 4, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA5, sigma, sigmasigma);  pEvsCHist->SetBinContent( 5, sigma); pEvsCHist->SetBinError( 5, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA6, sigma, sigmasigma);  pEvsCHist->SetBinContent( 6, sigma); pEvsCHist->SetBinError( 6, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA7, sigma, sigmasigma);  pEvsCHist->SetBinContent( 7, sigma); pEvsCHist->SetBinError( 7, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA8, sigma, sigmasigma);  pEvsCHist->SetBinContent( 8, sigma); pEvsCHist->SetBinError( 8, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA9, sigma, sigmasigma);  pEvsCHist->SetBinContent( 9, sigma); pEvsCHist->SetBinError( 9, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA11, sigma, sigmasigma); pEvsCHist->SetBinContent(10, sigma); pEvsCHist->SetBinError(10, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA12, sigma, sigmasigma); pEvsCHist->SetBinContent(11, sigma); pEvsCHist->SetBinError(11, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA13, sigma, sigmasigma); pEvsCHist->SetBinContent(12, sigma); pEvsCHist->SetBinError(12, sigmasigma);
    AnalysisHelper::CalculatePerformance(pPFA14, sigma, sigmasigma); pEvsCHist->SetBinContent(13, sigma); pEvsCHist->SetBinError(13, sigmasigma);

    AnalysisHelper::CalculatePerformance(pPFAudsHM20);
    AnalysisHelper::CalculatePerformance(pPFAudsHM10);
    AnalysisHelper::CalculatePerformance(pPFAuds);
    AnalysisHelper::CalculatePerformance(pPFAudsHP10);
    AnalysisHelper::CalculatePerformance(pPFAudsHP20);

    TH1F *pPFA1Clone = static_cast<TH1F*>(pPFA1->Clone());
    pPFA1Clone->Add(pPFA2);
    pPFA1Clone->Add(pPFA3);
    pPFA1Clone->Add(pPFA4);
    pPFA1Clone->Add(pPFA5);  std::cout << " < 0.5 : ";   AnalysisHelper::CalculatePerformance(pPFA1Clone);
    pPFA1Clone->Add(pPFA6);  std::cout << " < 0.6 : ";   AnalysisHelper::CalculatePerformance(pPFA1Clone);
    pPFA1Clone->Add(pPFA7);  std::cout << " < 0.7 : ";   AnalysisHelper::CalculatePerformance(pPFA1Clone);
    pPFA1Clone->Add(pPFA8);  std::cout << " < 0.8 : ";   AnalysisHelper::CalculatePerformance(pPFA1Clone);
    pPFA1Clone->Add(pPFA9);  std::cout << " < 0.9 : ";   AnalysisHelper::CalculatePerformance(pPFA1Clone);
    pPFA1Clone->Add(pPFA11); std::cout << " < 0.925 : "; AnalysisHelper::CalculatePerformance(pPFA1Clone);
    pPFA1Clone->Add(pPFA12); std::cout << " < 0.95 : ";  AnalysisHelper::CalculatePerformance(pPFA1Clone);
    pPFA1Clone->Add(pPFA13); std::cout << " < 0.975 : "; AnalysisHelper::CalculatePerformance(pPFA1Clone);
    pPFA1Clone->Add(pPFA14); std::cout << " < 1.0 : ";   AnalysisHelper::CalculatePerformance(pPFA1Clone);

    TH1F *pPFA9Clone = static_cast<TH1F*>(pPFA9->Clone());
    pPFA9Clone->Add(pPFA11);
    pPFA9Clone->Add(pPFA12); std::cout << " 0.8-0.95 : "; AnalysisHelper::CalculatePerformance(pPFA9Clone);
    std::cout << " < 0.7 A : ";    AnalysisHelper::CalculatePerformance(pPFAL7A);
    std::cout << " < 0.7 A ud : "; AnalysisHelper::CalculatePerformance(pPFAL7Aud);
    std::cout << " < 0.7 A s : ";  AnalysisHelper::CalculatePerformance(pPFAL7As);
    std::cout << " < 0.7 B : ";    AnalysisHelper::CalculatePerformance(pPFAL7B);

    AnalysisHelper::CalculatePerformance(pPFAMZ);
    AnalysisHelper::CalculatePerformance(pPFAMW);
    AnalysisHelper::CalculatePerformance(pPFAMZa);
    AnalysisHelper::CalculatePerformance(pPFAMWa);
    AnalysisHelper::CalculatePerformance(pPFAnu);
    AnalysisHelper::CalculatePerformance(pPFAnufwd);
    AnalysisHelper::CalculatePerformance(pPFAudscb);
    AnalysisHelper::CalculatePerformance(pPFAcb);
    AnalysisHelper::CalculatePerformance(pPFAQQ, false);
    AnalysisHelper::CalculatePerformance(pPFAQQ8, false);
    AnalysisHelper::CalculatePerformance(pPFADMZQQ8, false);
    AnalysisHelper::CalculatePerformance(pPFADMZOMZQQ8, false);
    AnalysisHelper::CalculatePerformance(pPFADMZ, false);
    AnalysisHelper::CalculatePerformance(pPFADMZOMZ, false);
    AnalysisHelper::CalculatePerformance(pPFADMZ8, false);
    AnalysisHelper::CalculatePerformance(pPFADMZP8, false);

    // Write out histograms
    if (!outputRootFileName.empty())
    {
        std::cout << "Will write histograms to file : " << outputRootFileName << std::endl;
        TFile *pTOutputFile = new TFile(outputRootFileName.c_str(), "RECREATE");
        pTOutputFile->cd();
        pEvsCHist->Write();
        pNPFO->Write();
        pPFAMZ->Write();
        pPFAMW->Write();
        pPFAMZa->Write();
        pPFAMWa->Write();
        pPFA->Write();
        pPFAnu->Write();
        pPFAnufwd->Write();
        pPFAudscb->Write();
        pPFAuds->Write();
        pPFAudsHM20->Write();
        pPFAudsHM10->Write();
        pPFAudsHP10->Write();
        pPFAudsHP20->Write();
        pPFAFudsHM20->Write();
        pPFAFudsHM10->Write();
        pPFAFudsHP10->Write();
        pPFAFudsHP20->Write();
        pPFAcb->Write();
        pPFA1->Write();
        pPFA2->Write();
        pPFA3->Write();
        pPFA4->Write();
        pPFA5->Write();
        pPFA6->Write();
        pPFA7->Write();
        pPFAL7A->Write();
        pPFAL7Aud->Write();
        pPFAL7As->Write();
        pPFAL7Ac->Write();
        pPFAL7Ab->Write();
        pPFAL7B->Write();
        pPFAFL7A->Write();
        pPFAFL7Aud->Write();
        pPFAFL7As->Write();
        pPFAFL7Ac->Write();
        pPFAFL7Ab->Write();
        pPFAQQ->Write();
        pPFAQQ8->Write();
        pPFA8->Write();
        pPFA9->Write();
        pPFA10->Write();
        pPFA11->Write();
        pPFA12->Write();
        pPFA13->Write();
        pPFA14->Write();
        pPFADMZ->Write();
        pPFADMZ8->Write();
        pPFADMZQQ8->Write();
        pPFADMZP8->Write();
        pPFADMZOMZ->Write();
        pPFADMZOMZQQ8->Write();
        pTOutputFile->Close();
        delete pTOutputFile;
    }

    // Tidy up
    delete pNPFO;
    delete pPFAMZ;
    delete pPFAMW;
    delete pPFAMZa;
    delete pPFAMWa;
    delete pPFA;
    delete pPFAnu;
    delete pPFAnufwd;
    delete pPFAudscb;
    delete pPFAuds;
    delete pPFAudsHM20;
    delete pPFAudsHM10;
    delete pPFAudsHP10;
    delete pPFAudsHP20;
    delete pPFAFudsHM20;
    delete pPFAFudsHM10;
    delete pPFAFudsHP10;
    delete pPFAFudsHP20;
    delete pPFAcb;
    delete pPFA1;
    delete pPFA2;
    delete pPFA3;
    delete pPFA4;
    delete pPFA5;
    delete pPFA6;
    delete pPFA7;
    delete pPFAL7A;
    delete pPFAL7Aud;
    delete pPFAL7As;
    delete pPFAL7Ac;
    delete pPFAL7Ab;
    delete pPFAL7B;
    delete pPFAFL7A;
    delete pPFAFL7Aud;
    delete pPFAFL7As;
    delete pPFAFL7Ac;
    delete pPFAFL7Ab;
    delete pPFAQQ;
    delete pPFAQQ8;
    delete pPFA8;
    delete pPFA9;
    delete pPFA10;
    delete pPFA11;
    delete pPFA12;
    delete pPFA13;
    delete pPFA14;
    delete pPFADMZ;
    delete pPFADMZ8;
    delete pPFADMZQQ8;
    delete pPFADMZP8;
    delete pPFADMZOMZ;
    delete pPFADMZOMZQQ8;
    delete pPFA1Clone;
    delete pPFA9Clone;
    delete pEvsCHist;
}
