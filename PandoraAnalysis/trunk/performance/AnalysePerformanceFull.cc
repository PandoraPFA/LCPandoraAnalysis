/**
 *  @file   PandoraAnalysis/performance/AnalysePerformanceFull.cc
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
    try
    {
        const int nArgs(argc - 1);

        if ((nArgs < 1) || (nArgs > 2))
        {
            std::cout << std::endl
                      << "Usage: ./AnalysePerformanceFull inputFileName [outputFileName]" << std::endl << std::endl
                      << "  inputFileName  : file containing pandora pfo analysis tree" << std::endl
                      << "  outputFileName : optional output root file, for histogram output" << std::endl << std::endl;
            return 1;
        }

        const std::string inputFileName(argv[1]);
        const std::string outputRootFileName((nArgs == 2) ? argv[2] : "");

        TFile *pTFile = new TFile(inputFileName.c_str(), "READ");

        if (pTFile->IsZombie())
        {
            std::cout << "Error opening file " << inputFileName << std::endl;
            delete pTFile;
            return 1;
        }

        AnalysePerformance(pTFile, outputRootFileName);
        pTFile->Close();
    }
    catch (std::exception &exception)
    {
        std::cout << "Exception caught " << exception.what() << std::endl;
        return 1;
    }

    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysePerformance(TFile *pTFile, const std::string &outputRootFileName)
{
    // Pfo analysis tree details
    TTree *pTTree = NULL;
    pTFile->GetObject("PfoAnalysisTree", pTTree);

    if (!pTTree)
    {
        std::cout << "Error opening root tree: PfoAnalysisTree " << std::endl;
        return;
    }

    const unsigned int nTreeEntries(pTTree->GetEntries());

    int nPfosTotal(0), nPfosNeutralHadrons(0), nPfosPhotons(0), nPfosTracks(0), qPdg(0);
    float pfoEnergyTotal(0.f), pfoEnergyNeutralHadrons(0.f), pfoEnergyPhotons(0.f), pfoEnergyTracks(0.f), pfoMassTotal(0.f),
        mcEnergyENu(0.f), mcEnergyFwd(0.f), eQQ(0.f), eQ1(0.f), eQ2(0.f), costQQ(0.f), costQ1(0.f), costQ2(0.f), mQQ(0.f), thrust(0.f);

    pTTree->SetBranchAddress("nPfosTotal", &nPfosTotal);
    pTTree->SetBranchAddress("nPfosNeutralHadrons", &nPfosNeutralHadrons);
    pTTree->SetBranchAddress("nPfosPhotons", &nPfosPhotons);
    pTTree->SetBranchAddress("nPfosTracks", &nPfosTracks);
    pTTree->SetBranchAddress("pfoEnergyTotal", &pfoEnergyTotal);
    pTTree->SetBranchAddress("pfoEnergyNeutralHadrons", &pfoEnergyNeutralHadrons);
    pTTree->SetBranchAddress("pfoEnergyPhotons", &pfoEnergyPhotons);
    pTTree->SetBranchAddress("pfoEnergyTracks", &pfoEnergyTracks);
    pTTree->SetBranchAddress("pfoMassTotal", &pfoMassTotal);
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

    // Book histograms
    TH1F *pNPFO = new TH1F("fNPFO", "number of pfos ", 200, 0., 200.);
    TH1F *pPFA_MZ = new TH1F("fPFA_MZ", "vector boson mass", 200, 50., 150.);
    TH1F *pPFA_MW = new TH1F("fPFA_MW", "vector boson mass", 200, 50., 150.);
    TH1F *pPFA_MZa = new TH1F("fPFA_MZa", "vector boson mass", 200, 50., 150.);
    TH1F *pPFA_MWa = new TH1F("fPFA_MWa", "vector boson mass", 200, 50., 150.);
    TH1F *pPFA = new TH1F("fPFA", "total PFA energy", 10000, 0., 5000.);
    TH1F *pPFA_nu = new TH1F("fPFA_nu", "total energy + nu", 10000, 0., 5000.);
    TH1F *pPFA_nufwd = new TH1F("fPFA_nufwd", "total energy + nu + fwd", 10000, 0., 5000.);
    TH1F *pPFA_udscb = new TH1F("fPFA_udscb", "total energy", 10000, 0., 5000.);
    TH1F *pPFA_uds = new TH1F("fPFA_uds", "total energy", 10000, 0., 5000.);
    TH1F *pPFA_udsHM20 = new TH1F("fPFA_udsHM20", "total energy", 10000, 0., 5000.);
    TH1F *pPFA_udsHM10 = new TH1F("fPFA_udsHM10", "total energy", 10000, 0., 5000.);
    TH1F *pPFA_udsHP10 = new TH1F("fPFA_udsHP10", "total energy", 10000, 0., 5000.);
    TH1F *pPFA_udsHP20 = new TH1F("fPFA_udsHP20", "total energy", 10000, 0., 5000.);
    TH1F *pPFA_cb = new TH1F("fPFA_cb", "total energy", 10000, 0., 5000.);
    TH1F *pPFA_1 = new TH1F("fPFA_1", "total energy 0.0-0.1", 10000, 0., 5000.);
    TH1F *pPFA_2 = new TH1F("fPFA_2", "total energy 0.1-0.2", 10000, 0., 5000.);
    TH1F *pPFA_3 = new TH1F("fPFA_3", "total energy 0.2-0.3", 10000, 0., 5000.);
    TH1F *pPFA_4 = new TH1F("fPFA_4", "total energy 0.3-0.4", 10000, 0., 5000.);
    TH1F *pPFA_5 = new TH1F("fPFA_5", "total energy 0.4-0.5", 10000, 0., 5000.);
    TH1F *pPFA_6 = new TH1F("fPFA_6", "total energy 0.5-0.6", 10000, 0., 5000.);
    TH1F *pPFA_7 = new TH1F("fPFA_7", "total energy 0.6-0.7", 10000, 0., 5000.);
    TH1F *pPFA_L7A = new TH1F("fPFA_L7A", "total energy <0.7", 10000, 0., 5000.);
    TH1F *pPFA_L7Aud= new TH1F("fPFA_L7Aud", "total energy <0.7 ud", 10000, 0., 5000.);
    TH1F *pPFA_L7As = new TH1F("fPFA_L7As", "total energy <0.7 s", 10000, 0., 5000.);
    TH1F *pPFA_L7Ac = new TH1F("fPFA_L7Ac", "total energy <0.7 c", 10000, 0., 5000.);
    TH1F *pPFA_L7Ab = new TH1F("fPFA_L7Ab", "total energy <0.7 b", 10000, 0., 5000.);
    TH1F *pPFA_L7ANoFwd = new TH1F("fPFA_L7ANoFwd", "total energy <0.7 - bad", 10000, 0., 5000.);
    TH1F *pPFA_QQ   = new TH1F("fPFA_QQ", "total energy - true E", 2000, -500., 500.);
    TH1F *pPFA_QQ8  = new TH1F("fPFA_QQ8", "total energy - true E 8", 2000, -500., 500.);
    TH1F *pPFA_8  = new TH1F("fPFA_8", "total energy 0.7-0.8", 10000, 0., 5000.);
    TH1F *pPFA_9  = new TH1F("fPFA_9", "total energy 0.8-0.9", 10000, 0., 5000.);
    TH1F *pPFA_10 = new TH1F("fPFA_10","total energy 0.9-1.0", 10000, 0., 5000.);
    TH1F *pPFA_11 = new TH1F("fPFA_11","total energy 0.9-0.925", 10000, 0., 5000.);
    TH1F *pPFA_12 = new TH1F("fPFA_12","total energy 0.925-0.95", 10000, 0., 5000.);
    TH1F *pPFA_13 = new TH1F("fPFA_13","total energy 0.95-0.975", 10000, 0., 5000.);
    TH1F *pPFA_14 = new TH1F("fPFA_14","total energy 0.975-1.0", 10000, 0., 5000.);
    TH1F *pPFA_DMZ = new TH1F("fPFA_DMZ", "Delta Mz", 200, -50., 50.);
    TH1F *pPFA_DMZ8 = new TH1F("fPFA_DMZ8", "Delta Mz", 200, -50., 50.);
    TH1F *pPFA_DMZQQ8 = new TH1F("fPFA_DMZQQ8", "Delta Mz", 200, -50., 50.);
    TH1F *pPFA_DMZP8 = new TH1F("fPFA_DMZP8", "Delta Mz", 200, -50., 50.);
    TH1F *pPFA_DMZOMZ = new TH1F("fPFA_DMZOMZ" , "Delta Mz / Mz", 200, -25., 25.);
    TH1F *pPFA_DMZOMZQQ8 = new TH1F("fPFA_DMZOMZQQ8" , "Delta Mz / Mz", 200, -25., 25.);

    // Fill histograms with event information
    for (unsigned int iTree = 0; iTree < nTreeEntries; ++iTree)
    {
        pTTree->GetEntry(iTree);

        pNPFO->Fill(nPfosTotal);
        pPFA->Fill(pfoEnergyTotal);
        pPFA_nu->Fill(pfoEnergyTotal + mcEnergyENu);
        pPFA_nufwd->Fill(pfoEnergyTotal + mcEnergyENu + mcEnergyFwd);

        if (std::fabs(costQQ) < 0.9f)
        {
            pPFA_MWa->Fill(pfoMassTotal);
            pPFA_MZa->Fill(pfoMassTotal);

            if (std::fabs(mQQ - 80.3) < 4.f)
                pPFA_MW->Fill(pfoMassTotal);

            if (std::fabs(mQQ - 91.2f) < 4.f)
                pPFA_MZ->Fill(pfoMassTotal);

            if (mQQ > 0.1f)
            {
                pPFA_QQ->Fill(pfoEnergyTotal - eQQ);
                pPFA_DMZ->Fill(pfoMassTotal - mQQ);

                if (mQQ > 75.f)
                    pPFA_DMZOMZ->Fill(100.f * (pfoMassTotal - mQQ) / pfoMassTotal);

                if (std::fabs(costQQ) < 0.8f)
                    pPFA_DMZ8->Fill(pfoMassTotal - mQQ);

                if (std::fabs(costQQ) > 0.8f)
                    pPFA_DMZP8->Fill(pfoMassTotal - mQQ);

                if (std::fabs(costQ1) < 0.8f && std::fabs(costQ2) < 0.8f)
                {
                    pPFA_DMZQQ8->Fill(pfoMassTotal - mQQ);
                    pPFA_DMZOMZQQ8->Fill(100.f * (pfoMassTotal - mQQ) / pfoMassTotal);
                    pPFA_QQ8->Fill(pfoEnergyTotal - eQQ);
                }
            }
        }

        if (qPdg >= 1 && qPdg <= 3)
            pPFA_uds->Fill(pfoEnergyTotal);

        if (qPdg >= 4 && qPdg <= 5)
            pPFA_cb->Fill(pfoEnergyTotal);

        if (qPdg >= 1 && qPdg <= 5)
            pPFA_udscb->Fill(pfoEnergyTotal);

        if (thrust <= 0.7f)
        {
            if (qPdg <= 2)
                pPFA_L7Aud->Fill(pfoEnergyTotal + mcEnergyENu);

            if (qPdg == 3)
                pPFA_L7As->Fill(pfoEnergyTotal + mcEnergyENu);

            if (qPdg == 4)
                pPFA_L7Ac->Fill(pfoEnergyTotal + mcEnergyENu);

            if (qPdg == 5)
                pPFA_L7Ab->Fill(pfoEnergyTotal + mcEnergyENu);
        }

        if (qPdg >= 1 && qPdg <= 3)
        {
            if (thrust <= 0.1f)
                pPFA_1->Fill(pfoEnergyTotal);

            if (thrust > 0.1f && thrust <= 0.2f)
                pPFA_2->Fill(pfoEnergyTotal);

            if (thrust > 0.2f && thrust <= 0.3f)
                pPFA_3->Fill(pfoEnergyTotal);

            if (thrust > 0.3f && thrust <= 0.4f)
                pPFA_4->Fill(pfoEnergyTotal);

            if (thrust > 0.4f && thrust <= 0.5f)
                pPFA_5->Fill(pfoEnergyTotal);

            if (thrust > 0.5f && thrust <= 0.6f)
                pPFA_6->Fill(pfoEnergyTotal);

            if (thrust <= 0.7f)
            {
                pPFA_L7A->Fill(pfoEnergyTotal + mcEnergyENu);

                if (mcEnergyFwd / pfoEnergyTotal < 0.01f)
                    pPFA_L7ANoFwd->Fill(pfoEnergyTotal + mcEnergyENu);

                pPFA_udsHM20->Fill(pfoEnergyTotal + mcEnergyENu - pfoEnergyNeutralHadrons * 0.1f);
                pPFA_udsHM10->Fill(pfoEnergyTotal + mcEnergyENu - pfoEnergyNeutralHadrons * 0.05f);
                pPFA_udsHP10->Fill(pfoEnergyTotal + mcEnergyENu + pfoEnergyNeutralHadrons * 0.05f);
                pPFA_udsHP20->Fill(pfoEnergyTotal + mcEnergyENu + pfoEnergyNeutralHadrons * 0.1f);
            }

            if (thrust > 0.6f && thrust <= 0.7f)
                pPFA_7->Fill(pfoEnergyTotal);

            if (thrust > 0.7f && thrust <= 0.8f)
                pPFA_8->Fill(pfoEnergyTotal);

            if (thrust > 0.8f && thrust <= 0.9f)
                pPFA_9->Fill(pfoEnergyTotal);

            if (thrust > 0.9f)
                pPFA_10->Fill(pfoEnergyTotal);

            if (thrust > 0.9f  && thrust <= 0.925f)
                pPFA_11->Fill(pfoEnergyTotal);

            if (thrust > 0.925f && thrust <= 0.95f)
                pPFA_12->Fill(pfoEnergyTotal);

            if (thrust > 0.95f && thrust <= 0.975f)
                pPFA_13->Fill(pfoEnergyTotal);

            if (thrust > 0.975f)
                pPFA_14->Fill(pfoEnergyTotal);
        }
    }

    // Examine histograms and investigate variation of resolution with cos(theta)
    float cosThetaBins[14] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 1.0};
    TH1F *pResVsCosThetaHist = new TH1F("ResVsCosTheta", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", 13, cosThetaBins);
    pResVsCosThetaHist->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
    pResVsCosThetaHist->SetXTitle("|cos(#theta)|");

    float resolution(0.f), resolutionError(0.f);
    AnalysisHelper::CalculatePerformance(pPFA,    resolution, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_1,  resolution, resolutionError); pResVsCosThetaHist->SetBinContent( 1, resolution); pResVsCosThetaHist->SetBinError( 1, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_2,  resolution, resolutionError); pResVsCosThetaHist->SetBinContent( 2, resolution); pResVsCosThetaHist->SetBinError( 2, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_3,  resolution, resolutionError); pResVsCosThetaHist->SetBinContent( 3, resolution); pResVsCosThetaHist->SetBinError( 3, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_4,  resolution, resolutionError); pResVsCosThetaHist->SetBinContent( 4, resolution); pResVsCosThetaHist->SetBinError( 4, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_5,  resolution, resolutionError); pResVsCosThetaHist->SetBinContent( 5, resolution); pResVsCosThetaHist->SetBinError( 5, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_6,  resolution, resolutionError); pResVsCosThetaHist->SetBinContent( 6, resolution); pResVsCosThetaHist->SetBinError( 6, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_7,  resolution, resolutionError); pResVsCosThetaHist->SetBinContent( 7, resolution); pResVsCosThetaHist->SetBinError( 7, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_8,  resolution, resolutionError); pResVsCosThetaHist->SetBinContent( 8, resolution); pResVsCosThetaHist->SetBinError( 8, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_9,  resolution, resolutionError); pResVsCosThetaHist->SetBinContent( 9, resolution); pResVsCosThetaHist->SetBinError( 9, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_11, resolution, resolutionError); pResVsCosThetaHist->SetBinContent(10, resolution); pResVsCosThetaHist->SetBinError(10, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_12, resolution, resolutionError); pResVsCosThetaHist->SetBinContent(11, resolution); pResVsCosThetaHist->SetBinError(11, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_13, resolution, resolutionError); pResVsCosThetaHist->SetBinContent(12, resolution); pResVsCosThetaHist->SetBinError(12, resolutionError);
    AnalysisHelper::CalculatePerformance(pPFA_14, resolution, resolutionError); pResVsCosThetaHist->SetBinContent(13, resolution); pResVsCosThetaHist->SetBinError(13, resolutionError);

    AnalysisHelper::CalculatePerformance(pPFA_udsHM20);
    AnalysisHelper::CalculatePerformance(pPFA_udsHM10);
    AnalysisHelper::CalculatePerformance(pPFA_uds);
    AnalysisHelper::CalculatePerformance(pPFA_udsHP10);
    AnalysisHelper::CalculatePerformance(pPFA_udsHP20);

    TH1F *pPFA_1Clone = static_cast<TH1F*>(pPFA_1->Clone());
    pPFA_1Clone->Add(pPFA_2);
    pPFA_1Clone->Add(pPFA_3);
    pPFA_1Clone->Add(pPFA_4);
    pPFA_1Clone->Add(pPFA_5);  std::cout << " < 0.5 : ";   AnalysisHelper::CalculatePerformance(pPFA_1Clone);
    pPFA_1Clone->Add(pPFA_6);  std::cout << " < 0.6 : ";   AnalysisHelper::CalculatePerformance(pPFA_1Clone);
    pPFA_1Clone->Add(pPFA_7);  std::cout << " < 0.7 : ";   AnalysisHelper::CalculatePerformance(pPFA_1Clone);
    pPFA_1Clone->Add(pPFA_8);  std::cout << " < 0.8 : ";   AnalysisHelper::CalculatePerformance(pPFA_1Clone);
    pPFA_1Clone->Add(pPFA_9);  std::cout << " < 0.9 : ";   AnalysisHelper::CalculatePerformance(pPFA_1Clone);
    pPFA_1Clone->Add(pPFA_11); std::cout << " < 0.925 : "; AnalysisHelper::CalculatePerformance(pPFA_1Clone);
    pPFA_1Clone->Add(pPFA_12); std::cout << " < 0.95 : ";  AnalysisHelper::CalculatePerformance(pPFA_1Clone);
    pPFA_1Clone->Add(pPFA_13); std::cout << " < 0.975 : "; AnalysisHelper::CalculatePerformance(pPFA_1Clone);
    pPFA_1Clone->Add(pPFA_14); std::cout << " < 1.0 : ";   AnalysisHelper::CalculatePerformance(pPFA_1Clone);

    TH1F *pPFA_9Clone = static_cast<TH1F*>(pPFA_9->Clone());
    pPFA_9Clone->Add(pPFA_11);
    pPFA_9Clone->Add(pPFA_12); std::cout << " 0.8-0.95 : "; AnalysisHelper::CalculatePerformance(pPFA_9Clone);
    std::cout << " < 0.7 A : ";    AnalysisHelper::CalculatePerformance(pPFA_L7A);
    std::cout << " < 0.7 A ud : "; AnalysisHelper::CalculatePerformance(pPFA_L7Aud);
    std::cout << " < 0.7 A s : ";  AnalysisHelper::CalculatePerformance(pPFA_L7As);
    std::cout << " < 0.7 A NoFwd : ";    AnalysisHelper::CalculatePerformance(pPFA_L7ANoFwd);

    AnalysisHelper::CalculatePerformance(pPFA_MZ);
    AnalysisHelper::CalculatePerformance(pPFA_MW);
    AnalysisHelper::CalculatePerformance(pPFA_MZa);
    AnalysisHelper::CalculatePerformance(pPFA_MWa);
    AnalysisHelper::CalculatePerformance(pPFA_nu);
    AnalysisHelper::CalculatePerformance(pPFA_nufwd);
    AnalysisHelper::CalculatePerformance(pPFA_udscb);
    AnalysisHelper::CalculatePerformance(pPFA_cb);
    AnalysisHelper::CalculatePerformance(pPFA_QQ, false);
    AnalysisHelper::CalculatePerformance(pPFA_QQ8, false);
    AnalysisHelper::CalculatePerformance(pPFA_DMZQQ8, false);
    AnalysisHelper::CalculatePerformance(pPFA_DMZOMZQQ8, false);
    AnalysisHelper::CalculatePerformance(pPFA_DMZ, false);
    AnalysisHelper::CalculatePerformance(pPFA_DMZOMZ, false);
    AnalysisHelper::CalculatePerformance(pPFA_DMZ8, false);
    AnalysisHelper::CalculatePerformance(pPFA_DMZP8, false);

    // Write out histograms
    if (!outputRootFileName.empty())
    {
        std::cout << "Will write histograms to file : " << outputRootFileName << std::endl;
        TFile *pTOutputFile = new TFile(outputRootFileName.c_str(), "RECREATE");
        pTOutputFile->cd();
        pResVsCosThetaHist->Write();
        pNPFO->Write();
        pPFA_MZ->Write();
        pPFA_MW->Write();
        pPFA_MZa->Write();
        pPFA_MWa->Write();
        pPFA->Write();
        pPFA_nu->Write();
        pPFA_nufwd->Write();
        pPFA_udscb->Write();
        pPFA_uds->Write();
        pPFA_udsHM20->Write();
        pPFA_udsHM10->Write();
        pPFA_udsHP10->Write();
        pPFA_udsHP20->Write();
        pPFA_cb->Write();
        pPFA_1->Write();
        pPFA_2->Write();
        pPFA_3->Write();
        pPFA_4->Write();
        pPFA_5->Write();
        pPFA_6->Write();
        pPFA_7->Write();
        pPFA_L7A->Write();
        pPFA_L7Aud->Write();
        pPFA_L7As->Write();
        pPFA_L7Ac->Write();
        pPFA_L7Ab->Write();
        pPFA_L7ANoFwd->Write();
        pPFA_QQ->Write();
        pPFA_QQ8->Write();
        pPFA_8->Write();
        pPFA_9->Write();
        pPFA_10->Write();
        pPFA_11->Write();
        pPFA_12->Write();
        pPFA_13->Write();
        pPFA_14->Write();
        pPFA_DMZ->Write();
        pPFA_DMZ8->Write();
        pPFA_DMZQQ8->Write();
        pPFA_DMZP8->Write();
        pPFA_DMZOMZ->Write();
        pPFA_DMZOMZQQ8->Write();
        pTOutputFile->Close();
        delete pTOutputFile;
    }

    // Tidy up
    delete pNPFO;
    delete pPFA_MZ;
    delete pPFA_MW;
    delete pPFA_MZa;
    delete pPFA_MWa;
    delete pPFA;
    delete pPFA_nu;
    delete pPFA_nufwd;
    delete pPFA_udscb;
    delete pPFA_uds;
    delete pPFA_udsHM20;
    delete pPFA_udsHM10;
    delete pPFA_udsHP10;
    delete pPFA_udsHP20;
    delete pPFA_cb;
    delete pPFA_1;
    delete pPFA_2;
    delete pPFA_3;
    delete pPFA_4;
    delete pPFA_5;
    delete pPFA_6;
    delete pPFA_7;
    delete pPFA_L7A;
    delete pPFA_L7Aud;
    delete pPFA_L7As;
    delete pPFA_L7Ac;
    delete pPFA_L7Ab;
    delete pPFA_L7ANoFwd;
    delete pPFA_QQ;
    delete pPFA_QQ8;
    delete pPFA_8;
    delete pPFA_9;
    delete pPFA_10;
    delete pPFA_11;
    delete pPFA_12;
    delete pPFA_13;
    delete pPFA_14;
    delete pPFA_DMZ;
    delete pPFA_DMZ8;
    delete pPFA_DMZQQ8;
    delete pPFA_DMZP8;
    delete pPFA_DMZOMZ;
    delete pPFA_DMZOMZQQ8;
    delete pPFA_1Clone;
    delete pPFA_9Clone;
    delete pResVsCosThetaHist;
}
