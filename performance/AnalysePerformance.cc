/**
 *  @file   PandoraAnalysis/performance/AnalysePerformance.cc
 * 
 *  @brief  Implementation of pandora analyse performance binary.
 * 
 *  $Log: $
 */

#include "TFile.h"
#include "TH1F.h"
#include "TChain.h"

#include "AnalysisHelper.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

using namespace pandora_analysis;

void AnalysePerformance(TChain *pTChain, const std::string &outputRootFileName);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    try
    {
        const int nArgs(argc - 1);

        if ((nArgs < 1) || (nArgs > 2))
        {
            std::cout << std::endl
                      << "Usage: ./AnalysePerformance inputFileName [outputFileName]" << std::endl << std::endl
                      << "  inputFileName  : file containing pandora pfo analysis tree" << std::endl
                      << "  outputFileName : optional output root file, for histogram output" << std::endl << std::endl;
            return 1;
        }

        const std::string inputFileName(argv[1]);
        const std::string outputRootFileName((nArgs == 2) ? argv[2] : "");

        TChain *pTChain = new TChain("PfoAnalysisTree");
        pTChain->Add(inputFileName.c_str());

        if (0 == pTChain->GetEntries())
        {
            std::cout << "Error opening PfoAnalysisTree " << std::endl;
            return 1;
        }

        AnalysePerformance(pTChain, outputRootFileName);
        delete pTChain;
    }
    catch (std::exception &exception)
    {
        std::cout << "Exception caught " << exception.what() << std::endl;
        return 1;
    }

    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysePerformance(TChain *pTChain, const std::string &outputRootFileName)
{
    // Define regions for which to create and investigate pfo energy spectra
    const unsigned int nRegionBins(13);
    float pRegionBinEdges[nRegionBins + 1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 1.0};

    TH1F *pResVsCosThetaHist = new TH1F("ResVsCosTheta", "RMS_{90}(E_{j}) / Mean_{90}(E_{j}) vs |cos(#theta)|", nRegionBins, pRegionBinEdges);
    pResVsCosThetaHist->SetYTitle("RMS_{90}(E_{j}) / Mean_{90}(E_{j}) [%]");
    pResVsCosThetaHist->SetXTitle("|cos(#theta)|");

    // Book histograms
    TH1F **pRegionHistograms = new TH1F*[nRegionBins];
    TH1F *pPFAL7A = new TH1F("fPFA_L7A", "TotalEnergy<0.7A", 100000, 0., 5000.);

    for (unsigned int i = 0; i < nRegionBins; ++i)
    {
        std::ostringstream name, title;
        name << "fPFA_" << i;
        title << "TotalEnergy_" << pRegionBinEdges[i] << "-" << pRegionBinEdges[i + 1];
        pRegionHistograms[i] = new TH1F(name.str().c_str(), title.str().c_str(), 100000, 0., 5000.);
    }

    // Loop over entries in pfo analysis tree, creating energy spectra for specified regions
    const unsigned int nTreeEntries(pTChain->GetEntries());

    int qPdg(0);
    float pfoEnergyTotal(0.f), mcEnergyENu(0.f), thrust(0.f);
    pTChain->SetBranchAddress("pfoEnergyTotal", &pfoEnergyTotal);
    pTChain->SetBranchAddress("mcEnergyENu", &mcEnergyENu);
    pTChain->SetBranchAddress("qPdg", &qPdg);
    pTChain->SetBranchAddress("thrust", &thrust);

    for (unsigned int iTree = 0; iTree < nTreeEntries; ++iTree)
    {
        pTChain->GetEntry(iTree);

        if ((qPdg < 1) || (qPdg > 3))
            continue;

        if (thrust <= 0.7f)
            pPFAL7A->Fill(pfoEnergyTotal + mcEnergyENu, 1.);

        for (unsigned int i = 0; i < nRegionBins; ++i)
        {
            if ((thrust >= pRegionBinEdges[i]) && (thrust < pRegionBinEdges[i + 1]))
                pRegionHistograms[i]->Fill(pfoEnergyTotal, 1.);
        }
    }

    // Extract performance figures from energy spectra histograms
    float resolution(0.f), resolutionError(0.f);
    AnalysisHelper::CalculatePerformance(pPFAL7A, resolution, resolutionError);

    for (unsigned int i = 0; i < nRegionBins; ++i)
    {
        resolution = 0.f; resolutionError = 0.f;
        AnalysisHelper::CalculatePerformance(pRegionHistograms[i], resolution, resolutionError);
        pResVsCosThetaHist->SetBinContent(i + 1, resolution); pResVsCosThetaHist->SetBinError(i + 1, resolutionError);
    }

    // Write histograms to output file, if specified
    if (!outputRootFileName.empty())
    {
        std::cout << "Will write histograms to file : " << outputRootFileName << std::endl;
        TFile *pTOutputFile = new TFile(outputRootFileName.c_str(), "RECREATE");
        pTOutputFile->cd();
        pResVsCosThetaHist->GetYaxis()->SetRangeUser(0.f, 10.f);
        pResVsCosThetaHist->Write();
        pPFAL7A->Write();

        for (unsigned int i = 0; i < nRegionBins; ++i)
            pRegionHistograms[i]->Write();

        pTOutputFile->Close();
        delete pTOutputFile;
    }

    // Tidy up
    for (unsigned int i = 0; i < nRegionBins; ++i)
        delete pRegionHistograms[i];

    delete pPFAL7A;
    delete pResVsCosThetaHist;
    delete [] pRegionHistograms;
}
