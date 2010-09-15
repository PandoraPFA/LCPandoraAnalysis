/**
 *  @file   PandoraAnalysis/src/AnalysePerformance.cc
 * 
 *  @brief  Implementation of pandora analyse performance binary.
 * 
 *  $Log: $
 */

#include "TFile.h"
#include "TH1F.h"
#include "TRandom.h"

#include "AnalysePerformance.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>

int main(int argc, char **argv)
{
    if ((argc < 2) || (argc > 3))
    {
        std::cout << std::endl
                  << "Usage: ./analysePerformance inputFileList outputRootFile" << std::endl << std::endl
                  << "  inputFileList  : file listing input root files, one per line, list terminated by \"END\" " << std::endl
                  << "  outputRootFile : optional output root file name, for output of analysis histograms" << std::endl << std::endl;
        exit(1) ;
    }

    m_pRandom = new TRandom();
    const std::string inputFileList(argv[1]);
    std::string outputRootFileName((argc == 3) ? argv[2] : "");

    // Read list of input files
    ifstream infile(inputFileList.c_str());

    if (!infile)
    {
        std::cerr << std::endl << "ATTN: Can't open input file list : " << inputFileList.c_str() << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<std::string>files;

    while (getline(infile, line, '\n'))
    {
        const unsigned int p(line.find("END"));

        if ((p != std::string::npos) && (p < 2))
            break;

        if (line[0] != 'C')
        {
            std::cout << "Adding file to list : " << line << std::endl;
            files.push_back(line);
        }
    }

    if (files.empty())
    {
        std::cout << std::endl << "ATTN: Input file list is empty " << std::endl;
    }

    if ((files.size() > 1) && (!outputRootFileName.empty()))
    {
        std::cout << std::endl << "ATTN: Histograms will only be written for first file in list " << std::endl;
    }

    // Process the files in the list
    std::string printString = "P";

    for (unsigned int i = 0; i < files.size(); ++i)
    {
        if (files[i].at(0) == printString)
        {
            std::cout << files[i] << std::endl;
        }
        else
        {
            std::cout << std::endl << "Processing file: "<< files[i] << std::endl;

            if (!ifstream(files[i].c_str()))
            {
                std::cerr << "Can't open file : " << files[i] << std::endl;
            }
            else
            {
                TFile *pTFile = new TFile(files[i].c_str(), "READ");
                AnalyseHistograms(pTFile, outputRootFileName);
                outputRootFileName.clear();
                pTFile->Close();
                delete pTFile;
            }
        }
    }

    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalyseHistograms(TFile* pTFile, const std::string &outputRootFileName)
{
    float xbins[14] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95, 0.975, 1.0};
    TH1F *pEvsCHist = new TH1F("EvsC", "sigmaE/E vs costheta", 13, xbins);

    CalcRms((TH1F*)pTFile->Get("fPFA"));
    CalcRms((TH1F*)pTFile->Get("fPFA1"));  pEvsCHist->SetBinContent( 1, m_sigma); pEvsCHist->SetBinError( 1, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA2"));  pEvsCHist->SetBinContent( 2, m_sigma); pEvsCHist->SetBinError( 2, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA3"));  pEvsCHist->SetBinContent( 3, m_sigma); pEvsCHist->SetBinError( 3, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA4"));  pEvsCHist->SetBinContent( 4, m_sigma); pEvsCHist->SetBinError( 4, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA5"));  pEvsCHist->SetBinContent( 5, m_sigma); pEvsCHist->SetBinError( 5, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA6"));  pEvsCHist->SetBinContent( 6, m_sigma); pEvsCHist->SetBinError( 6, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA7"));  pEvsCHist->SetBinContent( 7, m_sigma); pEvsCHist->SetBinError( 7, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA8"));  pEvsCHist->SetBinContent( 8, m_sigma); pEvsCHist->SetBinError( 8, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA9"));  pEvsCHist->SetBinContent( 9, m_sigma); pEvsCHist->SetBinError( 9, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA11")); pEvsCHist->SetBinContent(10, m_sigma); pEvsCHist->SetBinError(10, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA12")); pEvsCHist->SetBinContent(11, m_sigma); pEvsCHist->SetBinError(11, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA13")); pEvsCHist->SetBinContent(12, m_sigma); pEvsCHist->SetBinError(12, m_sigmasigma);
    CalcRms((TH1F*)pTFile->Get("fPFA14")); pEvsCHist->SetBinContent(13, m_sigma); pEvsCHist->SetBinError(13, m_sigmasigma);

    CalcRms((TH1F*)pTFile->Get("fPFAudsHM20"));
    CalcRms((TH1F*)pTFile->Get("fPFAudsHM10"));
    CalcRms((TH1F*)pTFile->Get("fPFAuds"));
    CalcRms((TH1F*)pTFile->Get("fPFAudsHP10"));
    CalcRms((TH1F*)pTFile->Get("fPFAudsHP20"));

    TH1F *pTH1F = (TH1F*)pTFile->Get("fPFA1");
    pTH1F->Add((TH1F*)pTFile->Get("fPFA2"));
    pTH1F->Add((TH1F*)pTFile->Get("fPFA3"));
    pTH1F->Add((TH1F*)pTFile->Get("fPFA4"));
    pTH1F->Add((TH1F*)pTFile->Get("fPFA5"));  std::cout << " < 0.5 : ";    CalcRms(pTH1F);
    pTH1F->Add((TH1F*)pTFile->Get("fPFA6"));  std::cout << " < 0.6 : ";    CalcRms(pTH1F);
    pTH1F->Add((TH1F*)pTFile->Get("fPFA7"));  std::cout << " < 0.7 : ";    CalcRms(pTH1F);
    pTH1F->Add((TH1F*)pTFile->Get("fPFA8"));  std::cout << " < 0.8 : ";    CalcRms(pTH1F);
    pTH1F->Add((TH1F*)pTFile->Get("fPFA9"));  std::cout << " < 0.9 : ";    CalcRms(pTH1F);
    pTH1F->Add((TH1F*)pTFile->Get("fPFA11")); std::cout << " < 0.925 : ";  CalcRms(pTH1F);
    pTH1F->Add((TH1F*)pTFile->Get("fPFA12")); std::cout << " < 0.95 : ";   CalcRms(pTH1F);
    pTH1F->Add((TH1F*)pTFile->Get("fPFA13")); std::cout << " < 0.975 : ";  CalcRms(pTH1F);
    pTH1F->Add((TH1F*)pTFile->Get("fPFA14")); std::cout << " < 1.0 : ";    CalcRms(pTH1F);

    pTH1F = (TH1F*)pTFile->Get("fPFA9");
    pTH1F->Add((TH1F*)pTFile->Get("fPFA11"));
    pTH1F->Add((TH1F*)pTFile->Get("fPFA12")); std::cout << " 0.8-0.95 : "; CalcRms(pTH1F);

    std::cout << " < 0.7 A : ";    CalcRms((TH1F*)pTFile->Get("fPFAL7A"));
    std::cout << " < 0.7 A ud : "; CalcRms((TH1F*)pTFile->Get("fPFAL7Aud"));
    std::cout << " < 0.7 A s : ";  CalcRms((TH1F*)pTFile->Get("fPFAL7As"));
    std::cout << " < 0.7 B : ";    CalcRms((TH1F*)pTFile->Get("fPFAL7B"));

    CalcRms((TH1F*)pTFile->Get("fPFAMZ"));
    CalcRms((TH1F*)pTFile->Get("fPFAMW"));
    CalcRms((TH1F*)pTFile->Get("fPFAMZa"));
    CalcRms((TH1F*)pTFile->Get("fPFAMWa"));
    CalcRms((TH1F*)pTFile->Get("fPFAnu"));
    CalcRms((TH1F*)pTFile->Get("fPFAnufwd"));
    CalcRms((TH1F*)pTFile->Get("fPFAudscb"));
    CalcRms((TH1F*)pTFile->Get("fPFAcb"));

    const bool doNotFixCentre(false);
    CalcRms((TH1F*)pTFile->Get("fPFAQQ"), doNotFixCentre);
    CalcRms((TH1F*)pTFile->Get("fPFAQQ8"), doNotFixCentre);
    CalcRms((TH1F*)pTFile->Get("fPFADMZQQ8"), doNotFixCentre);
    CalcRms((TH1F*)pTFile->Get("fPFADMZOMZQQ8"), doNotFixCentre);
    CalcRms((TH1F*)pTFile->Get("fPFADMZ"), doNotFixCentre);
    CalcRms((TH1F*)pTFile->Get("fPFADMZOMZ"), doNotFixCentre);
    CalcRms((TH1F*)pTFile->Get("fPFADMZ8"), doNotFixCentre);
    CalcRms((TH1F*)pTFile->Get("fPFADMZP8"), doNotFixCentre);

    if (!outputRootFileName.empty())
    {
        std::cout << "Will write histograms to file : " << outputRootFileName << std::endl;
        TFile *pTOutputFile = new TFile(outputRootFileName.c_str(), "RECREATE");
        pTOutputFile->cd();
        pEvsCHist->Write();
        pTOutputFile->Close();
        delete pTOutputFile;
    }

    delete pEvsCHist;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalcRms(TH1 *pTH1F, bool fixDistributionCentre)
{
    if (NULL == pTH1F)
        return;

    // Calculate raw properties of distribution
    float sum = 0., total = 0.;
    double sx = 0., sxx = 0.;
    const unsigned int nbins(pTH1F->GetNbinsX());

    for (unsigned int i = 0; i <= nbins; ++i)
    {
        const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
        const float yi(pTH1F->GetBinContent(i));
        sx += yi * binx;
        sxx += yi * binx * binx;
        total += yi;
    }

    const float rawMean(sx / total);
    const float rawMeanSquared(sxx / total);
    const float rawRms(std::sqrt(rawMeanSquared - rawMean * rawMean));

    sum = 0.;
    unsigned int is0 = 0;

    for (unsigned int i = 0; (i <= nbins) && (sum < total / 10.); ++i)
    {
        sum += pTH1F->GetBinContent(i);
        is0 = i;
    }

    // Calculate truncated properties
    static const float FLOAT_MAX(std::numeric_limits<float>::max());

    float rmsmin(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX);
    m_mean = m_low = m_rms = m_sigma = m_sigmasigma = FLOAT_MAX;
    m_high = 0;

    for (unsigned int istart = 0; istart <= is0; ++istart)
    {
        double sumn = 0.;
        double csum = 0.;
        double sumx = 0.;
        double sumxx = 0.;
        unsigned int iend = 0;

        for (unsigned int i = istart; (i <= nbins) && (csum < 0.9 * total); ++i)
        {
            const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
            const float yi(pTH1F->GetBinContent(i));
            csum += yi;

            if (sumn < 0.9 * total)
            {
                sumn += yi;
                sumx += yi * binx;
                sumxx+= yi * binx * binx;
                iend = i;
            }
        }

        const float mean(sumx / sumn);
        const float meanSquared(sumxx / sumn);
        const float rms(std::sqrt(meanSquared - mean * mean));

        if (rms < rmsmin)
        {
            m_mean = mean;
            m_rms = rms;
            m_low = pTH1F->GetBinLowEdge(istart);
            m_high = pTH1F->GetBinLowEdge(iend);
            rmsmin = rms;

            if (fixDistributionCentre)
            {
                float centre = 91.2;

                if ((mean > 2500.) && (mean < 3500.))
                    centre = 3000.;

                if ((mean > 1500.) && (mean < 2500.))
                    centre = 2000.;

                if ((mean > 700.) && (mean < 1500.))
                    centre = 1000.;

                if ((mean > 400.) && (mean < 700.))
                    centre = 500.;

                if ((mean > 250.) && (mean < 400.))
                    centre = 360.;

                if ((mean > 150.) && (mean < 250.))
                    centre = 200.;

                m_sigma = rms / mean * sqrt(centre);
                m_sigmasigma = m_sigma / std::sqrt(total);

                frac = rms / mean * std::sqrt(2) * 100.;
                efrac = frac / std::sqrt(total);
            }
            else
            {
                m_sigma = rms;
                m_sigmasigma = m_sigma / std::sqrt(total);
            }
        }
    }

    std::cout << " " << pTH1F->GetName() << ", rawrms " << rawRms <<  ", rms90 : " << rmsmin
              << ", " << m_low << "-" << m_high << ", mean : " << m_mean
              << ", sigma : " << m_sigma <<  "+-" << m_sigmasigma;

    (fixDistributionCentre) ? (std::cout << ", sE/E : " << frac << "+-" << efrac << std::endl) : (std::cout << std::endl);
}
