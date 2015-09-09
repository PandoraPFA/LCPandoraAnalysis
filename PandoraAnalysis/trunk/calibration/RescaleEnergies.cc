/**
 *  @file   PandoraAnalysis/calibration/RescaleEnergies.cc
 * 
 *  @brief  Implementation of rescale energies binary.
 * 
 *  $Log: $
 */

#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TTree.h"

#include "AnalysisHelper.h"

#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <iostream>
#include <string>

using namespace pandora_analysis;

/**
 *  @brief  Parameters class
 */
class Parameters
{
public:
    /**
     *  @brief Default constructor
     */
    Parameters();

    bool            m_shouldPlotResults;            ///< Whether to plot resulting energy distribution
    std::string     m_inputFileNames;               ///< The name of the file(s) containing the Pandora analysis tree(s)
    double          m_eCalToEm;                     ///< ECalToEm calibration factor
    double          m_hCalToEm;                     ///< HCalToEm calibration factor
    double          m_eCalToHad;                    ///< ECalToHad calibration factor
    double          m_hCalToHad;                    ///< HCalToHad calibration factor
    double          m_muon;                         ///< Muon calibration factor
    TChain         *m_pTChain;                      ///< The chain of Pandora analysis trees
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get sigma
 * 
 *  @param  parameters the parameters
 */
void GetSigma(const Parameters &parameters);

/**
 *  @brief  Pause
 */
void Pause();

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 * 
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  parameters to receive the application parameters
 * 
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], Parameters &parameters);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    try
    {
        Parameters parameters;

        if (!ParseCommandLine(argc, argv, parameters))
            return 1;

        int myargc = 0;
        char* myargv = (char *)"";
        TApplication *pApplication = new TApplication("PandoraMonitoring", &myargc, &myargv);
        pApplication->SetReturnFromRun(kTRUE);

        parameters.m_pTChain->Add(parameters.m_inputFileNames.c_str());

        if (0 == parameters.m_pTChain->GetEntries())
        {
            std::cout << "Error opening PfoAnalysisTree " << std::endl;
            return 1;
        }

        GetSigma(parameters);

        delete parameters.m_pTChain;
        parameters.m_pTChain = NULL;
        delete pApplication;
    }
    catch (std::exception &exception)
    {
        std::cout << "Exception caught " << exception.what() << std::endl;
        return 1;
    }

    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GetSigma(const Parameters &parameters)
{
    const unsigned int nEntries(parameters.m_pTChain->GetEntries());

    int qPdg(0);
    float pfoEnergyTracks(0.f), pfoECalToEmEnergy(0.f), pfoHCalToEmEnergy(0.f), pfoECalToHadEnergy(0.f), pfoHCalToHadEnergy(0.f),
        pfoOtherEnergy(0.f), pfoMuonToEnergy(0.f), mcEnergyENu(0.f), thrust(0.f);

    parameters.m_pTChain->SetBranchAddress("pfoEnergyTracks", &pfoEnergyTracks);
    parameters.m_pTChain->SetBranchAddress("pfoECalToEmEnergy", &pfoECalToEmEnergy);
    parameters.m_pTChain->SetBranchAddress("pfoHCalToEmEnergy", &pfoHCalToEmEnergy);
    parameters.m_pTChain->SetBranchAddress("pfoECalToHadEnergy", &pfoECalToHadEnergy);
    parameters.m_pTChain->SetBranchAddress("pfoHCalToHadEnergy", &pfoHCalToHadEnergy);
    parameters.m_pTChain->SetBranchAddress("pfoOtherEnergy", &pfoOtherEnergy);
    parameters.m_pTChain->SetBranchAddress("pfoMuonToEnergy", &pfoMuonToEnergy);
    parameters.m_pTChain->SetBranchAddress("mcEnergyENu", &mcEnergyENu);
    parameters.m_pTChain->SetBranchAddress("qPdg", &qPdg);
    parameters.m_pTChain->SetBranchAddress("thrust", &thrust);

    TH1F *pPFAL7A = new TH1F("fPFA_L7A", "TotalEnergy<0.7A", 10000, 0., 5000.);

    for (unsigned int iTree = 0; iTree < nEntries; ++iTree)
    {
        parameters.m_pTChain->GetEntry(iTree);

        //if ((qPdg < 1) || (qPdg > 3))
        //    continue;

        if (thrust <= 0.7f)
        {
            float totalRecoEnergy(pfoEnergyTracks);
            totalRecoEnergy += parameters.m_eCalToEm  * static_cast<double>(pfoECalToEmEnergy);
            totalRecoEnergy += parameters.m_hCalToEm  * static_cast<double>(pfoHCalToEmEnergy);
            totalRecoEnergy += parameters.m_eCalToHad * static_cast<double>(pfoECalToHadEnergy);
            totalRecoEnergy += parameters.m_hCalToHad * static_cast<double>(pfoHCalToHadEnergy);
            totalRecoEnergy += parameters.m_muon      * static_cast<double>(pfoMuonToEnergy);
            totalRecoEnergy += pfoOtherEnergy;

            pPFAL7A->Fill(totalRecoEnergy + mcEnergyENu, 1.);
        }
    }

    // Extract performance figures from energy spectra histogram
    float sigma(0.f), sigmasigma(0.f);
    AnalysisHelper::CalculatePerformance(pPFAL7A, sigma, sigmasigma, true, true);

    std::cout << "   sigma: " << sigma << " (ECalToEm: " << parameters.m_eCalToEm << ", HCalToEm: " << parameters.m_hCalToEm
              << ", ECalToHad: " << parameters.m_eCalToHad << ", HCalToHad: " << parameters.m_hCalToHad << ", Muon: " << parameters.m_muon << ") "
              << std::endl << std::endl;

    if (parameters.m_shouldPlotResults)
    {
        TCanvas *pCanvas = new TCanvas();
        pCanvas->SetFillColor(kWhite);
        pCanvas->SetHighLightColor(kWhite);
        pCanvas->Draw();
        pPFAL7A->Draw();
        Pause();
        delete pCanvas;
    }

    delete pPFAL7A;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Pause()
{
#ifdef __unix__
    std::cout << "Press return to continue ..." << std::endl;
    int flag = fcntl(1, F_GETFL, 0);

    int key = 0;
    while(true)
    {
        gSystem->ProcessEvents();
        (void) fcntl(1, F_SETFL, flag | O_NONBLOCK);
        key = getchar();

        if((key == '\n') || (key == '\r'))
            break;

        usleep(1000);
    }

    (void) fcntl(1, F_SETFL, flag);
#else
    std::cout << "Pause() is only implemented for unix operating systems." << std::endl;
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], Parameters &parameters)
{
    int c(0);

    while (((c = getopt(argc, argv, "i:p::a:b:c:d:e:h")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'i':
            parameters.m_inputFileNames = optarg;
            break;
        case 'p':
            parameters.m_shouldPlotResults = true;
            break;
        case 'a':
            parameters.m_eCalToEm = atof(optarg);
            break;
        case 'b':
            parameters.m_hCalToEm = atof(optarg);
            break;
        case 'c':
            parameters.m_eCalToHad = atof(optarg);
            break;
        case 'd':
            parameters.m_hCalToHad = atof(optarg);
            break;
        case 'e':
            parameters.m_muon = atof(optarg);
            break;
        case 'h':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -i name   (mandatory, input file name(s), can include wildcards if string is in quotes)" << std::endl
                      << "    -p        (optional, whether to plot results,      default 0)" << std::endl
                      << "    -a value  (optional, ECalToEM calibration factor,  default 1)" << std::endl
                      << "    -b value  (optional, HCalToEM calibration factor,  default 1)" << std::endl
                      << "    -c value  (optional, ECalToHad calibration factor, default 1)" << std::endl
                      << "    -d value  (optional, HCalToHad calibration factor, default 1)" << std::endl
                      << "    -e value  (optional, Muon calibration factor,      default 1)" << std::endl
                      << std::endl;
            return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

Parameters::Parameters() :
    m_shouldPlotResults(false),
    m_inputFileNames(""),
    m_eCalToEm(1.),
    m_hCalToEm(1.),
    m_eCalToHad(1.),
    m_hCalToHad(1.),
    m_muon(1.),
    m_pTChain(new TChain("PfoAnalysisTree"))
{
}
