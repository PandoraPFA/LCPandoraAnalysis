/**
 *  @file   PandoraAnalysis/calibration/PandoraPFACalibrate_HadronicScale_ChiSquareMethod.cc
 * 
 *  @brief  Calculate E/HCalToHadGeVCalibration using the chi squared method on kaonL events
 * 
 *  $Log: $
 */

#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TTree.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <stdexcept>

/**
 *  @brief  ChiSquaredMethod class
 */
class ChiSquaredMethod 
{
public:
    /**
     *  @brief  Constructor
    */
    ChiSquaredMethod();

    /**
     *  @brief  Destructor
     */
    ~ChiSquaredMethod();

    /**
     *  @brief  Perform the chi squared method for hadronic energy scale calibration in PandoraPFA
    */
    void Process();

// Outputs
    float               m_eCalToHadInterceptMinChi2;    ///< Best fit ECalToHad intercept
    float               m_hCalToHadInterceptMinChi2;    ///< Best fit HCalToHad intercept
    float               m_minChi2;                      ///< Best fit measure, chi squared of fit
    int                 m_nEntriesHist;                 ///< Number of entries in histogram

// Non trivial setting on initialisation
    float               m_trueEnergy;                   ///< Total energy of calibration particle
    std::string         m_inputKaonLRootFiles;          ///< Input root files for E/HCalToHad Calibration
    std::string         m_outputPath;                   ///< Output path to send plots to
    int                 m_numberHCalLayers;             ///< Number of layers in the HCal

private:
    /**
     *  @brief  Create the histogram needed for the chi squared method
    */
    void CreateHistogram();

    /**
     *  @brief  Fill the histogram needed for the chi squared method
    */
    void FillHistogram();

    /**
     *  @brief  Perform the chi squared method.
    */
    void CSM();

    /**
     *  @brief  Does the point defined by m_pfoHCalToHadEnergy and m_pfoECalToHadEnergy fall within 3 sigme of the ideal distribution
    */
    bool ThreeSigmaCut();

    /**
     *  @brief  Plot m_histogram.
    */
    void Plot();

// Non trivial setting on initialisation
    float               m_kineticEnergy;                ///< Kinetic energy of calibration particle, from m_trueEnergy
    int                 m_binNumber;                    ///< Number of bins (x/y symmetric), taken as 150
    float               m_maxHistogramEnergy;           ///< Max energy on histogram (x/y symmetric), taken as 1.5 * m_trueEnergy
    float               m_hCalToHadResolutionConstant;  ///< HCal coefficient for stoichastic (sigma_E/E ~ coefficeint /sqrt(E)) process default 0.55
    float               m_eCalToHadResolutionConstant;  ///< ECal coefficient for stoichastic (sigma_E/E ~ coefficeint /sqrt(E)) process default 0.55

// Trivial setting on initialisation
    TChain             *m_pTChain;                      ///< Chain of root files
    TH2F               *m_histogram;                    ///< 2D histogram of ECalToHad vs HCalToHad energy
    float               m_pfoHCalToHadEnergy;           ///< Error from stoichastic term from single particle resolution : sigmaE/E ~ ResConstant / sqrt(E)
    float               m_pfoECalToHadEnergy;           ///< Error from stoichastic term from single particle resolution : sigmaE/E ~ ResConstant / sqrt(E)

typedef std::vector<float> FloatVector;
    FloatVector         m_eCalEnergyAccepted;           ///< Vector of ECalToHad energies passing three sigma cuts
    FloatVector         m_hCalEnergyAccepted;           ///< Vector of HCalToHad energies passing three sigma cuts
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 * 
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  eCalDigitisation to receive the application parameters
 * 
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], ChiSquaredMethod &chiSquaredMethod);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();

    try
    {
        ChiSquaredMethod chiSquaredMethod;

        if (!ParseCommandLine(argc, argv, chiSquaredMethod))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        chiSquaredMethod.Process();

        std::string dataFileName(chiSquaredMethod.m_outputPath + "Calibration.txt");
        ofstream data_file(dataFileName.c_str(), std::ios_base::app);
        data_file << "_____________________________________________________________________________________" << std::endl;
        data_file << "Hadronic Energy Scale PandoraPFA Calibration performed using the Chi Squared Method." << std::endl << std::endl;
        data_file << "For kaon energy                                    : " << chiSquaredMethod.m_trueEnergy << " : " <<std::endl;
        data_file << "m_eCalToHadInterceptMinChi2                        : " << chiSquaredMethod.m_eCalToHadInterceptMinChi2 << " : " <<std::endl;
        data_file << "m_hCalToHadInterceptMinChi2                        : " << chiSquaredMethod.m_hCalToHadInterceptMinChi2 << " : " <<std::endl;
        data_file << "m_minChi2                                          : " << chiSquaredMethod.m_minChi2 << " : " <<std::endl;
        data_file << "The total number of events considered was          : " << chiSquaredMethod.m_nEntriesHist << " : " <<std::endl<<std::endl;

        std::cout << "_____________________________________________________________________________________" << std::endl;
        std::cout << "For kaon energy                                    : " << chiSquaredMethod.m_trueEnergy << " : " <<std::endl;
        std::cout << "m_eCalToHadInterceptMinChi2                        : " << chiSquaredMethod.m_eCalToHadInterceptMinChi2 << " : " <<std::endl;
        std::cout << "m_hCalToHadInterceptMinChi2                        : " << chiSquaredMethod.m_hCalToHadInterceptMinChi2 << " : " <<std::endl;
        std::cout << "m_minChi2                                          : " << chiSquaredMethod.m_minChi2 << " : " <<std::endl;
        std::cout << "The total number of events considered was          : " << chiSquaredMethod.m_nEntriesHist << " : " <<std::endl<< std::endl;

        data_file.close();
        delete pTApplication;
    }

    catch (std::exception &exception)
    {
        std::cout << "Exception caught " << exception.what() << std::endl;
        delete pTApplication;
        return 1;
    }
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ChiSquaredMethod::ChiSquaredMethod() :
    m_eCalToHadInterceptMinChi2(std::numeric_limits<float>::max()),
    m_hCalToHadInterceptMinChi2(std::numeric_limits<float>::max()),
    m_minChi2(std::numeric_limits<float>::max()),
    m_nEntriesHist(0),
    m_trueEnergy(std::numeric_limits<float>::max()),
    m_inputKaonLRootFiles(""),
    m_outputPath(""),
    m_numberHCalLayers(48),
    m_kineticEnergy(std::numeric_limits<float>::max()),
    m_binNumber(150),
    m_maxHistogramEnergy(std::numeric_limits<float>::max()),
    m_hCalToHadResolutionConstant(0.55), // Resolution constants from single particle calibration
    m_eCalToHadResolutionConstant(0.55), // Resolution constants from single particle calibration
    m_pTChain(NULL),
    m_histogram(NULL),
    m_pfoHCalToHadEnergy(std::numeric_limits<float>::max()),
    m_pfoECalToHadEnergy(std::numeric_limits<float>::max()),
    m_eCalEnergyAccepted(NULL),
    m_hCalEnergyAccepted(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ChiSquaredMethod::~ChiSquaredMethod()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::Process()
{
    m_maxHistogramEnergy = 1.5 * m_trueEnergy;
    m_kineticEnergy = m_trueEnergy - 0.497672;
    this->CreateHistogram();
    this->FillHistogram();
    this->CSM();
    this->Plot();
}
//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::CreateHistogram()
{
    std::string Name = "histCSM";
    std::string Title = "Hadronic Energy Scale PandoraPFA Calibration, Chi Squared Method";
    m_histogram = new TH2F(Name.c_str(),Title.c_str(),m_binNumber,0,m_maxHistogramEnergy,m_binNumber,0,m_maxHistogramEnergy);
    m_histogram->GetYaxis()->SetTitle("Hadronic Energy Measured in the ECal / GeV");
    m_histogram->GetXaxis()->SetTitle("Hadronic Energy Measured in the HCal / GeV");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::FillHistogram()
{
    int nPfoTargetsTotal;
    int nPfoTargetsNeutralHadrons;
    int nPfosTotal;
    int nPfosNeutralHadrons;
    int pfoMinHCalLayerToEdge;

    m_pTChain = new TChain("PfoAnalysisTree");
    m_pTChain->Add(m_inputKaonLRootFiles.c_str());

    unsigned int nEntriesHistogram(m_pTChain->GetEntries());

    m_pTChain->SetBranchAddress("pfoECalToHadEnergy",&m_pfoECalToHadEnergy);
    m_pTChain->SetBranchAddress("pfoHCalToHadEnergy",&m_pfoHCalToHadEnergy);
    m_pTChain->SetBranchAddress("nPfoTargetsTotal",&nPfoTargetsTotal);
    m_pTChain->SetBranchAddress("nPfoTargetsNeutralHadrons",&nPfoTargetsNeutralHadrons);
    m_pTChain->SetBranchAddress("nPfosTotal",&nPfosTotal);
    m_pTChain->SetBranchAddress("nPfosNeutralHadrons",&nPfosNeutralHadrons);
    m_pTChain->SetBranchAddress("pfoMinHCalLayerToEdge",&pfoMinHCalLayerToEdge);

    int containedLayerExclusion = ceil(m_numberHCalLayers * 0.1);

    for (unsigned int i = 0; i < nEntriesHistogram ; i++) 
    {
        m_pTChain->GetEntry(i);

        if (ThreeSigmaCut() && nPfoTargetsTotal == 1 && nPfoTargetsNeutralHadrons == 1 && nPfosTotal == 1 && nPfosNeutralHadrons == 1 && 
            pfoMinHCalLayerToEdge > containedLayerExclusion)
        {
            m_eCalEnergyAccepted.push_back(m_pfoECalToHadEnergy);
            m_hCalEnergyAccepted.push_back(m_pfoHCalToHadEnergy);
            m_nEntriesHist++;
            m_histogram->Fill(m_pfoHCalToHadEnergy,m_pfoECalToHadEnergy);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::CSM()
{
    float perpDist = std::numeric_limits<float>::max();
    float perpDistSigma = std::numeric_limits<float>::max();
    float eCalToHadIntercept = std::numeric_limits<float>::max();
    float hCalToHadIntercept = std::numeric_limits<float>::max();
    float Chi2 = std::numeric_limits<float>::max();

    for (unsigned int p = 0; p < 100; p++)
    {
        eCalToHadIntercept = (m_kineticEnergy * 0.75) + (m_kineticEnergy * 0.01 * p * 0.5);

        for (unsigned int q = 0; q < 100 ; q++)
        {
            Chi2 = 0;

            hCalToHadIntercept = (m_kineticEnergy * 0.75) + (m_kineticEnergy * 0.01 * q * 0.5);

            for (unsigned int j = 0; j < m_eCalEnergyAccepted.size() ; j++)
            {
                float eCalSigmaAccepted = m_eCalToHadResolutionConstant * sqrt(m_eCalEnergyAccepted[j]);
                float hCalSigmaAccepted = m_hCalToHadResolutionConstant * sqrt(m_hCalEnergyAccepted[j]);
                
                perpDist = abs( ( (m_hCalEnergyAccepted[j] * eCalToHadIntercept) + (m_eCalEnergyAccepted[j] * hCalToHadIntercept) - (hCalToHadIntercept * eCalToHadIntercept) ) / sqrt( (hCalToHadIntercept * hCalToHadIntercept) + (eCalToHadIntercept * eCalToHadIntercept) ) );
                perpDistSigma = sqrt( ( pow (hCalSigmaAccepted * eCalToHadIntercept, 2) + pow(eCalSigmaAccepted * hCalToHadIntercept, 2) ) / ( pow(hCalToHadIntercept,2) + pow(eCalToHadIntercept,2) ) );

                Chi2 = Chi2 + pow(perpDist/perpDistSigma ,2);
            }

            if (Chi2 < m_minChi2)
            {
                m_eCalToHadInterceptMinChi2 = eCalToHadIntercept;
                m_hCalToHadInterceptMinChi2 = hCalToHadIntercept;
                m_minChi2 = Chi2;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ChiSquaredMethod::ThreeSigmaCut()
{
    float idealHCalToHadIntercept = m_kineticEnergy;
    float idealECalToHadIntercept = m_kineticEnergy;

    float pfoECalToHadEnergySigma = m_eCalToHadResolutionConstant * sqrt(m_pfoECalToHadEnergy);
    float pfoHCalToHadEnergySigma = m_hCalToHadResolutionConstant * sqrt(m_pfoHCalToHadEnergy);

    // Formula for perpendicular distance of point m_pfoE/HCalToHad to ideal distribution
    float perpDist = abs( ( (m_pfoHCalToHadEnergy * idealECalToHadIntercept) + (m_pfoECalToHadEnergy * idealHCalToHadIntercept) - (idealHCalToHadIntercept * idealECalToHadIntercept) ) / sqrt( (idealHCalToHadIntercept * idealHCalToHadIntercept) + (idealECalToHadIntercept * idealECalToHadIntercept) ) );
    float perpDistSigma = sqrt( ( pow (pfoHCalToHadEnergySigma * idealECalToHadIntercept, 2) + pow(pfoECalToHadEnergySigma * idealHCalToHadIntercept, 2) ) / ( pow(idealHCalToHadIntercept,2) + pow(idealECalToHadIntercept,2) ) );

    if ( perpDist <= 3 * perpDistSigma)
        return true;

    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::Plot()
{
    TCanvas *pCanvas = new TCanvas("2D Historgram","2D Histogram",200,10,1400,1400);
    pCanvas->cd();
    m_histogram->Draw("COLZ");

    pCanvas->SaveAs(m_outputPath + "PandoraPFA_Calibration_Hadronic_Energy_Scale_Chi_Sqaured_Method_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_KaonL.png");
    pCanvas->SaveAs(m_outputPath + "PandoraPFA_Calibration_Hadronic_Energy_Scale_Chi_Sqaured_Method_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_KaonL.C");
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], ChiSquaredMethod &chiSquaredMethod)
{
    int c(0);

    while (((c = getopt(argc, argv, "e:i:g:o:f")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'e':
            chiSquaredMethod.m_trueEnergy = atof(optarg);
            break;
        case 'i':
            chiSquaredMethod.m_inputKaonLRootFiles = optarg;
            break;
        case 'g':
            chiSquaredMethod.m_numberHCalLayers = atoi(optarg);
            break;
        case 'o':
            chiSquaredMethod.m_outputPath = optarg;
            break;
        case 'f':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -e value  (mandatory, true energy of KaonL being used for calibration)" << std::endl
                      << "    -i        (mandatory, input file name(s), can include wildcards if string is in quotes)" << std::endl
                      << "    -g value  (mandatory, number of HCal layers in simulation, default 48)" << std::endl
                      << "    -o value  (mandatory, output path to send results to)" << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}