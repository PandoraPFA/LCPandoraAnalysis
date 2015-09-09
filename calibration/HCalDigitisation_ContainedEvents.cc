/**
 *  @file   PandoraAnalysis/calibration/HCalDigitisation_ContainedEvents.cc
 * 
 *  @brief  Mean calo hit energy for kaonL events contained in HCal.  Used for setting digitisation constant in 
 *          HCal (CalibrHCALBarrel or CalibrHCALEndcap).
 * 
 *  $Log: $
 */

#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TROOT.h"

#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <limits>
#include <stdexcept>

/**
 *  @brief  HCalDigitisation class
 */
class HCalDigitisation 
{
public:
    /**
     *  @brief  Constructor
    */
    HCalDigitisation();

    /**
     *  @brief  Destructor
     */
    ~HCalDigitisation();

    /**
     *  @brief  Find the limits for the percentage fit and perform Gaussian fit over this range
    */
    void Process();

// Inputs Set By Parsing Command Line
    std::string     m_inputKaonLRootFiles;  ///< Input root files - KaonL
    float           m_trueEnergy;           ///< True energy (opposed to kinetic) of particle being simulated
    float           m_calibrationAccuracy;  ///< Fractional accuracy target for reconstructed energy
    std::string     m_outputPath;           ///< Output path to send results
    float           m_fitPercentage;        ///< Percentage of continuous data with narrowest range for Gaussian fit
    int             m_numberHCalLayers;     ///< Number of layers in the HCal
    std::string     m_element;              ///< Detector component (Barrel, EndCap or Other)
    float           m_lowerCosThetaCut;     ///< Lower cosTheta cut defining m_element
    float           m_upperCosTheraCut;     ///< Lower cosTheta cut defining m_element

// Outputs
    float           m_amplitude;            ///< Amplitude of Gaussian fit
    float           m_mean;                 ///< Mean of Gaussian fit
    float           m_stdDev;               ///< Standard deviation of Gaussian fit
    int             m_nEventsHCalHist;      ///< Number of events in m_histogram

private:
    /**
     *  @brief  Set m_maxHistogramEnergy
    */
    void PrepareHistogram();

    /**
     *  @brief  Create HCal calo hit energy histrogram (m_histogram)
    */
    void CreateHistogram();

    /**
     *  @brief  Fill HCal calo hit energy histrogram 
    */
    void FillHistogram();

    /**
     *  @brief  Set fitRangeLow, fitRangeHigh and RMSFitPercentage
    */
    void RMSFitPercentageRange();

    /**
     *  @brief  Perform Gaussian fit to histogram
    */
    void Fit();

typedef std::vector<float> FloatVector;

// Trivial setting on initialisation
    TChain         *m_pTChain;              ///< Chain of root files
    TH1F           *m_histogram;            ///< Histogram
    float           m_fitRangeLow;          ///< Low fit edge
    float           m_fitRangeHigh;         ///< High fit edge
    float           m_rMSFitRange;          ///< RMS of the fitRange data
    float           m_maxHistogramEnergy;   ///< Max energy on m_histogram
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 * 
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  hCalDigitisation to receive the application parameters
 * 
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], HCalDigitisation &hCalDigitisation);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();

    try
    {
        HCalDigitisation hCalDigitisation;

        if (!ParseCommandLine(argc, argv, hCalDigitisation))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        hCalDigitisation.Process();

        std::string dataFileName(hCalDigitisation.m_outputPath + "Calibration.txt");
        ofstream    data_file(dataFileName.c_str(), std::ios_base::app);
        data_file << "_____________________________________________________________________________________" << std::endl;
        std::cout << "_____________________________________________________________________________________" << std::endl;

        data_file << "Digitisation of the HCal " + hCalDigitisation.m_element + "                    : " << std::endl;
        data_file << "Element                                            : " << hCalDigitisation.m_element << std::endl;
        data_file << "Lower abs(cosTheta) cut defining element           : " << hCalDigitisation.m_lowerCosThetaCut << std::endl;
        data_file << "Upper abs(cosTheta) cut defining element           : " << hCalDigitisation.m_upperCosTheraCut << std::endl;
        data_file << "Continuous Fit percentage of data (narrowest RMS)  : " << hCalDigitisation.m_fitPercentage << std::endl;
        data_file << "Gaussian Fit Information:                          : " << std::endl;
        data_file << "Amplitude                                          : " << hCalDigitisation.m_amplitude << " : " <<std::endl;
        data_file << "HCal " + hCalDigitisation.m_element + " Digi Mean                              : " << hCalDigitisation.m_mean << " : " <<std::endl;
        data_file << "Standard Deviation                                 : " << hCalDigitisation.m_stdDev << " : " <<std::endl;
        data_file << "The total number of events considered was          : " << hCalDigitisation.m_nEventsHCalHist << " : " <<std::endl<<std::endl;

        std::cout << "Digitisation of the HCal " + hCalDigitisation.m_element + "                    : " << std::endl;
        std::cout << "Element                                            : " << hCalDigitisation.m_element << std::endl;
        std::cout << "Lower abs(cosTheta) cut defining element           : " << hCalDigitisation.m_lowerCosThetaCut << std::endl;
        std::cout << "Upper abs(cosTheta) cut defining element           : " << hCalDigitisation.m_upperCosTheraCut << std::endl;
        std::cout << "Continuous Fit percentage of data (narrowest RMS)  : " << hCalDigitisation.m_fitPercentage << std::endl;
        std::cout << "Gaussian Fit Information:                          : " << std::endl;
        std::cout << "Amplitude                                          : " << hCalDigitisation.m_amplitude << " : " <<std::endl;
        std::cout << "HCal " + hCalDigitisation.m_element + " Digi Mean                              : " << hCalDigitisation.m_mean << " : " <<std::endl;
        std::cout << "Standard Deviation                                 : " << hCalDigitisation.m_stdDev << " : " <<std::endl;
        std::cout << "The total number of events considered was          : " << hCalDigitisation.m_nEventsHCalHist << " : " <<std::endl<<std::endl;

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

HCalDigitisation::HCalDigitisation() :
    m_inputKaonLRootFiles(""),
    m_trueEnergy(0.f),
    m_calibrationAccuracy(0.02),
    m_outputPath(""),
    m_fitPercentage(90.f),
    m_numberHCalLayers(48),
    m_element(""),
    m_lowerCosThetaCut(0.f),
    m_upperCosTheraCut(1.f),
    m_amplitude(std::numeric_limits<float>::max()),
    m_mean(std::numeric_limits<float>::max()),
    m_stdDev(std::numeric_limits<float>::max()),
    m_nEventsHCalHist(0),
    m_pTChain(NULL),
    m_histogram(NULL),
    m_fitRangeLow(std::numeric_limits<float>::max()),
    m_fitRangeHigh(std::numeric_limits<float>::max()),
    m_rMSFitRange(std::numeric_limits<float>::max()),
    m_maxHistogramEnergy(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HCalDigitisation::~HCalDigitisation()
{
}

//--------------------------------------------------------------------------------------------------------------------------------------

void HCalDigitisation::Process()
{
    m_pTChain = new TChain("PfoAnalysisTree");
    m_pTChain->Add(m_inputKaonLRootFiles.c_str());
    this->PrepareHistogram();
    this->CreateHistogram();
    this->FillHistogram();
    this->RMSFitPercentageRange();
    this->Fit();
}

//--------------------------------------------------------------------------------------------------------------------------------------

void HCalDigitisation::PrepareHistogram()
{
    FloatVector *pfoTargetCosTheta(NULL);
    float totalCaloHitEnergy, hCalTotalCaloHitEnergy;
    int nPfoTargetsTotal, nPfoTargetsNeutralHadrons, pfoMinHCalLayerToEdge;

    m_pTChain->SetBranchAddress("TotalCaloHitEnergy",&totalCaloHitEnergy);
    m_pTChain->SetBranchAddress("HCalTotalCaloHitEnergy",&hCalTotalCaloHitEnergy);
    m_pTChain->SetBranchAddress("nPfoTargetsTotal",&nPfoTargetsTotal);
    m_pTChain->SetBranchAddress("nPfoTargetsNeutralHadrons",&nPfoTargetsNeutralHadrons);
    m_pTChain->SetBranchAddress("pfoTargetCosTheta", &pfoTargetCosTheta);
    m_pTChain->SetBranchAddress("pfoMinHCalLayerToEdge",&pfoMinHCalLayerToEdge);

    int containedLayerExclusion = ceil(m_numberHCalLayers * 0.1);

    for (unsigned int i = 0; i < m_pTChain->GetEntries(); i++) 
    {
        m_pTChain->GetEntry(i);

        bool isContained(containedLayerExclusion < pfoMinHCalLayerToEdge && (totalCaloHitEnergy-hCalTotalCaloHitEnergy) < (m_trueEnergy*0.05));

        if (1 == nPfoTargetsTotal && 1 == nPfoTargetsNeutralHadrons && isContained)
        {
            float cosTheta = std::fabs(pfoTargetCosTheta->at(0));

            if (m_lowerCosThetaCut < cosTheta && cosTheta < m_upperCosTheraCut)
            {
                if (hCalTotalCaloHitEnergy > m_maxHistogramEnergy)
                {
                    m_maxHistogramEnergy = hCalTotalCaloHitEnergy;
                }
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

void HCalDigitisation::CreateHistogram()
{
    std::string Name = "CaloHitEnergyHCal" + m_element;
    std::string Title = "Calorimeter Hit Energy HCal " + m_element + " (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)";
    int binNumber = static_cast<int>( m_maxHistogramEnergy / ( m_calibrationAccuracy * m_trueEnergy ) );
    m_histogram = new TH1F(Name.c_str(), Title.c_str(), binNumber, 0., m_maxHistogramEnergy);
    m_histogram->GetXaxis()->SetTitle("Calorimeter Hit Energy / GeV");
    m_histogram->GetYaxis()->SetTitle("Entries");
}

//--------------------------------------------------------------------------------------------------------------------------------------

void HCalDigitisation::FillHistogram()
{
    FloatVector *pfoTargetCosTheta(NULL);
    float totalCaloHitEnergy, hCalTotalCaloHitEnergy;
    int nPfoTargetsTotal, nPfoTargetsNeutralHadrons, pfoMinHCalLayerToEdge;

    m_pTChain->SetBranchAddress("TotalCaloHitEnergy",&totalCaloHitEnergy);
    m_pTChain->SetBranchAddress("HCalTotalCaloHitEnergy",&hCalTotalCaloHitEnergy);
    m_pTChain->SetBranchAddress("nPfoTargetsTotal",&nPfoTargetsTotal);
    m_pTChain->SetBranchAddress("nPfoTargetsNeutralHadrons",&nPfoTargetsNeutralHadrons);
    m_pTChain->SetBranchAddress("pfoTargetCosTheta", &pfoTargetCosTheta);
    m_pTChain->SetBranchAddress("pfoMinHCalLayerToEdge",&pfoMinHCalLayerToEdge);

    int containedLayerExclusion = ceil(m_numberHCalLayers * 0.1);

    for (unsigned int i = 0; i < m_pTChain->GetEntries(); i++) 
    {
        m_pTChain->GetEntry(i);

        bool isContained(containedLayerExclusion < pfoMinHCalLayerToEdge && (totalCaloHitEnergy-hCalTotalCaloHitEnergy) < (m_trueEnergy*0.05));

        if (1 == nPfoTargetsTotal && 1 == nPfoTargetsNeutralHadrons && isContained)
        {
            float cosTheta = std::fabs(pfoTargetCosTheta->at(0));

            if (m_lowerCosThetaCut < cosTheta && cosTheta < m_upperCosTheraCut)
            {
                m_histogram->Fill(hCalTotalCaloHitEnergy);
                m_nEventsHCalHist++;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HCalDigitisation::RMSFitPercentageRange()
{
    static const float FLOAT_MAX(std::numeric_limits<float>::max());

    if (NULL == m_histogram)
        return;

    if (5 > m_histogram->GetEntries())
    {
        std::cout << m_histogram->GetName() << " (" << m_histogram->GetEntries() << " entries) - skipped" << std::endl;
        return;
    }

    // Calculate raw properties of distribution (ie rms100)
    float sum = 0., total = 0.;
    double sx = 0., sxx = 0.;
    const unsigned int nbins(m_histogram->GetNbinsX());

    for (unsigned int i = 0; i <= nbins; ++i)
    {
        const float binx(m_histogram->GetBinLowEdge(i) + (0.5 * m_histogram->GetBinWidth(i)));
        const float yi(m_histogram->GetBinContent(i));
        sx += yi * binx;
        sxx += yi * binx * binx;
        total += yi;
    }

    const float rawMean(sx / total);
    const float rawMeanSquared(sxx / total);
    const float rawRms(std::sqrt(rawMeanSquared - rawMean * rawMean));

    sum = 0.;
    unsigned int is0 = 0;

    //  The /10 comes from the fact that for rms 90 the start point for the fit must occur in the first 10% of the data.
    float frac = (1 - (m_fitPercentage/100.0));
    for (unsigned int i = 0; (i <= nbins) && (sum < total * frac); ++i)
    {
        sum += m_histogram->GetBinContent(i);
        is0 = i;
    }

    // Calculate truncated properties
    float rmsmin(FLOAT_MAX), mean(FLOAT_MAX), low(FLOAT_MAX);
    float high(0.f);

    for (unsigned int istart = 0; istart <= is0; ++istart)
    {
        double sumn = 0.;
        double csum = 0.;
        double sumx = 0.;
        double sumxx = 0.;
        unsigned int iend = 0;

        for (unsigned int i = istart; (i <= nbins) && (csum < (m_fitPercentage/100) * total); ++i)
        {
            const float binx(m_histogram->GetBinLowEdge(i) + (0.5 * m_histogram->GetBinWidth(i)));
            const float yi(m_histogram->GetBinContent(i));
            //csum is the sum of yi from istart and is used to stop the sum when this exceeds X% of data.
            csum += yi;

            if (sumn < (m_fitPercentage/100) * total)
            {
                // These variables define the final sums required once we have considered X% of data, anything else is 
                // continuously overwritten.
                sumn += yi;
                sumx += yi * binx;
                sumxx+= yi * binx * binx;
                iend = i;
            }
        }

        const float localMean(sumx / sumn);
        const float localMeanSquared(sumxx / sumn);
        // Standard deviation formula
        const float localRms(std::sqrt(localMeanSquared - localMean * localMean));

        if (localRms < rmsmin)
        {
            mean = localMean;
            if (istart==0)
            {
                low = 0;
                m_fitRangeLow = 0;
            }
            else
            {
                low = m_histogram->GetBinLowEdge(istart);
                m_fitRangeLow=m_histogram->GetBinLowEdge(istart) + (0.5 * m_histogram->GetBinWidth(istart));
            }
            
            high = m_histogram->GetBinLowEdge(iend);
            rmsmin = localRms;
            m_fitRangeHigh=m_histogram->GetBinLowEdge(iend) + (0.5 * m_histogram->GetBinWidth(iend));
        }
    }
    
    m_rMSFitRange = rmsmin;
    
    std::cout << m_histogram->GetName() << " (" << m_histogram->GetEntries() << " entries), rawrms: " << rawRms << ", rmsx: " << rmsmin
              << " (" << low << "-" << high << "), low_fit and high_fit " << " (" << m_fitRangeLow << "-" << m_fitRangeHigh 
              << "), << mean: " << mean << std::endl;
}

//--------------------------------------------------------------------------------------------------------------------------------------

void HCalDigitisation::Fit()
{
    try
    {
        if (m_histogram->GetEntries() != 0)
        {
            // Fit Function
            TF1 *Gaussian_Fit_Func = new TF1("Gaussian_Fit_Func","[0] * TMath::Exp( -0.5 * [2] * TMath::Power(x-[1],2) )");

            // Initial Param Guess
            float Hist_Amp = m_histogram->GetBinContent(m_histogram->GetMaximumBin());
            float Hist_Mean = m_histogram->GetMean();
            float Hist_RMS = std::pow(m_histogram->GetRMS(),-2);

            Gaussian_Fit_Func->SetParameters(0, Hist_Amp);
            Gaussian_Fit_Func->SetParameters(1, Hist_Mean);
            Gaussian_Fit_Func->SetParameters(2, Hist_RMS);
            Gaussian_Fit_Func->SetParLimits(2,0,100);

            // Perform Fit and Get Results
            TFitResultPtr pTFitResultPtr = m_histogram->Fit(Gaussian_Fit_Func, "SM", "", m_fitRangeLow, m_fitRangeHigh);

            bool IsValidFit = pTFitResultPtr->IsValid();
            int FitQuality = pTFitResultPtr->CovMatrixStatus();
            int count = 0;

            // While loop until fit converges
            while (
            IsValidFit != 1 ||
            FitQuality != 3 ||
            Gaussian_Fit_Func->GetParameter(1) < m_fitRangeLow ||
            Gaussian_Fit_Func->GetParameter(1) > m_fitRangeHigh
            )
            {
                const float sign1(((static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX)) > 0.5f) ? 1.f : -1.f);
                const float sign2(((static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX)) > 0.5f) ? 1.f : -1.f);
                const float sign3(((static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX)) > 0.5f) ? 1.f : -1.f);

                float RescaleFactor1 = 1 + (0.5 * sign1 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float RescaleFactor2 = 1 + (0.5 * sign2 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float RescaleFactor3 = 1 + (0.5 * sign3 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));

                float LowerFraction1 = 1 - std::numeric_limits<float>::epsilon() - (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float LowerFraction2 = 1 - std::numeric_limits<float>::epsilon() - (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float LowerFraction3 = 1 - std::numeric_limits<float>::epsilon() - (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));

                float UpperFraction1 = 1 + std::numeric_limits<float>::epsilon() + (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float UpperFraction2 = 1 + std::numeric_limits<float>::epsilon() + (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float UpperFraction3 = 1 + std::numeric_limits<float>::epsilon() + (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));

                Gaussian_Fit_Func->SetParameters(0, RescaleFactor1 * Hist_Amp);
                Gaussian_Fit_Func->SetParameters(1, RescaleFactor2 * Hist_Mean);
                Gaussian_Fit_Func->SetParameters(2, RescaleFactor3 * Hist_RMS);

                Gaussian_Fit_Func->SetParLimits(0, LowerFraction1 * RescaleFactor1 * Hist_Amp, UpperFraction1 * RescaleFactor1 * Hist_Amp);
                Gaussian_Fit_Func->SetParLimits(1, LowerFraction2 * RescaleFactor2 * Hist_Mean, UpperFraction2 * RescaleFactor2 * Hist_Mean);
                Gaussian_Fit_Func->SetParLimits(2, LowerFraction3 * RescaleFactor3 * Hist_RMS, UpperFraction3 * RescaleFactor3 * Hist_RMS);

                pTFitResultPtr = m_histogram->Fit(Gaussian_Fit_Func, "SM", "", m_fitRangeLow, m_fitRangeHigh);

                IsValidFit = pTFitResultPtr->IsValid();
                FitQuality = pTFitResultPtr->CovMatrixStatus();
                count++;

                std::cout << "IsValidFit          :" << IsValidFit << std::endl;
                std::cout << "FitQuality          :" << FitQuality << std::endl;
                std::cout << "Count               :" << count << std::endl;

                if (count > 1000)
                    throw std::runtime_error("Fitting attempts exceeding maximum limit (1000).");
            }

            m_amplitude = Gaussian_Fit_Func->GetParameter(0);
            m_mean = Gaussian_Fit_Func->GetParameter(1);
            m_stdDev = std::pow(Gaussian_Fit_Func->GetParameter(2),-0.5);

            std::string canvasName = m_element + "Digitisation";
            std::string canvasTitle = m_element + " Digitisation";
            TCanvas *pCanvas = new TCanvas(canvasName.c_str(),canvasTitle.c_str(),200,10,600,500);
            pCanvas->cd();
            m_histogram->GetYaxis()->SetTitle("Entries");
            std::string xAxisTitle = "Calorimeter Hit Energy in HCal " + m_element + "/ GeV";
            m_histogram->GetXaxis()->SetTitle(xAxisTitle.c_str());
            m_histogram->Draw("");
            Gaussian_Fit_Func->Draw("same");

            std::string pngOutputFilename = m_outputPath + "Calorimeter_Hit_Energies_HCal_" + m_element + "_Digitisation.png";
            std::string dotCOutputFilename = m_outputPath + "Calorimeter_Hit_Energies_HCal_" + m_element + "_Digitisation.C";

            pCanvas->SaveAs(pngOutputFilename.c_str());
            pCanvas->SaveAs(dotCOutputFilename.c_str());
            delete pCanvas;
        }
        
        else
        {
            std::cout << "Histogram is empty, no fit possible." << std::endl;
        }
    }
    
    catch (const std::runtime_error& e)
    {
        std::cout << "An exception occurred. " << e.what() << '\n';
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], HCalDigitisation &hCalDigitisation)
{
    int c(0);

    while (((c = getopt(argc, argv, "a:b:c:d:e:f:g:i:j:k")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'a':
            hCalDigitisation.m_inputKaonLRootFiles = optarg;
            break;
        case 'b':
            hCalDigitisation.m_trueEnergy = atof(optarg);
            break;
        case 'c':
            hCalDigitisation.m_calibrationAccuracy = atof(optarg);
            break;
        case 'd':
            hCalDigitisation.m_outputPath = optarg;
            break;
        case 'e':
            hCalDigitisation.m_fitPercentage = atof(optarg);
            break;
        case 'f':
            hCalDigitisation.m_numberHCalLayers = atoi(optarg);
            break;
        case 'g':
            hCalDigitisation.m_element = optarg;
            break;
        case 'i':
            hCalDigitisation.m_lowerCosThetaCut = atof(optarg);
            break;
        case 'j':
            hCalDigitisation.m_upperCosTheraCut = atof(optarg);
            break;
        case 'k':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -a        (mandatory, input file name(s), can include wildcards if string is in quotes)           " << std::endl
                      << "    -b value  (mandatory, true energy of kaonL being used for calibration)                            " << std::endl
                      << "    -c value  (optional, fractional accuracy to calibrate CalibrHCAL to, default 0.02)                " << std::endl
                      << "    -d        (mandatory, output path to send results to)                                             " << std::endl
                      << "    -e value  (optional, fit percentage used for calibration, default 90% of data with narrowest rms) " << std::endl
                      << "    -f value  (optional, number of HCal layers in simulation, default 48)                             " << std::endl
                      << "    -g        (mandatory, element of the detector being calibrated (Barrel or EndCap))                " << std::endl
                      << "    -i value  (mandatory, lower abs cos theta cut defining element, default 0)                        " << std::endl
                      << "    -j value  (mandatory, upper abs cos theta cut defining element, default 1)                        " << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}
