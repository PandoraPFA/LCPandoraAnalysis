/**
 *  @file   PandoraAnalysis/calibration/PandoraPFA_SoftwareCompensation.cc
 * 
 *  @brief  Implementation of software compensation binary. Used for extracting software compensation weights with a minimizer
 * 
 *  $Log: $
 */

#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TTree.h"
#include "TStyle.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Fit/ParameterSettings.h"

#include "AnalysisHelper.h"

#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

using namespace pandora_analysis;

bool devDebug = true;

// Global typedefs
typedef std::vector<float>        FloatVector;
typedef std::vector<double>       DoubleVector;
typedef std::vector<int>          IntVector;
typedef std::vector<std::string>  StringVector;
typedef std::vector<TH1F*>        TH1FVector;
typedef std::vector<TProfile*>    TProfileVector;

class SoftwareCompensation;

/**
 *  @brief  Event class
 */
class Event
{
public:
    /**
     *  @brief  Constructor
     *  
     *  @param softwareCompensation the software compensation instance from which to read inputs
     */
    Event(const SoftwareCompensation &softwareCompensation);
  
public:
    int            m_firstPseudoLayer;        ///< Pseudolayer at the IP
    float          m_trueEnergy;              ///< The true (mc or beam) kaon0L energy (unit GeV)
    float          m_reconstructedEnergy;     ///< The kaon0L reconstructed energy (unit GeV)
    float          m_eCalEnergy;              ///< Sum of the kaon0L ecal hit energies 
    FloatVector    m_hCalBinEnergies;         ///< HCal hit energies
    FloatVector    m_hCalHitEnergies;         ///< HCal hit energy bins
    IntVector      m_caloHitPseudoLayers;     ///< Pseudolayers for the calorimeter hits
    IntVector      m_caloHitIsIsolated;       ///< Is calorimeter hit isolated
    IntVector      m_caloHitType;             ///< Type of calorimeter hits
    FloatVector    m_caloHitEnergies;         ///< Energy of the calorimeter hits
};

typedef std::vector<Event>       EventVector;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SoftwareCompensation class
 */
class SoftwareCompensation
{
public:
  /**
   *  @brief  Constructor
   */
    SoftwareCompensation();
    
    /**
     *  @brief  Destructor
     */
    ~SoftwareCompensation();
    
    /**
     *  @brief  Set true energies to read from input files
     *  
     *  @param  energiesStr energy list as string, comma separated (i.e '10:20:30:40')
     */
    void SetTrueEnergies(const std::string &energiesStr);
    
    /**
     *  @brief  Process the minimization
     */
    void Process();
    
    /**
     *  @brief  Get the calo hit energy bin given a hit energy
     *  
     *  @param  hitEnergy  the calo hit energy (unit GeV)
     *  @param  cellVolume the cell volume (unit mm)
     *  
     *  @return            [description]
     */
    int FindCaloHitBin(const float hitEnergy, const float cellVolume) const;
    
    /**
     *  @brief  Get the number calo hit density bins
     *  
     *  @return The number calo hit density bins
     */
    int GetNDensityBins() const;
    
    /**
     *  @brief  The minuit function used for minimizing the chi2 for software compensation
     * 
     *  @param  params the software compensation weights (usually 9)
     * 
     *  @return  chi2
     */
    double MinuitChi2(const double *params);
    
    /**
     *  @brief  Whether to skip the current event while reading the ROOT TChain
     *  
     *  @return skip or not
     */
    bool SkipCurrentEvent() const;
    
    /**
     *  @brief  Get the compensated energy for a particular event
     *  
     *  @param  event     the event to consider
     *  @param  params    the software compensation weights
     *  @param  pfoEnergy the pfo energy to consider. Might be the true beam energy or first estimate of hadronic shower energy
     *  
     *  @return           the compensated energy
     */
    double GetCompensatedEnergy(const Event &event, const double *params, double pfoEnergy);

    /**
     *  @brief Apply the CleanClusters logic to the cluster energy as it appears in the LCContent library
     *
     *  @param  event                 the event to consider
     *  @param  clusterHadronicEnergy the hadronic energy of the cluster
     */
    void CleanCluster(const Event &event, double &correctedHadronicEnergy);

    /**
     *  @brief Find the hadronic energy in a given pseudolayer
     *
     *  @param  event       the event to consider
     *  @param  pseudoLayer the pseudoLayer to consider
     */
   float GetHadronicEnergyInLayer(const Event &event, const unsigned int pseudoLayer) const;

    /**
     *  @brief  Draw the list of plots on a canvas
     *  
     *  @param  plots       the list of plots (must be TH1* like)
     *  @param  pTLegend    the preconfigured legend (without entries) to draw
     *  @param  canvasTitle the canvas title
     *  @param  xTitle      the X axis title
     *  
     *  @return             the canvas on which to plots are drawn
     */
    template <typename T>
    TCanvas *DrawPlots(const std::vector<T*> &plots, TLegend *pTLegend, const std::string &canvasTitle, const std::string &xTitle) const;

    
    // Inputs Set By Parsing Command Line
    std::string                  m_filePattern;                ///< The global file pattern for reading input data (root file) - Kaon0L
    std::vector<std::string>     m_trueEnergies;               ///< The list of true (mc or beam) energies
    std::string                  m_trueEnergiesStr;            ///< The list of true (mc or beam) energies - as string
    std::string                  m_outputPath;                 ///< Output path to send results
    std::string                  m_treeName;                   ///< The root tree name

    int                          m_firstPseudoLayer;           ///< Pseudolayers at the IP
    FloatVector                  m_densityBinEdges;            ///< List of energy bin edges (size N)
    FloatVector                  m_densityBins;                ///< List of energy bins (size N-1) 
    int                          m_nMinuitParameters;          ///< The number of minuit parameters (9)
    bool                         m_minimizeUsingTrueEnergy;    ///< Whether to minimize the chi2 using true energy or reconstructed energy
    DoubleVector                 m_softCompGuessParameters;    ///< The guess parameters for software compensation provided to the ROOT minimizer (9)
    DoubleVector                 m_softCompFinalParameters;    ///< The final software compensation parameters after minimization (9)
    EventVector                  m_eventVector;                ///< The event list after reading input data

private:
    /**
     *  @brief  Set the TChain branch addresses
     */
    void SetBranchAddresses();
    
    /**
     *  @brief  Read the root input files - Kaon0L
     */
    void ReadInputFiles();
    
    /**
     *  @brief  Create the monitoring histograms - pre minimization step
     */
    void CreateHistograms();
    
    /**
     *  @brief  Perform the minimization
     */
    void PerformMinimization();
    
    /**
     *  @brief  Fill the histograms - post minimization step
     */
    void FillHistograms();
    
    /**
     *  @brief  Save the filled histograms in a root file
     */
    void SaveHistograms();
      
    /**
     *  @brief  Factory method to create a canvas
     *  
     *  @param  title the canvas title
     *  
     *  @return canvas
     */
    TCanvas *CreateCanvas(const std::string &title) const;
    
private:
    TChain         *m_pTChain;                ///< Chain of root files
    float           m_rawClusterEnergy;       ///< Raw cluster energy
    FloatVector    *m_hitEnergies;            ///< All hit energies (ECal + HCal + Muon)
    FloatVector    *m_cellSize0;              ///< Cell sizes 0
    FloatVector    *m_cellSize1;              ///< Cell sizes 1
    FloatVector    *m_cellThickness;          ///< Cell thicknesses
    IntVector      *m_hitType;                ///< Hit types (ECal, HCal or Muon)
    IntVector      *m_caloHitPseudoLayers;    ///< Pseudolayers of the calorimeter hits
    IntVector      *m_caloHitIsIsolated;      ///< Is calorimeter hit isolated

    float           m_trueEnergy;             ///< The true energy while reading input files 
    float           m_minCleanHitEnergy;          ///< Min calo hit hadronic energy to consider cleaning hit/cluster
    float           m_minCleanHitEnergyFraction;  ///< Min fraction of cluster energy represented by hit to consider cleaning
    float           m_minCleanCorrectedHitEnergy; ///< Min value of new hit hadronic energy estimate after cleaning

    TGraph         *m_gChi2Evolution;         ///< The chi2 evolution graph
    TH1FVector      m_hECalEnergyVector;      ///< The list of ECal total hit energy plots (1 per true energy)
    TProfileVector  m_hHCalEnergyVector;      ///< The list of HCal hit energy sum plots (1 per true energy)
    TH1FVector      m_hHCalHitEnergyVector;   ///< The list of HCal hit energy plots (1 per true energy)
    TH1FVector      m_hTotalEnergyVector;     ///< The list of total pfo energy plots (1 per true energy)
    TH1FVector      m_hSoftCompEnergyVector;  ///< The list of software compensated energy plots (1 per true energy)
    
    SoftwareCompensation(const SoftwareCompensation &) = delete;
    SoftwareCompensation &operator=(const SoftwareCompensation &) = delete;
    friend class Event;
};


//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 * 
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  softwareCompensation to receive the application parameters
 * 
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], SoftwareCompensation &softwareCompensation);
void Tokenize(const std::string &inputString, std::vector<std::string> &tokens, const std::string &delimiter = ":");
std::string &Replace(std::string &str, const std::string &variable, const std::string &value);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    try
    {
        SoftwareCompensation softwareCompensation;

        if (!ParseCommandLine(argc, argv, softwareCompensation))
            return 1;

        int myargc = 0;
        char* myargv = (char *)"";
        TApplication *pApplication = new TApplication("PandoraMonitoring", &myargc, &myargv);
        pApplication->SetReturnFromRun(kTRUE);
        gStyle->SetOptStat(0);

        softwareCompensation.Process();
        
        std::string dataFileName(softwareCompensation.m_outputPath + "Calibration.txt");
        std::ofstream data_file(dataFileName.c_str(), std::ios_base::app);

        data_file << "_____________________________________________________________________________________" << std::endl;
        std::cout << "_____________________________________________________________________________________" << std::endl;
        
        data_file << "Software compensation weight determination         : " << std::endl;
        data_file << "File pattern                                       : " << softwareCompensation.m_filePattern << std::endl;
        data_file << "Input energies                                     : " << softwareCompensation.m_trueEnergiesStr << std::endl;
        data_file << "True energy used                                   : " << (softwareCompensation.m_minimizeUsingTrueEnergy ? "true" : "false") << std::endl;
        data_file << "Minimizer                                          : Minuit (Migrad)" << std::endl;
        data_file << "Software compensation parameters                   : " << std::endl;
        for (unsigned int p = 0 ; p < softwareCompensation.m_softCompFinalParameters.size() ; p++)
            data_file << "Parameter " << p << "                                        : " << softwareCompensation.m_softCompFinalParameters.at(p) << std::endl;
        data_file << "The total number of events considered was          : " << softwareCompensation.m_eventVector.size() << std::endl;
        
        std::cout << "Software compensation weight determination         : " << std::endl;
        std::cout << "File pattern                                       : " << softwareCompensation.m_filePattern << std::endl;
        std::cout << "Input energies                                     : " << softwareCompensation.m_trueEnergiesStr << std::endl;
        std::cout << "True energy used                                   : " << (softwareCompensation.m_minimizeUsingTrueEnergy ? "true" : "false") << std::endl;
        std::cout << "Minimizer                                          : Minuit (Migrad)" << std::endl;
        std::cout << "Software compensation parameters                   : " << std::endl;
        for (unsigned int p = 0 ; p < softwareCompensation.m_softCompFinalParameters.size() ; p++)
            std::cout << "Parameter " << p << "                                        : " << softwareCompensation.m_softCompFinalParameters.at(p) << std::endl;
        std::cout << "The total number of events considered was          : " << softwareCompensation.m_eventVector.size() << std::endl;
        
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
//------------------------------------------------------------------------------------------------------------------------------------------

Event::Event(const SoftwareCompensation &softwareCompensation) :
    m_firstPseudoLayer(0),
    m_trueEnergy(softwareCompensation.m_trueEnergy),
    m_reconstructedEnergy(softwareCompensation.m_rawClusterEnergy),
    m_eCalEnergy(0.f),
    m_hCalBinEnergies(),
    m_hCalHitEnergies(),
    m_caloHitPseudoLayers(),
    m_caloHitIsIsolated(),
    m_caloHitType(),
    m_caloHitEnergies()
{
    const unsigned int nHits(softwareCompensation.m_hitEnergies->size());
    
    m_hCalBinEnergies.clear();
    m_hCalBinEnergies.resize(softwareCompensation.GetNDensityBins(), 0.f);
    m_hCalHitEnergies.reserve(nHits);
    
    for (unsigned int i = 0 ; i<nHits ; i++)
    {
        const float hitEnergy(softwareCompensation.m_hitEnergies->at(i));
        const int hitType(softwareCompensation.m_hitType->at(i));
        const bool invalidEnergy(hitEnergy <= 0 || hitEnergy >= 10000);

        m_caloHitPseudoLayers.push_back(softwareCompensation.m_caloHitPseudoLayers->at(i));
        m_caloHitIsIsolated.push_back(softwareCompensation.m_caloHitIsIsolated->at(i));
        m_caloHitEnergies.push_back(hitEnergy);
        m_caloHitType.push_back(hitType);

        if (invalidEnergy)
            continue;
        
        if(2 == hitType) // hCal hits
        {
            m_hCalHitEnergies.push_back(hitEnergy);
            
            const float cellSize0(softwareCompensation.m_cellSize0->at(i));
            const float cellSize1(softwareCompensation.m_cellSize1->at(i));
            const float cellThickness(softwareCompensation.m_cellThickness->at(i));
            const float cellVolume((cellSize0*cellSize1*cellThickness) / 1000000.f);
            
            const int bin(softwareCompensation.FindCaloHitBin(hitEnergy, cellVolume));
            m_hCalBinEnergies[bin] += hitEnergy;
        }
        else if(1 == hitType)
        {
            m_eCalEnergy += hitEnergy;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

SoftwareCompensation::SoftwareCompensation() :
    m_filePattern(""),
    m_trueEnergies(),
    m_trueEnergiesStr(""),
    m_outputPath(""),
    m_treeName("SoftwareCompensationTrainingTree"),
    m_firstPseudoLayer(0),
    m_densityBinEdges(),
    m_densityBins(),
    m_nMinuitParameters(0),
    m_minimizeUsingTrueEnergy(false),
    m_softCompGuessParameters(),
    m_softCompFinalParameters(),
    m_eventVector(),
    m_pTChain(NULL),
    m_rawClusterEnergy(0.f),
    m_hitEnergies(NULL),
    m_cellSize0(NULL),
    m_cellSize1(NULL),
    m_cellThickness(NULL),
    m_hitType(NULL),
    m_caloHitPseudoLayers(NULL),
    m_caloHitIsIsolated(NULL),
    m_trueEnergy(0.f),
    m_minCleanHitEnergy(0.5f),
    m_minCleanHitEnergyFraction(0.01f),
    m_minCleanCorrectedHitEnergy(0.1f),
    m_gChi2Evolution(NULL),
    m_hECalEnergyVector(),
    m_hHCalEnergyVector(),
    m_hHCalHitEnergyVector(),
    m_hTotalEnergyVector(),
    m_hSoftCompEnergyVector()
{
    m_densityBinEdges = {0,  2.,  5., 7.5, 9.5,  13., 16.,  20., 23.5,  28., 1e6};
    m_densityBins.resize(m_densityBinEdges.size()-1);
    m_softCompGuessParameters = {2.4, -0.06, 0.0008, -0.09, -0.004, -0.00008, 0.05, 0.07, -0.1};
    m_nMinuitParameters = m_softCompGuessParameters.size();
    m_softCompFinalParameters.resize(m_nMinuitParameters);
    
    for (unsigned int i = 0 ; i<m_densityBinEdges.size()-1 ; i++)
    {
        m_densityBins[i] = (i == m_densityBinEdges.size()-2) ? 30 : (m_densityBinEdges[i] + m_densityBinEdges[i+1])/2.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

SoftwareCompensation::~SoftwareCompensation()
{
    if (NULL != m_gChi2Evolution) delete m_gChi2Evolution;
    
    for (unsigned int i = 0 ; i < m_hECalEnergyVector.size() ; i++) delete m_hECalEnergyVector.at(i);
    for (unsigned int i = 0 ; i < m_hHCalEnergyVector.size() ; i++) delete m_hHCalEnergyVector.at(i);
    for (unsigned int i = 0 ; i < m_hHCalHitEnergyVector.size() ; i++) delete m_hHCalHitEnergyVector.at(i);
    for (unsigned int i = 0 ; i < m_hTotalEnergyVector.size() ; i++) delete m_hTotalEnergyVector.at(i);
    for (unsigned int i = 0 ; i < m_hSoftCompEnergyVector.size() ; i++) delete m_hSoftCompEnergyVector.at(i);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SoftwareCompensation::SetTrueEnergies(const std::string &energiesStr)
{
    m_trueEnergiesStr = energiesStr;
    m_trueEnergies.clear();
    Tokenize(energiesStr, m_trueEnergies, ":");
}

//------------------------------------------------------------------------------------------------------------------------------------------

int SoftwareCompensation::GetNDensityBins() const
{ 
    return m_densityBins.size(); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SoftwareCompensation::SkipCurrentEvent() const
{
    if(m_hitEnergies->empty())
        return true;
        
    const bool invalidSingleParticle(m_rawClusterEnergy <= 0);
    
    if(invalidSingleParticle)
        return true;
    
    bool hasMuonHits(false);
    bool hasHCalHits(false);
    
    for (unsigned int t = 0 ; t < m_hitType->size() ; t++)
    {
        if(3 == m_hitType->at(t))
        {
            hasMuonHits = true;
            break;
        }
        
        if(2 == m_hitType->at(t))
        {
            hasHCalHits = true;
        }
    }
    
    if (hasMuonHits || !hasHCalHits)
    {    
        return true;
    }
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SoftwareCompensation::Process()
{
    this->ReadInputFiles();
    this->CreateHistograms();
    this->PerformMinimization();
    this->FillHistograms();
    this->SaveHistograms();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SoftwareCompensation::ReadInputFiles()
{
    m_eventVector.clear();
    
    for (unsigned int e = 0 ; e < m_trueEnergies.size() ; e++)
    {
        const std::string energyStr(m_trueEnergies.at(e));
        std::stringstream ss(energyStr);
        m_trueEnergy = 0.f;
        std::string fileName(m_filePattern);
        Replace(fileName, "energy", energyStr);
        
        if( (ss >> m_trueEnergy).fail() )
            throw std::runtime_error("SoftwareCompensation::ReadInputFiles: Couldn't convert energy string to float !");
        
        if (devDebug)
            std::cout << "File name : " << fileName << " for energy " << m_trueEnergy << " GeV" << std::endl;
        
        m_pTChain = new TChain(m_treeName.c_str());
        m_pTChain->Add(fileName.c_str());
        this->SetBranchAddresses();
        int nEventsForThisEnergy(0);
        
        for (unsigned int i = 0; i < m_pTChain->GetEntries(); i++)
        {
            m_pTChain->GetEntry(i);
            
            if(this->SkipCurrentEvent())
                continue;

            m_eventVector.push_back(Event(*this));
            nEventsForThisEnergy++;
        }
        
        if (devDebug)
            std::cout << "Got " << nEventsForThisEnergy << " events ..." << std::endl;
        
        delete m_pTChain;
    }
    
    if (devDebug)
        std::cout << "Got " << m_eventVector.size() << " events over " << m_trueEnergies.size() << " energies" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SoftwareCompensation::PerformMinimization()
{
    ROOT::Math::Minimizer* pMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    pMinimizer->SetMaxFunctionCalls(10000000);
    pMinimizer->SetMaxIterations(1000000);
    pMinimizer->SetTolerance(0.001);
    pMinimizer->SetPrintLevel(1);

    ROOT::Math::Functor minimizerFunctor(this, &SoftwareCompensation::MinuitChi2, m_nMinuitParameters);
    pMinimizer->SetFunction(minimizerFunctor);
    
    for (unsigned int i = 0 ; i < m_nMinuitParameters ; i++)
    {
        std::stringstream name; name << "p" << i;
        pMinimizer->SetVariable(i, name.str().c_str(), m_softCompGuessParameters[i], 0.01);
        std::cout << "Initial parameter '" << name.str() << "' value is " << m_softCompGuessParameters[i] << std::endl;
    }
    
    pMinimizer->Minimize();
    const double *finalParameters(pMinimizer->X());
    std::cout << "Reached minimum with Chi2 : " << pMinimizer->MinValue() << std::endl;
    
    for (unsigned int i = 0 ; i < m_nMinuitParameters ; i++)
    {
        ROOT::Fit::ParameterSettings parameterSettings;
        pMinimizer->GetVariableSettings(i, parameterSettings);
        std::cout << "  => Parameter " << parameterSettings.Name() << ", final value : " << finalParameters[i] << std::endl;
        m_softCompFinalParameters[i] = finalParameters[i];
    }
    
    delete pMinimizer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SoftwareCompensation::CreateHistograms()
{
    m_gChi2Evolution = new TGraph();
    m_gChi2Evolution->SetTitle("Chi2 Evolution");
    
    for (unsigned int e = 0 ; e < m_trueEnergies.size() ; e++)
    {
        const std::string energyStr(m_trueEnergies.at(e));
        std::stringstream ss(energyStr);
        std::string hName, hTitle;
        const int nDensityBins(this->GetNDensityBins());
        
        hName = "ECalEnergy" + energyStr + "GeV";
        hTitle = "ECal energy";
        TH1F *pECalEnergyHistogram = new TH1F(hName.c_str(), hTitle.c_str(), 500, 0, 100);
        
        hName = "HCalEnergy" + energyStr + "GeV";
        hTitle = "HCal energy";
        TProfile *pHCalEnergyHistogram = new TProfile(hName.c_str(), hTitle.c_str(), nDensityBins, 0, nDensityBins-1);
        
        hName = "HCalHitEnergy" + energyStr + "GeV";
        hTitle = "HCal hit energy";
        TH1F *pHCalHitEnergyHistogram = new TH1F(hName.c_str(), hTitle.c_str(), 600, 0, 3);
        
        hName = "TotalEnergy" + energyStr + "GeV";
        hTitle = "Total energy";
        TH1F *pTotalEnergyHistogram = new TH1F(hName.c_str(), hTitle.c_str(), 600, 0, 120);
        
        hName = "CompensatedEnergy" + energyStr + "GeV";
        hTitle = "Compensated energy";
        TH1F *pSoftCompEnergyHistogram = new TH1F(hName.c_str(), hTitle.c_str(), 600, 0, 120);
        
        pECalEnergyHistogram->SetDirectory(0); m_hECalEnergyVector.push_back(pECalEnergyHistogram);
        pHCalEnergyHistogram->SetDirectory(0); m_hHCalEnergyVector.push_back(pHCalEnergyHistogram);
        pHCalHitEnergyHistogram->SetDirectory(0); m_hHCalHitEnergyVector.push_back(pHCalHitEnergyHistogram);
        pTotalEnergyHistogram->SetDirectory(0); m_hTotalEnergyVector.push_back(pTotalEnergyHistogram);
        pSoftCompEnergyHistogram->SetDirectory(0); m_hSoftCompEnergyVector.push_back(pSoftCompEnergyHistogram);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SoftwareCompensation::FillHistograms()
{
    for (unsigned int e = 0 ; e < m_trueEnergies.size() ; e++)
    {
        const std::string energyStr(m_trueEnergies.at(e));
        std::stringstream ss(energyStr);
        float trueEnergy = 0.f;
        std::string fileName(m_filePattern);
        Replace(fileName, "energy", energyStr);
        
        if( (ss >> trueEnergy).fail() )
            throw std::runtime_error("SoftwareCompensation::FillHistograms: Couldn't convert energy string to float !");
        
        m_pTChain = new TChain(m_treeName.c_str());
        m_pTChain->Add(fileName.c_str());
        this->SetBranchAddresses();
        
        TH1F     *pECalEnergyHistogram(m_hECalEnergyVector.at(e));
        TProfile *pHCalEnergyHistogram(m_hHCalEnergyVector.at(e));
        TH1F     *pHCalHitEnergyHistogram(m_hHCalHitEnergyVector.at(e));
        TH1F     *pTotalEnergyHistogram(m_hTotalEnergyVector.at(e));
        TH1F     *pSoftCompEnergyHistogram(m_hSoftCompEnergyVector.at(e));
        
        for (unsigned int i = 0; i < m_pTChain->GetEntries(); i++)
        {
            m_pTChain->GetEntry(i);
            
            if(this->SkipCurrentEvent())
                continue;

            Event event(*this);
            float totalEnergy(0.f);
            
            if(event.m_eCalEnergy > 0.f)
                pECalEnergyHistogram->Fill(event.m_eCalEnergy);
              
            totalEnergy += event.m_eCalEnergy;
            
            for (unsigned int h = 0 ; h < event.m_hCalHitEnergies.size() ; h++)
            {
                const float hitEnergy(event.m_hCalHitEnergies.at(h));
                
                if(hitEnergy > 0.f)
                    pHCalHitEnergyHistogram->Fill(hitEnergy); // FIXME original code if Fill(hitEnergy/0.025). Why ???
            }
            
            for (unsigned int bin = 0 ; bin < event.m_hCalBinEnergies.size() ; bin++)
            {
                const float binEnergy(event.m_hCalBinEnergies.at(bin));
                pHCalEnergyHistogram->Fill(bin, binEnergy);
                totalEnergy += binEnergy;
            }
            
            pTotalEnergyHistogram->Fill(totalEnergy);
            
            const double *params(&m_softCompFinalParameters[0]);
            const double softCompEnergy(this->GetCompensatedEnergy(event, params, event.m_reconstructedEnergy));
            pSoftCompEnergyHistogram->Fill(softCompEnergy);
        }
        
        delete m_pTChain;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SoftwareCompensation::SaveHistograms()
{    
    TCanvas *pChi2Canvas = this->CreateCanvas("Chi2 evolution");
    pChi2Canvas->SetTitle("Chi2 evolution");
    pChi2Canvas->cd();
    m_gChi2Evolution->Draw("al");
    m_gChi2Evolution->GetXaxis()->SetTitle("Call");
    pChi2Canvas->Update();

    TLegend *pTLegend = new TLegend(0.6, 0.4, 0.9, 0.85);    
    pTLegend->SetLineStyle(0);
    pTLegend->SetLineWidth(0);
    pTLegend->SetLineColor(0);
    
    TCanvas *pECalCanvas           = this->DrawPlots(m_hECalEnergyVector,     pTLegend, "ECal energy",                 "ECal energy [GeV]");
    TCanvas *pHCalCanvas           = this->DrawPlots(m_hHCalEnergyVector,     pTLegend, "HCal bin energy",             "Bin energy [GeV]");
    TCanvas *pHCalHitCanvas        = this->DrawPlots(m_hHCalHitEnergyVector,  pTLegend, "HCal hit energy",             "HCal hit energy [GeV]");
    TCanvas *pTotalEnergyCanvas    = this->DrawPlots(m_hTotalEnergyVector,    pTLegend, "Total energy",                "Reconstructed energy [GeV]");
    TCanvas *pSoftCompEnergyCanvas = this->DrawPlots(m_hSoftCompEnergyVector, pTLegend, "Software compensated energy", "Software compensated energy [GeV]");
    
    //-------------
    // Write output
    TFile *pTFile = new TFile("SoftwareCompensationMonitoring.root", "RECREATE");
    pTFile->cd();
    pChi2Canvas->Write();
    pECalCanvas->Write();
    pHCalCanvas->Write();
    pHCalHitCanvas->Write();
    pTotalEnergyCanvas->Write();
    pSoftCompEnergyCanvas->Write();
    pTFile->Close();
    
    delete pTLegend;
    delete pChi2Canvas;
    delete pECalCanvas;
    delete pHCalCanvas;
    delete pHCalHitCanvas;
    delete pTotalEnergyCanvas;
    delete pSoftCompEnergyCanvas;
    delete pTFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TCanvas *SoftwareCompensation::DrawPlots(const std::vector<T*> &plots, TLegend *pTLegend, const std::string &canvasTitle, const std::string &xTitle) const
{
    assert(plots.size() == m_trueEnergies.size());
    double maximum(std::numeric_limits<double>::min());
    TCanvas *pTCanvas = this->CreateCanvas(canvasTitle);
    pTCanvas->cd();
    pTLegend->Clear();
    
    for (unsigned int p = 0 ; p < plots.size() ; p++)
    {   
        TH1 *pHistogram(plots.at(p));
        pHistogram->SetLineWidth(2);
        pHistogram->SetLineColor(p+1);
        const std::string legendLabel(m_trueEnergies.at(p) + " GeV");
        const std::string drawOption(p == 0 ? "" : "same");
        pTLegend->AddEntry(pHistogram, legendLabel.c_str(), "l");
        pHistogram->GetXaxis()->SetTitle(xTitle.c_str());
        pHistogram->DrawNormalized(drawOption.c_str());        
        const double localMaximum(pHistogram->GetMaximum()/pHistogram->Integral());
        
        if(localMaximum > maximum)
            maximum = localMaximum;
    }
    
    for (unsigned int p = 0 ; p < plots.size() ; p++)
    {
        TH1 *pHistogram(plots.at(p));
        pHistogram->GetYaxis()->SetRangeUser(0, maximum*1.1);
    }
    
    pTLegend->Draw();
    pTCanvas->Update();
    return pTCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SoftwareCompensation::SetBranchAddresses()
{
    m_pTChain->SetBranchAddress("FirstPseudoLayer",&m_firstPseudoLayer);
    m_pTChain->SetBranchAddress("RawEnergyOfCluster",&m_rawClusterEnergy);
    m_pTChain->SetBranchAddress("HitEnergies",&m_hitEnergies);
    m_pTChain->SetBranchAddress("CellSize0",&m_cellSize0);
    m_pTChain->SetBranchAddress("CellSize1", &m_cellSize1);
    m_pTChain->SetBranchAddress("CellThickness", &m_cellThickness);
    m_pTChain->SetBranchAddress("HitType", &m_hitType);
    m_pTChain->SetBranchAddress("PseudoLayer", &m_caloHitPseudoLayers);
    m_pTChain->SetBranchAddress("IsIsolated", &m_caloHitIsIsolated);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int SoftwareCompensation::FindCaloHitBin(const float hitEnergy, const float cellVolume) const
{
    const float hitEnergyDensity(hitEnergy/cellVolume);
    
    for (unsigned int i = 0 ; i<this->GetNDensityBins() ; i++)
    {
        if ( (m_densityBinEdges[i] <= hitEnergyDensity) && (hitEnergyDensity < m_densityBinEdges[i+1]) )
            return i;
    }
    
    throw std::runtime_error("SoftwareCompensation::FindCaloHitBin: Bin not found !");
}

//------------------------------------------------------------------------------------------------------------------------------------------

double SoftwareCompensation::MinuitChi2(const double *params)
{
    static int nChi2Calls(0);
    const unsigned int nEvents(m_eventVector.size());
    double chi2Total(0.f);
        
    for (unsigned int i = 0 ; i < nEvents ; i++)
    {
        try
        {
          Event &event(m_eventVector.at(i));
          const double pfoEnergy(m_minimizeUsingTrueEnergy ? event.m_trueEnergy : event.m_reconstructedEnergy);
          const double softCompTotalEnergy(this->GetCompensatedEnergy(event, params, pfoEnergy));
          const double trueEnergy(static_cast<double>(event.m_trueEnergy));          
          const double chi2Event((softCompTotalEnergy-trueEnergy)*(softCompTotalEnergy-trueEnergy)/(0.5*trueEnergy));
          
          if (std::numeric_limits<double>::max() - chi2Total < chi2Event)
              throw std::runtime_error("Beyond numeric limits !");
          
          chi2Total += chi2Event;
        }
        catch(std::runtime_error & e)
        {
          return std::numeric_limits<double>::max();
        }
    }
    
    if (NULL != m_gChi2Evolution) 
        m_gChi2Evolution->SetPoint(nChi2Calls, nChi2Calls, chi2Total/nEvents);
    
    nChi2Calls++;
    
    if (nChi2Calls % 100 == 0)
        std::cout << "Chi2 (call no " << nChi2Calls << "): total = " << chi2Total << " , chi2/nEvents = " << chi2Total/nEvents  << std::endl;
    
    return chi2Total/nEvents;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double SoftwareCompensation::GetCompensatedEnergy(const Event &event, const double *params, double pfoEnergy)
{
    static const double MAX_EXP_ARG(std::log(std::numeric_limits<double>::max()));
    double compensatedEnergy(0);
    compensatedEnergy += static_cast<double>(event.m_eCalEnergy);

    if(params[8]*pfoEnergy >= MAX_EXP_ARG)
        throw std::runtime_error("Beyond numeric limits !");
    
    for (unsigned int bin = 0 ; bin < event.m_hCalBinEnergies.size() ; bin++)
    {
        const double binEnergy(static_cast<double>(event.m_hCalBinEnergies.at(bin)));
        const double binDensity(static_cast<double>(m_densityBins[bin]));
        const double softCompWeight1(params[0] + params[1]*pfoEnergy + params[2]*pfoEnergy*pfoEnergy);
        const double softCompWeight2(params[3] + params[4]*pfoEnergy + params[5]*pfoEnergy*pfoEnergy);
        const double softCompWeight3(params[6] / (params[7] + std::exp(params[8]*pfoEnergy)));
        
        if(softCompWeight2*binDensity >= MAX_EXP_ARG)
            throw std::runtime_error("Beyond numeric limits !");
        
        const double softCompWeightTotal(softCompWeight1*std::exp(softCompWeight2*binDensity) + softCompWeight3);
        const double compensatedHCalEnergy(softCompWeightTotal*binEnergy);
        compensatedEnergy += compensatedHCalEnergy;    
    }
    
    this->CleanCluster(event, compensatedEnergy);
    return compensatedEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SoftwareCompensation::CleanCluster(const Event &event, double &correctedHadronicEnergy)
{
    const int firstPseudoLayer(event.m_firstPseudoLayer);
    const float clusterHadronicEnergy(event.m_reconstructedEnergy);
    
    for (unsigned int hit = 0; hit < event.m_caloHitEnergies.size(); hit++)
    {
        const unsigned int pseudoLayer(event.m_caloHitPseudoLayers.at(hit));
        const int isIsolated(event.m_caloHitIsIsolated.at(hit));
        const int hitType(event.m_caloHitType.at(hit));

        if (hitType != 1 || isIsolated == 1) continue; 

        const float hitHadronicEnergy(event.m_caloHitEnergies.at(hit));

        if ((hitHadronicEnergy > m_minCleanHitEnergy) && (hitHadronicEnergy / clusterHadronicEnergy > m_minCleanHitEnergyFraction))
        {
            float energyInPreviousLayer(0.f);

            if (pseudoLayer > firstPseudoLayer)
            {
                energyInPreviousLayer = this->GetHadronicEnergyInLayer(event, pseudoLayer - 1);
            }

            float energyInNextLayer(0.f);

            if (pseudoLayer < std::numeric_limits<unsigned int>::max())
            {
                energyInNextLayer = this->GetHadronicEnergyInLayer(event, pseudoLayer + 1);
            }

            const float energyInCurrentLayer = this->GetHadronicEnergyInLayer(event, pseudoLayer);
            float energyInAdjacentLayers(energyInPreviousLayer + energyInNextLayer);

            if (pseudoLayer > firstPseudoLayer)
                energyInAdjacentLayers /= 2.f;

            float newHitHadronicEnergy(energyInAdjacentLayers - energyInCurrentLayer + hitHadronicEnergy);
            newHitHadronicEnergy = std::max(newHitHadronicEnergy, m_minCleanCorrectedHitEnergy);

            if (newHitHadronicEnergy < hitHadronicEnergy)
                correctedHadronicEnergy += newHitHadronicEnergy - hitHadronicEnergy;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float SoftwareCompensation::GetHadronicEnergyInLayer(const Event &event, const unsigned int pseudoLayer) const
{
    float hadronicEnergy(0.f);

    for (unsigned int hit = 0; hit < event.m_caloHitEnergies.size(); hit++)
    {
        const unsigned int activePseudoLayer(event.m_caloHitPseudoLayers.at(hit));
        const float hitHadronicEnergy(event.m_caloHitEnergies.at(hit));
        const int isIsolated(event.m_caloHitIsIsolated.at(hit));

        if (pseudoLayer == activePseudoLayer && isIsolated == 0) 
        {
            hadronicEnergy += hitHadronicEnergy;
        }
    }
    return hadronicEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TCanvas *SoftwareCompensation::CreateCanvas(const std::string &title) const
{
  return new TCanvas(title.c_str(), title.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], SoftwareCompensation &softwareCompensation)
{
    int c(0);

    while (((c = getopt(argc, argv, "d:e:f:gt:")) != -1) || (argc == 1))
    {
        switch (c)
        {

        case 'd':
            softwareCompensation.m_outputPath = optarg;
            break;
        case 'e':
            softwareCompensation.SetTrueEnergies(optarg);
            break;
        case 'f':
            softwareCompensation.m_filePattern = optarg;
            break;
        case 'g':
            softwareCompensation.m_minimizeUsingTrueEnergy = true;
            break;
        case 't':
            softwareCompensation.m_treeName = optarg;
            break;
        case 'h':
        default:
            std::cout << std::endl << "PandoraPFA_SoftwareCompensation " << std::endl
                      << "    -d path      (mandatory, output path to send results to)                                                    " << std::endl
                      << "    -e energies  (mandatory, input energies comma separated, i.e '10:20:30:40:50:60')                           " << std::endl
                      << "    -f pattern   (mandatory, input file pattern. Should contains ${energy} key, i.e ./File_${energy}.root)      " << std::endl
                      << "    -g           (if specified, the minimization will use the true (mc or beam) energy instead of reconstructed)" << std::endl
                      << "    -t treename  (mandatory, the tree name to look in root files, usually SoftwareCompensationTrainingTree)     " << std::endl
                      << std::endl;
            return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Tokenize(const std::string &inputString, std::vector<std::string> &tokens, const std::string &delimiter)
{
	std::string::size_type lastPos = inputString.find_first_not_of(delimiter, 0);
	std::string::size_type pos     = inputString.find_first_of(delimiter, lastPos);

	while ((std::string::npos != pos) || (std::string::npos != lastPos))
	{
		tokens.push_back(inputString.substr(lastPos, pos - lastPos));
		lastPos = inputString.find_first_not_of(delimiter, pos);
		pos = inputString.find_first_of(delimiter, lastPos);
	}
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string &Replace(std::string &str, const std::string &variable, const std::string &value)
{
	std::string replaceVar = "%{" + variable + "}";
	size_t pos = str.find(replaceVar);

	while( pos != std::string::npos )
	{
		str.replace( pos , replaceVar.size() , value );
		pos = str.find(replaceVar);
	}

	return str;
}
