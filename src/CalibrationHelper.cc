/**
 *  @file   PandoraAnalysis/src/CalibrationHelper.cc
 * 
 *  @brief  Implementation of the calibration helper class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/SimCalorimeterHit.h"

#include "CalorimeterHitType.h"

#include "marlin/Global.h"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "CalibrationHelper.h"

using dd4hep::DetType;

namespace pandora_analysis
{

/**
* @brief CaloHitException class
*/
class CaloHitException : public std::exception
{
public:
    const char *what() const throw();
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
* @brief DivisionByZeroException class
*/
class DivisionByZeroException : public std::exception
{
public:
    const char *what() const throw();
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
* @brief HCalHitPositionException class
*/
class HCalHitPositionException : public std::exception
{
public:
    const char *what() const throw();
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CalibrationHelper::Settings::Settings()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CalibrationHelper::CalibrationHelper(const Settings &settings) :
    m_settings(settings),
    m_pfoMinHCalLayerToEdge(0),
    m_hCalCaloHitEnergyFromEdge(0.f),
    m_totalCaloHitEnergy(0.f),
    m_eCalTotalCaloHitEnergy(0.f),
    m_hCalTotalCaloHitEnergy(0.f),
    m_muonTotalCaloHitEnergy(0.f),
    m_bCalTotalCaloHitEnergy(0.f),
    m_lHCalTotalCaloHitEnergy(0.f),
    m_lCalTotalCaloHitEnergy(0.f),
    m_totalSimCaloHitEnergy(0.f),
    m_eCalTotalSimCaloHitEnergy(0.f),
    m_hCalTotalSimCaloHitEnergy(0.f),
    m_muonTotalSimCaloHitEnergy(0.f),
    m_bCalTotalSimCaloHitEnergy(0.f),
    m_lHCalTotalSimCaloHitEnergy(0.f),
    m_lCalTotalSimCaloHitEnergy(0.f),
    m_hECalDirectionCorrectedCaloHitEnergy(NULL),
    m_hHCalDirectionCorrectedCaloHitEnergy(NULL),
    m_hMuonDirectionCorrectedCaloHitEnergy(NULL),
    m_hHCalBarrelDirectionCorrectedSimCaloHit(NULL),
    m_hHCalEndCapDirectionCorrectedSimCaloHit(NULL),
    m_hHCalOtherDirectionCorrectedSimCaloHit(NULL),
    m_hECalDirectionCorrectedSimCaloHit(NULL),
    m_hHCalBarrelDirectionCorrectionSimCaloHit(NULL),
    m_hHCalEndCapDirectionCorrectionSimCaloHit(NULL),
    m_hHCalOtherDirectionCorrectionSimCaloHit(NULL),
    m_hECalBarrelDirectionCorrectionSimCaloHit(NULL),
    m_hECalEndCapDirectionCorrectionSimCaloHit(NULL),
    m_hECalOtherDirectionCorrectionSimCaloHit(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CalibrationHelper::~CalibrationHelper()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::SetBranchAddresses(TTree *pTTree)
{
    pTTree->Branch("pfoMinHCalLayerToEdge", &m_pfoMinHCalLayerToEdge, "pfoMinHCalLayerToEdge/I");
    pTTree->Branch("hCalCaloHitEnergyFromEdge", &m_hCalCaloHitEnergyFromEdge, "hCalCaloHitEnergyFromEdge/F");
    pTTree->Branch("TotalCaloHitEnergy", &m_totalCaloHitEnergy, "TotalCaloHitEnergy/F");
    pTTree->Branch("ECalTotalCaloHitEnergy", &m_eCalTotalCaloHitEnergy, "ECalTotalCaloHitEnergy/F");
    pTTree->Branch("HCalTotalCaloHitEnergy", &m_hCalTotalCaloHitEnergy, "HCalTotalCaloHitEnergy/F");
    pTTree->Branch("MuonTotalCaloHitEnergy", &m_muonTotalCaloHitEnergy, "MuonTotalCaloHitEnergy/F");
    pTTree->Branch("BCalTotalCaloHitEnergy", &m_bCalTotalCaloHitEnergy, "BCalTotalCaloHitEnergy/F");
    pTTree->Branch("LHCalTotalCaloHitEnergy", &m_lHCalTotalCaloHitEnergy, "LHCalTotalCaloHitEnergy/F");
    pTTree->Branch("LCalTotalCaloHitEnergy", &m_lCalTotalCaloHitEnergy, "LCalTotalCaloHitEnergy/F");
    pTTree->Branch("TotalSimCaloHitEnergy", &m_totalSimCaloHitEnergy, "TotalSimCaloHitEnergy/F");
    pTTree->Branch("ECalTotalSimCaloHitEnergy", &m_eCalTotalSimCaloHitEnergy, "ECalTotalSimCaloHitEnergy/F");
    pTTree->Branch("HCalTotalSimCaloHitEnergy", &m_hCalTotalSimCaloHitEnergy, "HCalTotalSimCaloHitEnergy/F");
    pTTree->Branch("MuonTotalSimCaloHitEnergy", &m_muonTotalSimCaloHitEnergy, "MuonTotalSimCaloHitEnergy/F");
    pTTree->Branch("BCalTotalSimCaloHitEnergy", &m_bCalTotalSimCaloHitEnergy, "BCalTotalSimCaloHitEnergy/F");
    pTTree->Branch("LHCalTotalSimCaloHitEnergy", &m_lHCalTotalSimCaloHitEnergy, "LHCalTotalSimCaloHitEnergy/F");
    pTTree->Branch("LCalTotalSimCaloHitEnergy", &m_lCalTotalSimCaloHitEnergy, "LCalTotalSimCaloHitEnergy/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::CreateHistograms()
{
    m_hHCalBarrelDirectionCorrectionSimCaloHit = new TH1F("HCalBarrelDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the HCal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalBarrelDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Direction Correction -> (DirCorr * SimCaloHit)");
    m_hHCalBarrelDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hHCalEndCapDirectionCorrectionSimCaloHit = new TH1F("HCalEndCapDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the HCal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalEndCapDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Direction Correction -> (DirCorr * SimCaloHit)");
    m_hHCalEndCapDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hHCalOtherDirectionCorrectionSimCaloHit = new TH1F("HCalOtherDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the HCal Other (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalOtherDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Direction Correction -> (DirCorr * SimCaloHit)");
    m_hHCalOtherDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hECalBarrelDirectionCorrectionSimCaloHit = new TH1F("ECalBarrelDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the ECal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsPhotons && Contained in ECal)", 200, 0., 1.0);
    m_hECalBarrelDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Direction Correction -> (DirCorr * SimCaloHit)");
    m_hECalBarrelDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hECalEndCapDirectionCorrectionSimCaloHit = new TH1F("ECalEndCapDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the ECal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsPhotons && Contained in ECal)", 200, 0., 1.0);
    m_hECalEndCapDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Direction Correction -> (DirCorr * SimCaloHit)");
    m_hECalEndCapDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hECalOtherDirectionCorrectionSimCaloHit = new TH1F("ECalOtherDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the ECal Other (1==nPfoTargetsTotal && 1==nPfoTargetsPhotons && Contained in ECal)", 200, 0., 1.0);
    m_hECalOtherDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Direction Correction -> (DirCorr * SimCaloHit)");
    m_hECalOtherDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hHCalBarrelDirectionCorrectedSimCaloHit = new TH1F("HCalDirectionCorrectedSimCaloHitBarrel", "Distribution of Direction Corrected SimCaloHits in the HCal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hHCalBarrelDirectionCorrectedSimCaloHit->GetXaxis()->SetTitle("Direction Corrected SimCaloHit Measurement");
    m_hHCalBarrelDirectionCorrectedSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hHCalEndCapDirectionCorrectedSimCaloHit = new TH1F("HCalDirectionCorrectedSimCaloHitEndCap", "Distribution of Direction Corrected SimCaloHits in the HCal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hHCalEndCapDirectionCorrectedSimCaloHit->GetXaxis()->SetTitle("Direction Corrected SimCaloHit Measurement");
    m_hHCalEndCapDirectionCorrectedSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hHCalOtherDirectionCorrectedSimCaloHit = new TH1F("HCalDirectionCorrectedSimCaloHitOther", "Distribution of Direction Corrected SimCaloHits in the HCal Other (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.005);
    m_hHCalOtherDirectionCorrectedSimCaloHit->GetXaxis()->SetTitle("Direction Corrected SimCaloHit Measurement");
    m_hHCalOtherDirectionCorrectedSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hECalDirectionCorrectedSimCaloHit = new TH1F("ECalDirectionCorrectedSimCaloHit", "Distribution of Direction Corrected SimCaloHits in the ECal (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hECalDirectionCorrectedSimCaloHit->GetXaxis()->SetTitle("Direction Corrected SimCaloHit Measurement");
    m_hECalDirectionCorrectedSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hECalDirectionCorrectedCaloHitEnergy = new TH1F("ECalDirectionCorrectedCaloHitEnergy", "1==nPfoTargetsTotal && 1==nPfoTargetsTracks", 500, 0., 0.1);
    m_hECalDirectionCorrectedCaloHitEnergy->GetXaxis()->SetTitle("Direction Corrected Calo Hit Energy in ECal");
    m_hECalDirectionCorrectedCaloHitEnergy->GetYaxis()->SetTitle("Entries");

    m_hHCalDirectionCorrectedCaloHitEnergy = new TH1F("HCalDirectionCorrectedCaloHitEnergy", "1==nPfoTargetsTotal && 1==nPfoTargetsTracks", 500, 0., 0.1);
    m_hHCalDirectionCorrectedCaloHitEnergy->GetXaxis()->SetTitle("Direction Corrected Calo Hit Energy in HCal");
    m_hHCalDirectionCorrectedCaloHitEnergy->GetYaxis()->SetTitle("Entries");

    m_hMuonDirectionCorrectedCaloHitEnergy = new TH1F("MuonDirectionCorrectedCaloHitEnergy", "1==nPfoTargetsTotal && 1==nPfoTargetsTracks", 500, 0., 1.0);
    m_hMuonDirectionCorrectedCaloHitEnergy->GetXaxis()->SetTitle("Direction Corrected Calo Hit Energy in Muon Chamber");
    m_hMuonDirectionCorrectedCaloHitEnergy->GetYaxis()->SetTitle("Entries");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::SetHistogramDirectories(TFile *pTFile)
{
    m_hHCalBarrelDirectionCorrectionSimCaloHit->SetDirectory(pTFile);
    m_hHCalEndCapDirectionCorrectionSimCaloHit->SetDirectory(pTFile);
    m_hHCalOtherDirectionCorrectionSimCaloHit->SetDirectory(pTFile);
    m_hECalBarrelDirectionCorrectionSimCaloHit->SetDirectory(pTFile);
    m_hECalEndCapDirectionCorrectionSimCaloHit->SetDirectory(pTFile);
    m_hECalOtherDirectionCorrectionSimCaloHit->SetDirectory(pTFile);
    m_hHCalBarrelDirectionCorrectedSimCaloHit->SetDirectory(pTFile);
    m_hHCalEndCapDirectionCorrectedSimCaloHit->SetDirectory(pTFile);
    m_hHCalOtherDirectionCorrectedSimCaloHit->SetDirectory(pTFile);
    m_hECalDirectionCorrectedSimCaloHit->SetDirectory(pTFile);
    m_hECalDirectionCorrectedCaloHitEnergy->SetDirectory(pTFile);
    m_hHCalDirectionCorrectedCaloHitEnergy->SetDirectory(pTFile);
    m_hMuonDirectionCorrectedCaloHitEnergy->SetDirectory(pTFile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::WriteHistograms()
{
    m_hHCalBarrelDirectionCorrectionSimCaloHit->Write();
    m_hHCalEndCapDirectionCorrectionSimCaloHit->Write();
    m_hHCalOtherDirectionCorrectionSimCaloHit->Write();
    m_hECalBarrelDirectionCorrectionSimCaloHit->Write();
    m_hECalEndCapDirectionCorrectionSimCaloHit->Write();
    m_hECalOtherDirectionCorrectionSimCaloHit->Write();
    m_hHCalBarrelDirectionCorrectedSimCaloHit->Write();
    m_hHCalEndCapDirectionCorrectedSimCaloHit->Write();
    m_hHCalOtherDirectionCorrectedSimCaloHit->Write();
    m_hECalDirectionCorrectedSimCaloHit->Write();
    m_hECalDirectionCorrectedCaloHitEnergy->Write();
    m_hHCalDirectionCorrectedCaloHitEnergy->Write();
    m_hMuonDirectionCorrectedCaloHitEnergy->Write();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::Clear()
{
    m_pfoMinHCalLayerToEdge = 0;
    m_hCalCaloHitEnergyFromEdge = 0.f;
    m_totalCaloHitEnergy = 0.f;
    m_eCalTotalCaloHitEnergy = 0.f;
    m_hCalTotalCaloHitEnergy = 0.f;
    m_muonTotalCaloHitEnergy = 0.f;
    m_bCalTotalCaloHitEnergy = 0.f;
    m_lHCalTotalCaloHitEnergy = 0.f;
    m_lCalTotalCaloHitEnergy = 0.f;
    m_totalSimCaloHitEnergy = 0.f;
    m_eCalTotalSimCaloHitEnergy = 0.f;
    m_hCalTotalSimCaloHitEnergy = 0.f;
    m_muonTotalSimCaloHitEnergy = 0.f;
    m_bCalTotalSimCaloHitEnergy = 0.f;
    m_lHCalTotalSimCaloHitEnergy = 0.f;
    m_lCalTotalSimCaloHitEnergy = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::Calibrate(const EVENT::LCEvent *pLCEvent, const ParticleVector &particleVector,
    const int nPfoTargetsTotal, const int nPfoTargetsTracks, const int nPfoTargetsNeutralHadrons, const int nPfoTargetsPhotons,
    const float pfoTargetsEnergyTotal)
{
    m_pfoMinHCalLayerToEdge = this->GetMinNHCalLayersFromEdge(particleVector);
    m_hCalCaloHitEnergyFromEdge = this->GetHCalCaloHitEnergyFromEdge(pLCEvent, m_settings.m_hCalCollections);

    this->ReadSimCaloHitEnergies(pLCEvent, m_settings.m_eCalCollectionsSimCaloHit, m_eCalTotalSimCaloHitEnergy);
    this->ReadSimCaloHitEnergies(pLCEvent, m_settings.m_hCalBarrelCollectionsSimCaloHit, m_hCalTotalSimCaloHitEnergy);
    this->ReadSimCaloHitEnergies(pLCEvent, m_settings.m_hCalEndCapCollectionsSimCaloHit, m_hCalTotalSimCaloHitEnergy);
    this->ReadSimCaloHitEnergies(pLCEvent, m_settings.m_hCalOtherCollectionsSimCaloHit,  m_hCalTotalSimCaloHitEnergy);
    this->ReadSimCaloHitEnergies(pLCEvent, m_settings.m_muonCollectionsSimCaloHit, m_muonTotalSimCaloHitEnergy);
    this->ReadSimCaloHitEnergies(pLCEvent, m_settings.m_bCalCollectionsSimCaloHit, m_bCalTotalSimCaloHitEnergy);
    this->ReadSimCaloHitEnergies(pLCEvent, m_settings.m_lHCalCollectionsSimCaloHit, m_lHCalTotalSimCaloHitEnergy);
    this->ReadSimCaloHitEnergies(pLCEvent, m_settings.m_lCalCollectionsSimCaloHit, m_lCalTotalSimCaloHitEnergy);

    m_totalSimCaloHitEnergy = m_eCalTotalSimCaloHitEnergy + m_hCalTotalSimCaloHitEnergy + m_muonTotalSimCaloHitEnergy + m_bCalTotalSimCaloHitEnergy + m_lHCalTotalSimCaloHitEnergy + m_lCalTotalSimCaloHitEnergy;

    this->ReadCaloHitEnergies(pLCEvent, m_settings.m_eCalCollections, m_eCalTotalCaloHitEnergy);
    this->ReadCaloHitEnergies(pLCEvent, m_settings.m_hCalCollections, m_hCalTotalCaloHitEnergy);
    this->ReadCaloHitEnergies(pLCEvent, m_settings.m_muonCollections, m_muonTotalCaloHitEnergy);
    this->ReadCaloHitEnergies(pLCEvent, m_settings.m_bCalCollections, m_bCalTotalCaloHitEnergy);
    this->ReadCaloHitEnergies(pLCEvent, m_settings.m_lHCalCollections, m_lHCalTotalCaloHitEnergy);
    this->ReadCaloHitEnergies(pLCEvent, m_settings.m_lCalCollections, m_lCalTotalCaloHitEnergy);

    m_totalCaloHitEnergy = m_eCalTotalCaloHitEnergy + m_hCalTotalCaloHitEnergy + m_muonTotalCaloHitEnergy + m_bCalTotalCaloHitEnergy +
        m_lHCalTotalCaloHitEnergy + m_lCalTotalCaloHitEnergy;

    if (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)
    {
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_hCalBarrelCollectionsSimCaloHit, 1, m_hHCalBarrelDirectionCorrectedSimCaloHit);
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_hCalEndCapCollectionsSimCaloHit, 1, m_hHCalEndCapDirectionCorrectedSimCaloHit);
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_hCalOtherCollectionsSimCaloHit, 1, m_hHCalOtherDirectionCorrectedSimCaloHit);
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_eCalCollectionsSimCaloHit, 1, m_hECalDirectionCorrectedSimCaloHit);

        this->AddDirectionCorrectedCaloHitEntries(pLCEvent, m_settings.m_eCalCollections, m_hECalDirectionCorrectedCaloHitEnergy);
        this->AddDirectionCorrectedCaloHitEntries(pLCEvent, m_settings.m_hCalCollections, m_hHCalDirectionCorrectedCaloHitEnergy);
        this->AddDirectionCorrectedCaloHitEntries(pLCEvent, m_settings.m_muonCollections, m_hMuonDirectionCorrectedCaloHitEnergy);
    }

    if (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && ((m_totalCaloHitEnergy-m_hCalTotalCaloHitEnergy) < (0.01*pfoTargetsEnergyTotal)) && m_pfoMinHCalLayerToEdge > 5 )
    {
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_hCalBarrelCollectionsSimCaloHit, 0, m_hHCalBarrelDirectionCorrectionSimCaloHit);
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_hCalEndCapCollectionsSimCaloHit, 0, m_hHCalEndCapDirectionCorrectionSimCaloHit);
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_hCalOtherCollectionsSimCaloHit, 0, m_hHCalOtherDirectionCorrectionSimCaloHit);
    }

    if(1==nPfoTargetsTotal && 1==nPfoTargetsPhotons && ((m_totalCaloHitEnergy-m_eCalTotalCaloHitEnergy) < (0.01*pfoTargetsEnergyTotal)))
    {
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_eCalBarrelCollectionsSimCaloHit, 0, m_hECalBarrelDirectionCorrectionSimCaloHit);
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_eCalEndCapCollectionsSimCaloHit, 0, m_hECalEndCapDirectionCorrectionSimCaloHit);
        this->AddSimCaloHitEntries(pLCEvent, m_settings.m_eCalOtherCollectionsSimCaloHit, 0, m_hECalOtherDirectionCorrectionSimCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int CalibrationHelper::GetMinNHCalLayersFromEdge(const ParticleVector &particleVector) const
{
    int minNHCalLayersFromEdge = std::numeric_limits<int>::max();

    for (ParticleVector::const_iterator pIter = particleVector.begin(), pIterEnd = particleVector.end(); pIter != pIterEnd; ++pIter)
    {
        const EVENT::ReconstructedParticle *pPfo = *pIter;
        const ClusterVec &clusterVec = pPfo->getClusters();

        for (ClusterVec::const_iterator cIter = clusterVec.begin(), cIterEnd = clusterVec.end(); cIter != cIterEnd; ++cIter)
        {
            const CalorimeterHitVec &calorimeterHitVec = (*cIter)->getCalorimeterHits();

            for (CalorimeterHitVec::const_iterator hIter = calorimeterHitVec.begin(), hIterEnd = calorimeterHitVec.end(); hIter != hIterEnd; ++hIter)
            {
                const EVENT::CalorimeterHit *pCalorimeterHit = *hIter;
                const CHT cht(pCalorimeterHit->getType());

                if (cht.is(CHT::hcal))
                {
                    try
                    {
                        const int trialMinHCALLayerToEdge = this->GetNHCalLayersFromEdge(pCalorimeterHit);

                        if (minNHCalLayersFromEdge > trialMinHCALLayerToEdge)
                        {
                            minNHCalLayersFromEdge = trialMinHCALLayerToEdge;
                        }
                    }
                    catch (std::out_of_range &e)
                    {
                        std::cout << "CalibrationHelper::GetMinNHCalLayersFromEdge: out of range error: " << e.what() <<std::endl;
                    }
                    catch (...)
                    {
                        std::cout << "CalibrationHelper::GetMinNHCalLayersFromEdge: Unknown exception." <<std::endl;
                    }
                }
            }
        }
    }
    return minNHCalLayersFromEdge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CalibrationHelper::GetHCalCaloHitEnergyFromEdge(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames) const
{
    float hCalCaloHitEnergyFromEdge(0.f);
    
    for (EVENT::LCStrVec::const_iterator iter = collectionNames.begin(), iterEnd = collectionNames.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            if (NULL != pLCCollection)
            {
                const int nElements(pLCCollection->getNumberOfElements());

                for (int iHit = 0; iHit < nElements; ++iHit)
                {
                    const CalorimeterHit *pCalorimeterHit = dynamic_cast<const CalorimeterHit*>(pLCCollection->getElementAt(iHit));
                    
                    if (NULL == pCalorimeterHit)
                    {
                        streamlog_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type CalorimeterHit " << std::endl;
                        throw;
                    }
                    
                    const CHT cht(pCalorimeterHit->getType());

                    if (cht.is(CHT::hcal))
                    {
                        try
                        {
                            const int trialMinHCALLayerToEdge = this->GetNHCalLayersFromEdge(pCalorimeterHit);
                            
                            if(trialMinHCALLayerToEdge < 5)
                            {
                                const float hitEnergy(pCalorimeterHit->getEnergy());
                                hCalCaloHitEnergyFromEdge += hitEnergy;                              
                            }
                        }
                        catch (std::out_of_range &e)
                        {
                            std::cout << "CalibrationHelper::GetHCalCaloHitEnergyFromEdge: out of range error: " << e.what() <<std::endl;
                        }
                        catch (...)
                        {
                            std::cout << "CalibrationHelper::GetHCalCaloHitEnergyFromEdge: Unknown exception." <<std::endl;
                        }
                    }
                }
            }
        }
        catch (DataNotAvailableException &)
        {
            streamlog_out(DEBUG) << "No Collection : " << (*iter) << std::endl;
        }
    }
    
    return hCalCaloHitEnergyFromEdge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int CalibrationHelper::GetNHCalLayersFromEdge(const EVENT::CalorimeterHit *const pCaloHit) const
{
    /* TODO maybe check hit is HCAL */

    CHT cht(pCaloHit->getType());

  
    if (cht.is(CHT::hcal) != true)
    {
        throw CaloHitException();
    }

    // Geometry information from dd4hep
    // Get Ring HCal Extension
    const dd4hep::rec::LayeredCalorimeterData *pHCalRingExtension(this->GetExtension( ( DetType::CALORIMETER | DetType::HADRONIC | DetType::ENDCAP | DetType::AUXILIARY ),( DetType::FORWARD )));
    const unsigned int hCalRingOuterSymmetryOrder(pHCalRingExtension->outer_symmetry);
    const float hCalRingOuterPhi0(pHCalRingExtension->outer_phi0);
    const float hCalRingOuterR(pHCalRingExtension->extent[1]/dd4hep::mm);
    const float hCalRingInnerZ(pHCalRingExtension->extent[2]/dd4hep::mm);
    const float hCalRingOuterZ(pHCalRingExtension->extent[3]/dd4hep::mm);
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &hCalRingLayers(this->GetExtension(( DetType::CALORIMETER | DetType::HADRONIC | DetType::ENDCAP | DetType::AUXILIARY ), ( DetType::FORWARD ))->layers);
    const float hCalRingLayerThickness((hCalRingLayers.back().inner_thickness+hCalRingLayers.back().outer_thickness)/dd4hep::mm);

    // Get HCal Barrel Extension
    const dd4hep::rec::LayeredCalorimeterData *pHCalBarrelExtension(this->GetExtension( ( DetType::CALORIMETER | DetType::HADRONIC | DetType::BARREL), (DetType::AUXILIARY |  DetType::FORWARD ) ));
    const unsigned int hCalBarrelOuterSymmetry(pHCalBarrelExtension->outer_symmetry);
    const float hCalBarrelOuterPhi0(pHCalBarrelExtension->outer_phi0);
    const float hCalBarrelOuterR(pHCalBarrelExtension->extent[1]/dd4hep::mm);
    const float hCalBarrelOuterZ(pHCalBarrelExtension->extent[3]/dd4hep::mm);
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &hCalBarrelLayers(this->GetExtension(( DetType::CALORIMETER | DetType::HADRONIC | DetType::BARREL), ( DetType::AUXILIARY | DetType::FORWARD ))->layers);
    const float hCalBarrelLayerThickness((hCalBarrelLayers.back().inner_thickness+hCalBarrelLayers.back().outer_thickness)/dd4hep::mm);

    // Get Ring HCal Endcap Extension
    const dd4hep::rec::LayeredCalorimeterData *pHCalEndcapExtension(this->GetExtension(( DetType::CALORIMETER | DetType::HADRONIC | DetType::ENDCAP),( DetType::AUXILIARY | DetType::FORWARD )));
    const unsigned int hCalEndCapOuterSymmetryOrder(pHCalEndcapExtension->outer_symmetry);
    const float hCalEndCapOuterPhiCoordinate(pHCalEndcapExtension->outer_phi0);
    const float hCalEndCapOuterR(pHCalEndcapExtension->extent[1]/dd4hep::mm);
    const float hCalEndCapInnerZ(pHCalEndcapExtension->extent[2]/dd4hep::mm);
    const float hCalEndCapOuterZ(pHCalEndcapExtension->extent[3]/dd4hep::mm);
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer> &hCalEndCapLayers(this->GetExtension(( DetType::CALORIMETER | DetType::HADRONIC | DetType::ENDCAP), ( DetType::AUXILIARY | DetType::FORWARD ))->layers);
    const float hCalEndCapLayerThickness((hCalEndCapLayers.back().inner_thickness+hCalEndCapLayers.back().outer_thickness)/dd4hep::mm);

    /* TODO maybe sanity-check any variables here*/

    if (std::fabs(hCalRingLayerThickness) < std::numeric_limits<float>::epsilon() ||
        std::fabs(hCalBarrelLayerThickness) < std::numeric_limits<float>::epsilon() ||
        std::fabs(hCalEndCapLayerThickness) < std::numeric_limits<float>::epsilon())
    {
        throw DivisionByZeroException();
    }

    // Calo hit coordinate calculations
    const float hCalBarrelMaximumRadius(this->GetMaximumRadius(pCaloHit, hCalBarrelOuterSymmetry, hCalBarrelOuterPhi0));
    const float hCalRingMaximumRadius(this->GetMaximumRadius(pCaloHit, hCalRingOuterSymmetryOrder, hCalRingOuterPhi0));
    const float hCalEndcapMaximumRadius(this->GetMaximumRadius(pCaloHit, hCalEndCapOuterSymmetryOrder, hCalEndCapOuterPhiCoordinate));
    const float caloHitAbsZ(std::fabs(pCaloHit->getPosition()[2]));

    float radialDistanceToEdge(std::numeric_limits<float>::max());
    float rearDistanceToEdge(std::numeric_limits<float>::max());

    if (caloHitAbsZ < hCalBarrelOuterZ)
    {
        radialDistanceToEdge = (hCalBarrelOuterR - hCalBarrelMaximumRadius) / hCalBarrelLayerThickness;

        rearDistanceToEdge = (hCalBarrelOuterZ-caloHitAbsZ) / hCalBarrelLayerThickness
                           + (hCalRingOuterZ-hCalRingInnerZ) / hCalRingLayerThickness
                           + (hCalEndCapOuterZ-hCalEndCapInnerZ) / hCalEndCapLayerThickness;
    }
    else if ( hCalRingInnerZ < caloHitAbsZ && caloHitAbsZ < hCalRingOuterZ )
    {
        radialDistanceToEdge = ( hCalRingOuterR - hCalRingMaximumRadius ) / hCalRingLayerThickness;

        rearDistanceToEdge = (hCalRingOuterZ-caloHitAbsZ) / hCalRingLayerThickness
                           + (hCalEndCapOuterZ-hCalEndCapInnerZ) / hCalEndCapLayerThickness;
    }
    else if (  hCalEndCapInnerZ < caloHitAbsZ && caloHitAbsZ < hCalEndCapOuterZ )
    {
        radialDistanceToEdge = (hCalEndCapOuterR - hCalEndcapMaximumRadius) / hCalEndCapLayerThickness;

        rearDistanceToEdge = (hCalEndCapOuterZ-caloHitAbsZ) / hCalEndCapLayerThickness;
    }
    else
    {
        /* Could throw an exception here*/
        throw HCalHitPositionException();
    }

    return static_cast<int>(std::min(radialDistanceToEdge, rearDistanceToEdge));

}

//------------------------------------------------------------------------------------------------------------------------------------------

float CalibrationHelper::GetMaximumRadius(const EVENT::CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const
{
    const float *pCaloHitPosition(pCaloHit->getPosition());

    if (symmetryOrder <= 2)
        return std::sqrt((pCaloHitPosition[0] * pCaloHitPosition[0]) + (pCaloHitPosition[1] * pCaloHitPosition[1]));

    float maximumRadius(0.f);
    static const float twoPi = static_cast<float>(2. * std::acos(-1.));

    for (unsigned int i = 0; i < symmetryOrder; ++i)
    {
        const float phi = phi0 + i * twoPi / static_cast<float>(symmetryOrder);
        float radius = pCaloHitPosition[0] * std::cos(phi) + pCaloHitPosition[1] * std::sin(phi);

        if (radius > maximumRadius)
            maximumRadius = radius;
    }

    return maximumRadius;
}

//------------------------------------------------------------------------------------------------------------------------------------------

dd4hep::rec::LayeredCalorimeterData *CalibrationHelper::GetExtension(unsigned int includeFlag, unsigned int excludeFlag) const
{
    dd4hep::rec::LayeredCalorimeterData *pExtension(NULL);
    dd4hep::Detector &theDetector = dd4hep::Detector::getInstance();
    const std::vector< dd4hep::DetElement> &theDetectors(dd4hep::DetectorSelector(theDetector).detectors(  includeFlag, excludeFlag ));
    streamlog_out( DEBUG2 ) << " GetExtension :  includeFlag: " << DetType( includeFlag ) << " excludeFlag: " << DetType( excludeFlag )
                            << "  found : " << theDetectors.size() << "  - first det: " << theDetectors.at(0).name() << std::endl ;

    if(theDetectors.size()  != 1)
    {
        std::stringstream es;
        es << " GetExtension: selection is not unique (or empty)  includeFlag: " << DetType( includeFlag ) << " excludeFlag: " << DetType( excludeFlag ) << " --- found detectors : " ;

        for(unsigned int i = 0, N = theDetectors.size(); i<N; ++i)
        {
            es << theDetectors.at(i).name() << ", " ;
        }
        throw std::runtime_error(es.str());
    }
    pExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();
    return pExtension;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::ReadCaloHitEnergies(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, float &hitEnergySum) const
{
    for (EVENT::LCStrVec::const_iterator iter = collectionNames.begin(), iterEnd = collectionNames.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            if (NULL != pLCCollection)
            {
                const int nElements(pLCCollection->getNumberOfElements());

                for (int iHit = 0; iHit < nElements; ++iHit)
                {
                    const CalorimeterHit *pCalorimeterHit = dynamic_cast<const CalorimeterHit*>(pLCCollection->getElementAt(iHit));

                    if (NULL == pCalorimeterHit)
                    {
                        streamlog_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type CalorimeterHit " << std::endl;
                        throw;
                    }

                    const float hitEnergy(pCalorimeterHit->getEnergy());
                    hitEnergySum += hitEnergy;
                }
            }
        }
        catch (DataNotAvailableException &)
        {
            streamlog_out(DEBUG) << "No Collection : " << (*iter) << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::ReadSimCaloHitEnergies(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, float &hitEnergySum) const
{
    for (EVENT::LCStrVec::const_iterator iter = collectionNames.begin(), iterEnd = collectionNames.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            if (NULL != pLCCollection)
            {
                const int nElements(pLCCollection->getNumberOfElements());

                for (int iHit = 0; iHit < nElements; ++iHit)
                {
                    const SimCalorimeterHit *pSimCalorimeterHit = dynamic_cast<const SimCalorimeterHit*>( pLCCollection->getElementAt(iHit));

                    if (NULL == pSimCalorimeterHit)
                    {
                        streamlog_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type CalorimeterHit " << std::endl;
                        throw;
                    }

                    const float hitEnergy(pSimCalorimeterHit->getEnergy());
                    hitEnergySum += hitEnergy;
                }
            }
        }
        catch (DataNotAvailableException &)
        {
            streamlog_out(DEBUG) << "No Collection : " << (*iter) << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::AddSimCaloHitEntries(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, const unsigned int SimCaloHit_DC, TH1F *pTH1F) const
{
    const dd4hep::rec::LayeredCalorimeterData *pECalEndcapExtension(this->GetExtension((DetType::CALORIMETER | DetType::ELECTROMAGNETIC | DetType::ENDCAP), ( DetType::AUXILIARY | DetType::FORWARD)));
    const float zOfEndCap(static_cast<float>(pECalEndcapExtension->extent[2]/dd4hep::mm));

    for (EVENT::LCStrVec::const_iterator iter = collectionNames.begin(), iterEnd = collectionNames.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            if (NULL != pLCCollection)
            {
                const int nElements = pLCCollection->getNumberOfElements();

                for (int iHit = 0; iHit < nElements; ++iHit)
                {
                    const SimCalorimeterHit *pSimCalorimeterHit = dynamic_cast<const SimCalorimeterHit*>( pLCCollection->getElementAt(iHit));

                    if (NULL == pSimCalorimeterHit)
                    {
                        streamlog_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type SimCalorimeterHit " << std::endl;
                        throw;
                    }

                    const float SimCaloHitEnergy = pSimCalorimeterHit->getEnergy();
                    const float x(pSimCalorimeterHit->getPosition()[0]);
                    const float y(pSimCalorimeterHit->getPosition()[1]);
                    const float z(pSimCalorimeterHit->getPosition()[2]);
                    const float r(std::sqrt(x * x + y * y + z * z));
                    const float correction((std::fabs(z) < zOfEndCap) ? r / std::sqrt(x * x + y * y) : r / std::fabs(z));

                    if (correction < std::numeric_limits<float>::epsilon())
                    {
                        /*/ TODO*/
                        throw DivisionByZeroException();
                    }

                    if(SimCaloHit_DC)
                        pTH1F->Fill(SimCaloHitEnergy / correction);

                    else
                        pTH1F->Fill(1.f / correction);

                }
            }
        }
        catch (DataNotAvailableException &)
        {
            streamlog_out(DEBUG) << "No Collection : " << *iter << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::AddDirectionCorrectedCaloHitEntries(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, TH1F *pTH1F) const
{
    const dd4hep::rec::LayeredCalorimeterData *pECalEndcapExtension(this->GetExtension((DetType::CALORIMETER | DetType::ELECTROMAGNETIC | DetType::ENDCAP), ( DetType::AUXILIARY | DetType::FORWARD)));
    const float zOfEndCap(static_cast<float>(pECalEndcapExtension->extent[2]/dd4hep::mm));
   
    for (EVENT::LCStrVec::const_iterator iter = collectionNames.begin(), iterEnd = collectionNames.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(iter->c_str());

            if (NULL != pLCCollection)
            {
                const int nElements(pLCCollection->getNumberOfElements());

                for (int iHit = 0; iHit < nElements; ++iHit)
                {
                    const CalorimeterHit *pCalorimeterHit = dynamic_cast<const CalorimeterHit*>(pLCCollection->getElementAt(iHit));

                    if (NULL == pCalorimeterHit)
                    {
                        streamlog_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type CalorimeterHit " << std::endl;
                        throw;
                    }

                    const float hitEnergy(pCalorimeterHit->getEnergy());
                    const float x(pCalorimeterHit->getPosition()[0]);
                    const float y(pCalorimeterHit->getPosition()[1]);
                    const float z(pCalorimeterHit->getPosition()[2]);
                    const float r(std::sqrt(x * x + y * y + z * z));
                    const float correction((std::fabs(z) < zOfEndCap) ? r / std::sqrt(x * x + y * y) : r / std::fabs(z));

                    if (correction < std::numeric_limits<float>::epsilon())
                    {
                        /*/ TODO*/
                        throw DivisionByZeroException();
                    }

                    const float Direction_Corrected_hitEnergy = hitEnergy / correction;
                    pTH1F->Fill(Direction_Corrected_hitEnergy);
                }
            }
        }
        catch (DataNotAvailableException &)
        {
            streamlog_out(DEBUG) << "No Collection : " << (*iter) << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const char *CaloHitException::what() const throw()
{
    return "CalibrationHelper::GetNHCalLayersFromEdge - Calorimeter hit type is not in HCal. ";
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const char *DivisionByZeroException::what() const throw()
{
    return "CalibrationHelper::GetNHCalLayersFromEdge,AddSimCaloHitEntries,AddDirectionCorrectedCaloHitEntries - Possibility of division by zero. ";
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const char *HCalHitPositionException::what() const throw()
{
    return "CalibrationHelper::GetNHCalLayersFromEdge - HCal hit does not appear to be in the HCal. ";
}

} // namespace pandora_analysis
