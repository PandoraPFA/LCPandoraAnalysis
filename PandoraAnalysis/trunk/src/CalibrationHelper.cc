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

#include "gear/LayerLayout.h"
#include "gear/CalorimeterParameters.h"

#include "marlin/Global.h"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "CalibrationHelper.h"

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

CalibrationHelper::Settings::Settings() :
    m_hCalRingOuterSymmetryOrder(8),
    m_hCalRingOuterPhi0(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CalibrationHelper::CalibrationHelper(const Settings &settings) :
    m_settings(settings),
    m_pfoMinHCalLayerToEdge(0),
    m_totalCaloHitEnergy(0.f),
    m_eCalTotalCaloHitEnergy(0.f),
    m_hCalTotalCaloHitEnergy(0.f),
    m_muonTotalCaloHitEnergy(0.f),
    m_bCalTotalCaloHitEnergy(0.f),
    m_lHCalTotalCaloHitEnergy(0.f),
    m_lCalTotalCaloHitEnergy(0.f),
    m_hECalDirectionCorrectedCaloHitEnergy(NULL),
    m_hHCalDirectionCorrectedCaloHitEnergy(NULL),
    m_hMuonDirectionCorrectedCaloHitEnergy(NULL),
    m_hHCalBarrelDirectionCorrectedADC(NULL),
    m_hHCalEndCapDirectionCorrectedADC(NULL),
    m_hHCalOtherDirectionCorrectedADC(NULL),
    m_hECalDirectionCorrectedADC(NULL),
    m_hHCalBarrelDirectionCorrectionADC(NULL),
    m_hHCalEndCapDirectionCorrectionADC(NULL),
    m_hHCalOtherDirectionCorrectionADC(NULL)
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
    pTTree->Branch("TotalCaloHitEnergy", &m_totalCaloHitEnergy, "TotalCaloHitEnergy/F");
    pTTree->Branch("ECalTotalCaloHitEnergy", &m_eCalTotalCaloHitEnergy, "ECalTotalCaloHitEnergy/F");
    pTTree->Branch("HCalTotalCaloHitEnergy", &m_hCalTotalCaloHitEnergy, "HCalTotalCaloHitEnergy/F");
    pTTree->Branch("MuonTotalCaloHitEnergy", &m_muonTotalCaloHitEnergy, "MuonTotalCaloHitEnergy/F");
    pTTree->Branch("BCalTotalCaloHitEnergy", &m_bCalTotalCaloHitEnergy, "BCalTotalCaloHitEnergy/F");
    pTTree->Branch("LHCalTotalCaloHitEnergy", &m_lHCalTotalCaloHitEnergy, "LHCalTotalCaloHitEnergy/F");
    pTTree->Branch("LCalTotalCaloHitEnergy", &m_lCalTotalCaloHitEnergy, "LCalTotalCaloHitEnergy/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::CreateHistograms()
{
    m_hHCalBarrelDirectionCorrectionADC = new TH1F("HCalBarrelDirectionCorrectionADC", "Distribution of Direction Corrections for ADCs in the HCal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalBarrelDirectionCorrectionADC->GetXaxis()->SetTitle("Direction Correction -> (DirCorr * ADC)");
    m_hHCalBarrelDirectionCorrectionADC->GetYaxis()->SetTitle("Entries");

    m_hHCalEndCapDirectionCorrectionADC = new TH1F("HCalEndCapDirectionCorrectionADC", "Distribution of Direction Corrections for ADCs in the HCal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalEndCapDirectionCorrectionADC->GetXaxis()->SetTitle("Direction Correction -> (DirCorr * ADC)");
    m_hHCalEndCapDirectionCorrectionADC->GetYaxis()->SetTitle("Entries");

    m_hHCalOtherDirectionCorrectionADC = new TH1F("HCalOtherDirectionCorrectionADC", "Distribution of Direction Corrections for ADCs in the HCal Other (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalOtherDirectionCorrectionADC->GetXaxis()->SetTitle("Direction Correction -> (DirCorr * ADC)");
    m_hHCalOtherDirectionCorrectionADC->GetYaxis()->SetTitle("Entries");

    m_hHCalBarrelDirectionCorrectedADC = new TH1F("HCalDirectionCorrectedADCBarrel", "Distribution of Direction Corrected ADCs in the HCal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hHCalBarrelDirectionCorrectedADC->GetXaxis()->SetTitle("Direction Corrected ADC Measurement");
    m_hHCalBarrelDirectionCorrectedADC->GetYaxis()->SetTitle("Entries");

    m_hHCalEndCapDirectionCorrectedADC = new TH1F("HCalDirectionCorrectedADCEndCap", "Distribution of Direction Corrected ADCs in the HCal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hHCalEndCapDirectionCorrectedADC->GetXaxis()->SetTitle("Direction Corrected ADC Measurement");
    m_hHCalEndCapDirectionCorrectedADC->GetYaxis()->SetTitle("Entries");

    m_hHCalOtherDirectionCorrectedADC = new TH1F("HCalDirectionCorrectedADCOther", "Distribution of Direction Corrected ADCs in the HCal Other (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.005);
    m_hHCalOtherDirectionCorrectedADC->GetXaxis()->SetTitle("Direction Corrected ADC Measurement");
    m_hHCalOtherDirectionCorrectedADC->GetYaxis()->SetTitle("Entries");

    m_hECalDirectionCorrectedADC = new TH1F("ECalDirectionCorrectedADC", "Distribution of Direction Corrected ADCs in the ECal (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hECalDirectionCorrectedADC->GetXaxis()->SetTitle("Direction Corrected ADC Measurement");
    m_hECalDirectionCorrectedADC->GetYaxis()->SetTitle("Entries");

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
    m_hHCalBarrelDirectionCorrectionADC->SetDirectory(pTFile);
    m_hHCalEndCapDirectionCorrectionADC->SetDirectory(pTFile);
    m_hHCalOtherDirectionCorrectionADC->SetDirectory(pTFile);
    m_hHCalBarrelDirectionCorrectedADC->SetDirectory(pTFile);
    m_hHCalEndCapDirectionCorrectedADC->SetDirectory(pTFile);
    m_hHCalOtherDirectionCorrectedADC->SetDirectory(pTFile);
    m_hECalDirectionCorrectedADC->SetDirectory(pTFile);
    m_hECalDirectionCorrectedCaloHitEnergy->SetDirectory(pTFile);
    m_hHCalDirectionCorrectedCaloHitEnergy->SetDirectory(pTFile);
    m_hMuonDirectionCorrectedCaloHitEnergy->SetDirectory(pTFile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::WriteHistograms()
{
    m_hHCalBarrelDirectionCorrectionADC->Write();
    m_hHCalEndCapDirectionCorrectionADC->Write();
    m_hHCalOtherDirectionCorrectionADC->Write();
    m_hHCalBarrelDirectionCorrectedADC->Write();
    m_hHCalEndCapDirectionCorrectedADC->Write();
    m_hHCalOtherDirectionCorrectedADC->Write();
    m_hECalDirectionCorrectedADC->Write();
    m_hECalDirectionCorrectedCaloHitEnergy->Write();
    m_hHCalDirectionCorrectedCaloHitEnergy->Write();
    m_hMuonDirectionCorrectedCaloHitEnergy->Write();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::Clear()
{
    m_pfoMinHCalLayerToEdge = 0;
    m_totalCaloHitEnergy = 0.f;
    m_eCalTotalCaloHitEnergy = 0.f;
    m_hCalTotalCaloHitEnergy = 0.f;
    m_muonTotalCaloHitEnergy = 0.f;
    m_bCalTotalCaloHitEnergy = 0.f;
    m_lHCalTotalCaloHitEnergy = 0.f;
    m_lCalTotalCaloHitEnergy = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::Calibrate(const EVENT::LCEvent *pLCEvent, const ParticleVector &particleVector, 
    const int nPfoTargetsTotal, const int nPfoTargetsTracks, const int nPfoTargetsNeutralHadrons, const float pfoTargetsEnergyTotal)
{
    m_pfoMinHCalLayerToEdge = this->GetMinNHCalLayersFromEdge(particleVector, m_settings.m_hCalRingOuterSymmetryOrder, m_settings.m_hCalRingOuterPhi0);

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
        this->AddADCEntries(pLCEvent, m_settings.m_hCalBarrelCollectionsADC, 1, m_hHCalBarrelDirectionCorrectedADC);
        this->AddADCEntries(pLCEvent, m_settings.m_hCalEndCapCollectionsADC, 1, m_hHCalEndCapDirectionCorrectedADC);
        this->AddADCEntries(pLCEvent, m_settings.m_hCalOtherCollectionsADC, 1, m_hHCalOtherDirectionCorrectedADC);
        this->AddADCEntries(pLCEvent, m_settings.m_eCalCollectionsADC, 1, m_hECalDirectionCorrectedADC);

        this->AddDirectionCorrectedCaloHitEntries(pLCEvent, m_settings.m_eCalCollections, m_hECalDirectionCorrectedCaloHitEnergy);
        this->AddDirectionCorrectedCaloHitEntries(pLCEvent, m_settings.m_hCalCollections, m_hHCalDirectionCorrectedCaloHitEnergy);
        this->AddDirectionCorrectedCaloHitEntries(pLCEvent, m_settings.m_muonCollections, m_hMuonDirectionCorrectedCaloHitEnergy);
    }

    if (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && ((m_totalCaloHitEnergy-m_hCalTotalCaloHitEnergy) < (0.01*pfoTargetsEnergyTotal)) && m_pfoMinHCalLayerToEdge > 5 )
    {
        this->AddADCEntries(pLCEvent, m_settings.m_hCalBarrelCollectionsADC, 0, m_hHCalBarrelDirectionCorrectionADC);
        this->AddADCEntries(pLCEvent, m_settings.m_hCalEndCapCollectionsADC, 0, m_hHCalEndCapDirectionCorrectionADC);
        this->AddADCEntries(pLCEvent, m_settings.m_hCalOtherCollectionsADC, 0, m_hHCalOtherDirectionCorrectionADC);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int CalibrationHelper::GetMinNHCalLayersFromEdge(const ParticleVector &particleVector, const int hCalRingOuterSymmetryOrder, const float hCalRingOuterPhi0) const
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
                        const int trialMinHCALLayerToEdge = this->GetNHCalLayersFromEdge(pCalorimeterHit,hCalRingOuterSymmetryOrder,hCalRingOuterPhi0);

                        if (minNHCalLayersFromEdge > trialMinHCALLayerToEdge)
                            minNHCalLayersFromEdge = trialMinHCALLayerToEdge;
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

int CalibrationHelper::GetNHCalLayersFromEdge(const EVENT::CalorimeterHit *const pCaloHit, const int hCalRingOuterSymmetryOrder, const float hCalRingOuterPhi0) const
{
    /* TODO maybe check hit is HCAL */

    CHT cht(pCaloHit->getType());

    if (cht.is(CHT::hcal) != true)
    {
        throw CaloHitException();
    }

    /* TODO FIX WHERE THESE PARAMETERS ARE SET => m_hCalRingOuterSymmetryOrder and m_hCalRingOuterPhi0*/
    // Geometry information from gear - additional parameters read-in by calling processor
    const float hCalRingOuterR(marlin::Global::GEAR->getHcalRingParameters().getExtent()[1]);
    const float hCalRingInnerZ(marlin::Global::GEAR->getHcalRingParameters().getExtent()[2]);
    const float hCalRingOuterZ(marlin::Global::GEAR->getHcalRingParameters().getExtent()[3]);
    const gear::LayerLayout &hCalRingLayerLayout(marlin::Global::GEAR->getHcalRingParameters().getLayerLayout());
    const float hCalRingLayerThickness(hCalRingLayerLayout.getThickness(hCalRingLayerLayout.getNLayers() - 1));

    const float hCalBarrelOuterR(marlin::Global::GEAR->getHcalBarrelParameters().getExtent()[1]);
    const float hCalBarrelOuterZ(marlin::Global::GEAR->getHcalBarrelParameters().getExtent()[3]);
    const float hCalBarrelOuterPhi0((std::find(marlin::Global::GEAR->getHcalBarrelParameters().getIntKeys().begin(),
        marlin::Global::GEAR->getHcalBarrelParameters().getIntKeys().end(),
        "Hcal_outer_polygon_phi0") != marlin::Global::GEAR->getHcalBarrelParameters().getIntKeys().end() ?
        marlin::Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_phi0") : 0));
    const unsigned int hCalBarrelOuterSymmetry((std::find(marlin::Global::GEAR->getHcalBarrelParameters().getIntKeys().begin(),
        marlin::Global::GEAR->getHcalBarrelParameters().getIntKeys().end(),
        "Hcal_outer_polygon_order") != marlin::Global::GEAR->getHcalBarrelParameters().getIntKeys().end() ?
        marlin::Global::GEAR->getHcalBarrelParameters().getIntVal("Hcal_outer_polygon_order") : 0));
    const gear::LayerLayout &hCalBarrelLayerLayout(marlin::Global::GEAR->getHcalBarrelParameters().getLayerLayout()); 
    const float hCalBarrelLayerThickness(hCalBarrelLayerLayout.getThickness(hCalBarrelLayerLayout.getNLayers() - 1));

    const float hCalEndCapOuterR(marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[1]);
    const float hCalEndCapInnerZ(marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[2]);
    const float hCalEndCapOuterZ(marlin::Global::GEAR->getHcalEndcapParameters().getExtent()[3]);
    const gear::LayerLayout &hCalEndCapLayerLayout(marlin::Global::GEAR->getHcalEndcapParameters().getLayerLayout());
    const float hCalEndCapLayerThickness(hCalEndCapLayerLayout.getThickness(hCalEndCapLayerLayout.getNLayers() - 1));
    const int hCalEndCapOuterSymmetryOrder(hCalBarrelOuterSymmetry);
    const float hCalEndCapOuterPhiCoordinate(hCalBarrelOuterPhi0);

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
                        m_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type CalorimeterHit " << std::endl;
                        throw;
                    }

                    const float hitEnergy(pCalorimeterHit->getEnergy());
                    hitEnergySum += hitEnergy;
                }
            }
        }
        catch (DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " << (*iter) << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::AddADCEntries(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, const unsigned int ADC_DC, TH1F *pTH1F) const
{
    const gear::CalorimeterParameters &ecalEndCapParameters(marlin::Global::GEAR->getEcalEndcapParameters());
    const float zOfEndCap = static_cast<float>(ecalEndCapParameters.getExtent()[2]);

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
                        m_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type SimCalorimeterHit " << std::endl;
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

                    if(ADC_DC)
                        pTH1F->Fill(SimCaloHitEnergy / correction);

                    else
                        pTH1F->Fill(1.f / correction);

                }
            }
        }
        catch (DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " << *iter << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CalibrationHelper::AddDirectionCorrectedCaloHitEntries(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, TH1F *pTH1F) const
{
    const gear::CalorimeterParameters &ecalEndCapParameters(marlin::Global::GEAR->getEcalEndcapParameters());
    float zOfEndCap = static_cast<float>(ecalEndCapParameters.getExtent()[2]);

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
                        m_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type CalorimeterHit " << std::endl;
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
            m_out(DEBUG) << "No Collection : " << (*iter) << std::endl;
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
    return "CalibrationHelper::GetNHCalLayersFromEdge,AddADCEntries,AddDirectionCorrectedCaloHitEntries - Possibility of division by zero. ";
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const char *HCalHitPositionException::what() const throw()
{
    return "CalibrationHelper::GetNHCalLayersFromEdge - HCal hit does not appear to be in the HCal. ";
}

} // namespace pandora_analysis
