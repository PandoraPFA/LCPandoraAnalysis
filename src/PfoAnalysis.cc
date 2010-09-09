/**
 *  @file   PandoraAnalysis/src/PfoAnalysis.cc
 * 
 *  @brief  Implementation of the pfo analysis class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1.h"

#include "PfoAnalysis.h"

#include <cmath>

PfoAnalysis pfoAnalysis;

//------------------------------------------------------------------------------------------------------------------------------------------

PfoAnalysis::PfoAnalysis() :
    Processor("PfoAnalysis")
{
    _description = "PfoAnalysis analyses output of PandoraPFA" ;

    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
                             "InputQuarkParticleCollections", 
                             "Names of input quark particle collections",
                             m_inputQuarkParticleCollections,
                             StringVector());

    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
                             "InputMCParticleCollections", 
                             "Names of input mc particle collections",
                             m_inputMCParticleCollections,
                             StringVector());

    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
                             "InputParticleCollections", 
                             "Names of input reconstructed particle collections",
                             m_inputParticleCollections,
                             StringVector());

    std::string rootFile("MyPFOAnalysis.root");
    registerProcessorParameter( "RootFile",
                                "Name of the Track collection used for clustering",
                                m_rootFile,
                                rootFile );

    registerProcessorParameter( "Printing",
                                "Set the debug print level",
                                m_printing,
                                (int)0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::init()
{
    m_nRun = 0;
    m_nEvt = 0;
    this->Clear();

    fNPFO = new TH1F("fNPFO", "number of pfos ", 200, 0., 200.);
    fPFAMZ = new TH1F("fPFAMZ", "vector boson mass", 200, 50., 150.);
    fPFAMW = new TH1F("fPFAMW", "vector boson mass", 200, 50., 150.);
    fPFAMZa = new TH1F("fPFAMZa", "vector boson mass", 200, 50., 150.);
    fPFAMWa = new TH1F("fPFAMWa", "vector boson mass", 200, 50., 150.);
    fEq = new TH1F("fEq", "Eq", 200, 0., 1000.);
    fPFA = new TH1F("fPFA", "total PFA energy", 10000, 0., 5000.);
    fPFAnu = new TH1F("fPFAnu", "total energy + nu", 10000, 0., 5000.);
    fPFAnufwd = new TH1F("fPFAnufwd", "total energy + nu + fwd", 10000, 0., 5000.);
    fPFAudscb = new TH1F("fPFAudscb", "total energy", 10000, 0., 5000.);
    fPFAuds = new TH1F("fPFAuds", "total energy", 10000, 0., 5000.);
    fPFAudsHM20 = new TH1F("fPFAudsHM20", "total energy", 10000, 0., 5000.);
    fPFAudsHM10 = new TH1F("fPFAudsHM10", "total energy", 10000, 0., 5000.);
    fPFAudsHP10 = new TH1F("fPFAudsHP10", "total energy", 10000, 0., 5000.);
    fPFAudsHP20 = new TH1F("fPFAudsHP20", "total energy", 10000, 0., 5000.);
    fPFAFudsHM20 = new TH1F("fPFAFudsHM20", "total energy",5000, 0., 250.);
    fPFAFudsHM10 = new TH1F("fPFAFudsHM10", "total energy",5000, 0., 250.);
    fPFAFudsHP10 = new TH1F("fPFAFudsHP10", "total energy",5000, 0., 250.);
    fPFAFudsHP20 = new TH1F("fPFAFudsHP20", "total energy",5000, 0., 250.);
    fPFAcb = new TH1F("fPFAcb", "total energy", 10000, 0., 5000.);
    fPFA1 = new TH1F("fPFA1", "total energy 0.-0.1", 10000, 0., 5000.);
    fPFA2 = new TH1F("fPFA2", "total energy 0.1-0.2", 10000, 0., 5000.);
    fPFA3 = new TH1F("fPFA3", "total energy 0.2-0.3", 10000, 0., 5000.);
    fPFA4 = new TH1F("fPFA4", "total energy 0.3-0.4", 10000, 0., 5000.);
    fPFA5 = new TH1F("fPFA5", "total energy 0.4-0.5", 10000, 0., 5000.);
    fPFA6 = new TH1F("fPFA6", "total energy 0.5-0.6", 10000, 0., 5000.);
    fPFA7 = new TH1F("fPFA7", "total energy 0.6-0.7", 10000, 0., 5000.);
    fPFAL7A = new TH1F("fPFAL7A", "total energy <0.7", 10000, 0., 5000.);
    fPFAL7Aud= new TH1F("fPFAL7Aud", "total energy <0.7 ud", 10000, 0., 5000.);
    fPFAL7As = new TH1F("fPFAL7As", "total energy <0.7 s", 10000, 0., 5000.);
    fPFAL7Ac = new TH1F("fPFAL7Ac", "total energy <0.7 c", 10000, 0., 5000.);
    fPFAL7Ab = new TH1F("fPFAL7Ab", "total energy <0.7 b", 10000, 0., 5000.);
    fPFAL7B = new TH1F("fPFAL7B", "total energy <0.7 - bad", 10000, 0., 5000.);
    fPFAFL7A = new TH1F("fPFAFL7A", "total energy <0.7",5000, 0., 250.);
    fPFAFL7Aud = new TH1F("fPFAFL7Aud", "total energy <0.7 ud",5000, 0., 250.);
    fPFAFL7As = new TH1F("fPFAFL7As", "total energy <0.7 s",5000, 0., 250.);
    fPFAFL7Ac = new TH1F("fPFAFL7Ac", "total energy <0.7 c",5000, 0., 250.);
    fPFAFL7Ab = new TH1F("fPFAFL7Ab", "total energy <0.7 b",5000, 0., 250.);
    fPFAQQ   = new TH1F("fPFAQQ", "total energy - true E", 2000, -500., 500.);
    fPFAQQ8  = new TH1F("fPFAQQ8", "total energy - true E 8", 2000, -500., 500.);
    fPFA8  = new TH1F("fPFA8", "total energy 0.7-0.8", 10000, 0., 5000.);
    fPFA9  = new TH1F("fPFA9", "total energy 0.8-0.9", 10000, 0., 5000.);
    fPFA10 = new TH1F("fPFA10","total energy 0.9-1.", 10000, 0., 5000.);
    fPFA11 = new TH1F("fPFA11","total energy 0.9-0.925", 10000, 0., 5000.);
    fPFA12 = new TH1F("fPFA12","total energy 0.925-0.95", 10000, 0., 5000.);
    fPFA13 = new TH1F("fPFA13","total energy 0.95-0.975", 10000, 0., 5000.);
    fPFA14 = new TH1F("fPFA14","total energy 0.975-1.", 10000, 0., 5000.);
    fPFADMZ = new TH1F("fPFADMZ", "Delta Mz", 200, -50., 50.);
    fPFADMZ8 = new TH1F("fPFADMZ8", "Delta Mz", 200, -50., 50.);
    fPFADMZQQ8 = new TH1F("fPFADMZQQ8", "Delta Mz", 200, -50., 50.);
    fPFADMZP8 = new TH1F("fPFADMZP8", "Delta Mz", 200, -50., 50.);
    fPFADMZOMZ = new TH1F("fPFADMZOMZ" , "Delta Mz / Mz", 200, -25., 25.);
    fPFADMZOMZQQ8 = new TH1F("fPFADMZOMZQQ8" , "Delta Mz / Mz", 200, -25., 25.);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::processRunHeader(lcio::LCRunHeader *pLCRunHeader)
{
    m_nRun++;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::processEvent(lcio::LCEvent *pLCEvent)
{
    m_nEvt++;
    this->Clear();

    // Extract quark particle collection
    for (StringVector::const_iterator iter = m_inputQuarkParticleCollections.begin(), iterEnd = m_inputQuarkParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
            {
                ReconstructedParticle *pReconstructedParticle = dynamic_cast<ReconstructedParticle*>(pLCCollection->getElementAt(i));
                m_quarkpfovec.push_back(pReconstructedParticle);
            }

            this->MakeQuarkVariables();
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Could not extract input quark particle collection: " << *iter << std::endl;
        }
    }

    // Extract mc particle collection
    for (StringVector::const_iterator iter = m_inputMCParticleCollections.begin(), iterEnd = m_inputMCParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
            {
                ReconstructedParticle *pReconstructedParticle = dynamic_cast<ReconstructedParticle*>(pLCCollection->getElementAt(i));
                m_mcpfovec.push_back(pReconstructedParticle);
            }
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Could not extract input mc particle collection: " << *iter << std::endl;
        }
    }

    // Extract reconstructed particle collection
    for (StringVector::const_iterator iter = m_inputParticleCollections.begin(), iterEnd = m_inputParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
            {
                ReconstructedParticle *pReconstructedParticle = dynamic_cast<ReconstructedParticle*>(pLCCollection->getElementAt(i));
                m_pfovec.push_back(pReconstructedParticle);
            }
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Could not extract input particle collection: " << *iter << std::endl;
        }
    }

    this->AnalysePFAPerformance();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::check(lcio::LCEvent *pLCEvent)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::end()
{
    if (m_printing > -1)
    {
        std::cout << "PfoAnalysis::end()  " << name()
                  << " processed " << m_nEvt << " events in " << m_nRun << " runs " << std::endl
                  << "Rootfile: " << m_rootFile.c_str() << std::endl;
    }

    TFile *hfile = new TFile(m_rootFile.c_str(), "recreate");

    fEq->TH1F::Write();
    fPFA->TH1F::Write();
    fPFAMW->TH1F::Write();
    fPFAMZ->TH1F::Write();
    fPFAMWa->TH1F::Write();
    fPFAMZa->TH1F::Write();
    fPFAnu->TH1F::Write();
    fPFAnufwd->TH1F::Write();
    fPFAuds->TH1F::Write();
    fPFAudsHP20->TH1F::Write();
    fPFAudsHP10->TH1F::Write();
    fPFAudsHM10->TH1F::Write();
    fPFAudsHM20->TH1F::Write();
    fPFAFudsHP20->TH1F::Write();
    fPFAFudsHP10->TH1F::Write();
    fPFAFudsHM10->TH1F::Write();
    fPFAFudsHM20->TH1F::Write();
    fPFAcb->TH1F::Write();
    fPFAudscb->TH1F::Write();
    fPFA1->TH1F::Write();
    fPFA2->TH1F::Write();
    fPFA3->TH1F::Write();
    fPFA4->TH1F::Write();
    fPFA5->TH1F::Write();
    fPFA6->TH1F::Write();
    fPFA7->TH1F::Write();
    fPFAQQ->TH1F::Write();
    fPFAQQ8->TH1F::Write();
    fPFAL7A->TH1F::Write();
    fPFAL7Aud->TH1F::Write();
    fPFAL7As->TH1F::Write();
    fPFAL7Ac->TH1F::Write();
    fPFAL7Ab->TH1F::Write();
    fPFAL7B->TH1F::Write();
    fPFAFL7A->TH1F::Write();
    fPFAFL7Aud->TH1F::Write();
    fPFAFL7As->TH1F::Write();
    fPFAFL7Ac->TH1F::Write();
    fPFAFL7Ab->TH1F::Write();
    fPFA8->TH1F::Write();
    fPFA9->TH1F::Write();
    fPFA10->TH1F::Write();
    fPFA11->TH1F::Write();
    fPFA12->TH1F::Write();
    fPFA13->TH1F::Write();
    fPFA14->TH1F::Write();
    fPFADMZ->TH1F::Write();
    fPFADMZOMZ->TH1F::Write();
    fPFADMZQQ8->TH1F::Write();
    fPFADMZOMZQQ8->TH1F::Write();
    fPFADMZ8->TH1F::Write();
    fPFADMZP8->TH1F::Write();

    hfile->Close();
    delete hfile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::AnalysePFAPerformance()
{
    // Loop over reconstructed PFOs
    float totenergy = 0, totgammaenergy = 0, tothadenergy = 0;
    float totmom[3] = {0., 0., 0.};
    unsigned int ntrack = 0, nclust = 0, npfo = 0;

    for (unsigned int i = 0; i < m_pfovec.size(); ++i)
    {
        npfo++;
        const TrackVec &tracks = m_pfovec[i]->getTracks();
        const ClusterVec &clusters = m_pfovec[i]->getClusters();

        if (!tracks.empty())
        {
            ntrack++;
        }
        else
        {
            nclust++;
        }

        totenergy += m_pfovec[i]->getEnergy();
        totmom[0] += m_pfovec[i]->getMomentum()[0];
        totmom[1] += m_pfovec[i]->getMomentum()[1];
        totmom[2] += m_pfovec[i]->getMomentum()[2];

        if (m_printing > 0)
        {
            std::cout << " RECOPFO : " << i << " " << m_pfovec[i]->getType() << " " << m_pfovec[i]->getEnergy()
                      << "  " << tracks.size() << ":" << clusters.size() << " charge : " << m_pfovec[i]->getCharge() << std::endl;
        }

        if ((m_pfovec[i]->getType() != 22) && tracks.empty())
        {
            tothadenergy += m_pfovec[i]->getEnergy();
        }

        if ((m_pfovec[i]->getType() == 22) &&tracks.empty())
        {
            totgammaenergy += m_pfovec[i]->getEnergy();
        }
    }

    const float totmass(std::sqrt(totenergy * totenergy - totmom[0] * totmom[0] - totmom[1] * totmom[1] - totmom[2] * totmom[2]));

    // Loop over MC PFOs and find energy in "primary" neutrinos
    float enu = 0, fwdenergy = 0;

    for (unsigned int imc = 0; imc < m_mcpfovec.size(); ++imc)
    {
        const int pdgCode(m_mcpfovec[imc]->getType());

        if ((std::abs(pdgCode) == 12) || (std::abs(pdgCode) == 14) || (std::abs(pdgCode) == 16))
        {
            enu += m_mcpfovec[imc]->getEnergy();
        }

        // Find polar angle of MC PFO
        const float px(m_mcpfovec[imc]->getMomentum()[0]);
        const float py(m_mcpfovec[imc]->getMomentum()[1]);
        const float pz(m_mcpfovec[imc]->getMomentum()[2]);

        if (std::fabs(pz) / std::sqrt(px * px + py * py + pz * pz) > 0.98)
        {
            fwdenergy += m_mcpfovec[imc]->getEnergy();
        }
    }

    if (m_printing > -1)
    {
        std::cout << " EVENT                : " << m_nEvt << std::endl;
        std::cout << " NPFOs                : " << ntrack + nclust << " ( " << ntrack << " + " << nclust << " )"  << std::endl; 
        std::cout << " RECONSTRUCTED ENERGY : " << totenergy << std::endl; 
        std::cout << " RECO ENERGY + ENU    : " << totenergy + enu << std::endl; 
    }

    // Fill histograms with event information
    fNPFO->Fill(npfo);

    if (std::fabs(m_costz) < 0.90)
        fPFAMWa->Fill(totmass);

    if (std::fabs(m_mz - 80.3) < 4 && std::fabs(m_costz) < 0.90)
        fPFAMW->Fill(totmass);

    if (std::fabs(m_costz) < 0.90)
        fPFAMZa->Fill(totmass);

    if (std::fabs(m_mz - 91.2) < 4 && std::fabs(m_costz) < 0.90)
        fPFAMZ->Fill(totmass);

    if (m_mz > 0.1 && std::fabs(m_costz) < 0.90)
    {
        fPFAQQ->Fill(totenergy - m_ez);
        fPFADMZ->Fill(totmass - m_mz);

        if (m_mz > 75)
            fPFADMZOMZ->Fill(100 * (totmass - m_mz) / totmass);

        if (std::fabs(m_costz) < 0.8)
            fPFADMZ8->Fill(totmass - m_mz);

        if (std::fabs(m_costz) > 0.8)
            fPFADMZP8->Fill(totmass - m_mz);

        if (std::fabs(m_costq1) < 0.8 && std::fabs(m_costq2) < 0.8)
        {
            fPFADMZQQ8->Fill(totmass - m_mz);
            fPFADMZOMZQQ8->Fill(100 * (totmass - m_mz) / totmass);
            fPFAQQ8->Fill(totenergy - m_ez);
        }
    }

    fPFA->Fill(totenergy, 1.);
    fPFAnu->Fill(totenergy + enu, 1.);
    fPFAnufwd->Fill(totenergy + enu + fwdenergy, 1.);

    if (m_qpdg >= 1 && m_qpdg <= 3)
        fPFAuds->Fill(totenergy, 1.);

    if (m_qpdg >= 4 && m_qpdg <= 5)
        fPFAcb->Fill(totenergy, 1.);

    if (m_qpdg >= 1 && m_qpdg <= 5)
        fPFAudscb->Fill(totenergy, 1.);

    if (m_qpdg >= 1 && m_qpdg <= 3)
    {
        std::cout << " Thrust : " << m_thrust << std::endl;

        if (m_thrust <= 0.1)
            fPFA1->Fill(totenergy, 1.);

        if (m_thrust > 0.1 && m_thrust <= 0.2)
            fPFA2->Fill(totenergy, 1.);

        if (m_thrust > 0.2 && m_thrust <= 0.3)
            fPFA3->Fill(totenergy, 1.);

        if (m_thrust > 0.3 && m_thrust <= 0.4)
            fPFA4->Fill(totenergy, 1.);

        if (m_thrust > 0.4 && m_thrust <= 0.5)
            fPFA5->Fill(totenergy, 1.);

        if (m_thrust > 0.5 && m_thrust <= 0.6)
            fPFA6->Fill(totenergy, 1.);

        if (m_thrust <= 0.7)
        {
            std::cout << "Filling L7A : " << totenergy+enu << std::endl; 

            fPFAL7A->Fill(totenergy + enu, 1.);
            fPFAFL7A->Fill(totenergy + enu, 1.);

            if (m_qpdg <= 2)
                fPFAL7Aud->Fill(totenergy + enu, 1.);

            if (m_qpdg == 3)
                fPFAL7As->Fill(totenergy + enu, 1.);

            if (m_qpdg == 4)
                fPFAL7Ac->Fill(totenergy + enu, 1.);

            if (m_qpdg == 5)
                fPFAL7Ab->Fill(totenergy + enu, 1.);

            if (m_qpdg <= 2)
                fPFAFL7Aud->Fill(totenergy + enu, 1.);

            if (m_qpdg == 3)
                fPFAFL7As->Fill(totenergy + enu, 1.);

            if (m_qpdg == 4)
                fPFAFL7Ac->Fill(totenergy + enu, 1.);

            if (m_qpdg == 5)
                fPFAFL7Ab->Fill(totenergy + enu, 1.);

            fPFAudsHM20->Fill(totenergy + enu - tothadenergy * 0.1, 1.);
            fPFAudsHM10->Fill(totenergy + enu - tothadenergy * 0.05, 1.);
            fPFAudsHP10->Fill(totenergy + enu + tothadenergy * 0.05, 1.);
            fPFAudsHP20->Fill(totenergy + enu + tothadenergy * 0.1, 1.);

            fPFAFudsHM20->Fill(totenergy + enu  -tothadenergy * 0.1, 1.);
            fPFAFudsHM10->Fill(totenergy + enu - tothadenergy * 0.05, 1.);
            fPFAFudsHP10->Fill(totenergy + enu + tothadenergy * 0.05, 1.);
            fPFAFudsHP20->Fill(totenergy + enu + tothadenergy * 0.1, 1.);
        }

        if (m_thrust <= 0.7 && fwdenergy / totenergy < 0.01)
            fPFAL7B->Fill(totenergy + enu, 1.);

        if (m_thrust > 0.6 && m_thrust <= 0.7)
            fPFA7->Fill(totenergy, 1.);

        if (m_thrust > 0.7 && m_thrust <= 0.8)
            fPFA8->Fill(totenergy, 1.);

        if (m_thrust > 0.8 && m_thrust <= 0.9)
            fPFA9->Fill(totenergy, 1.);

        if (m_thrust > 0.9)
            fPFA10->Fill(totenergy, 1.);

        if (m_thrust > 0.9  && m_thrust <= 0.925)
            fPFA11->Fill(totenergy, 1.);

        if (m_thrust > 0.925 && m_thrust <= 0.95)
            fPFA12->Fill(totenergy, 1.);

        if (m_thrust > 0.95 && m_thrust <= 0.975)
            fPFA13->Fill(totenergy, 1.);

        if (m_thrust > 0.975)
            fPFA14->Fill(totenergy, 1.);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::Clear()
{
    m_pfovec.clear();
    m_mcpfovec.clear();
    m_quarkpfovec.clear();

    m_ez = 0.;
    m_costq1 = 0.;
    m_costq2 = 0.;
    m_costz = 0.;
    m_mz = 0.;
    m_thrust = 0;
    m_qpdg = 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::MakeQuarkVariables()
{
    if (!m_quarkpfovec.empty())
    {
        m_qpdg = std::abs(m_quarkpfovec[0]->getType());
        float energytot = 0;
        float costtot = 0;

        for (unsigned int i = 0; i < m_quarkpfovec.size(); ++i)
        {
            const float px(m_quarkpfovec[i]->getMomentum()[0]);
            const float py(m_quarkpfovec[i]->getMomentum()[1]);
            const float pz(m_quarkpfovec[i]->getMomentum()[2]);
            const float energy(m_quarkpfovec[i]->getEnergy());
            const float p(std::sqrt(px * px + py * py + pz * pz));
            const float cost(std::fabs(pz) / p);
            energytot += energy;
            costtot += cost * energy;
        }

        m_thrust = costtot / energytot;
    }

    if (m_quarkpfovec.size() == 2)
    {
        TLorentzVector q1(m_quarkpfovec[0]->getMomentum()[0], m_quarkpfovec[0]->getMomentum()[1], m_quarkpfovec[0]->getMomentum()[2], m_quarkpfovec[0]->getEnergy());
        TLorentzVector q2(m_quarkpfovec[1]->getMomentum()[0], m_quarkpfovec[1]->getMomentum()[1], m_quarkpfovec[1]->getMomentum()[2], m_quarkpfovec[1]->getEnergy());
        TLorentzVector z = q1 + q2;
        m_mz = z.M();
        m_ez = m_quarkpfovec[0]->getEnergy() + m_quarkpfovec[1]->getEnergy(); 

        float pq1[3];
        pq1[0] = m_quarkpfovec[0]->getMomentum()[0];
        pq1[1] = m_quarkpfovec[0]->getMomentum()[1];
        pq1[2] = m_quarkpfovec[0]->getMomentum()[2];

        float pq2[3];
        pq2[0] = m_quarkpfovec[1]->getMomentum()[0];
        pq2[1] = m_quarkpfovec[1]->getMomentum()[1];
        pq2[2] = m_quarkpfovec[1]->getMomentum()[2];

        float pz[3];
        pz[0] = m_quarkpfovec[0]->getMomentum()[0] + m_quarkpfovec[1]->getMomentum()[0];
        pz[1] = m_quarkpfovec[0]->getMomentum()[1] + m_quarkpfovec[1]->getMomentum()[1];
        pz[2] = m_quarkpfovec[0]->getMomentum()[2] + m_quarkpfovec[1]->getMomentum()[2];

        const float pq1tot(std::sqrt(pq1[0] * pq1[0] + pq1[1] * pq1[1]  +pq1[2] * pq1[2]));
        const float pq2tot(std::sqrt(pq2[0] * pq2[0] + pq2[1] * pq2[1] + pq2[2] * pq2[2]));
        fEq->Fill(pq1tot);
        fEq->Fill(pq2tot);

        const float pztot(std::sqrt(pz[0] * pz[0] + pz[1] * pz[1] + pz[2] * pz[2]));

        if (0. != pztot)
            m_costz  = pz[2] / pztot;

        if (0. != pq1tot)
            m_costq1 = pq1[2] / pq1tot;

        if (0. != pq2tot)
            m_costq2 = pq2[2] / pq2tot;

        m_eq1 = pq1tot;
        m_eq2 = pq2tot;

        std::cout << " Eqq    = " << m_ez << std::endl
                  << "  q1    = " << pq1tot << std::endl
                  << "  q2    = " << pq2tot << std::endl
                  << " Costq1 = " << m_costq1 << std::endl
                  << " Costq2 = " << m_costq2 << std::endl
                  << " Costz  = " << m_costz << std::endl
                  << " Mqq    = " << m_mz  << std::endl
                  << " Thrust = " << m_thrust  << std::endl
                  << " QPDG   = " << m_qpdg  << std::endl;
    }
}
