// -*- C++ -*-
//
// Package:    HCALTest/HCAL_MET_Ana
// Class:      HCAL_MET_Ana
//
/**\class HCAL_MET_Ana HCAL_MET_Ana.cc HCALTest/HCAL_MET_Ana/plugins/HCAL_MET_Ana.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//DataFormats and Geometry
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

//ROOT includes
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"

//STL headers
#include <vector>
#include <string>
#include <iostream>
#include <numeric>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class HCAL_MET_Ana : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit HCAL_MET_Ana(const edm::ParameterSet&);
        ~HCAL_MET_Ana();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        bool print_channel;

        edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitsToken_;
        edm::EDGetTokenT<std::vector<reco::CaloMET>> CaloMETToken_;
        edm::EDGetTokenT<std::vector<reco::CaloMET>> CaloMETBEToken_;
        edm::EDGetTokenT<std::vector<reco::Muon>> MuonToken_;
        edm::EDGetTokenT<std::vector<reco::Vertex>> VertexToken_;
                
        TH1F * PU_h;
        TH1F * CaloMET_h, * CaloMET_phi_h;
        TH1F * CaloMETBE_h, * CaloMETBE_phi_h, * CaloMETBE20_phi_h;
        TH1F * myCaloMETBE_h, * myCaloMETBE_phi_h;
        TH1F * myCaloMETBE_Muon_h, * myCaloMETBE_Muon_phi_h;
        TH1F * myCaloMETBE1_Muon_h, * myCaloMETBE1_Muon_phi_h;
        TH1F * myCaloETBE_h;
        TH2F * CaloMETBE_vs_PU_h, * CaloMETBE_phi_vs_PU_h;
        TH2F * myCaloMETBE_Muon_vs_PU_h, * myCaloMETBE_Muon_phi_vs_PU_h;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HCAL_MET_Ana::HCAL_MET_Ana(const edm::ParameterSet& iConfig):
    print_channel(iConfig.getUntrackedParameter<bool>("print_channel"))
{
    //now do what ever initialization is needed
    HBHERecHitsToken_ = consumes<HBHERecHitCollection>(edm::InputTag("hbhereco"));
    CaloMETToken_ = consumes<std::vector<reco::CaloMET>>(edm::InputTag("caloMet"));
    CaloMETBEToken_ = consumes<std::vector<reco::CaloMET>>(edm::InputTag("caloMetBE"));
    MuonToken_ = consumes<std::vector<reco::Muon>>(edm::InputTag("muons"));
    VertexToken_ = consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));
        
    edm::Service<TFileService> TFS;
    PU_h = TFS->make<TH1F>("PU_h", "PU_h", 100, 0.0, 100.0);
    CaloMET_h = TFS->make<TH1F>("CaloMET_h", "CaloMET_h", 100, 0.0, 200.0);
    CaloMET_phi_h = TFS->make<TH1F>("CaloMET_phi_h", "CaloMET_phi_h", 100, -3.2, 3.2);
    CaloMETBE_h = TFS->make<TH1F>("CaloMETBE_h", "CaloMETBE_h", 100, 0.0, 200.0);
    CaloMETBE_phi_h = TFS->make<TH1F>("CaloMETBE_phi_h", "CaloMETBE_phi_h", 100, -3.2, 3.2);
    CaloMETBE20_phi_h = TFS->make<TH1F>("CaloMETBE20_phi_h", "CaloMETBE20_phi_h", 100, -3.2, 3.2);

    myCaloMETBE_h = TFS->make<TH1F>("myCaloMETBE_h", "myCaloMETBE_h", 100, 0.0, 200.0);
    myCaloMETBE_phi_h = TFS->make<TH1F>("myCaloMETBE_phi_h", "myCaloMETBE_phi_h", 100, -3.2, 3.2);
    myCaloMETBE_Muon_h = TFS->make<TH1F>("myCaloMETBE_Muon_h", "myCaloMETBE_Muon_h", 100, 0.0, 200.0);
    myCaloMETBE_Muon_phi_h = TFS->make<TH1F>("myCaloMETBE_Muon_phi_h", "myCaloMETBE_Muon_phi_h", 100, -3.2, 3.2);
    myCaloMETBE1_Muon_h = TFS->make<TH1F>("myCaloMETBE1_Muon_h", "myCaloMETBE1_Muon_h", 100, 0.0, 200.0);
    myCaloMETBE1_Muon_phi_h = TFS->make<TH1F>("myCaloMETBE1_Muon_phi_h", "myCaloMETBE1_Muon_phi_h", 100, -3.2, 3.2);

    myCaloETBE_h = TFS->make<TH1F>("myCaloETBE_h", "myCaloETBE_h", 200, 0.0, 1000.0);

    CaloMETBE_vs_PU_h = TFS->make<TH2F>("CaloMETBE_vs_PU_h", "CaloMETBE_vs_PU_h", 100, 0.0, 100.0, 100, 0.0, 200.0);
    CaloMETBE_phi_vs_PU_h = TFS->make<TH2F>("CaloMETBE_phi_vs_PU_h", "CaloMETBE_phi_vs_PU_h", 100, 0.0, 100.0, 100, -3.2, 3.2);
    myCaloMETBE_Muon_vs_PU_h = TFS->make<TH2F>("myCaloMETBE_Muon_vs_PU_h", "myCaloMETBE_Muon_vs_PU_h", 100, 0.0, 100.0, 100, 0.0, 200.0);
    myCaloMETBE_Muon_phi_vs_PU_h = TFS->make<TH2F>("myCaloMETBE_Muon_phi_vs_PU_h", "myCaloMETBE_Muon_phi_vs_PU_h", 100, 0.0, 100.0, 100, -3.2, 3.2);
}


HCAL_MET_Ana::~HCAL_MET_Ana()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------

void HCAL_MET_Ana::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<std::vector<reco::Muon>> MuonHandle;
    iEvent.getByToken(MuonToken_, MuonHandle);
    auto Muons = MuonHandle.product();

    TLorentzVector MuonTotTLV; 
    for(auto Muon : *Muons)
    {
        TLorentzVector MuonTLV;
        MuonTLV.SetPtEtaPhiM(Muon.p4().Pt(), Muon.p4().Eta(), Muon.p4().Phi(), Muon.p4().M());
        MuonTotTLV += MuonTLV;
    }

    edm::Handle<std::vector<reco::Vertex>> VertexHandle;
    iEvent.getByToken(VertexToken_, VertexHandle);
    auto Vertex = VertexHandle.product();
    int nPU = Vertex->size();
    PU_h->Fill(nPU);

    edm::ESHandle<CaloGeometry> CaloGeoHandle;
    iSetup.get<CaloGeometryRecord>().get(CaloGeoHandle);
    CaloGeometry CaloGeo = *CaloGeoHandle;

    edm::Handle<HBHERecHitCollection> HBHERecHitsHandle;
    iEvent.getByToken(HBHERecHitsToken_, HBHERecHitsHandle);
    auto HBHERecHits = HBHERecHitsHandle.product();

    TLorentzVector HBHETotTLV, HBHE1TotTLV;
    float HBHEET = 0; 
    for(auto HBHERecHit : *HBHERecHits)
    {
        auto Hid = HBHERecHit.id();
        auto RawId = Hid.rawId();
        auto SubDet = Hid.subdet();
        auto Depth = Hid.depth();
        auto Ieta = Hid.ieta();
        auto Iphi = Hid.iphi();
        auto Energy = HBHERecHit.energy();

        auto Eta = CaloGeo.getPosition(Hid).eta();
        auto Phi = CaloGeo.getPosition(Hid).phi();

        auto Theta = 2 * TMath::ATan(exp(-1 * Eta));
        auto ET = Energy * TMath::Sin(Theta);

        if(Energy > 0)
        {
            TLorentzVector HBHERecHitTLV;
            HBHERecHitTLV.SetPtEtaPhiM(ET, Eta, Phi, 0);
            HBHETotTLV += HBHERecHitTLV;
            if(ET > 1) {HBHE1TotTLV += HBHERecHitTLV;}
            HBHEET += ET;
        } 

        if(print_channel) std::cout << Hid << ", " << RawId << ", " << SubDet << ", " << Depth << ", " << Ieta << ", " << Eta << ", " << Iphi << ", " << Phi << ", " << Energy << std::endl;
    }
    //std::cout << iEvent.id().event() << ", " << HBHETotTLV.Pt() << ", " << HBHETotTLV.Eta() << ", " << HBHETotTLV.Phi() << std::endl;

    TLorentzVector myCaloMETBE = - HBHETotTLV;
    myCaloMETBE_h->Fill(myCaloMETBE.Pt());
    myCaloMETBE_phi_h->Fill(myCaloMETBE.Phi());
    myCaloETBE_h->Fill(HBHEET);

    TLorentzVector myCaloMETBE_Muon = - HBHETotTLV - MuonTotTLV;
    auto myCaloMETBE_Muon_pt = myCaloMETBE_Muon.Pt();
    auto myCaloMETBE_Muon_phi = myCaloMETBE_Muon.Phi();

    myCaloMETBE_Muon_h->Fill(myCaloMETBE_Muon_pt);
    myCaloMETBE_Muon_phi_h->Fill(myCaloMETBE_Muon_phi);
    myCaloMETBE_Muon_vs_PU_h->Fill(nPU, myCaloMETBE_Muon_pt);
    myCaloMETBE_Muon_phi_vs_PU_h->Fill(nPU, myCaloMETBE_Muon_phi);

    TLorentzVector myCaloMETBE1_Muon = - HBHE1TotTLV - MuonTotTLV;
    myCaloMETBE1_Muon_h->Fill(myCaloMETBE1_Muon.Pt());
    myCaloMETBE1_Muon_phi_h->Fill(myCaloMETBE1_Muon.Phi());

    edm::Handle<std::vector<reco::CaloMET>> CaloMETHandle;
    iEvent.getByToken(CaloMETToken_, CaloMETHandle);
    auto CaloMET = CaloMETHandle.product();

    if(CaloMET->size() == 1)
    {
        CaloMET_h->Fill(CaloMET->at(0).p4().Pt());
        CaloMET_phi_h->Fill(CaloMET->at(0).p4().Phi());
    }

    edm::Handle<std::vector<reco::CaloMET>> CaloMETBEHandle;
    iEvent.getByToken(CaloMETBEToken_, CaloMETBEHandle);
    auto CaloMETBE = CaloMETBEHandle.product();

    if(CaloMETBE->size() == 1)
    {
        auto CaloMETBE_pt = CaloMETBE->at(0).p4().Pt();
        auto CaloMETBE_phi = CaloMETBE->at(0).p4().Phi();

        CaloMETBE_h->Fill(CaloMETBE_pt);
        CaloMETBE_phi_h->Fill(CaloMETBE_phi);
        if(CaloMETBE_pt > 20){CaloMETBE20_phi_h->Fill(CaloMETBE_phi);}
        CaloMETBE_vs_PU_h->Fill(nPU, CaloMETBE_pt);
        CaloMETBE_phi_vs_PU_h->Fill(nPU, CaloMETBE_phi);
    }
}


// ------------ method called once each job just before starting event loop  ------------
void HCAL_MET_Ana::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HCAL_MET_Ana::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HCAL_MET_Ana::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HCAL_MET_Ana);
