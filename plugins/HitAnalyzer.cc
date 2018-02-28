
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// For event setup
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

// ROOT includes
#include "TLorentzVector.h" 
#include "TMath.h" 
#include "TTree.h" 

// For output
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// For btag
#include "DataFormats/BTauReco/interface/JetTag.h"


// For jet
#include "DataFormats/JetReco/interface/PFJet.h" 
#include "DataFormats/JetReco/interface/PFJetCollection.h" 

// #include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
// #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
// #include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

// For pixel clusters
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"

// Pixel topology
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 


#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"


//  For gen particle
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// For vertices
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"



class HitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HitAnalyzer(const edm::ParameterSet&);
      ~HitAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void reset( void );
      
      edm::EDGetTokenT< reco::GenParticleCollection>          genPtoken;
      edm::EDGetTokenT< reco::PFJetCollection >               ak4CHStoken;
      edm::EDGetTokenT< reco::PFJetCollection >               ak8CHStoken;
      edm::EDGetTokenT< edmNew::DetSetVector<SiPixelCluster>> clusterToken;
      edm::EDGetTokenT< reco::VertexCollection >              svToken;
      edm::EDGetTokenT< reco::TrackCollection >               trackToken;
      edm::EDGetTokenT< reco::JetTagCollection >              csv2Token;
     
           //
      // vector<reco::Vertex>                  "inclusiveSecondaryVertices"   ""                "RECO"
      //   vector<reco::TrackExtra>              "generalTracks"             ""                "RECO"
      //     vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >    "pfGhostTrackVertexTagInfos"   ""                "RECO"
      //     vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >    "pfInclusiveSecondaryVertexFinderCvsLTagInfos"   ""                "RECO"
      //     vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >    "pfInclusiveSecondaryVertexFinderTagInfos"   ""                "RECO"
      //     vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >    "pfSecondaryVertexTagInfos"   ""                "RECO"
      //       vector<reco::PFCandidate>             "particleFlow"              ""                "RECO"
      //         edm::ValueMap<float>                  "offlinePrimaryVertices"    ""                "RECO"
      //           dm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfChargeBJetTags"          ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedCvsBJetTags"     ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedCvsLJetTags"     ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedInclusiveSecondaryVertexV2BJetTags"   ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedMVAV2BJetTags"   ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfCombinedSecondaryVertexV2BJetTags"   ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfGhostTrackBJetTags"      ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfJetBProbabilityBJetTags"   ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfJetProbabilityBJetTags"   ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfSimpleInclusiveSecondaryVertexHighEffBJetTags"   ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfSimpleSecondaryVertexHighEffBJetTags"   ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfTrackCountingHighEffBJetTags"   ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "softPFElectronBJetTags"    ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "softPFMuonBJetTags"        ""                "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfDeepCMVAJetTags"         "probb"           "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfDeepCSVJetTags"          "probb"           "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfDeepCMVAJetTags"         "probbb"          "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfDeepCSVJetTags"          "probbb"          "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfDeepCMVAJetTags"         "probc"           "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfDeepCSVJetTags"          "probc"           "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfDeepCMVAJetTags"         "probcc"          "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfDeepCMVAJetTags"         "probudsg"        "RECO"
      //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfDeepCSVJetTags"          "probudsg"        "RECO"

      // ----------member data ---------------------------
      
       int                  nJets;
       std::vector<double>  jet_pt;
       std::vector<double>  jet_eta;
       std::vector<double>  jet_phi;
       std::vector<double>  jet_mass;
       std::vector<int>     jet_pdgId;
       std::vector<double>  jet_bTag;
       std::vector<double>  dr_jetGen;
       
       int                  nGenParticles;
       std::vector<double>  genParticle_pt;
       std::vector<double>  genParticle_eta;
       std::vector<double>  genParticle_phi;
       std::vector<double>  genParticle_mass;
       std::vector<int>     genParticle_pdgId;
       std::vector<int>     genParticle_status;
       
       int                               nDetUnits;
       std::vector<unsigned int >        detUnit_detType;
       std::vector<unsigned int >        detUnit_subdetId;
       
       std::vector<int>                  nClusters;
       std::vector<std::vector<double>>  cluster_x;
       std::vector<std::vector<double>>  cluster_y;
       std::vector<std::vector<double>>  cluster_globalz;
       std::vector<std::vector<double>>  cluster_globalx;
       std::vector<std::vector<double>>  cluster_globaly;
       std::vector<std::vector<double>>  cluster_globalPhi;
       std::vector<std::vector<double>>  cluster_globalR;
       
       
       TTree *tree;
       
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
HitAnalyzer::HitAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>( "tree", "tree" );
  
  tree->Branch( "nJets"             , &nJets );
  tree->Branch( "jet_pt"            , &jet_pt );
  tree->Branch( "jet_eta"           , &jet_eta );
  tree->Branch( "jet_phi"           , &jet_phi );
  tree->Branch( "jet_mass"          , &jet_mass );
  tree->Branch( "jet_pdgId"         , &jet_pdgId );
  tree->Branch( "jet_bTag"          , &jet_bTag );
  
  tree->Branch( "dr_jetGen"         , &dr_jetGen );
  
  tree->Branch( "nGenParticles"     , &nGenParticles );
  tree->Branch( "genParticle_pt"    , &genParticle_pt );
  tree->Branch( "genParticle_eta"   , &genParticle_eta );
  tree->Branch( "genParticle_phi"   , &genParticle_phi );
  tree->Branch( "genParticle_mass"  , &genParticle_mass );
  tree->Branch( "genParticle_pdgId" , &genParticle_pdgId );
  tree->Branch( "genParticle_status", &genParticle_status );
  
  tree->Branch( "nDetUnits"         ,&nDetUnits           );
  tree->Branch( "detUnit_detType"   ,&detUnit_detType     );
  tree->Branch( "detUnit_subdetId"  ,&detUnit_subdetId    );
  
  tree->Branch( "nClusters"         ,&nClusters           );
  tree->Branch( "cluster_x"         ,&cluster_x           );
  tree->Branch( "cluster_y"         ,&cluster_y           );
  tree->Branch( "cluster_globalz"   ,&cluster_globalz     );
  tree->Branch( "cluster_globalx"   ,&cluster_globalx     );
  tree->Branch( "cluster_globaly"   ,&cluster_globaly     );
  tree->Branch( "cluster_globalPhi" ,&cluster_globalPhi   );
  tree->Branch( "cluster_globalR"   ,&cluster_globalR     );
  
  
    
   
   std::string labelgenP("genParticles");
   std::string labelAK8s("ak8PFJetsCHS");
   std::string labelAK4s("ak4PFJetsCHS");
   std::string labelClusters("siPixelClusters");
   std::string labelSVs("inclusiveSecondaryVertices");
   std::string labelTracks("generalTracks");
   std::string labelCSV("pfCombinedSecondaryVertexV2BJetTags");
   genPtoken      = consumes< reco::GenParticleCollection         >(edm::InputTag(labelgenP));
   ak8CHStoken    = consumes< reco::PFJetCollection               >(edm::InputTag(labelAK8s));
   ak4CHStoken    = consumes< reco::PFJetCollection               >(edm::InputTag(labelAK4s));
   clusterToken   = consumes< edmNew::DetSetVector<SiPixelCluster>>(edm::InputTag(labelClusters));
   svToken        = consumes< reco::VertexCollection              >(edm::InputTag(labelSVs));
   trackToken     = consumes< reco::TrackCollection               >(edm::InputTag(labelTracks));
   csv2Token      = consumes< reco::JetTagCollection              >(edm::InputTag(labelCSV));

}


HitAnalyzer::~HitAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
  HitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // std::cout<< " -- NEW EVENT --" << std::endl;
    
  using namespace edm;
   
  // Make sure vectors are clear
  reset();
   
  // // Get event setup
  // edm::ESHandle<TrackerGeometry> geom;
  // iSetup.get<TrackerDigiGeometryRecord>().get( geom );
  // const TrackerGeometry& theTracker(*geom);
   
  // Get handles
  Handle<edmNew::DetSetVector<SiPixelCluster> > clusters ; iEvent.getByToken( clusterToken , clusters );
  Handle<reco::PFJetCollection                > ak8CHS   ; iEvent.getByToken( ak8CHStoken  , ak8CHS   );
  Handle<reco::GenParticleCollection          > genPs    ; iEvent.getByToken( genPtoken    , genPs    );
  Handle<reco::JetTagCollection               > CSVs     ; iEvent.getByToken( csv2Token    , CSVs    );
  const reco::JetTagCollection & bTags = *(CSVs.product());
  
   
  // Loop over jets
  
  for ( reco::PFJetCollection::const_iterator jet = ak8CHS->begin(); jet != ak8CHS->end(); ++jet ) {
    if (jet->pt()<200.) continue;
    TLorentzVector fatJet;
    fatJet.SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(),jet->mass());
    
    double minDR = 99.;
    int pdgId   = 0; 
    int status  = 0; 
    TLorentzVector selectedGenP ;
    for ( reco::GenParticleCollection::const_iterator gp = genPs->begin(); gp != genPs->end(); ++gp ) {
      if (gp->status()<20 || gp->status()>29) continue;
      TLorentzVector genP;
      genP.SetPtEtaPhiM(gp->pt(),gp->eta(),gp->phi(),gp->mass());
      float dR = fatJet.DeltaR(genP);
      if (dR < minDR ){
        minDR = dR;
        pdgId = gp->pdgId();
        status = gp->status();
        selectedGenP = genP;
      }
    }
    if (pdgId==0 || minDR > 0.8) continue;
    
    jet_pt.push_back(jet->pt());
    jet_eta.push_back(jet->eta());
    jet_phi.push_back(jet->phi());
    jet_mass.push_back(jet->mass());
    jet_pdgId.push_back(pdgId);
    dr_jetGen.push_back(minDR);
  
    genParticle_pt  .push_back(selectedGenP.Pt());
    genParticle_eta .push_back(selectedGenP.Eta());
    genParticle_phi .push_back(selectedGenP.Phi());
    genParticle_mass.push_back(selectedGenP.M());
    genParticle_pdgId.push_back(pdgId);
    genParticle_status.push_back(status);
    
    // b-tag infos
    double match = 0.8;
    double csv2 = -99.;
    for (unsigned int i = 0; i != bTags.size(); ++i) {
      if (bTags[i].first->pt()<170.) continue;
      TLorentzVector bTagJet;
      bTagJet.SetPtEtaPhiM(bTags[i].first->pt(),bTags[i].first->eta(),bTags[i].first->phi(),bTags[i].first->mass());
      float dR = fatJet.DeltaR(bTagJet);
      if (dR > match ) continue;
        match = dR;
        csv2 = bTags[i].second; 
    }
    jet_bTag.push_back(csv2);
  }
  nGenParticles = genParticle_pt.size();
  nJets         = jet_pt.size();
  

  
   
  // Get vector of detunit ids and loop over
  const edmNew::DetSetVector<SiPixelCluster>& input = *clusters;
  
  int numberOfDetUnits = 0;
  for ( edmNew::DetSetVector<SiPixelCluster>::const_iterator detUnit = input.begin(); detUnit != input.end(); ++detUnit ) {
    unsigned int detid = detUnit->detId();
    DetId detId = DetId(detid);       // Get the Detid object
    unsigned int detType=detId.det(); // det type, pixel=1
    if(detType!=1) continue; // look only at pixels
    numberOfDetUnits++;
    unsigned int subid=detId.subdetId(); //subdetector type, barrel=1
    detUnit_detType.push_back(detType);
    detUnit_subdetId.push_back(subid);
    
    // // // Get the geom-detector
    // const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detId) );
    // const PixelTopology * topol = &(theGeomDet->specificTopology());
    int numberOfClusters = 0;
    std::vector<double>  _cluster_x;
    std::vector<double>  _cluster_y;
    std::vector<double>  _cluster_globalz;
    std::vector<double>  _cluster_globalx;
    std::vector<double>  _cluster_globaly;
    std::vector<double>  _cluster_globalPhi;
    std::vector<double>  _cluster_globalR;
    
    for ( edmNew::DetSet<SiPixelCluster>::const_iterator clustIt = detUnit->begin(); clustIt != detUnit->end(); ++clustIt ) {
      numberOfClusters++;
      
      _cluster_x.push_back(clustIt->x());
      _cluster_y.push_back(clustIt->x());
      _cluster_globalz.push_back(clustIt->x());
      _cluster_globalx.push_back(clustIt->x());
      _cluster_globaly.push_back(clustIt->x());
      _cluster_globalPhi.push_back(clustIt->x());
      _cluster_globalR.push_back(clustIt->x());
   
    //   // get global position of the cluster
    //   LocalPoint lp = topol->localPosition(MeasurementPoint(clustIt->x(),clustIt->y()));
    //   // float lx = lp.x(); // local cluster position in cm
    //   // float ly = lp.y();
    //
    //   GlobalPoint clustgp = theGeomDet->surface().toGlobal( lp );
    //   double gZ = clustgp.z();  // global z
    //   double gX = clustgp.x();
    //   double gY = clustgp.y();
    //   TVector3 v(gX,gY,gZ);
    //   float gPhi = v.Phi(); // phi of the hit
    //   float gR = v.Perp(); // r of the hit
    //
    }
  nClusters.push_back(numberOfClusters);
  cluster_x        .push_back(_cluster_x);
  cluster_y        .push_back(_cluster_y);
  cluster_globalz  .push_back(_cluster_globalz);
  cluster_globalx  .push_back(_cluster_globalx);
  cluster_globaly  .push_back(_cluster_globaly);
  cluster_globalPhi.push_back(_cluster_globalPhi);  
  cluster_globalR  .push_back(_cluster_globalR);
    
    
  }
  nDetUnits = numberOfDetUnits;
  
  tree->Fill();
}

// Private methods
void HitAnalyzer::reset( void ){

  nJets = 0;
  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_mass.clear();
  jet_pdgId.clear();
  jet_bTag.clear();
  dr_jetGen.clear();
  
  nGenParticles = 0;
  genParticle_pt.clear();
  genParticle_eta.clear();
  genParticle_phi.clear();
  genParticle_mass.clear();
  genParticle_pdgId.clear();
  genParticle_status.clear();
  
  nClusters.clear();
  cluster_x.clear();
  cluster_y.clear();
  cluster_globalz.clear();
  cluster_globalx.clear();
  cluster_globaly.clear();
  cluster_globalPhi.clear();
  cluster_globalR.clear();
  
  nDetUnits = 0.;
  detUnit_detType.clear();
  detUnit_subdetId.clear();
  
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
HitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HitAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HitAnalyzer);
