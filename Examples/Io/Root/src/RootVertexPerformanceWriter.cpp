// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootVertexPerformanceWriter.hpp"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootVertexPerformanceWriter::RootVertexPerformanceWriter(
    const ActsExamples::RootVertexPerformanceWriter::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputVertices, "RootVertexPerformanceWriter", lvl),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  if (cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.outputTreename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  if (m_cfg.inputAllTruthParticles.empty()) {
    throw std::invalid_argument("Collection with all truth particles missing");
  }
  if (m_cfg.inputSelectedTruthParticles.empty()) {
    throw std::invalid_argument(
        "Collection with selected truth particles missing");
  }
  if (m_cfg.inputAssociatedTruthParticles.empty()) {
    throw std::invalid_argument(
        "Collection with track-associated truth particles missing");
  }
  if (m_cfg.inputFittedTracks.empty()) {
    throw std::invalid_argument(
        "Collection with all fitted track parameters missing");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
    m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path);
    }
  }
  m_outputFile->cd();
  m_outputTree =
      new TTree(m_cfg.outputTreename.c_str(), m_cfg.outputTreename.c_str());
  m_outputTree_Truth = 
      new TTree(m_cfg.outputTreename_Truth.c_str(), m_cfg.outputTreename_Truth.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("diffx", &m_diffx);
    m_outputTree->Branch("diffy", &m_diffy);
    m_outputTree->Branch("diffz", &m_diffz);
    m_outputTree->Branch("nRecoVtx", &m_nrecoVtx);
    m_outputTree->Branch("nTrueVtx", &m_ntrueVtx);
    m_outputTree->Branch("nVtxDetectorAcceptance", &m_nVtxDetAcceptance);
    m_outputTree->Branch("nVtxReconstructable", &m_nVtxReconstructable);

    m_outputTree_Truth->Branch("truth_vtx_vx", &m_truth_vtx_vx);
    m_outputTree_Truth->Branch("truth_vtx_vy", &m_truth_vtx_vy);
    m_outputTree_Truth->Branch("truth_vtx_vz", &m_truth_vtx_vz);

    m_outputTree_Truth->Branch("truth_vtx_trk_d0", &m_truth_vtx_trk_d0);
    m_outputTree_Truth->Branch("truth_vtx_trk_z0", &m_truth_vtx_trk_z0);
    m_outputTree_Truth->Branch("truth_vtx_trk_phi", &m_truth_vtx_trk_phi);
    m_outputTree_Truth->Branch("truth_vtx_trk_theta", &m_truth_vtx_trk_theta);
    m_outputTree_Truth->Branch("truth_vtx_trk_qp", &m_truth_vtx_trk_qp);
    m_outputTree_Truth->Branch("truth_vtx_trk_time", &m_truth_vtx_trk_time);

    m_outputTree_Truth->Branch("truth_vtx_trk_err_d0", &m_truth_vtx_trk_err_d0);
    m_outputTree_Truth->Branch("truth_vtx_trk_err_z0", &m_truth_vtx_trk_err_z0);
    m_outputTree_Truth->Branch("truth_vtx_trk_err_phi", &m_truth_vtx_trk_err_phi);
    m_outputTree_Truth->Branch("truth_vtx_trk_err_theta", &m_truth_vtx_trk_err_theta);
    m_outputTree_Truth->Branch("truth_vtx_trk_err_qp", &m_truth_vtx_trk_err_qp);
    m_outputTree_Truth->Branch("truth_vtx_trk_err_time", &m_truth_vtx_trk_err_time);
    m_outputTree_Truth->Branch("truth_vtx_trk_vtxID", &m_truth_vtx_trk_vtxID);
  }
}

ActsExamples::RootVertexPerformanceWriter::~RootVertexPerformanceWriter() {
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootVertexPerformanceWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    m_outputTree_Truth->Write();
  }
  return ProcessCode::SUCCESS;
}

int ActsExamples::RootVertexPerformanceWriter::
    getNumberOfReconstructableVertices(
        const SimParticleContainer& collection) const {
  // map for finding frequency
  std::map<int, int> fmap;

  std::vector<int> reconstructableTruthVertices;

  // traverse the array for frequency
  for (const auto& p : collection) {
    int secVtxId = p.particleId().vertexSecondary();
    if (secVtxId != 0) {
      // truthparticle from secondary vtx
      continue;
    }
    int priVtxId = p.particleId().vertexPrimary();
    fmap[priVtxId]++;
  }

  // iterate over the map
  for (auto it : fmap) {
    // Require at least 2 tracks
    if (it.second > 1) {
      reconstructableTruthVertices.push_back(it.first);
    }
  }

  return reconstructableTruthVertices.size();
}

std::vector<PV > ActsExamples::RootVertexPerformanceWriter::
    getTruthVerticesVec(
        const SimParticleContainer& collection) {
  // map for finding frequency
  std::vector<PV > PV_list;
  
  long unsigned int r = (*collection.rbegin()).particleId().vertexPrimary();
  PV_list.resize(r); 
  
  auto it = collection.cbegin();
  int idx_particle = 0;

  for (long unsigned int i = 0; i < r; i++)
  {
    PV_list[i].PV_loc[0] = (*it).position()[0];
    PV_list[i].PV_loc[1] = (*it).position()[1];
    PV_list[i].PV_loc[2] = (*it).position()[2];

    while ((*it).particleId().vertexPrimary() == i+1 )
    {
      PV_list[i].track_ID.push_back(idx_particle);
      ++it;
      ++idx_particle;
    }
    
  }
  
  return PV_list;
}

void ActsExamples::RootVertexPerformanceWriter::writeTruthInfo(
    const TrackParametersContainer& inputFittedTracks,
    std::vector<PV> PV_list) {
  for (size_t i = 0; i < PV_list.size(); ++i) {
    if (PV_list[i].track_ID.size() > 1) {
      m_truth_vtx_vx.push_back(PV_list[i].PV_loc[0]);
      m_truth_vtx_vy.push_back(PV_list[i].PV_loc[1]);
      m_truth_vtx_vz.push_back(PV_list[i].PV_loc[2]);

      for (size_t j = 0; j < PV_list[i].track_ID.size(); ++j) {
        const auto& boundParam = inputFittedTracks[PV_list[i].track_ID[j]];
        const auto& parameter = boundParam.parameters();

        m_truth_vtx_trk_d0.push_back(parameter[Acts::eBoundLoc0]);
        m_truth_vtx_trk_z0.push_back(parameter[Acts::eBoundLoc1]);
        m_truth_vtx_trk_phi.push_back(parameter[Acts::eBoundPhi]);
        m_truth_vtx_trk_theta.push_back(parameter[Acts::eBoundTheta]);
        m_truth_vtx_trk_qp.push_back(parameter[Acts::eBoundQOverP]);
        m_truth_vtx_trk_time.push_back(parameter[Acts::eBoundTime]);


        if (boundParam.covariance().has_value()) {
          const auto& covariance = *boundParam.covariance();
          m_truth_vtx_trk_err_d0.push_back(
              sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
          m_truth_vtx_trk_err_z0.push_back(
              sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
          m_truth_vtx_trk_err_phi.push_back(
              sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
          m_truth_vtx_trk_err_theta.push_back(
              sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
          m_truth_vtx_trk_err_qp.push_back(
              sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
          m_truth_vtx_trk_err_time.push_back(
              sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));
        } else {
          m_truth_vtx_trk_err_d0.push_back(NaNfloat);
          m_truth_vtx_trk_err_z0.push_back(NaNfloat);
          m_truth_vtx_trk_err_phi.push_back(NaNfloat);
          m_truth_vtx_trk_err_theta.push_back(NaNfloat);
          m_truth_vtx_trk_err_qp.push_back(NaNfloat);
          m_truth_vtx_trk_err_time.push_back(NaNfloat);
        }
        m_truth_vtx_trk_vtxID.push_back(m_truth_vtx_vx.size() - 1);
      }
    }
  }
}

int ActsExamples::RootVertexPerformanceWriter::getNumberOfTruePriVertices(
    const SimParticleContainer& collection) const {
  // Vector to store indices of all primary vertices
  std::set<int> allPriVtxIds;
  for (const auto& p : collection) {
    int priVtxId = p.particleId().vertexPrimary();
    int secVtxId = p.particleId().vertexSecondary();
    if (secVtxId != 0) {
      // truthparticle from secondary vtx
      continue;
    }
    // Insert to set, removing duplicates
    allPriVtxIds.insert(priVtxId);
  }
  // Size of set corresponds to total number of primary vertices
  return allPriVtxIds.size();
}

ActsExamples::ProcessCode ActsExamples::RootVertexPerformanceWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);
  // Event Number 
  m_eventNr = ctx.eventNumber;
  m_nrecoVtx = vertices.size();

  ACTS_DEBUG("Number of reco vertices in event: " << m_nrecoVtx);
  if (m_outputFile == nullptr)
    return ProcessCode::SUCCESS;

  // Read truth particle input collection
  const auto& allTruthParticles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputAllTruthParticles);
  // Get number of generated true primary vertices
  m_ntrueVtx = getNumberOfTruePriVertices(allTruthParticles);

  ACTS_INFO("Total number of generated truth particles in event : "
            << allTruthParticles.size());
  ACTS_INFO(
      "Total number of generated truth primary vertices : " << m_ntrueVtx);

  // Read selected truth particle input collection
  const auto& selectedTruthParticles = ctx.eventStore.get<SimParticleContainer>(
      m_cfg.inputSelectedTruthParticles);
  // Get number of detector-accepted true primary vertices
  m_nVtxDetAcceptance = getNumberOfTruePriVertices(selectedTruthParticles);

  ACTS_INFO("Total number of selected truth particles in event : "
            << selectedTruthParticles.size());
  ACTS_INFO("Total number of detector-accepted truth primary vertices : "
            << m_nVtxDetAcceptance);

  // Read track-associated truth particle input collection
  const auto& associatedTruthParticles =
      ctx.eventStore.get<SimParticleContainer>(
          m_cfg.inputAssociatedTruthParticles);//
  // Get number of track-associated true primary vertices
  m_nVtxReconstructable =
      getNumberOfReconstructableVertices(associatedTruthParticles);

  std::vector<PV > PV_list = getTruthVerticesVec(associatedTruthParticles);


  ACTS_INFO("Total number of reco track-associated truth particles in event : "
            << associatedTruthParticles.size());
  ACTS_INFO("Total number of reco track-associated truth primary vertices : "
            << m_nVtxReconstructable);

  /*****************  Start x,y,z resolution plots here *****************/
  // Matching tracks at vertex to fitted tracks that are in turn matched
  // to truth particles. Match reco and true vtx if >50% of tracks match

  const auto& inputFittedTracks =
      ctx.eventStore.get<std::vector<Acts::BoundTrackParameters>>(
          m_cfg.inputFittedTracks); //

  if (associatedTruthParticles.size() != inputFittedTracks.size()) {
    ACTS_WARNING(
        "Number of fitted tracks and associated truth particles do not match. "
        "Not able to match fitted tracks at reconstructed vertex to truth "
        "vertex.");
  } else {
    
    writeTruthInfo(inputFittedTracks, PV_list);
    // Loop over all reco vertices and find associated truth particles
    std::vector<SimParticleContainer> truthParticlesAtVtxContainer;
    for (const auto& vtx : vertices) {
      const auto tracks = vtx.tracks();
      // Store all associated truth particles to current vtx
      SimParticleContainer particleAtVtx;

      std::vector<int> contributingTruthVertices;

      for (const auto& trk : tracks) {
        Acts::BoundTrackParameters origTrack = *(trk.originalParams);

        // Find associated truth particle now
        int idx = 0;
        for (const auto& particle : associatedTruthParticles) {
          if (origTrack.parameters() == inputFittedTracks[idx].parameters()) {
            particleAtVtx.insert(particleAtVtx.end(), particle);

            int priVtxId = particle.particleId().vertexPrimary();
            contributingTruthVertices.push_back(priVtxId);
          }
          idx++;
        }
      }  // end loop tracks

      // Now find true vtx with most matching tracks at reco vtx
      // and check if it contributes more than 50 of all tracks
      std::map<int, int> fmap;
      for (int priVtxId : contributingTruthVertices) {
        fmap[priVtxId]++;
      }
      int maxOccurrenceId = -1;
      int maxOccurence = -1;
      for (auto it : fmap) {
        if (it.second > maxOccurence) {
          maxOccurence = it.second;
          maxOccurrenceId = it.first;
        }
      }

      // Match reco to truth vertex if at least 50% of tracks match
      if ((double)fmap[maxOccurrenceId] / tracks.size() >
          m_cfg.minTrackVtxMatchFraction) {
        for (const auto& particle : associatedTruthParticles) { // full loop of 
          int priVtxId = particle.particleId().vertexPrimary();
          int secVtxId = particle.particleId().vertexSecondary();

          if (secVtxId != 0) {
            // truthparticle from secondary vtx
            continue;
          }

          if (priVtxId == maxOccurrenceId) {
            // Vertex found, fill varibles
            const auto& truePos = particle.position();

            m_diffx.push_back(vtx.position()[0] - truePos[0]);
            m_diffy.push_back(vtx.position()[1] - truePos[1]);
            m_diffz.push_back(vtx.position()[2] - truePos[2]);
            // Next vertex now
            break;
          }
        }
      }
    }  // end loop vertices
  }

  // fill the variables
  m_outputTree->Fill();
  m_outputTree_Truth->Fill();

  m_diffx.clear();
  m_diffy.clear();
  m_diffz.clear();

  m_truth_vtx_vx.clear();
  m_truth_vtx_vy.clear();
  m_truth_vtx_vz.clear();

  m_truth_vtx_trk_d0.clear(); 
  m_truth_vtx_trk_z0.clear(); 
  m_truth_vtx_trk_phi.clear(); 
  m_truth_vtx_trk_theta.clear(); 
  m_truth_vtx_trk_qp.clear(); 
  m_truth_vtx_trk_time.clear(); 
  m_truth_vtx_trk_vtxID.clear(); 

  m_truth_vtx_trk_err_d0.clear(); 
  m_truth_vtx_trk_err_z0.clear(); 
  m_truth_vtx_trk_err_phi.clear(); 
  m_truth_vtx_trk_err_theta.clear(); 
  m_truth_vtx_trk_err_qp.clear(); 
  m_truth_vtx_trk_err_time.clear(); 

  return ProcessCode::SUCCESS;
}
