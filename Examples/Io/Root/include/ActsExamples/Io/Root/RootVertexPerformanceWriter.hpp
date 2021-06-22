// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <mutex>
#include <vector>

class TFile;
class TTree;

struct PV {
  double PV_loc[3];
  std::vector<int> track_ID;
};

namespace ActsExamples {

/// @class RootVertexPerformanceWriter
///
/// Writes out the number of reconstructed primary vertices along with
/// the number of primary vertices in detector acceptance as well as
/// reconstructable primary vertices after track fitting.
/// Additionally it matches the reco vertices to their truth vertices
/// and write out the difference in x,y and z position.
class RootVertexPerformanceWriter final
    : public WriterT<std::vector<Acts::Vertex<Acts::BoundTrackParameters>>> {
 public:
  struct Config {
    /// All input truth particle collection
    std::string inputAllTruthParticles;
    /// Selected input truth particle collection
    std::string inputSelectedTruthParticles;
    /// Truth particles associated to fitted tracks
    std::string inputAssociatedTruthParticles;
    /// All event fitted tracks
    std::string inputFittedTracks;
    /// Input vertex collection.
    std::string inputVertices;
    /// Input reconstruction time.
    std::string inputTime;
    /// Output directory.
    std::string outputDir;
    /// Output filename.
    std::string outputFilename = "vertexingperformance.root";
    /// Name of the output tree.
    std::string outputTreename = "vertextree";
    /// Name of the output tree of truth info.
    std::string outputTreename_Truth = "Truth_Vertex";
    /// Name of the output tree of reco info.
    std::string outputTreename_Reco = "Reco_Vertex";
    /// File access mode.
    std::string fileMode = "RECREATE";
    /// Common root file.
    TFile* rootFile = nullptr;
    /// Minimum fraction of tracks matched between truth
    /// and reco vertices to be matched for resolution plots
    double minTrackVtxMatchFraction = 0.5;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  RootVertexPerformanceWriter(const Config& cfg, Acts::Logging::Level lvl);
  ~RootVertexPerformanceWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  ProcessCode writeT(
      const AlgorithmContext& ctx,
      const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices)
      final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  TTree* m_outputTree_Truth{nullptr};  ///< The output tree of truth info
  TTree* m_outputTree_Reco{nullptr};  ///< The output tree of reco info 

  int m_eventId{0}; 

  /// For truth vtx information 
  std::vector<double> m_truth_vtx_vx;
  std::vector<double> m_truth_vtx_vy;
  std::vector<double> m_truth_vtx_vz;

  /// Event-unique particle identifier a.k.a barcode.
  std::vector<uint64_t> m_truth_particle_Id;
  /// Particle type a.k.a. PDG particle number
  std::vector<int32_t> m_truth_particle_Type;
  /// Production process type, i.e. what generated the particle.
  std::vector<uint32_t> m_truth_particle_process;
  /// Production position components in mm.
  std::vector<double> m_truth_particle_vx;
  std::vector<double> m_truth_particle_vy;
  std::vector<double> m_truth_particle_vz;
  // Production time in ns.
  std::vector<double> m_truth_particle_vt;
  /// Total momentum in GeV
  std::vector<double> m_truth_particle_p;
  /// Momentum components in GeV.
  std::vector<double> m_truth_particle_px;
  std::vector<double> m_truth_particle_py;
  std::vector<double> m_truth_particle_pz;
  /// Mass in GeV.
  std::vector<double> m_truth_particle_m;
  /// Charge in e.
  std::vector<double> m_truth_particle_q;
  // Derived kinematic quantities
  /// Direction pseudo-rapidity.
  std::vector<double> m_truth_particle_eta;
  /// Direction angle in the transverse plane.
  std::vector<double> m_truth_particle_phi;
  /// Transverse momentum in GeV.
  std::vector<double> m_truth_particle_pt;
  // Decoded particle identifier; see Barcode definition for details.
  std::vector<uint32_t> m_truth_particle_vertexPrimary;
  std::vector<uint32_t> m_truth_particle_vertexSecondary;
  std::vector<uint32_t> m_truth_particle_particle;
  std::vector<uint32_t> m_truth_particle_generation;
  std::vector<uint32_t> m_truth_particle_subParticle;


  /// The track parameter associated to the truth vtx 
  std::vector<double> m_truth_vtx_fitted_trk_d0;
  std::vector<double> m_truth_vtx_fitted_trk_z0;
  std::vector<double> m_truth_vtx_fitted_trk_phi;
  std::vector<double> m_truth_vtx_fitted_trk_theta;
  std::vector<double> m_truth_vtx_fitted_trk_qp;
  std::vector<double> m_truth_vtx_fitted_trk_time;

  std::vector<double> m_truth_vtx_fitted_trk_err_d0;
  std::vector<double> m_truth_vtx_fitted_trk_err_z0;
  std::vector<double> m_truth_vtx_fitted_trk_err_phi;
  std::vector<double> m_truth_vtx_fitted_trk_err_theta;
  std::vector<double> m_truth_vtx_fitted_trk_err_qp;
  std::vector<double> m_truth_vtx_fitted_trk_err_time;

  std::vector<int> m_truth_vtx_fitted_trk_vtxID;


  /// For reco vtx information 
  std::vector<double> m_reco_vtx_vx;
  std::vector<double> m_reco_vtx_vy;
  std::vector<double> m_reco_vtx_vz;

  std::vector<double> m_reco_vtx_fitquality_chiSquared;
  std::vector<double> m_reco_vtx_fitquality_nDoF;

  std::vector<double> m_reco_vtx_err_vx_vx;
  std::vector<double> m_reco_vtx_err_vx_vy;
  std::vector<double> m_reco_vtx_err_vx_vz;
  std::vector<double> m_reco_vtx_err_vy_vy;
  std::vector<double> m_reco_vtx_err_vy_vz;
  std::vector<double> m_reco_vtx_err_vz_vz;

  /// The track parameter associated to the reco vtx 
  std::vector<double> m_reco_vtx_fitted_trk_d0;
  std::vector<double> m_reco_vtx_fitted_trk_z0;
  std::vector<double> m_reco_vtx_fitted_trk_phi;
  std::vector<double> m_reco_vtx_fitted_trk_theta;
  std::vector<double> m_reco_vtx_fitted_trk_qp;
  std::vector<double> m_reco_vtx_fitted_trk_time; 

  std::vector<double> m_reco_vtx_fitted_trk_err_d0;
  std::vector<double> m_reco_vtx_fitted_trk_err_z0;
  std::vector<double> m_reco_vtx_fitted_trk_err_phi;
  std::vector<double> m_reco_vtx_fitted_trk_err_theta;
  std::vector<double> m_reco_vtx_fitted_trk_err_qp;
  std::vector<double> m_reco_vtx_fitted_trk_err_time;

  std::vector<double> m_reco_vtx_fitted_trk_chi2Track;
  std::vector<double> m_reco_vtx_fitted_trk_ndf;
  std::vector<double> m_reco_vtx_fitted_trk_vertexCompatibility;
  std::vector<double> m_reco_vtx_fitted_trk_trackWeight;

  std::vector<int> m_reco_vtx_fitted_trk_vtxID;

  std::vector<float>
      m_diffx;  ///< Difference in x positon between reco and true vtx
  std::vector<float>
      m_diffy;  ///< Difference in y positon between reco and true vtx
  std::vector<float>
      m_diffz;  ///< Difference in z positon between reco and true vtx

  int m_nrecoVtx = -1;           ///< Number of reconstructed vertices
  int m_ntrueVtx = -1;           ///< Number of true vertices
  int m_nVtxDetAcceptance = -1;  ///< Number of vertices in detector acceptance
  int m_nVtxReconstructable =
      -1;  ///< Max. number of reconstructable vertices (detector acceptance +
           ///< tracking efficiency)
  int m_timeMS = -1;  ///< Reconstruction time in ms

  int getNumberOfReconstructableVertices(
      const SimParticleContainer& collection) const;

  int getNumberOfTruePriVertices(const SimParticleContainer& collection) const;

  std::vector<PV> getTruthVerticesVec(
      const SimParticleContainer& collection);

  void writeTruthInfo(std::vector<PV> PV_list, 
  const SimParticleContainer& collection, 
  const TrackParametersContainer& inputFittedTracks);

  void writeRecoInfo(const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices);


};

}  // namespace ActsExamples
