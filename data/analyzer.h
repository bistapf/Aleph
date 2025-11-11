#ifndef SELECTIONUTILS_H
#define SELECTIONUTILS_H

/*
  Selection utilities for filtering particles and events.
  Fully RDataFrame compatible (using ROOT::VecOps::RVec).

  Includes:
    - sel_charged: selects reconstructed particles by absolute charge.
    - sel_class_filter: filters events based on their class bit.
    - sel_runs_filter: filters events by allowed run numbers.
    - get_isEl / get_isMu / get_isChargedHad / get_isNeutralHad / get_isGamma:
      classify jet constituents by particle type.

  Example usage in RDataFrame:

    df = df.Define("charged_particles",
                   "FCCAnalyses::AlephSelection::sel_charged(1)(ReconstructedParticles)")
           .Filter("FCCAnalyses::AlephSelection::sel_class_filter(16)(ClassBitset)")
           .Filter("FCCAnalyses::AlephSelection::sel_runs_filter(allowedRuns)(EventHeader)");

    df = df.Define("isMu", "FCCAnalyses::AlephSelection::get_isMu(JetConstituents)")
           .Define("n_muons_per_jet", "Sum(isMu)");
*/
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include <set>
#include <bitset>
#include <cmath>
#include <ROOT/RVec.hxx>

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"
#include "FCCAnalyses/ReconstructedParticle2MC.h"
#include "FCCAnalyses/JetClusteringUtils.h"
// #include "FCCAnalyses/ExternalRecombiner.h"
#include "FCCAnalyses/MCParticle.h"

#include "edm4hep/MCParticleData.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterData.h"
#include "edm4hep/CalorimeterHitData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/EDM4hepVersion.h"
#include "edm4hep/RecDqdx.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"

#include <iostream>
#include <algorithm>



namespace FCCAnalyses { namespace AlephSelection {

namespace rv = ROOT::VecOps;

// -----------------------------------
// Type aliases for jet constituent data
// -----------------------------------

using FCCAnalysesJetConstituents = rv::RVec<edm4hep::ReconstructedParticleData>;
using FCCAnalysesJetConstituentsData = rv::RVec<float>;

// ----------------------------
// Event and particle selectors
// ----------------------------

//////////////////////////////////////////////////////////////////////////////////////////
// -----------------------------------
// The following will first filter the 
// event to check if it's a qq event
// then it return the PID of qq 
// -----------------------------------

float getJetPID(const ROOT::VecOps::RVec<uint32_t>& ClassBit,
                   const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles) {
    // Check if bit 15 (16-1) is set in the first element of ClassBit
    if (ClassBit.empty() || !std::bitset<32>(ClassBit[0])[15]) {
        return -1.0f; // Not selected
    }

    // Bit 15 is true — find the first quark (|PDG| in 1..5)
    float result = -1.0f;
    for (const auto& particle : particles) {
        if (std::abs(particle.PDG) > 0 && std::abs(particle.PDG) < 6) {
            result = static_cast<float>(std::abs(particle.PDG));
            break;
        }
    }

    return result;
}

/// Selects charged particles based on their absolute charge.
struct sel_charged {
  const int m_charge;
  sel_charged(int arg_charge) : m_charge(arg_charge) {};

  edm4hep::ReconstructedParticleCollection
  operator()(const edm4hep::ReconstructedParticleCollection& in_coll) const {
    edm4hep::ReconstructedParticleCollection result;
    result.setSubsetCollection();

    for (const auto& i : in_coll) {
      if (std::abs(i.getCharge()) == m_charge) {
        result.push_back(i);
      }
    }
    return result;
  }
};

/// Filters events based on their class bit (RVec-compatible)
struct sel_class_filter {
  const int m_class;
  sel_class_filter(int arg_class) : m_class(arg_class) {};

  bool operator()(const ROOT::VecOps::RVec<uint32_t>& bitset_coll) const {
    if (bitset_coll.empty()) return false;
    std::bitset<32> bits(bitset_coll[0]);
    return bits[m_class - 1];
  }
};

/// Filters events by run number (RVec-compatible)
struct sel_runs_filter {
  const std::set<int>& m_runs_set;
  sel_runs_filter(const std::set<int>& arg_runs_set) : m_runs_set(arg_runs_set) {};

  bool operator()(const ROOT::VecOps::RVec<edm4hep::EventHeader>& event_header) const {
    if (event_header.empty()) return false;
    return m_runs_set.count(event_header[0].getRunNumber()) > 0;
  }
};

// --------------------------------------
// Jet constituent particle identification
// --------------------------------------

rv::RVec<FCCAnalysesJetConstituentsData>
get_isEl(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
      mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.000511) < 1e-5) ? 1.f : 0.f);
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isMu(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
      mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.105658) < 1e-3) ? 1.f : 0.f);
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isChargedHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
      mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.13957) < 1e-3) ? 1.f : 0.f);
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isNeutralHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
      mask.push_back((c.PDG == 130) ? 1.f : 0.f);
#else
      mask.push_back((c.type == 130) ? 1.f : 0.f);
#endif
    out.push_back(std::move(mask));
  }
  return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isGamma(const rv::RVec<FCCAnalysesJetConstituents>& jcs) {
  rv::RVec<FCCAnalysesJetConstituentsData> out;
  out.reserve(jcs.size());
  for (const auto& jet : jcs) {
    FCCAnalysesJetConstituentsData mask;
    mask.reserve(jet.size());
    for (const auto& c : jet)
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
      mask.push_back((c.PDG == 22) ? 1.f : 0.f);
#else
      mask.push_back((c.type == 22) ? 1.f : 0.f);
#endif
    out.push_back(std::move(mask));
  }
  return out;
}
// --------------------------------------
// Event primary vertex reconstruction
// --------------------------------------

struct get_EventPrimaryVertexP4 {
  int m_genstatus = 21; // default generator status for incoming hard subprocess
  get_EventPrimaryVertexP4() {};

  TLorentzVector operator()(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in) const {
    TLorentzVector result(-1e12, -1e12, -1e12, -1e12);
    bool found_py8 = false;

    // First, look for generatorStatus == m_genstatus (e.g., 21 for Pythia8 hard process incoming)
    for (const auto& p : in) {
      //if (p.generatorStatus == m_genstatus) {
      if (1) {
        // vertex.time is in seconds, convert to mm
        TLorentzVector res(p.vertex.x,
                           p.vertex.y,
                           p.vertex.z,
                           p.time * 1.0e3 * 2.99792458e+8);
        result = res;
        found_py8 = true;
        break;
      }
    }

    // Fallback: look for genStatus == 2 with non-zero z vertex
    if (!found_py8) {
      for (const auto& p : in) {
        if (p.generatorStatus == 2 && std::abs(p.vertex.z) > 1.e-12) {
          TLorentzVector res(p.vertex.x,
                             p.vertex.y,
                             p.vertex.z,
                             p.time * 1.0e3 * 2.99792458e+8);
          result = res;
          break;
        }
      }
    }

    return result;
  }
};

float get_EventType(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& in) {
    float result = -1;
    // Look for the first particle with non-zero type

    //std::cout<<"DEBUG: get_EventType called with "<<in.size()<<" particles\n"<<std::endl;

    // for (const auto& p : in) {
    //   std::cout<<"DEBUG: Particle PDG="<<p.PDG<<" px="<<p.momentum.x<<" py="<<p.momentum.y<<" pz="<<p.momentum.z<<" status="<<p.generatorStatus<<std::endl;
    // }

    // first particle encountered with pdg > 0 and < 6
    for (const auto& p : in) {
      if (std::abs(p.PDG) > 0 && std::abs(p.PDG) < 6) {
        result = static_cast<float>(abs(p.PDG));
        break;
      }
    }

    // std::cout<<""<<result<<"\n"<<std::endl;
    return result;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isType(const rv::RVec<FCCAnalysesJetConstituentsData>& jcs, float type) {
    rv::RVec<FCCAnalysesJetConstituentsData> out;
    out.reserve(jcs.size());

    for (const auto& jet : jcs) {
        FCCAnalysesJetConstituentsData mask;
        mask.reserve(jet.size());

        for (const auto& c : jet) {
            if (c == type)
                mask.push_back(1);
            else
                mask.push_back(0);
        }

        out.push_back(std::move(mask));
    }

    return out;
}

struct build_constituents_Types {
    rv::RVec<FCCAnalysesJetConstituentsData>
    operator()(const rv::RVec<edm4hep::ParticleIDData> &rpid,
               const std::vector<std::vector<int>> &indices) const
    { 
        rv::RVec<FCCAnalysesJetConstituentsData> jcs;
        for (const auto &jet_index : indices)
        {
            FCCAnalysesJetConstituentsData jc;
            for (const auto &const_index : jet_index)
            {
                jc.push_back(rpid.at(const_index).type);
            }
            jcs.push_back(jc);
        }
        return jcs;
    }
};

// struct build_constituents_Types {
//   // Make the operator static to match usage style
//   static rv::RVec<FCCAnalysesJetConstituentsData>
//   operator()(const rv::RVec<edm4hep::ParticleIDData> &rpid,
//              const std::vector<std::vector<int>> &indices) 
//   {
//     rv::RVec<FCCAnalysesJetConstituentsData> jcs;
//     for (const auto &jet_index : indices)
//     {
//       FCCAnalysesJetConstituentsData jc;
//       for (const auto &const_index : jet_index)
//       {
//         jc.push_back(rpid.at(const_index).type);
//       }
//       jcs.push_back(jc);
//     }
//     return jcs;
//   }
// };

// rv::RVec<FCCAnalysesJetConstituentsData>  build_constituents_Types(const rv::RVec<edm4hep::ParticleIDData> &rpid, const std::vector<std::vector<int>> &indices) 
// { 
//   rv::RVec<FCCAnalysesJetConstituentsData> jcs;
//    for (const auto &jet_index : indices)
//            { FCCAnalysesJetConstituentsData jc;
//               for (const auto &const_index : jet_index) 
//                   { 
//                     jc.push_back(rpid.at(const_index).type);
//                   } 
     
     
//      jcs.push_back(jc); } return jcs; }


struct build_constituents_dEdx{
    rv::RVec<rv::RVec<edm4hep::RecDqdxData>>
    operator()(const rv::RVec<edm4hep::ReconstructedParticleData> &recoParticles,
             const rv::RVec<int> &_recoParticlesIndices,
             const rv::RVec<edm4hep::RecDqdxData> &dEdxCollection,
             const rv::RVec<int> &_dEdxIndicesCollection, 
             const std::vector<std::vector<int>> &jet_indices) const
    { 
        rv::RVec<rv::RVec<edm4hep::RecDqdxData>> dedx_constituents;

        // The links dEdx -> Track and RecoPart -> Track are one-directional, we need a map to store
        // Track.index -> dEdx to not have to loop everytime 
        // in addition, the object itself is stored in <Collection>
        // while the relations (=indices we need for links) are on _<Collection>
        std::unordered_map<int, edm4hep::RecDqdxData> track_index_to_dEdx;
        for (size_t i = 0; i < _dEdxIndicesCollection.size(); ++i) {
          int track_index = _dEdxIndicesCollection[i];
          edm4hep::RecDqdxData dedx_value = dEdxCollection[i];
          track_index_to_dEdx[track_index] = dedx_value;
        }

        //now, for each jet loop over the indices of the jet constituents provided by teh JetClusteringUtils
        // retrieve the associated RecoParticle
        // from there get the link to the Track from the corresponding index collection
        for (const auto &jet_const_indices : jet_indices) { //loop over jets
          rv::RVec<edm4hep::RecDqdxData> jet_dEdx;

          for (int constituent_index : jet_const_indices) { // loop over jet constituents
            const auto &recoPart = recoParticles[constituent_index];

            //loop over tracks associated to the RecoPart (should always be one in Aleph data)
            for (int track = recoPart.tracks_begin; track < recoPart.tracks_end; ++track) {
                 int track_index = _recoParticlesIndices[track]; //this should be the same index used in the link from dEdx to track

                  //find the matching dEdx in the map
                  if (track_index_to_dEdx.count(track_index)) {
                    jet_dEdx.push_back(track_index_to_dEdx[track_index]);
                  }
            }
          }
          dedx_constituents.push_back(jet_dEdx); 
        }
        return dedx_constituents;
    }
};

//helpers to read the dEdx objects (to check if can reuse existing FCCAna functions instead):

rv::RVec<rv::RVec<float>> get_dEdx_type(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.type);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_dEdx_value(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.value);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_dEdx_error(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.error);
    }
    values.push_back(inner_values);
  }
  return values;
}





}} // namespace FCCAnalyses::AlephSelection

#endif

