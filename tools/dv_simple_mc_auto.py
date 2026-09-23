
import os
import sys

from ..tupling_maker_MC import *
from GaudiKernel.SystemOfUnits import MeV
from RecoConf.algorithms_thor import ParticleFilter, ParticleCombiner, ParticleContainersMerger
from DaVinciMCTools import MCTruthAndBkgCat
from PyConf.reading import get_pvs, get_odin, get_rec_summary
import Functors as F
from FunTuple import FunctorCollection
from FunTuple import FunTuple_Particles as Funtuple
from RecoConf.event_filters import require_pvs
from RecoConf.reconstruction_objects import make_pvs as _make_pvs
from RecoConf.standard_particles import standard_protoparticle_filter 
from PyConf.Algorithms import FunctionalParticleMaker
from RecoConf.standard_particles import get_all_track_selector
from RecoConf.mc_truth_matching import custom_truth_match_candidates
from RecoConf.protoparticles import make_charged_protoparticles
from PyConf.reading import reconstruction
from RecoConf.standard_particles import _make_particles
from RecoConf.standard_particles import  make_down_electrons , make_down_electrons_no_brem
_basic = "basic"
_composite = "composite"

N_BODY = 2

###########################################################
# Part 1: Tupling helpers (2-body only)
###########################################################

def generate_decay_descriptor(charge_combo_str):
    particle_list = ["e+" if c == "p" else "e-" for c in charge_combo_str]
    particle_list.sort(reverse=True)
    particle_string = " ".join(particle_list)

    combiner_descriptor = f"KS0 -> {particle_string}"

    fields = {
        "Photon": f"^(KS0 -> {particle_string})",
    }

    for i, word in enumerate(["One", "Two"]):
        marked = particle_list.copy()
        marked[i] = f"^{marked[i]}"
        fields[f"Track{word}"] = f"(KS0 -> {' '.join(marked)})"

    return combiner_descriptor, fields


def tupling_template(mytuple_name, decay_descriptor, photon_fields,
                     combination_cuts, composite_cuts, input_particles):
    

    KS0 = ParticleCombiner(        #gamma
        Inputs=[input_particles, input_particles],
        DecayDescriptor=decay_descriptor,
        name="make_vertex_{hash}",
        CombinationCut=combination_cuts,
        CompositeCut=composite_cuts,
    )

    MC_TRUTH = MCTruthAndBkgCat(
        KS0,
        name=f"MCTruthAndBkgCat_{mytuple_name}_{{hash}}"
    )


    truth_info_track = FunctorCollection({
        "ORIGIN_VTX_TYPE":  F.VALUE_OR(-1) @ MC_TRUTH(F.MC_VTX_TYPE @ F.MC_ORIGINVERTEX), #gives -1 if no match
        "TRUEKEY":          F.VALUE_OR(-1) @ MC_TRUTH(F.OBJECT_KEY),
        "TRUEFOURMOMENTUM": MC_TRUTH(F.FOURMOMENTUM),
        "TRUEORIGIN_VX":    MC_TRUTH(F.ORIGIN_VX),
        "TRUEORIGIN_VY":    MC_TRUTH(F.ORIGIN_VY),
        "TRUEORIGIN_VZ":    MC_TRUTH(F.ORIGIN_VZ),
    })
    truth_info_vertex = FunctorCollection({
        "TRUEKEY":          F.VALUE_OR(-1) @ MC_TRUTH(F.OBJECT_KEY),
        "TRUEM":            MC_TRUTH(F.MASS),
        "TRUEFOURMOMENTUM": MC_TRUTH(F.FOURMOMENTUM),
    })

    hlt1_lines = ['Hlt1GECPassThroughDecision']
    hlt2_lines = ['Hlt2QEE_DiElectronPrompt_PersistPhotons_Full', 'Hlt2QEE_DiElectronDisplaced_PersistPhotons_Full'] #they only check if the candidates passed, they don't select the candidates here
    pvs, odin, rec_sum = get_pvs(), get_odin(), get_rec_summary()

    rec_keys = [
        "nPVs", "nTracks", "nLongTracks", "nDownstreamTracks", "nUpstreamTracks",
        "nVeloTracks", "nBackTracks", "nGhosts", "nRich1Hits", "nRich2Hits",
        "nVPClusters", "nUTClusters", "nFTClusters", "nSPDhits",
        "eCalTot", "hCalTot", "nEcalClusters",
        "nMuonCoordsS0", "nMuonCoordsS1", "nMuonCoordsS2",
        "nMuonCoordsS3", "nMuonCoordsS4", "nMuonTracks",
    ] #list of variables to extract from the RecSummary, which contains summary information about the reconstructed event, such as the number of primary vertices, tracks, hits in different subdetectors, etc. 
    rec_info = {k: F.VALUE_OR(-1) @ F.RECSUMMARY_INFO(rec_sum, k) for k in rec_keys}
    rec_info.update({"ALLPVX": F.ALLPVX(pvs), "ALLPVY": F.ALLPVY(pvs), "ALLPVZ": F.ALLPVZ(pvs)})

    evt_vars = event_variables(pvs, odin, rec_sum, hlt1_lines, hlt2_lines)\
               + FunctorCollection(rec_info)

    variables = {
        "Photon": all_variables(pvs, MC_TRUTH, _composite) + truth_info_vertex,
        "TrackOne": all_variables(pvs, MC_TRUTH, _basic) + truth_info_track,
        "TrackTwo": all_variables(pvs, MC_TRUTH, _basic) + truth_info_track,
    }

    return Funtuple(
        name=mytuple_name, tuple_name="DecayTree",
        fields=photon_fields, variables=variables,
        event_variables=evt_vars, store_multiple_cand_info=True,
        inputs=KS0,
    )


###########################################################
# Part 2: Main Configuration
###########################################################

sys.path.append(os.path.join(os.environ['ANALYSIS_PRODUCTIONS_BASE'], 'achu_tomography_pair_downstream'))
from DaVinci import Options, make_config


def main(options: Options):
    composite_cuts =  F.require_all(
        F.CHI2DOF < 100 , 
        F.BPVIPCHI2() <5000
    )
    combination_cuts = F.require_all(
        F.MASS < 300 * MeV 
    )



    preselection = F.ALL

    loose_electrons = ParticleFilter(
        make_down_electrons_no_brem(),
        F.FILTER(F.ALL),
    )

    filtered_electrons = ParticleFilter(loose_electrons, F.FILTER(preselection))

    user_algorithms = {}
    for charge_combo in ['pp', 'pm', 'mm']:  #builds one algorithm per charge combination, with the corresponding decay descriptor and variable definitions
        decay_descriptor, photon_fields = generate_decay_descriptor(charge_combo)
        mytuple_name = f"Photon_{charge_combo}"

        my_tuple = tupling_template(
            mytuple_name=mytuple_name,
            decay_descriptor=decay_descriptor,
            combination_cuts=combination_cuts,
            composite_cuts=composite_cuts,
            input_particles=filtered_electrons,
            photon_fields=photon_fields
        )

        alg_name = f'Alg_tuple_{charge_combo}'
        user_algorithms[alg_name] = [require_pvs(_make_pvs()), my_tuple]

    print("\nGenerated the following algorithms:")
    for alg_name in user_algorithms:
        print(f"- {alg_name}")

    return make_config(options, user_algorithms)