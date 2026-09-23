import FunTuple.functorcollections as FC
import Functors as F
from FunTuple import FunctorCollection
from .tools.descriptor_writer import *
from DecayTreeFitter import DecayTreeFitter

_basic = "basic"
_composite = "composite"
_toplevel = "toplevel"

def all_variables(pvs, mctruth, ptype, candidates=None, ftAlg=None):

    if ptype not in [_basic, _composite]:
        Exception(f"I want {_basic} or {_composite}. Got {ptype}")
    all_vars = FunctorCollection({})

    comp = _composite == ptype or _toplevel == ptype  # is composite
    basic = _basic == ptype  # is not composite
    top = _toplevel == ptype  # the B
    
    ##top: B baryon in LHCb physics analysis
    ##
    ##comp: the composite particle in the decay channel. 
    ##      For example, KS0 and Jpsi in B2JpsiKS0 channel.
    ## 
    ##basic: some basic particles, including proton, keon, pion, electron(positron) and muon.  

    all_vars += FC.Kinematics()
    if basic:
        all_vars += FC.ParticleID(extra_info=True)
    all_vars += FC.MCHierarchy(mctruth_alg=mctruth)

    all_vars.update({"BPVIP": F.BPVIP(pvs)})
    all_vars.update({"BPVIPCHI2": F.BPVIPCHI2(pvs)})
    all_vars.update({"BPVX": F.BPVX(pvs)})
    all_vars.update({"BPVY": F.BPVY(pvs)})
    all_vars.update({"BPVZ": F.BPVZ(pvs)})

    if comp: 
        all_vars.update({"ALLPV_FD": F.ALLPV_FD(pvs)})
        all_vars.update({"ALLPV_IP": F.ALLPV_IP(pvs)})
        all_vars.update({"BPVVDX": F.BPVVDX(pvs)})
        all_vars.update({"BPVVDY": F.BPVVDY(pvs)})
        all_vars.update({"BPVVDZ": F.BPVVDZ(pvs)})

    all_vars.update({"CHI2": F.CHI2})
    all_vars.update({"CHI2DOF": F.CHI2DOF})

    if comp:
        all_vars.update({"END_VX": F.END_VX})
        all_vars.update({"END_VY": F.END_VY})
        all_vars.update({"END_VZ": F.END_VZ})
        all_vars.update({"OWNPVDLS": F.OWNPVDLS})

    all_vars.update({"ETA": F.ETA})

    if basic:
        all_vars.update({"GHOSTPROB": F.GHOSTPROB})
        all_vars.update({"ISMUON": F.ISMUON})
        all_vars.update({"INMUON": F.INMUON})
        all_vars.update({"INECAL": F.INECAL})
        all_vars.update({"HASBREM": F.HASBREM})
        all_vars.update({"BREMENERGY": F.BREMENERGY})
        all_vars.update({"BREMBENDCORR": F.BREMBENDCORR})
        all_vars.update({"CLUSTERMATCH": F.CLUSTERMATCH_CHI2})
        all_vars.update({"ELECTRONMATCH": F.ELECTRONMATCH_CHI2})

    if comp:
        all_vars.update({"MAXPT": F.MAX(F.PT)})
        all_vars.update({"MAXDOCA": F.MAXDOCA})
        all_vars.update({"MAXDOCACHI2": F.MAXDOCACHI2})
        all_vars.update({"MAXSDOCA": F.MAXSDOCA})
        all_vars.update({"MAXSDOCACHI2": F.MAXSDOCACHI2})

    if comp:
        all_vars.update({"MINPT": F.MIN(F.PT)})
    all_vars.update({"MINIP": F.MINIP(pvs)})
    all_vars.update({"MINIPCHI2": F.MINIPCHI2(pvs)})

    if basic:
        all_vars.update({"TRACKTYPE": F.VALUE_OR(-1) @ F.TRACKTYPE @ F.TRACK})
        all_vars.update({"TRACKHISTORY": F.VALUE_OR(-1) @ F.TRACKHISTORY @ F.TRACK})
        all_vars.update({"NDOF": F.VALUE_OR(-1) @ F.NDOF @ F.TRACK})
        all_vars.update({"NFTHITS": F.VALUE_OR(-1) @ F.NFTHITS @ F.TRACK})
        all_vars.update({"NHITS": F.VALUE_OR(-1) @ F.NHITS @ F.TRACK})
        all_vars.update({"NUTHITS": F.VALUE_OR(-1) @ F.NUTHITS @ F.TRACK})
        all_vars.update({"NVPHITS": F.VALUE_OR(-1) @ F.NVPHITS @ F.TRACK})

    all_vars.update({"OBJECT_KEY": F.OBJECT_KEY})
    all_vars.update({"PHI": F.PHI})
    return all_vars                                                                                                                      

def tistos_variables(Hlt1_decisions, Hlt2_decisions, data, isturbo, isHlt1=True):
    tistos_vars = FunctorCollection({})
    if isHlt1:
        tistos_vars += FC.HltTisTos( selection_type="Hlt1", trigger_lines=Hlt1_decisions, data=data)
    if not isturbo:
        tistos_vars += FC.HltTisTos( selection_type="Hlt2", trigger_lines=Hlt2_decisions, data=data)
    return tistos_vars

def event_variables(PVs, ODIN, rec_sum, hlt1_lines, sprucing_lines):
    evt_vars = FunctorCollection({
        "nPVs": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nPVs"),
        "nTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nTracks"),
        "nLongTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nLongTracks"),
        "nDownstreamTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nDownstreamTracks"),
        "nUpstreamTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nUpstreamTracks"),
        "nVeloTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nVeloTracks"),
        "nBackTracks": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nBackTracks"),
        "nGhosts": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nGhosts"),
        "nRich1Hits": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nRich1Hits"),
        "nRich2Hits": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nRich2Hits"),
        "nVPClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nVPClusters"),
        "nUTClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nUTClusters"),
        "nFTClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nFTClusters"),
        "eCalTot": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"eCalTot"),
        "hCalTot": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"hCalTot"),
        "nEcalClusters": F.VALUE_OR(-1) @F.RECSUMMARY_INFO(rec_sum,"nEcalClusters"),
        "ALLPVX": F.ALLPVX(PVs),
        "ALLPVY": F.ALLPVY(PVs),
        "ALLPVZ": F.ALLPVZ(PVs),
        })
    evt_vars += FC.EventInfo()
    evt_vars += FC.SelectionInfo( selection_type="Hlt1", trigger_lines=hlt1_lines)                                               
    evt_vars += FC.SelectionInfo( selection_type="Hlt2", trigger_lines=sprucing_lines)
    if ODIN:
        evt_vars.update({"EVENTTYPE": F.EVENTTYPE(ODIN)})
    evt_vars.update({"PV_SIZE": F.SIZE(PVs)})
    return evt_vars

def candidate_variables(pvs, particles, mctruth):

    abbrs = default_names(particles)
    names = all_particles(particles)
    variables_B = {abbr: all_variables(pvs, mctruth, particle_df["type"][name]) for abbr, name in zip(abbrs, names)}
    return variables_B

from DecayTreeFitter import DecayTreeFitter
def make_dtf_variables(pvs, data, particles, pv_constraint=False, mass_constraints=[], label_forLLDD=""):

    abbrs = default_names(particles)
    names = all_particles(particles)

    if pv_constraint: dtf_name = "DTF_PV_"
    else: dtf_name = "DTF_"
    dtf_name += "".join(f"{par}_" for par in particle_df["abbr"][mass_constraints])
    my_hash = '{hash}'

    DTF = DecayTreeFitter(
            name = f"{dtf_name}DecayTreeFitter_{my_hash}",
            input_particles = data,
            input_pvs = (pv_constraint and pvs),
            mass_constraints = mass_constraints,
            )

    dtf_dict = {}

    shared_variables = FunctorCollection(
                {   "ETA": F.ETA,
                    "PHI": F.PHI,
                    "BPVIPCHI2": F.BPVIPCHI2(pvs),
                    "BPVIP": F.BPVIP(pvs),
                    "BPVX": F.BPVX(pvs),
                    "BPVY": F.BPVY(pvs),
                    "BPVZ": F.BPVZ(pvs),
                    }
                ) + FC.Kinematics()

    dtf_quality_variables = FunctorCollection( {
            dtf_name+"_DTFCHI2": DTF.CHI2,
            dtf_name+"_DTFNDOF": DTF.NDOF,
            dtf_name+"_CTAU": DTF.CTAU,
            dtf_name+"_CTAUERR": DTF.CTAUERR,
            dtf_name+"_MERR": DTF.MASSERR,
            }
        )

    #make branches
    for abbr, name in zip(abbrs, names):
        is_basic = particle_df["type"][name]
        if is_basic:
            orig_variables = shared_variables + FunctorCollection(
                    {   "TX"          : F.TX,
                        "TY"          : F.TY,
                        "MINIPCHI2"   : F.MINIPCHI2(pvs),
                        "MINIP"       : F.MINIP(pvs),
                        "KEY"         : F.VALUE_OR(-1) @ F.OBJECT_KEY @ F.TRACK,
                        "TRGHOSTPROB": F.GHOSTPROB,
                        "TRACKPT": F.TRACK_PT,
                        "TRACKHISTORY": F.VALUE_OR(-1) @ F.TRACKHISTORY @ F.TRACK,
                        "QOVERP": F.QOVERP @ F.TRACK,
                        "TRCHI2DOF": F.CHI2DOF @ F.TRACK,
                        "NDOF": F.VALUE_OR(-1) @ F.NDOF @ F.TRACK,
                        }
                    )
        else:
            orig_variables = shared_variables + FunctorCollection(
                {   "MAXPT": F.MAX(F.PT),
                    "MINPT": F.MIN(F.PT),
                    "SUMPT": F.SUM(F.PT),
                    "MAXP": F.MAX(F.P),
                    "MINP": F.MIN(F.P),
                    "BPVDIRA": F.BPVDIRA(pvs),
                    "CHI2DOF": F.CHI2DOF, #CHI2VXNDOF
                    "BPVFDCHI2": F.BPVFDCHI2(pvs),
                    "BPVFD": F.BPVFD(pvs),
                    "BPVVDRHO": F.BPVVDRHO(pvs),
                    "BPVVDZ": F.BPVVDZ(pvs),
                    "BPVLTIME": F.BPVLTIME(pvs),
                    "END_VX": F.END_VX, #END_
                    "END_VY": F.END_VY,
                    "END_VZ": F.END_VZ,
                    }
                )
            orig_variables += dtf_quality_variables

        dtf_variables = FunctorCollection({ dtf_name + expr: DTF(func) for expr, func in orig_variables.get_thor_functors().items() })
        dtf_dict.update( {abbr: dtf_variables} )

    return dtf_dict

import re
def convert_cut(string_cut):
    cuts = string_cut.split("&")
    paras = ["BPV", "PID_", "MASS", "CHARGE", "CHI2DOF", "END", "IS_ID", "ABS_", "MAX", "SUM", "GHOSTPROB", "CHILD"]
    values = []
    for cut in cuts:
        value = cut
        for para in paras:
            if para in cut:
                value = value.replace(para, f"F.{para}")
        value = re.sub(r"\bPT\b", "F.PT", value)
        value = re.sub(r"\bP\b", "F.P", value)
        value = re.sub(r"\bCHI2\b", "F.CHI2", value)
        values.append(value)
    functor_cut = ",".join(values)
    return eval(f"F.require_all({functor_cut})")