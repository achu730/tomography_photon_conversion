import pandas as pd
particle_dict = {
                 "KS0":   {"abbr": "Photon", "type": "composite", "daughters": ["e+", "e-"]},
                 "e+":       {"abbr": "ep", "type": "basic", "daughters":[]}, 
                 "e-":       {"abbr": "em", "type": "basic", "daughters":[]}, 
            }
particle_df = pd.DataFrame(particle_dict).T #puts each particle name as a column, abbr/type/daughters as an index, and transpose that with T
#basic is reconstructed directly from a track, composite must be expanded into a sub-decay in the descriptor and reconstructed from the daughters

def descriptor(particles, cc=False):
    def expand(particle, top=True):
        if particle_df["type"][particle] == "basic":
            return particle

        body = " ".join(
            expand(d, top=False)
            for d in particle_df["daughters"][particle]
        )

        s = f"{particle} -> {body}"
        return s if top else f"({s})"

    d = f"({expand(particles[0])})"

    if cc:
        d = f"[{d}]cc"

    return d

descriptor(["KS0"])

def decay_branches(mother, daughters, decay_descriptor):
    pos = {}  # per-target search cursor, so repeated particles get successive carets
    if "cc" in decay_descriptor[-5:]:
        branch_descriptor = decay_descriptor[:-2] + "CC"
    else:
        branch_descriptor = decay_descriptor


    abbr_to_par = {particle_df["abbr"][p]: p for p in particle_df.index}

    def mark(target, key):
        start = pos.get(key, 0)
        i = branch_descriptor.find(target, start)
        if i == -1:
          raise ValueError(f"Could not find {target} in {branch_descriptor}")
        pos[key] = i + 1
        return branch_descriptor[:i] + "^" + branch_descriptor[i:]

    def branch_for(par):
        if particle_df["type"][par] == "composite":
            return mark(f"({par}", par)   # find "(KS0", caret before the "(" -> ^(KS0 -> ...)
        return mark(par, par)             # tracks: find "e+"/"e-" directly

    branches = {}
    mpar = abbr_to_par.get(mother)
    branches[mother] = branch_for(mpar) if mpar is not None else branch_descriptor

    for daughter in daughters:
        head = daughter.split("_")[0]     # "ep_Photon" -> "ep", "Photon" -> "Photon"
        par = abbr_to_par.get(head)
        if par is not None:
            branches[daughter] = branch_for(par)
    return branches
    
#{
#  "Photon":      "^(KS0 -> e+ e-)",
#  "ep-Photon": " (KS0 -> ^e+ e-)",
#  "em-Photon": "(KS0 -> e+ ^e-)",
#}

def all_particles(particles):
    all_pars = []               # start with an empty result list
    for particle in particles:  # for each particle you passed in
        all_pars += ( [ particle ] + particle_df["daughters"][particle]) # particle + its daughter particles, added to the string
    return all_pars
# e.g. if you pass  [ "KS0"], you will get   [ "KS0", "e+", "e-" ]

def default_names(particles):
    names = []
    for particle in particles:
        abbr = particle_df["abbr"][ particle ] # this particle's abbreviation 
        names.append( abbr )
        names += [f"{daughter_abbr}_{abbr}" for daughter_abbr in particle_df["abbr"][particle_df["daughters"][particle]]]
    return names

def default_branches(particles, cc=False):
    names = default_names(particles)   #builds the abbreviatede name list
    return decay_branches(names[0], names[1:], descriptor(particles, cc)) #mother, daughters, descriptor