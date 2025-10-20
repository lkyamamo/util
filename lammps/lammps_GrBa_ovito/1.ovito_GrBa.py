from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier, TimeAveragingModifier, CreateBondsModifier, BondAnalysisModifier
import numpy as np
import pandas as pd
import re

### CHANGE THESE PARAMETERS ###
frame_count = 1000
mapping = {1: "O", 2: "H"}
concentrations = {"O": 1/3, "H": 2/3}

gr_bins = 150
gr_cutoff = 12

ba_bins = 180
ba_cutoffs=  [[5.0, 1.15],
             [1.15, 1.75]]
###############################

DIFFRACTION_LENGTHS = {"Si": 4.1491/(10**5),
                       "O": 5.803/(10**5),
                       "H": 6.671/(10**5)
                       }

def create_table(tag, data):
    headers_name = []

    match tag:
        case "ba":
            table_name='bond-angle-distr[average]'
            headers_name.append("angle")
        case "gr":
            table_name='coordination-rdf[average]'
            headers_name.append("distance")

    table_data = data.tables[table_name].xy()
    headers_num = data.tables[table_name].y.component_names

    for name in headers_num:
        etypes = re.findall(r'\d+',name)
        etypes = [mapping[int(istring)] for istring in etypes]
        header_name = '-'.join(etypes)
        headers_name.append(header_name)

    dframe = pd.DataFrame(data=table_data,columns=headers_name)
    return dframe


filename = "../all_lammps.xyz"
pipeline = import_file(filename)

print("Number of MD frames:", pipeline.num_frames)
print("Number of MD frames used:", frame_count)

# partial gr 
gr_modifier = CoordinationAnalysisModifier(cutoff = gr_cutoff, number_of_bins = gr_bins, partial = True)
pipeline.modifiers.append(gr_modifier)
gr_averaging_modifier = TimeAveragingModifier(operate_on='table:coordination-rdf')
gr_averaging_modifier.interval = (pipeline.num_frames - frame_count - 1, pipeline.num_frames - 1)
pipeline.modifiers.append(gr_averaging_modifier)

# bond angle calculations
bond_modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
for i in range(len(ba_cutoffs)):
    for j in range(i,len(ba_cutoffs)):
        bond_modifier.set_pairwise_cutoff(type_a=i+1, type_b=j+1, cutoff=ba_cutoffs[i][j])

pipeline.modifiers.append(bond_modifier)
ba_modifier = BondAnalysisModifier(bins=ba_bins, partition=BondAnalysisModifier.Partition.ByParticleType)
pipeline.modifiers.append(ba_modifier)

ba_averaging_modifier = TimeAveragingModifier(operate_on='table:bond-angle-distr')
ba_averaging_modifier.interval = (pipeline.num_frames - frame_count - 1, pipeline.num_frames - 1)
pipeline.modifiers.append(ba_averaging_modifier)

data = pipeline.compute()

gr_df = create_table("gr", data)
ba_df = create_table("ba", data)

# normalizing ba distribution
names = [name for name in ba_df.columns if name != "angle"]
prefactors = (ba_df[names]).sum(axis=0)
for name, factor in zip(names, prefactors):
    ba_df[name] = ba_df[name]/factor

# neutron gr
c = 0
for value in mapping.values():
    c += DIFFRACTION_LENGTHS[value]*concentrations[value]
c2_in = 1/(c*c)

neutron_gr = pd.DataFrame({"neutron_gr":[0 for i in range(len(gr_df))]})

for name in gr_df.columns:
    if name == "distance":
        continue
    elements = name.split("-")
    coefficient = 1
    for e in elements:
        coefficient *= DIFFRACTION_LENGTHS[e]*concentrations[e]

    neutron_gr["neutron_gr"] = neutron_gr["neutron_gr"] + coefficient*gr_df[name]

neutron_gr = neutron_gr * c2_in
gr_df = pd.concat([gr_df, neutron_gr], axis=1)

pipeline2 = import_file(filename)

# total gr 
total_gr_modifier = CoordinationAnalysisModifier(cutoff = gr_cutoff, number_of_bins = gr_bins, partial=False)
pipeline2.modifiers.append(total_gr_modifier)
total_gr_averaging_modifier = TimeAveragingModifier(operate_on='table:coordination-rdf')
gr_averaging_modifier.interval = (pipeline2.num_frames - frame_count - 1, pipeline.num_frames - 1)
pipeline2.modifiers.append(gr_averaging_modifier)

data = pipeline2.compute()

table_data = data.tables['coordination-rdf[average]'].y
total_gr_df = pd.DataFrame(data=table_data[:],columns=["total_gr"])

gr_df = pd.concat([gr_df,total_gr_df], axis=1)

print(gr_df)
print(ba_df)

gr_df.to_csv("gr.csv")
ba_df.to_csv("ba.csv")



