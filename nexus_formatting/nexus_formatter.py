#this script translates target site information from Cassiopeia to a NEXUS format
from pickle import FALSE
import pandas as pd 
import numpy as np
import os 
import sys
import joblib
from joblib import Parallel, delayed
pd.options.display.max_colwidth = None

def df_to_nexus(dataf,output_file,no_chains,gamma,pf,it):
    ntax = dataf.shape[0]
    nchar = len(dataf["sequence"][0]) 
    with open(output_file, 'w') as f:
        f.write("#NEXUS\n")
        f.write("[ TITLE ]\n\n")
        f.write("BEGIN"+' '+"data;\n")
        f.write("\tDIMENSIONS"+" "+"NTAX="+str(ntax)+" "+"NCHAR="+str(nchar)+";\n")
        f.write("\tFORMAT"+" "+"MISSING=?"+" "+"GAP=-"+" "+"MATCHCHAR=."+" "+"datatype=nucleotide;\n")
        f.write("MATRIX\n\n")
        for i in range(0,ntax):
            f.write(dataf["combined"][i])
            f.write("\n")
        f.write(";\nEND;\n\n")
        f.write("Begin BayesPhylogenies;\n\n")
        f.write("\tautorun\n\n")  
        f.write("\tChains"+" "+str(no_chains)+";\n")
        f.write("\tCreatePart"+" "+"Def"+" "+"GTNR;\n")
        f.write("\thelp\n\n")  
        f.write("\tgamma"+" "+str(gamma)+"\n")
        f.write("\tSample"+" "+str(pf)+"\n")
        f.write("\titerations"+" "+str(it)+"\n")
        f.write("\tRoot RootSeq;\n")
        f.write("\trjpatterns\n\n")
        f.write("\tDebug;\n\n")
        f.write("end;")

def seq_format(df,lin_dir,lin_grp):
    seq_length = len(root_seq)
    intbc_list = df["intBC"].unique()

    cell_BC = []
    seq_l = []

    cell_BC.append("RootSeq")
    seq_l.append(root_seq * len(intbc_list))

    for i in df["cellBC"].unique():
        tmp = df[df["cellBC"] == i]
        tot_seq = ""
        for j in intbc_list:
            if (j in tmp["intBC"].unique()) == True:
                tmp_seq = tmp[tmp["intBC"] == j]["Seq"].values[0]
                tmp_seq = tmp_seq[intbc_cutoff:]
                tot_seq = tot_seq + tmp_seq
            if (j in tmp["intBC"].unique()) == False:
                tmp_seq = "?" * seq_length
                tot_seq = tot_seq + tmp_seq
        cell_BC.append(str(i))
        seq_l.append(tot_seq)

    datal  = list(zip(cell_BC,seq_l))
    df1 = pd.DataFrame(data= datal)
    df1.columns = ["cellBC","sequence"]
    df1.to_csv(lin_dir + "/lineagegrp_" + str(lin_grp) + "concate_seq.csv",index=False)
    return df1


def lineage_main(lin_id):
    lgrp = lin_id
    lineage_dir = output_dir + "/lineagegrp_" + str(lgrp) + "_tree"
    if os.path.exists(lineage_dir) == True:
        print("Lineage dir exists")
        sys.exit()
    if os.path.exists(lineage_dir) == False:
        os.mkdir(lineage_dir)
    
    input_data = data[data["LineageGroup_"+lineage_type] == lgrp]

    df1 = seq_format(input_data,lineage_dir,lgrp)

    spec_no_chains = 1
    spec_gamma = 4
    spec_pf = 100000
    spec_it = 200000000

    df1["combined"] = df1["cellBC"] +" "+ df1["sequence"]
    output_file = lineage_dir + "/lineagegrp_" + str(lgrp) + ".nexus"
    df_to_nexus(df1,output_file,spec_no_chains,spec_gamma,spec_pf,spec_it)

data = pd.read_csv("~/final_data_filtered_min10_avg2_seq.txt",sep='\t',low_memory=False)

#this only considers metastatic trees and non-metastatic trees within an animal in which a metastatic tree has been recovered
df = data.loc[data['LineageGroup_mets'] < 14,]
data = data.loc[data.animal_id.isin(df['animal_id'])]

#remove the intbc and all proceeding nucleotides
intbc_cutoff = 34
lineage_type = "mets"

output_dir = "~/trees_"+lineage_type


root_seq = open('~/root_seq.txt', "r").readline()
root_seq = root_seq[intbc_cutoff:]

lineage_groups = sorted(data["LineageGroup_"+lineage_type].unique())

num_cores = joblib.cpu_count()
Parallel(n_jobs =12, backend="loky", verbose=5)(delayed(lineage_main)(i) for i in lineage_groups)





