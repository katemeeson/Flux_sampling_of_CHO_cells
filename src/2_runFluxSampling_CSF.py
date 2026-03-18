import pandas as pd
import cobra
from re import sub
import numpy as np
from function_mapping_transcriptome_data_JW import map_transcriptome_data
from multiprocessing import Pool
import os
nslots = int(os.environ['NSLOTS'])
print(nslots)

transcriptomeData = pd.read_csv("DesEq2_counts_wDE_genes.csv",
                                    dtype = {'GeneID': 'str'})

model_path = "iCHO2441_221-107_producing.xml"
model = cobra.io.read_sbml_model(model_path)
model_genes = [x.id for x in model.genes]

transcriptomeData_sig = pd.concat([
    transcriptomeData["GeneID"],
    transcriptomeData.filter(regex=r"^B", axis=1)],
    axis = 1
)

transcriptomeData_sig_long = pd.melt(transcriptomeData_sig, id_vars = ["GeneID"], var_name = "day", value_name = "normCounts")
transcriptomeData_sig_long["day"] = [int(sub(r"[_D]", "", x[5:7])) for x in transcriptomeData_sig_long["day"]]
transcriptomeData_sig_long_mean = transcriptomeData_sig_long.groupby(['GeneID', 'day'], as_index=False).mean()
transcriptomeData_sig_mean = transcriptomeData_sig_long_mean.pivot_table(index='GeneID', columns='day', values='normCounts').reset_index()

def scaled_df(df, min_value, max_value):
    scaled = min_value + (df - df.min()) * (max_value - min_value) / (df.max() - df.min())
    return(scaled)

transcriptomeData_sig_scaled = pd.concat([transcriptomeData_sig_mean.GeneID,
                                          scaled_df(transcriptomeData_sig_mean.iloc[:,1:], 0, 1.5)],
                                          #zscore(transcriptomeData_sig_mean.iloc[:,1:], axis=0)],
                                          axis=1
                                        )
transcriptomeData_sig_scaled = transcriptomeData_sig_scaled[transcriptomeData_sig_scaled["GeneID"].isin(model_genes)]

def integrate_run_FluxSampling(phase, transcriptomeData_sig_scaled=transcriptomeData_sig_scaled, model_path = model_path, model_genes = model_genes):
    model = cobra.io.read_sbml_model(model_path)
    print("Constrain model for phase {0}".format(phase))
    dat = transcriptomeData_sig_scaled.loc[:, ["GeneID"] + phase]

    if len(dat.columns) > 2:
        dat["value"] = dat.loc[:,phase].mean(axis=1)
        dat = dat.loc[:, ["GeneID", "value"]]
    else:
        dat.rename(columns={phase[0]: 'value'}, inplace=True)

    transcripts = dat.set_index('GeneID')['value'].to_dict()

    for i in model_genes:
        if i not in transcripts.keys():
            transcripts[i] = 1000

    model_tr = map_transcriptome_data(model,transcripts)
    model_tr.write_sbml_model("iCHO2441_221-107_producing_{0}.xml".format(phase))

    # HERE : if reaction or gene in essential, bounds are open, else nothing

    with model_tr:
        # Run Flux Sampling.
        # From Strain 2023:
        # "The solution space was sampled 50,000,000 times, of which solutions
        #  were stored every 10,000 iterations, resulting in 5000 data‐points
        # per reaction."
        print("Run flux sampling")
        num_samples = 5000

        flux_samples = cobra.sampling.sample(model_tr, num_samples, thinning=10000, processes= nslots)

    flux_samples.to_csv("fluxSamples_2441_results_phase{0}.csv".format(phase))

timePhases = [[4],
              [5],
              [6,7,8],
              [11,12,14]]

pool = Pool(processes= nslots)

pool.map(integrate_run_FluxSampling, timePhases)
pool.close()

#phase = [6,7,8]
#solutions = {}
#for phase in timePhases:

    # print("Constrain model for phase {0}".format(phase))
    # dat = transcriptomeData_sig_scaled.loc[:, ["GeneID"] + phase]

    # if len(dat.columns) > 2:
    #     dat["value"] = dat.loc[:,phase].mean(axis=1)
    #     dat = dat.loc[:, ["GeneID", "value"]]
    # else:
    #     dat.rename(columns={phase[0]: 'value'}, inplace=True)

    # transcripts = dat.set_index('GeneID')['value'].to_dict()

    # for i in model_genes:
    #     if i not in transcripts.keys():
    #         transcripts[i] = 1000


    # model_tr = map_transcriptome_data(model,transcripts)

    # # HERE : if reaction or gene in essential, bounds are open, else nothing

    # with model_tr:
    #     # Run Flux Sampling.
    #     # From Strain 2023:
    #     # "The solution space was sampled 50,000,000 times, of which solutions
    #     #  were stored every 10,000 iterations, resulting in 5000 data‐points
    #     # per reaction."
    #     print("Run flux sampling")
    #     num_samples = 5000

    #    flux_samples = cobra.sampling.sample(model_tr, num_samples, thinning=10000, processes= nslots)

#    flux_samples.to_csv("fluxSamples_results_phase{0}.csv".format(phase))


