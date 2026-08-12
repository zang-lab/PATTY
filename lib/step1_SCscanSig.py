#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import time
import pandas as pd
import numpy as np
import copy,random
from joblib import load
# --------------------------
# custom package
# --------------------------
from PATTYpipe.Utility      import (sp, 
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   open_bed_file,
                                   add_coverPos,
                                   read_in_reads,
                                   read_in_reads_sc,
                                   scProcess)

# -------------------------- 
# main 
# --------------------------
def step1_SCscanSig(conf_dict,logfile):

    wlog('preprocess scCnT data with ArchR package',logfile)
#    scProcess(CnTdata,GENOME,readsCutoff)
    conf_dict['General']['scPackage'] = scProcess(conf_dict['General']['cuttag'],
                                                  conf_dict['General']['genome'],
                                                  conf_dict['options']['readcutoff'])
    if os.path.isfile("tmp_scProcess/tmp_cellDist.txt") and os.path.isfile("tmp_scProcess/tmp_highVarBin.bed"):
        pass
    else:
#    if conf_dict['General']['scPackage']  == "noPackage" and :
        ewlog("ArchR process was not done. PATTY cannot process scCnT data.",logfile)

    #if conf_dict['options']['binlist'] == "NA":
    rawbinFile = "tmp_scProcess/tmp_highVarBin.bed"
    #else:
    #    rawbinFile = conf_dict['options']['binlist']
    cmd = """zcat %s | %s intersect -a stdin -b %s -u > tmp_usebins.bed""" %(conf_dict['General']['genomebin'],conf_dict['General']['bedtools'],rawbinFile)
    tmplog = sp(cmd)#, logfile)
    binFile = "tmp_usebins.bed"

    # get cell-cell distance and target cell list :
    df = pd.read_csv("tmp_scProcess/tmp_cellDist.txt",header=0, index_col=0,sep="\t")
    distMat = df.values
    dist_cellnames = df.index.tolist()
    if len(conf_dict['options']["usecells"]) > 0:
        cell_list = list(set(dist_cellnames) & set(conf_dict['options']["usecells"])) 
    else:
        cell_list = dist_cellnames

    K=10

    binInfo = {}
    binSigATAC = {}
    binSigCnT = {}
    binextsize=1000
    with open_bed_file(binFile) as inf:
        for line in inf:
            ll = line.split()
            binname = ll[0]+":"+ll[1]+"-"+ll[2]
            binmid = int(ll[1])+ 100
            ext_left = binmid - int(binextsize/2)
            ext_right = binmid + int(binextsize/2)
            if ext_left >= 0:
                binInfo[binname] = ll
                for this_bin in range(ext_left, ext_right, 100):
                    relatedBinname = ll[0]+":"+str(this_bin)+"-"+str(this_bin+100)
                    binSigATAC[relatedBinname] = np.array([0]*100)
                    binSigCnT[relatedBinname] = np.zeros((len(cell_list),100))

    binList_real = list(sorted(binInfo.keys()))
    with open("tmp_binList.bed", "w") as f:
        for item in binList_real:
            f.write( "\t".join(binInfo[item])+"\n")

    wlog("read-in ATAC reads",logfile)
    [binSigATAC, ATAC_total] = read_in_reads(conf_dict['General']['atac'], binSigATAC)

    wlog("read-in sc CnT reads",logfile)
    CnT_total = {}
    for this_cell in cell_list:
        CnT_total[this_cell] = 0
    [binSigCnT, CnT_total] = read_in_reads_sc(conf_dict['General']['cuttag'], binSigCnT, CnT_total, cell_list)

    ### assign ATAC real bin value
    binSigATAC_real = np.zeros((len(binList_real), 100))
    for this_bin_idx in range(len(binList_real)):
        this_bin = binList_real[this_bin_idx]
        binChrm = this_bin.split(":")[0]#binInfo[]
        binStart = this_bin.split(":")[1].split("-")[0]
        binEnd = this_bin.split(":")[1].split("-")[1]
        binmid = int(binStart)+ 100
        ext_left = binmid - int(binextsize/2)
        ext_right = binmid + int(binextsize/2)
        outValue = []
        for related_bin in range(ext_left, ext_right, 100):
            relatedBinname = binChrm+":"+str(related_bin)+"-"+str(related_bin+100)
            this_array = binSigATAC[relatedBinname] / ATAC_total * 1e6
            ave_array = np.mean(this_array.reshape(10, 10), axis=1)
            outValue.extend(list(ave_array))
        binSigATAC_real[this_bin_idx,] = outValue

    ### assign CnT real bin value
    binSigCnT_real = {}
    for cellname in CnT_total.keys():
        binSigCnT_real[cellname] = np.zeros((len(binList_real), 100))

    split_bins = [b.split(":") for b in binList_real]
    binChrms = [split_bin[0] for split_bin in split_bins]
    starts_ends = [split_bin[1].split("-") for split_bin in split_bins]
    binStarts = np.array([int(se[0]) for se in starts_ends])
    binEnds = [int(se[1]) for se in starts_ends]
    binmids = binStarts + 100
    ext_lefts = binmids - int(binextsize / 2)
    ext_rights = binmids + int(binextsize / 2)
    
    # Precompute related bins for all indices
    all_related_bins = [np.arange(left, right, 100) for left, right in zip(ext_lefts, ext_rights)]
    all_relatedBinnames = [[binChrm + ":" + str(rb) + "-" + str(rb + 100) for rb in related_bins] for related_bins, binChrm in zip(all_related_bins, binChrms)]
    
    for this_bin_idx, relatedBinnames in enumerate(all_relatedBinnames):
        for cellindex, cellname in enumerate(cell_list):
            ave_values = []
            for relatedBinname in relatedBinnames:
                this_array = binSigCnT[relatedBinname][cellindex]
                if this_array.max() == 0:
                    ave_values.extend(np.zeros(10))
                else:
                    ave_values.extend(np.mean(this_array.reshape(10, 10), axis=1))  # Reshape and calculate mean
            binSigCnT_real[cellname][this_bin_idx, ] = ave_values

    wlog("sc correction",logfile)
    model = load(conf_dict['General']['model'])

    ATACdata_ext = np.expand_dims(binSigATAC_real, axis=2)
    metaCellPrediction = np.zeros((len(binList_real), len(cell_list)))
    for this_cell_idx in range(len(cell_list)):
        this_cell = cell_list[this_cell_idx]
        this_cell_index_distMat = dist_cellnames.index(this_cell)
        this_cell_row_distMat = distMat[this_cell_index_distMat,]
        top_nearest_index_distMat = np.argsort(this_cell_row_distMat)[1:(K+1)]
        top_nearest_cell_distMat = [dist_cellnames[i] for i in top_nearest_index_distMat]
        # make metaCell profile
        metaCellMat_CnT = np.zeros((len(binList_real),100))
        metaCellTotalReads_CnT = 0
        for cell_i in [this_cell]+top_nearest_cell_distMat:
            metaCellMat_CnT += binSigCnT_real[cell_i]
            metaCellTotalReads_CnT += CnT_total[cell_i]
        metaCellMatNorm_CnT = metaCellMat_CnT / metaCellTotalReads_CnT * 1e6
        # make data, prediction
        CnTdata_ext = np.expand_dims(metaCellMatNorm_CnT, axis=2)
        sigdata_ext = np.concatenate((CnTdata_ext, ATACdata_ext),axis=2)
        sigdata_flat = sigdata_ext.reshape(sigdata_ext.shape[0], -1)
        y_pred = model.predict_proba(sigdata_flat)[:, 1]
        metaCellPrediction[:,this_cell_idx] = y_pred
    np.savetxt("%s_correctMat.txt"%conf_dict['General']['outname'], metaCellPrediction, delimiter="\t",fmt="%.3e")
    return conf_dict


















