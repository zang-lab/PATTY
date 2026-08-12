#!/usr/bin/env python
 
# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
from string import *
import platform

# --------------------------
# custom package
# --------------------------

### tool function
import PATTYpipe
from PATTYpipe.Utility import (sp, 
                    raise_error,
                    wlog,
                    ewlog,
                    checkbedformat,
                    CMD
                    )

# --------------------------
# main 
# --------------------------

def step0_check_data(conf_dict,logfile):
    '''
    step0 integrate data 
    check and complement parameter
    '''
    ### check data path , format ,
    # input cuttag&atac dat
    if "/" in conf_dict['General']['cuttag']:
        if conf_dict['General']['cuttag'].startswith("/"):
            pass
        elif conf_dict['General']['cuttag'].startswith("~/"):
            homedir = os.path.expanduser("~")
            conf_dict['General']['cuttag'] = homedir +"/" + conf_dict['General']['cuttag'][1:]
        else:
            conf_dict['General']['cuttag'] = conf_dict['General']['startdir'] + conf_dict['General']['cuttag']
    else:
        conf_dict['General']['cuttag'] = conf_dict['General']['startdir'] + conf_dict['General']['cuttag']
    if not os.path.isfile(conf_dict['General']['cuttag']):
        ewlog("cuttag file %s not found"%(conf_dict['General']['cuttag']),logfile)
    if not conf_dict['General']['cuttag'].endswith('.bed') and not conf_dict['General']['cuttag'].endswith('.bed.gz'):
        ewlog('extenion of cuttag file is not .bed (nor .bed.gz)',logfile)
    checkbed = checkbedformat(conf_dict['General']['cuttag'])
    if checkbed == "fail":
        ewlog("cuttag file is not a PE/SE bed file",logfile)
    else:#elif checkbed in ["PE","SE"]:
        conf_dict['General']['format'] = checkbed
        wlog("detected bed file format is %s"%checkbed,logfile)

    if "/" in conf_dict['General']['atac']:
        if conf_dict['General']['atac'].startswith("/"):
            pass
        elif conf_dict['General']['atac'].startswith("~/"):
            homedir = os.path.expanduser("~")
            conf_dict['General']['atac'] = homedir +"/" + conf_dict['General']['atac'][1:]
        else:
            conf_dict['General']['atac'] = conf_dict['General']['startdir'] + conf_dict['General']['atac']
    else:
        conf_dict['General']['atac'] = conf_dict['General']['startdir'] + conf_dict['General']['atac']
    if not os.path.isfile(conf_dict['General']['atac']):
        ewlog("atac file %s not found"%(conf_dict['General']['atac']),logfile)
    if not conf_dict['General']['atac'].endswith('.bed') and not conf_dict['General']['atac'].endswith('.bed.gz'):
        ewlog('extenion of atac file is not .bed (nor .bed.gz)',logfile)
    checkbed = checkbedformat(conf_dict['General']['atac'])
    if checkbed == "fail":
        ewlog("atac file is not a PE/SE bed file",logfile)
    else:#elif checkbed in ["PE","SE"]:
        wlog("detected bed file format is %s"%checkbed,logfile)

    if "/" in conf_dict['General']['genomebin']:
        if conf_dict['General']['genomebin'].startswith("/"):
            pass
        elif conf_dict['General']['genomebin'].startswith("~/"):
            homedir = os.path.expanduser("~")
            conf_dict['General']['genomebin'] = homedir +"/" + conf_dict['General']['genomebin'][1:]
        else:
            conf_dict['General']['genomebin'] = conf_dict['General']['startdir'] + conf_dict['General']['genomebin']
    else:
        conf_dict['General']['genomebin'] = conf_dict['General']['startdir'] + conf_dict['General']['genomebin']
    if not os.path.isfile(conf_dict['General']['genomebin']):
        ewlog("genome-wide mappable bin file %s not found"%(conf_dict['General']['genomebin']),logfile)
    if not conf_dict['General']['genomebin'].endswith('.bed') and not conf_dict['General']['genomebin'].endswith('.bed.gz'):
        ewlog('extenion of genome-wide mappable bin file is not .bed or not .bed.gz',logfile)

    # --binlist
    #if conf_dict['options']['binlist']:
    #    if "/" in conf_dict['options']['binlist']:
    #        if conf_dict['options']['binlist'].startswith("/"):
    #            pass
    #        elif conf_dict['options']['binlist'].startswith("~/"):
    #            homedir = os.path.expanduser("~")
    #            conf_dict['options']['binlist'] = homedir +"/" + conf_dict['options']['binlist'][1:]
    #        else:
    #            conf_dict['options']['binlist'] = conf_dict['General']['startdir'] + conf_dict['options']['binlist']
    #    else:
    #        conf_dict['options']['binlist'] = conf_dict['General']['startdir'] + conf_dict['options']['binlist']
#
    #    if not os.path.isfile(conf_dict['options']['binlist']):
    #        wlog("external binlist file %s not found, PATTY use genome-wide mappable bins with sufficient reads"%(conf_dict['options']['binlist']),logfile)
    #        conf_dict['options']['binlist'] = "NA"
    #    checkbed = checkbedformat(conf_dict['options']['binlist'])
    #    if checkbed == "fail":
    #        wlog("external binlist file %s is not a bed file, PATTY use genome-wide mappable bins with sufficient reads"%(conf_dict['options']['binlist']),logfile)
    #        conf_dict['options']['binlist'] = "NA"
    #else:
    #    wlog("no external binlist file inputted, PATTY use genome-wide mappable bins file as target regions",logfile)
    #    conf_dict['options']['binlist'] = "NA"
        
    # check software
    OS = platform.system()
    # check system
    if OS == "Linux":
        bedtools_software = "bedtools_linux"
        bdg2bw_software = "bedGraphToBigWig_linux"
    elif OS == "Darwin":
        bedtools_software = "bedtools_mac"
        bdg2bw_software = "bedGraphToBigWig_mac"
    else:
        wlog("detected system is nither linux nor mac, try linux version",logfile)
        bedtools_software = "bedtools_linux"
        bdg2bw_software = "bedGraphToBigWig_linux"

    check_bedtools = sp("which bedtools")
    check_bdg2bw = sp("which bedGraphToBigWig")
    check_R = sp("which Rscript")

    if check_bedtools[0].decode("ascii") != "":
        wlog("bedtools installed",logfile)
        conf_dict['General']['bedtools'] = "bedtools"
    else:
        wlog("bedtools not installed in the default path, use PATTY-built-in bedtools",logfile)
        conf_dict['General']['bedtools'] = PATTYpipe.__path__[0]+"/external_script/%s"%bedtools_software

    if check_bdg2bw[0].decode("ascii") != "":
        wlog("bedGraphToBigWig(UCSCtools) installed",logfile)
        conf_dict['General']['bedGraphToBigWig'] = "bedGraphToBigWig"
    else:
        wlog("bedGraphToBigWig(UCSCtools) not installed in the default path, use PATTY-built-in bedGraphToBigWig",logfile)
        conf_dict['General']['bedGraphToBigWig'] = PATTYpipe.__path__[0]+"/external_script/%s"%bdg2bw_software

    if check_R[0].decode("ascii") != "":
        wlog("Rscript installed",logfile)
        conf_dict['General']['Rscript'] = "Rscript"
    else:
        ewlog("require Rscript installed in the default path",logfile)

    if conf_dict['General']['mode'] == "sc":

        # --cellnames
        if conf_dict['options']["cellnames"]:
            if "/" in conf_dict['options']["cellnames"]:
                if conf_dict['options']["cellnames"].startswith("/"):
                    pass
                elif conf_dict['options']["cellnames"].startswith("~/"):
                    homedir = os.path.expanduser("~")
                    conf_dict['options']["cellnames"] = homedir +"/" + conf_dict['options']["cellnames"][1:]
                else:
                    conf_dict['options']["cellnames"] = conf_dict['General']['startdir'] + conf_dict['options']["cellnames"]
            else:
                conf_dict['options']["cellnames"] = conf_dict['General']['startdir'] + conf_dict['options']["cellnames"]
            if os.path.isfile(conf_dict['options']["cellnames"]):
                conf_dict['options']["usecells"] = []
                wlog("readin used cellnames",logfile)
                inf = open(conf_dict['options']["cellnames"])
                for line in inf:
                    ll = line.strip().split("\t")
                    conf_dict['options']["usecells"].append(ll[0])
                inf.close()
                if len(conf_dict['options']["usecells"]) < 100:
                    wlog("less than 100 cells specified in cellnames file (%s cells). ignore --cellnames parameter"%(len(conf_dict['options']["usecells"])),logfile)
                    conf_dict['options']["usecells"] = []
            else:
                wlog("cellnames file %s not found, ignore --cellnames parameter"%(conf_dict['options']["cellnames"]),logfile)
                conf_dict['options']["usecells"] = []
        else:
            conf_dict['options']["usecells"] = []

    ### check trained model file
    conf_dict['General']['model'] = PATTYpipe.__path__[0]+"/refdata/%s_LR_CnTATAC_model.joblib"%(conf_dict['General']["factor"])

    return conf_dict

    
    
