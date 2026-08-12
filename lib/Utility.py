#!/usr/bin/env python
"""
Utilities used in PATTY pipeline
"""
import subprocess
import sys
import os
import math
import time
import random
import string
import gzip
import numpy

def CMD(cmd):
    os.system(cmd) 

def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
def sperr(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
   
def raise_error():
    '''
    Raise an error messgae and exit
    '''
    print('error occurs, check log file~!')
    sys.exit(1)


def detect_memory():
    meminfo={}#OrderedDict()
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                meminfo[line.split(':')[0].strip()] = line.split(':')[1].strip()
        totalM = meminfo['MemTotal'].split()
        #freeM = meminfo['MemFree'].split()
        if totalM[1].lower() == "kb":
            try:
                totalM_G = int(totalM[0])/1e6
                return totalM_G
            except:
                return 'NA'
        else:
            return 'NA'    
    except:
        return 'NA'

  
def wlog(message,logfile):
    '''
    print a message and write the message to logfile
    '''
    message = "### "+message
    print(message)
    os.system('echo "%s " >> %s'%(message,logfile))
    
def ewlog(message,logfile):
    '''
    print an error message and write the error message to logfile
    then exit Dr.seq
    error messages start with [ERROR]
    '''
    print("[ERROR] %s "%(message))
    os.system('echo "[ERROR] %s " >> %s'%(message,logfile))
    raise_error()
    
def rwlog(cmd,logfile) :
    '''
    print an (shell) command line and write the command line to logfile
    then conduct the command line
    command lines start with [CMD]
    '''
    print("[CMD] %s "%(cmd))
    os.system('echo "[CMD] %s " >> %s'%(cmd,logfile))
    CMD(cmd)

def rlogonly(cmd,logfile) :
    '''
    print an (shell) command line and write the command line to logfile
    then conduct the command line
    command lines start with [CMD]
    '''
    #print "[CMD] %s "%(cmd)
    os.system('echo "[CMD] %s " >> %s'%(cmd,logfile))
    CMD(cmd)

def checkbedformat(bedfile):
    if bedfile.endswith(".bed"):
        inf = open(bedfile)
        line = inf.readline()
        inf.close()
    else:
        inf = gzip.open(bedfile,'rb')
        line = inf.readline().decode("ascii")
        inf.close()

    ll = line.strip().split("\t")
    if len(ll) < 3:
        return "fail"
    try:
        a=int(ll[1])
        b=int(ll[2])
    except:
        return "fail"
    if len(ll) >= 6 and ll[5] in ["+","-"]:
        return "SE"
    elif len(ll)>=3:
        return "PE"
    else:
        return "fail"
    #peaknum = int(sp("wc -l %s"%(bedfile))[0].split()[0])
    #if peaknum < cutoff:
    #    return "lesspeak"
    #return peaknum#"pass"


### read in bed
def open_bed_file(filename):
    if filename.endswith(".bed.gz"):
        return gzip.open(filename, "rt")
    elif filename.endswith(".bed"): 
        return open(filename, "r")
    else:
        raise ValueError("Unsupported file type: " + filename)

## readin reads + map to bin

def add_coverPos(Read_cover_list,chrm,binSig):
    bin1 = int(Read_cover_list[0] / 100) * 100
    bin1_name = chrm+":"+str(bin1)+"-"+str(bin1+100)
    if bin1_name in binSig:
        binSig[bin1_name][(Read_cover_list[0] - bin1): min((Read_cover_list[-1] - bin1 + 1),100)] += 1
    bin2 = bin1+100
    bin2_name = chrm+":"+str(bin2)+"-"+str(bin2+100)
    if bin2_name in binSig:
        binSig[bin2_name][0: min((Read_cover_list[-1] - bin2 + 1),100)] += 1
    if (Read_cover_list[-1] - (bin2 + 100) ) >= 0:
        bin3 = bin2+100
        bin3_name = chrm+":"+str(bin3)+"-"+str(bin3+100)
        if bin3_name in binSig:
            binSig[bin3_name][0: min((Read_cover_list[-1] - bin3 + 1),100)] += 1
    return binSig

def read_in_reads(readFile,binSig):
    random.seed(1)
    total_reads_count= 0
    with open_bed_file(readFile) as inf:
        for line in inf:
            total_reads_count  += 1
            ll = line.split()
            reads_chrm = ll[0]
            reads_cover1=[]
            #reads_cover2=[]
            if len(ll) >= 6 and ll[5] in ["+","-"]:
                if ll[5] == "+":
                    plus = int(ll[1])
                    reads_cover1 = list(range(plus, plus+146))
                else:
                    minus = int(ll[2]) #- 1
                    reads_cover1 = list(range(minus-146,minus))
                binSig = add_coverPos(reads_cover1,reads_chrm,binSig)
            else:
                plus = int(ll[1])
                minus = int(ll[2])# - 1
                if random.randint(0,1) == 0:
                    reads_cover1 = list(range(plus, plus+146))
                else:
                    reads_cover1 = list(range(minus-146,minus))
                binSig = add_coverPos(reads_cover1,reads_chrm,binSig)
    return [binSig,total_reads_count]


def add_coverPos_sc(Read_cover_list,cellIdx,chrm,binSig):
    bin1 = int(Read_cover_list[0] / 100) * 100
    bin1_name = chrm+":"+str(bin1)+"-"+str(bin1+100)
    if bin1_name in binSig:
        binSig[bin1_name][cellIdx,(Read_cover_list[0] - bin1): min((Read_cover_list[-1] - bin1 + 1),100)] += 1
    bin2 = bin1+100
    bin2_name = chrm+":"+str(bin2)+"-"+str(bin2+100)
    if bin2_name in binSig:
        binSig[bin2_name][cellIdx,0: min((Read_cover_list[-1] - bin2 + 1),100)] += 1
    if (Read_cover_list[-1] - (bin2 + 100) ) >= 0:
        bin3 = bin2+100
        bin3_name = chrm+":"+str(bin3)+"-"+str(bin3+100)
        if bin3_name in binSig:
            binSig[bin3_name][cellIdx,0: min((Read_cover_list[-1] - bin3 + 1),100)] += 1
    return binSig

def read_in_reads_sc(readFile,binSig, total_reads_count, cell_list):
    random.seed(1)
    with open_bed_file(readFile) as inf:
        for line in inf:
            ll = line.split()
            reads_chrm = ll[0]
            reads_cell = ll[3]
            if not reads_cell in cell_list:
                continue
            cellindex = cell_list.index(reads_cell)
            if not reads_cell in total_reads_count.keys():
                continue
            total_reads_count[reads_cell]  += 1
            reads_cover1=[]
            #reads_cover2=[]
            plus = int(ll[1])
            minus = int(ll[2])# - 1
            if random.randint(0,1) == 0:
                reads_cover1 = list(range(plus, plus+146))
            else:
                reads_cover1 = list(range(minus-146,minus))
            binSig = add_coverPos_sc(reads_cover1,cellindex,reads_chrm,binSig)
    return [binSig,total_reads_count]


def prepare_chrom_size(genomeVer,outfile):
    mm10_sizes = {
    "chr1":195471971,"chr2":182113224,"chrX":171031299,"chr3":160039680,"chr4":156508116,
    "chr5":151834684,"chr6":149736546,"chr7":145441459,"chr10":130694993,"chr8":129401213,
    "chr14":124902244,"chr9":124595110,"chr11":122082543,"chr13":120421639,"chr12":120129022,
    "chr15":104043685,"chr16":98207768,"chr17":94987271,"chrY":91744698,"chr18":90702639,"chr19":61431566}

    hg38_sizes = {
    "chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,"chr5":181538259,
    "chr6":170805979,"chr7":159345973,"chrX":156040895,"chr8":145138636,"chr9":138394717,
    "chr11":135086622,"chr10":133797422,"chr12":133275309,"chr13":114364328,"chr14":107043718,
    "chr15":101991189,"chr16":90338345,"chr17":83257441,"chr18":80373285,"chr20":64444167,
    "chr19":58617616,"chrY":57227415,"chr22":50818468,"chr21":46709983}
    if genomeVer == "hg38":
        genome_dict = hg38_sizes
    else:
        genome_dict = mm10_sizes

    with open(outfile, "w") as f:
        for chrom, size in genome_dict.items():
            f.write(f"{chrom}\t{size}\n")


def scProcess(CnTdata,GENOME,readsCutoff):    
    scRscript="""
if(require("ArchR")){
  library(ArchR) 
  library(proxy)
  set.seed(1)
  
  CnTdata <- "%s"
  GENOME <- "%s"
  readsCutoff <- %s

  dir.create("tmp_scProcess")
  setwd("tmp_scProcess")
  set.seed(1)
  addArchRGenome(GENOME)
  addArchRThreads(threads = 1)
  
  inputFiles <- c(CnTdata)
  names(inputFiles)<-c("combine")
  
  ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    minTSS = 0, 
    minFrags = readsCutoff, 
    addTileMat = TRUE,
    force=TRUE,
    addGeneScoreMat = TRUE
  )
  
  projCnTpro <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = "tmp_ArchR",
    copyArrows = TRUE 
  )

  proj_CnTpro_3 <- addIterativeLSI(
      ArchRProj = projCnTpro,
      useMatrix = "TileMatrix", 
      name = "IterativeLSI", 
      iterations = 2, 
      clusterParams = list( #See Seurat::FindClusters
          resolution = c(0.2), 
          sampleCells = 10000, 
          n.start = 10
      ), 
      varFeatures = 25000, 
      dimsToUse = 1:30,
      seed=1,force=T
  )

  matSVD <- proj_CnTpro_3@reducedDims$IterativeLSI$matSVD
  pairwise_distances <- proxy::dist(as.matrix(matSVD), method = "Euclidean")
  distMat <- as.matrix(pairwise_distances)
  rownames(distMat) <- gsub("combine#", "", rownames(distMat))
  colnames(distMat) <- gsub("combine#", "", colnames(distMat))
  write.table(distMat, file="tmp_cellDist.txt",row.names=T,col.names=T,sep="\t",quote=F)
  
  options(scipen = 999)
  highVarBin_tmp <- cbind( as.vector(proj_CnTpro_3@reducedDims$IterativeLSI$LSIFeatures$seqnames),
                           as.numeric(proj_CnTpro_3@reducedDims$IterativeLSI$LSIFeatures$start),
                           as.numeric(proj_CnTpro_3@reducedDims$IterativeLSI$LSIFeatures$start) + 500)
  write.table(highVarBin_tmp, file="tmp_highVarBin.bed",row.names=F,col.names=F,sep="\t",quote=F)            
}else{
  simpleError("NoInstall")
}

"""%(CnTdata,GENOME,readsCutoff)
    outf = open("tmp_scRscript.r",'w')
    outf.write(scRscript)
    outf.close()

    tmplog = sp("Rscript tmp_scRscript.r")    
    if "simpleError" in tmplog[0].decode("ascii"):
        return("noPackage")
    else:
        return("yesPackage")

def split_chromosome_reads(bedfile,outname,scATAC10x,usechrom):
    if bedfile.endswith(".bed"):
        inf = open(bedfile)
    else:
        inf = gzip.open(bedfile,'rb')
    outf_chrM = open(outname + "_chrM.bed",'w')
    outf_chromatin = open(outname + "_chromatin.bed",'w')

    chrom_reads = {}
    for lineRaw in inf:
        if bedfile.endswith(".bed.gz"):
            line = lineRaw.decode("ascii")
        else:
            line = lineRaw
        ll = line.strip().split("\t")
        chrom=ll[0]
        if not chrom in usechrom:
            continue
        if not chrom in chrom_reads:
            chrom_reads[chrom] = 0
        chrom_reads[chrom] += 1
        if scATAC10x:
            if int(ll[1]) >= 4:
                newll = [ll[0], max(0,int(ll[1])-4), int(ll[2])+5, ll[3]  ]
                newline = "\t".join(map(str,newll))+"\n"
                if chrom == "chrM":
                    outf_chrM.write(newline)
                else:
                    outf_chromatin.write(newline)
        else:
            if chrom == "chrM":
                outf_chrM.write(line)
            else:
                outf_chromatin.write(line)
    inf.close()
    outf_chrM.close()
    outf_chromatin.close()
    return chrom_reads

def filter_highQcell_reads(outname,cutoff,usecells):

    cell_reads = {}
    allcells = []
    inf =open(outname + "_chromatin.bed")
    for line in inf:
        ll = line.strip().split("\t")
        chrom = ll[0]
        cellname = ll[3]
        if chrom == "chrM":
            continue
        if not cellname in cell_reads:
            cell_reads[cellname]=0
            allcells.append(cellname)
        cell_reads[cellname] += 1

    highQcells = []
    for cell in allcells:#cell_reads.keys():
        if cell_reads[cell] >= cutoff:
            highQcells.append(cell)

    inf.seek(0)

    if len(highQcells) < 100:
        return "fail"

    if len(usecells) == 0:
        usehighQcells = highQcells
    else:
#        usehighQcells = [value for value in highQcells if value in usecells]
        usehighQcells = [value for value in usecells if value in highQcells]
    if len(usehighQcells) < 100:
        finalcell = highQcells
        usetag = "highQ"
    else:
        finalcell = usehighQcells
        usetag = "highQuse"

    highQcellnum = len(finalcell)
    highQreadnum = 0
    outf = open(outname + "_highQcellReads.bed",'w')
    for line in inf:
        ll = line.strip().split("\t")
        chrom = ll[0]
        cellname = ll[3]
        if chrom == "chrM":
            continue
        if cellname in  finalcell:
            outf.write(line)
            highQreadnum += 1
    outf.close()
    inf.close()
    return [finalcell, len(usehighQcells), highQreadnum,len(cell_reads.keys()),usetag]


def fetchsignal_bw(bwsum,bwfile,chrm,start,end):
    cmd = '%s %s %s %s %s %s'%(bwsum,bwfile,chrm,start,end,end-start)
    rawsig = sp(cmd)[0].decode("ascii").strip().split()
    if len(rawsig) == (end-start):
        sig = [0 if "n" in x.lower() else int(x) for x in rawsig ]
        return sig
    else:
        return [0]*(end-start)


def bias_exp_cleavage_DNase(outname,peakfile,biasMat,kmer,bedtools,seq_dict,totalreads,dataformat):

    Cspan = 25
    kmer=int(kmer)
    flank = int(kmer/2)

    # extend peak to peak+Cspan, split to chromosome level
    chromosome_peak_dict = {}
    plus_cut_dict = {}
    minus_cut_dict = {}
    # split merge peak
    inf = open(peakfile)
    count = 0
    for line in inf:
        ll = line.split()
        chrom = ll[0]
        count += 1
        newll = [chrom, int(ll[1]) - Cspan, int(ll[2]) + Cspan, "mergePeak%s"%count]
        if not chrom in chromosome_peak_dict:
            chromosome_peak_dict[chrom] = open("%s_mergePeaks.bed"%(chrom),'w')
            plus_cut_dict[chrom] = open("%s_plusCuts.bed"%(chrom),'w')
            minus_cut_dict[chrom] = open("%s_minusCuts.bed"%(chrom),'w')
        chromosome_peak_dict[chrom].write("\t".join(map(str,newll))+"\n")
    inf.close()
    for chrom in chromosome_peak_dict.keys():
        chromosome_peak_dict[chrom].close()
    # split reads (to plus and minus)
    inf = open(totalreads)
    count = 0
    for line in inf:
        ll = line.split()
        chrom = ll[0]
        if not chrom in chromosome_peak_dict:
            continue
        if dataformat == "PE":
            count += 1
            newll = [chrom, ll[1] ,int(ll[1])+1, "c%s"%count,".","+"]
            plus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
            count += 1
            newll = [chrom, int(ll[2])-1 ,ll[2], "c%s"%count,".","-"]
            minus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
        else:
            count += 1
            if ll[5] == "+":
                newll = [chrom, ll[1] ,int(ll[1])+1, "c%s"%count,".","+"]
                plus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
            else:
                newll = [chrom, int(ll[2])-1 ,ll[2], "c%s"%count,".","-"]
                minus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
    inf.close()
    for chrom in chromosome_peak_dict.keys():
        plus_cut_dict[chrom].close()
        minus_cut_dict[chrom].close()


    outf_plus = open(outname + "_biasExpCuts_plus.bdg",'w')
    outf_minus = open(outname + "_biasExpCuts_minus.bdg",'w')
    outf_plusCuts = open(outname + "_cleavage_plus.bdg",'w')
    outf_minusCuts = open(outname + "_cleavage_minus.bdg",'w')

    ### for each chromosome, intersect and calculate cleavage pattern
    for chrom in chromosome_peak_dict.keys():
        OVcmd1 = """%s intersect -a %s -b %s -wao > %s """%(bedtools,"%s_mergePeaks.bed"%(chrom), "%s_plusCuts.bed"%(chrom), "%s_plusCutsOnPeak.bed"%(chrom))
        OVcmd2 = """%s intersect -a %s -b %s -wao > %s """%(bedtools,"%s_mergePeaks.bed"%(chrom), "%s_minusCuts.bed"%(chrom), "%s_minusCutsOnPeak.bed"%(chrom))
        os.system(OVcmd1)
        os.system(OVcmd2)

        #### readin bias vector
        mergePeaks = []
        mergePeak_dict = {}
        thisChrom_plus_bias_dict = {}
        thisChrom_minus_bias_dict = {}
        thisChrom_plus_cuts_dict = {}
        thisChrom_minus_cuts_dict = {}
        inf = open("%s_mergePeaks.bed"%(chrom))
        for line in inf:
            ll = line.strip().split("\t")
            chrm = ll[0]
            start = int(ll[1])
            end = int(ll[2])
            peakname = ll[3]
            plus_seq_all = seq_dict[chrm][(start-flank):(end+flank)]
            minus_seq_all = rev(plus_seq_all)
            plus_single_bias_enc_vector = []
            minus_single_bias_enc_vector = []
            for pos in range(end-start):
                plus_seq = plus_seq_all[pos:(pos+kmer)]
                minus_seq = minus_seq_all[pos:(pos+kmer)]
                if len(plus_seq) == kmer and not "N" in plus_seq:
                    plus_bias_enc = 2**biasMat[plus_seq]
                else:
                    plus_bias_enc = 1
    
                if len(minus_seq) == kmer and not "N" in minus_seq:
                    minus_bias_enc = 2**biasMat[minus_seq]
                else:
                    minus_bias_enc = 0
                plus_single_bias_enc_vector.append(plus_bias_enc)
                minus_single_bias_enc_vector.append(minus_bias_enc)
            Plus_Single_encBias = numpy.array(plus_single_bias_enc_vector)
            Minus_Single_encBias = numpy.array(minus_single_bias_enc_vector[::-1]) 
            thisChrom_plus_bias_dict[peakname] = Plus_Single_encBias
            thisChrom_minus_bias_dict[peakname] = Minus_Single_encBias
            thisChrom_plus_cuts_dict[peakname] = [0]*(end-start)#Plus_Single_encBias
            thisChrom_minus_cuts_dict[peakname] = [0]*(end-start)#Minus_Single_encBias
            mergePeaks.append(peakname)
            mergePeak_dict[peakname] = [chrm,start+Cspan,end-Cspan,peakname]
        inf.close()

        ### readin cuts 
        inf = open("%s_plusCutsOnPeak.bed"%(chrom))
        for line in inf:        
            ll = line.strip().split("\t")
            chrm = ll[0]
            peak_start = int(ll[1])#+Cspan
            peak_end = int(ll[2])#-Cspan
            peakname = ll[3]
            read_start = int(ll[5])
            if read_start != -1:
                relative_pos = read_start - peak_start
            thisChrom_plus_cuts_dict[peakname][relative_pos] += 1
        inf.close()

        inf = open("%s_minusCutsOnPeak.bed"%(chrom))
        for line in inf:        
            ll = line.strip().split("\t")
            chrm = ll[0]
            peak_start = int(ll[1])#+Cspan
            peak_end = int(ll[2])#-Cspan
            peakname = ll[3]
            read_start = int(ll[5])
            if read_start != -1:
                relative_pos = read_start - peak_start
            thisChrom_minus_cuts_dict[peakname][relative_pos] += 1
        inf.close()

        # calculate biasExpCuts
        for peakname in mergePeaks:
            peakinfo = mergePeak_dict[peakname]
            chrm = peakinfo[0]
            start = peakinfo[1]
            end = peakinfo[2]
            Plus_Single_encBias = thisChrom_plus_bias_dict[peakname]
            Minus_Single_encBias = thisChrom_minus_bias_dict[peakname]
            plus_vector = thisChrom_plus_cuts_dict[peakname]
            minus_vector = thisChrom_minus_cuts_dict[peakname]

            for outpos in range(Cspan,(end-start+Cspan)):
                this_plus_single_enc = Plus_Single_encBias[outpos]
                this_minus_single_enc = Minus_Single_encBias[outpos]
                this_plus_sum_enc = sum(Plus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
                this_minus_sum_enc = sum(Minus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
    
                this_plus = plus_vector[outpos]
                this_minus = minus_vector[outpos]
                this_plus_cuts_sum = sum(plus_vector[(outpos-Cspan):(outpos+Cspan)])
                this_minus_cuts_sum = sum(minus_vector[(outpos-Cspan):(outpos+Cspan)])
    
                out_chrm = chrm
                out_start = start + outpos - Cspan
                out_end = out_start+1
                
                expcut_plus_enc = this_plus_cuts_sum * (this_plus_single_enc/this_plus_sum_enc)
                expcut_minus_enc = this_minus_cuts_sum * (this_minus_single_enc/this_minus_sum_enc)
    
                outf_plus.write("\t".join( map(str, [out_chrm,out_start,out_end,round(expcut_plus_enc,6)] ))+"\n")
                outf_minus.write("\t".join( map(str, [out_chrm,out_start,out_end,round(expcut_minus_enc,6)] ))+"\n")
                outf_plusCuts.write("\t".join( map(str, [out_chrm,out_start,out_end,round(this_plus,6)] ))+"\n")
                outf_minusCuts.write("\t".join( map(str, [out_chrm,out_start,out_end,round(this_minus,6)] ))+"\n")
    outf_plus.close()
    outf_minus.close()
    outf_plusCuts.close()
    outf_minusCuts.close()
    inf.close()


def bias_exp_cleavage_ATAC(outname,peakfile,biasMat,kmer,bedtools,seq_dict,totalreads,dataformat):

    offset=9
    Cspan = 25
    kmer=int(kmer)
    flank = int(kmer/2)

    chromosome_peak_dict = {}
    plus_cut_dict = {}
    minus_cut_dict = {}
    # split merge peak
    inf = open(peakfile)
    count = 0
    for line in inf:
        ll = line.split()
        chrom = ll[0]
        count += 1
        newll = [chrom, int(ll[1]) - Cspan, int(ll[2]) + Cspan, "mergePeak%s"%count]
        if not chrom in chromosome_peak_dict:
            chromosome_peak_dict[chrom] = open("%s_mergePeaks.bed"%(chrom),'w')
            plus_cut_dict[chrom] = open("%s_plusCuts.bed"%(chrom),'w')
            minus_cut_dict[chrom] = open("%s_minusCuts.bed"%(chrom),'w')
        chromosome_peak_dict[chrom].write("\t".join(map(str,newll))+"\n")
    inf.close()
    for chrom in chromosome_peak_dict.keys():
        chromosome_peak_dict[chrom].close()
    # split reads (to plus and minus)
    inf = open(totalreads)
    count = 0
    for line in inf:
        ll = line.split()
        chrom = ll[0]
        if not chrom in chromosome_peak_dict:
            continue
        if dataformat == "PE":
            count += 1
            newll = [chrom, ll[1] ,int(ll[1])+1, "c%s"%count,".","+"]
            plus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
            count += 1
            newll = [chrom, int(ll[2])-1 ,ll[2], "c%s"%count,".","-"]
            minus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
        else:
            count += 1
            if ll[5] == "+":
                newll = [chrom, ll[1] ,int(ll[1])+1, "c%s"%count,".","+"]
                plus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
            else:
                newll = [chrom, int(ll[2])-1 ,ll[2], "c%s"%count,".","-"]
                minus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
    inf.close()
    for chrom in chromosome_peak_dict.keys():
        plus_cut_dict[chrom].close()
        minus_cut_dict[chrom].close()


    outf_plus = open(outname + "_biasExpCuts_plus.bdg",'w')
    outf_minus = open(outname + "_biasExpCuts_minus.bdg",'w')
    outf_plusCuts = open(outname + "_cleavage_plus.bdg",'w')
    outf_minusCuts = open(outname + "_cleavage_minus.bdg",'w')

    ### for each chromosome, intersect and calculate cleavage pattern
    for chrom in chromosome_peak_dict.keys():
        OVcmd1 = """%s intersect -a %s -b %s -wao > %s """%(bedtools,"%s_mergePeaks.bed"%(chrom), "%s_plusCuts.bed"%(chrom), "%s_plusCutsOnPeak.bed"%(chrom))
        OVcmd2 = """%s intersect -a %s -b %s -wao > %s """%(bedtools,"%s_mergePeaks.bed"%(chrom), "%s_minusCuts.bed"%(chrom), "%s_minusCutsOnPeak.bed"%(chrom))
        os.system(OVcmd1)
        os.system(OVcmd2)

        #### readin bias vector
        mergePeaks = []
        mergePeak_dict = {}
        thisChrom_plus_bias_dict = {}
        thisChrom_minus_bias_dict = {}
        thisChrom_plus_cuts_dict = {}
        thisChrom_minus_cuts_dict = {}
        inf = open("%s_mergePeaks.bed"%(chrom))
        for line in inf:
            ll = line.strip().split("\t")
            chrm = ll[0]
            start = int(ll[1])
            end = int(ll[2])
            peakname = ll[3]
            plus_single_bias_enc_vector = []
            minus_single_bias_enc_vector = []
            seqall = seq_dict[chrm][(start-flank-offset):(end+flank+offset)]#fetchseq_2bit(twoBitFa,seq2bit,chrm,start-Cspan-flank,end+Cspan+flank)
            for pos in range(offset,end-start+offset):
                plus_forward_seq = seqall[pos:(pos+kmer)]
                plus_reverse_seq = rev(seqall[(pos+offset):(pos+kmer+offset)])
                minus_forward_seq = rev(seqall[(pos+1):(pos+1+kmer)])
                minus_reverse_seq = seqall[(pos+1-offset):(pos+1+kmer-offset)]
                if len(plus_forward_seq) == kmer and not "N" in plus_forward_seq:
                    plus_bias = 2**biasMat[plus_forward_seq]
                else:
                    plus_bias = 1
                if len(minus_forward_seq) == kmer and not "N" in minus_forward_seq:
                    minus_bias = 2**biasMat[minus_forward_seq]
                else:
                    minus_bias = 1
                if len(plus_reverse_seq) == kmer and not "N" in plus_reverse_seq:
                    plus_reverse_bias = 2**biasMat[plus_reverse_seq]
                else:
                    plus_reverse_bias = 1
                if len(minus_reverse_seq) == kmer and not "N" in minus_reverse_seq:
                    minus_reverse_bias = 2**biasMat[minus_reverse_seq]
                else:
                    minus_reverse_bias = 1    
                plus_cb_bias = numpy.sqrt(plus_bias * plus_reverse_bias ) 
                minus_cb_bias = numpy.sqrt(minus_bias * minus_reverse_bias)
                plus_single_bias_enc_vector.append(plus_cb_bias)
                minus_single_bias_enc_vector.append(minus_cb_bias)
                
            Plus_Single_encBias = numpy.array(plus_single_bias_enc_vector)
            Minus_Single_encBias = numpy.array(minus_single_bias_enc_vector)
            thisChrom_plus_bias_dict[peakname] = Plus_Single_encBias
            thisChrom_minus_bias_dict[peakname] = Minus_Single_encBias
            thisChrom_plus_cuts_dict[peakname] = [0]*(end-start)#Plus_Single_encBias
            thisChrom_minus_cuts_dict[peakname] = [0]*(end-start)#Minus_Single_encBias
            mergePeaks.append(peakname)
            mergePeak_dict[peakname] = [chrm,start+Cspan,end-Cspan,peakname]
        inf.close()
        ### readin cuts 
        inf = open("%s_plusCutsOnPeak.bed"%(chrom))
        for line in inf:        
            ll = line.strip().split("\t")
            chrm = ll[0]
            peak_start = int(ll[1])#+Cspan
            peak_end = int(ll[2])#-Cspan
            peakname = ll[3]
            read_start = int(ll[5])
            if read_start != -1:
                relative_pos = read_start - peak_start
            thisChrom_plus_cuts_dict[peakname][relative_pos] += 1
        inf.close()

        inf = open("%s_minusCutsOnPeak.bed"%(chrom))
        for line in inf:        
            ll = line.strip().split("\t")
            chrm = ll[0]
            peak_start = int(ll[1])#+Cspan
            peak_end = int(ll[2])#-Cspan
            peakname = ll[3]
            read_start = int(ll[5])
            if read_start != -1:
                relative_pos = read_start - peak_start
            thisChrom_minus_cuts_dict[peakname][relative_pos] += 1
        inf.close()

        # calculate biasExpCuts
        for peakname in mergePeaks:
            peakinfo = mergePeak_dict[peakname]
            chrm = peakinfo[0]
            start = peakinfo[1]
            end = peakinfo[2]
            Plus_Single_encBias = thisChrom_plus_bias_dict[peakname]
            Minus_Single_encBias = thisChrom_minus_bias_dict[peakname]
            plus_vector = thisChrom_plus_cuts_dict[peakname]
            minus_vector = thisChrom_minus_cuts_dict[peakname]

            for outpos in range(Cspan,(end-start+Cspan)):
                this_plus_single_enc = Plus_Single_encBias[outpos]
                this_minus_single_enc = Minus_Single_encBias[outpos]
                this_plus_sum_enc = sum(Plus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
                this_minus_sum_enc = sum(Minus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
    
                this_plus = plus_vector[outpos]
                this_minus = minus_vector[outpos]
                this_plus_cuts_sum = sum(plus_vector[(outpos-Cspan):(outpos+Cspan)])
                this_minus_cuts_sum = sum(minus_vector[(outpos-Cspan):(outpos+Cspan)])
    
                out_chrm = chrm
                out_start = start + outpos - Cspan
                out_end = out_start+1
                
                expcut_plus_enc = this_plus_cuts_sum * (this_plus_single_enc/this_plus_sum_enc)
                expcut_minus_enc = this_minus_cuts_sum * (this_minus_single_enc/this_minus_sum_enc)
    
                outf_plus.write("\t".join( map(str, [out_chrm,out_start,out_end,round(expcut_plus_enc,6)] ))+"\n")
                outf_minus.write("\t".join( map(str, [out_chrm,out_start,out_end,round(expcut_minus_enc,6)] ))+"\n")
                outf_plusCuts.write("\t".join( map(str, [out_chrm,out_start,out_end,round(this_plus,6)] ))+"\n")
                outf_minusCuts.write("\t".join( map(str, [out_chrm,out_start,out_end,round(this_minus,6)] ))+"\n")
    outf_plus.close()
    outf_minus.close()
    outf_plusCuts.close()
    outf_minusCuts.close()
    inf.close()
  
  

def bias_peakXcell_mat(outname,bedtools,chrom_list, kmer, biasDict, seqDict, usecells, datatype,peakminreads):

    flank = int(int(kmer)/2)
    offset=9
    peakminreads = int(peakminreads)
    #if peakmaxreads == "NA":
    #    peakmaxreads = int(1e10)
    #else:
    #    peakmaxreads = int(peakmaxreads)

    peakfile = outname + "_summitEXT.bed"
    readfile = outname + "_highQcellReads.bed"
    peakFeatures = open(outname + "_peakFeatures.txt","w")
    peakXcellMat = open(outname + "_peakXcellMat.txt","w")

    newll = ['chrm','start','end','peakname','score','cutsSum','avebias']
    peakFeatures.write("\t".join(newll)+"\n")
    
    newll = ['peakname'] + usecells
    peakXcellMat.write("\t".join(newll)+"\n")

    for chrom in chrom_list:
        if chrom == "chrM":
            continue
        cmdpeak = """awk '{OFS="\\t";if($1=="%s") print $0}' %s > %s"""%(chrom,peakfile, outname+"_tmpSCpeaks.bed")
        tmplog = sp(cmdpeak)
        try:
            chrom_peaknum = int(sp('wc -l %s'%(outname+"_tmpSCpeaks.bed"))[0].decode("ascii").split()[0])
        except:
            #wlog('no peak detected for %s, skip %s'%(chrom,chrom),logfile)
            continue 
        if chrom_peaknum == 0:
            #wlog('no peak detected for %s, skip %s'%(chrom,chrom),logfile)
            continue
        cmdread = """awk '{OFS="\\t";if($1=="%s") print $1,$2,$2+1,$4,".","+\\n"$1,$3-1,$3,$4,".","-"}' %s > %s"""%(chrom,outname+"_highQcellReads.bed", outname+"_tmpSCreads.bed")
        tmplog = sp(cmdread)

        cmdassign = "%s intersect -a %s_tmpSCpeaks.bed -b %s_tmpSCreads.bed -wao | awk '{if($NF > 0 ) print $0}'|sort -k 4,4 > %s_scOVcleavage.bed"%(bedtools,outname,outname,outname)
        tmplog = sp(cmdassign)

        inf = open("%s_scOVcleavage.bed"%outname)
        this_peak = "NA"
        for line in inf:
            ll = line.strip().split("\t")
            if this_peak == "NA":
                this_peak = ll[3]
                this_loci = ll[:5]
                cell_count = [0]*len(usecells)
                cutsSum = 0
                biasSum = 0
                read_info = ll[5:11]
                if  read_info[3] in usecells:
                    this_bias = reads_level_bias(read_info,datatype,seqDict,biasDict,flank)
                    if this_bias != "NA":
                        cutsSum += 1
                        biasSum += this_bias
                        cell_count[usecells.index(read_info[3])] += 1
            elif this_peak != ll[3]:
                if cutsSum >= peakminreads :#and cutsSum <= peakmaxreads:
                    avebias = biasSum / cutsSum#,6)
                    newll =  this_loci + [cutsSum,avebias] 
                    peakFeatures.write("\t".join(map(str,newll))+"\n")
                    newll = [this_peak] + cell_count
                    peakXcellMat.write("\t".join(map(str,newll))+"\n")
                this_peak = ll[3]
                this_loci = ll[:5]
                cell_count = [0]*len(usecells)
                cutsSum = 0
                biasSum = 0
                read_info = ll[5:11]
                if  read_info[3] in usecells:
                    this_bias = reads_level_bias(read_info,datatype,seqDict,biasDict,flank)
                    if this_bias != "NA":
                        cutsSum += 1
                        biasSum += this_bias
                        cell_count[usecells.index(read_info[3])] += 1
            else:
                read_info = ll[5:11]
                if  read_info[3] in usecells:
                    this_bias = reads_level_bias(read_info,datatype,seqDict,biasDict,flank)
                    if this_bias != "NA":
                        cutsSum += 1
                        biasSum += this_bias
                        cell_count[usecells.index(read_info[3])] += 1
        
        if cutsSum >= peakminreads :#and cutsSum <= peakmaxreads:
            avebias = biasSum / cutsSum#,6)
            newll =  this_loci + [cutsSum,avebias] 
            peakFeatures.write("\t".join(map(str,newll))+"\n")
            newll = [this_peak] + cell_count
            peakXcellMat.write("\t".join(map(str,newll))+"\n")
        
        inf.close()

    peakFeatures.close()
    peakXcellMat.close()

def reads_level_bias(region, datatype, seqdict, biasDict, flank):
    flank = int(flank)
    offset=9
    chrm = region[0]
    start = int(region[1])#-100#upstream_ext
    end = int(region[2])#center+100#start + 200#fulllen
    strand = region[5]
    if datatype == "ATAC":
        if strand == "+":
            forward_seq = seqdict[chrm][(start-flank):(start+flank)].upper()
            reverse_seq = rev(seqdict[chrm][(start+offset-flank):(start+offset+flank)].upper())
        else:
            forward_seq = rev(seqdict[chrm][(end-flank):(end+flank)].upper())
            reverse_seq = seqdict[chrm][(end-offset-flank):(end-offset+flank)].upper()
    
        if forward_seq in biasDict and reverse_seq in biasDict:
            bias_score = (biasDict[forward_seq] + biasDict[reverse_seq])/2
        else:
            bias_score = "NA"
    else:

        if strand == "+":
            forward_seq = seqdict[chrm][(start-flank):(start+flank)].upper()
        else:
            forward_seq = rev(seqdict[chrm][(end-flank):(end+flank)].upper())
    
        if forward_seq in biasDict: #and biasDict.has_key(reverse_seq):
            bias_score = biasDict[forward_seq]# + biasDict[reverse_seq])/2
        else:
            bias_score = "NA"
            
    return bias_score


def bwsigAve(bwfile,chrm,start,end,software):
    cmd = '%s %s %s %s %s 1'%(software,bwfile,chrm,start,end)
#    print sp(cmd)
    bwS = sp(cmd)[0].strip()
    if len(bwS) == 0:#bwS == "":
        return 0
    else:
        sigAve = float( bwS)#numpy.array(map(float,sp(cmd)[0].strip().split("\t")))
        return sigAve#[CpGcount,aveME]

def readAnnotation(annotation):
    '''
    read full annotation file and output as a dictionary 
    file format is fixed to UCSC full annotation format
    '''
    inf = open(annotation)
    outdict = {}
    for line in inf:
        ll = line.split()
        outdict[ll[1]] = ll[12]
    return outdict
    
     
def createDIR(dirname):
    '''
    check dir name and create new dir
    '''
    if not os.path.isdir(dirname):
        os.system('mkdir %s'%(dirname))

       
