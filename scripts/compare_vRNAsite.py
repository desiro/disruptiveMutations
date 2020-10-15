#!/usr/bin/env python3
# script: compare.py
# author: Daniel Desiro'
# dependencies: numpy, RNAfold, VARNAv3-93.jar, inkscape
script_usage="""
usage
    compare.py -f <in_fasta> -p <out_prefix> [options]

version
    compare.py 0.0.1 (alpha)

dependencies
    x

description
    Validates significance of flexible regions in a sequence.

--prefix,-pfx
    output directory and prefix for result files

--vRNAsite,-vst
    vRNAsite table

--candidatePeak,-cdp
    define minimum peak MFE value (default: -10.0)

--candidateMFE,-cdm
    define minimum MFE value (default: -10.0)

--candidateOverlap,-cdo
    define the number of overlapping bases (default: 5)

--orientation,-ort
    only allow specific orientations of the vRNPs; A = arbitrary, U = upright,
    R = reverse (default: A) (choices: A,U,R,UR,AU,AR,AUR)

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--minlength,-mle
    minimum length for an interaction (default: 5)

--intra,-tra
    also do intra RNA combinations (default: False)

--reposition,-rep
    add reposition data if available (default: False)

reference
    Reference.
"""

import argparse as ap
import sys
import os
import re
import time
from operator import attrgetter, itemgetter
from itertools import combinations, product, combinations_with_replacement
from statistics import mean, median


### old


import random
try:
    from RNA import fold_compound, cvar, CONSTRAINT_DB, CONSTRAINT_DB_DEFAULT, bp_distance
except:
    print("Error: The mandatory library ViennaRNA 2.4 is missing, please read the README.md for more information!")
    exit()
from itertools import combinations, product, combinations_with_replacement
from multiprocessing import Process, Manager
from math import ceil, floor
from numpy import zeros
from pandas import DataFrame, set_option
# plotting
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from numpy import arange, mean, prod, array, matrix, zeros, nonzero
from numpy import sum as npsum
import pickle
from numpy import arange, mean, prod, array, asarray
import itertools
from subprocess import Popen, PIPE, call
from multiprocessing import Pool, Process, Manager, Lock
from copy import deepcopy
from itertools import permutations, product, combinations_with_replacement, combinations, chain
from collections import Counter

import networkx as nx
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt




################################################################################
## main
################################################################################

def main(opt):
    ############################################################################
    ## load fasta file
    time_s = getTime()
    print(f"Status: Read tables ...")
    tab_list = readTable(opt["var_vst"], **opt)
    time_s = getTime(time_s, f"Read table")
    ############################################################################
    ## create output folder
    opt = makeDir(**opt)
    ############################################################################
    ## compare data sets 103 69 215
    print(f"Status: Compare data sets ...")
    name = os.path.basename(opt["var_pfx"])
    compareData(tab_list, name, **opt)
    time_s = getTime(time_s, f"Compare data sets")
    return opt["var_pfx"]




################################################################################
## functions
################################################################################

def old_readTable(table, **opt):
    ## read table file
    ## TODO: update to vRNAsite with SPLASH data
    tab_list = list()
    with open(table, "r") as tabin:
        header = next(tabin)
        for i,line in enumerate(tabin):
            lnk = dict()
            line = line.strip().split()
            if len(line[0].split("-")) == 2: line[0] = f"{line[0]}-A"
            if len(line[3].split("-")) == 2: line[3] = f"{line[3]}-A"
            aSeq, aType, aOrt = line[0].split("-")
            bSeq, bType, bOrt = line[3].split("-")
            lnk["aOrt"], lnk["bOrt"] = aOrt, bOrt
            if opt["var_ort"]:
                if lnk["aOrt"] not in list(opt["var_ort"]) or lnk["bOrt"] not in list(opt["var_ort"]): continue
            # add to dictionary in correct order
            lnk["aSeq"], lnk["aType"], lnk["ai"], lnk["aj"] = aSeq, aType, line[1], line[2]
            lnk["bSeq"], lnk["bType"], lnk["bi"], lnk["bj"] = bSeq, bType, line[4], line[5]
            lnk["peak"], lnk["mfe"], lnk["rank"] = line[6], line[10], 0
            lnk["sp_mean"], lnk["sp_min"], lnk["sp_max"] = line[7], line[8], line[9]
            lnk["RNA"], lnk["structure"] = line[11], line[12]
            lnk["a_mfe"], lnk["a_structure"], lnk["b_mfe"], lnk["b_structure"], lnk["free_mfe"], lnk["free_structure"] = line[13:19]
            lk = links(**lnk)
            if lk.aj - lk.ai + 1 < opt["var_mle"]: continue
            if lk.bj - lk.bi + 1 < opt["var_mle"]: continue
            if lk.peak > opt["var_cdp"]: continue
            if lk.mfe > opt["var_cdm"]: continue
            print(lk.plot(" "))
            exit()
            tab_list.append(lk)
    tab_list = sorted(tab_list, key=attrgetter("peak"))
    for r,t in enumerate(tab_list, 1):
        t.rank = r
    return tab_list

def readTable(table, **opt):
    ## read table file
    ## TODO: update to vRNAsite with SPLASH data
    tab_list = list()
    with open(table, "r") as tabin:
        header = next(tabin)
        hlist = header.strip().split()
        for i,line in enumerate(tabin):
            lnk = dict()
            line = line.strip().split()
            if len(hlist) != len(line):
                print(f"Error: Not the same number of header and list elements!")
                sys.exit()
            for name,item in zip(hlist,line):
                lnk[name] = item
            ## check for orientation
            if len(lnk["aSeq"].split("-")) == 2: aSequence = f"{lnk['aSeq']}-A"
            if len(lnk["bSeq"].split("-")) == 2: bSequence = f"{lnk['bSeq']}-A"
            lnk["aSeq"], lnk["aType"], lnk["aOrt"] = aSequence.split("-")
            lnk["bSeq"], lnk["bType"], lnk["bOrt"] = bSequence.split("-")
            lnk["rank"] = 0
            if opt["var_ort"]:
                if lnk["aOrt"] not in list(opt["var_ort"]) or lnk["bOrt"] not in list(opt["var_ort"]): continue
            ## create class and remove items
            lk = links(**lnk)
            if lk.aj - lk.ai + 1 < opt["var_mle"]: continue
            if lk.bj - lk.bi + 1 < opt["var_mle"]: continue
            if lk.peak > opt["var_cdp"]: continue
            if lk.mfe > opt["var_cdm"]: continue
            tab_list.append(lk)
    tab_list = sorted(tab_list, key=attrgetter("peak"))
    for r,t in enumerate(tab_list, 1):
        t.rank = r
    return tab_list

class links(object):
    def __init__(self, **data):
        self.__dict__.update((k,trans(v)) for k,v in data.items())
    def trans(s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    def plot(self, sep):
        ldat = sep.join([f"{var}" for key,var in vars(self).items()])
        return ldat

def comp(tx, x):
    ## compare positions
    if min([t.aj for t in tx]) - max([t.ai for t in tx]) >= x and\
       min([t.bj for t in tx]) - max([t.bi for t in tx]) >= x: return True
    else: return False

def getMutants(tlist):
    ## get mutant names
    ta = [(tx.aSeq, tx.aType) for tx in tlist if tx.aType != "WT"]
    tb = [(tx.bSeq, tx.bType) for tx in tlist if tx.bType != "WT"]
    return sorted(list(set(ta+tb)))

def removeIntra(tlist):
    ## remove intra interactions
    return [tx for tx in tlist if tx.aSeq != tx.bSeq]

def separateInteractions(tlist):
    ## separate mutants and WT
    WT_list = [tx for tx in tlist if tx.aType == "WT" and tx.bType == "WT"]
    MT_list = [tx for tx in tlist if tx.aType != "WT"  or tx.bType != "WT"]
    return WT_list, MT_list

def plotList1(pList, pType, mType):
    pOut = ""
    for t in pList:
        line = f"{mType}\t{pType}\t{t.aSeq}-{t.aType}\t{t.ai}\t{t.aj}\t{t.bSeq}-{t.bType}\t{t.bi}\t{t.bj}\t{t.peak}\t{t.rank}\t{t.sp_mean}\t{t.sp_min}\t{t.sp_max}\t{t.mfe}\t{t.RNA}\t{t.structure}\t{t.aOrt}\t{t.bOrt}\t{t.a_mfe}\t{t.a_structure}\t{t.b_mfe}\t{t.b_structure}\t{t.free_structure}\t{t.free_mfe}"
        if opt["var_rep"]: line += f"\t{t.alen}\t{t.blen}\t{t.s_ai}\t{t.s_aj}\t{t.s_bi}\t{t.s_bj}\t{t.s_alen}\t{t.s_blen}"
        pOut += f"{line}\n"
    return pOut

def plotList2(pList, pType, mType):
    pOut = ""
    for (tw,tm) in pList:
        line = f"{mType}\t{pType}1\t{tw.aSeq}-{tw.aType}\t{tw.ai}\t{tw.aj}\t{tw.bSeq}-{tw.bType}\t{tw.bi}\t{tw.bj}\t{tw.peak}\t{tw.rank}\t{tw.sp_mean}\t{tw.sp_min}\t{tw.sp_max}\t{tw.mfe}\t{tw.RNA}\t{tw.structure}\t{tw.aOrt}\t{tw.bOrt}\t{tw.a_mfe}\t{tw.a_structure}\t{tw.b_mfe}\t{tw.b_structure}\t{tw.free_structure}\t{tw.free_mfe}"
        if opt["var_rep"]: line += f"\t{tw.alen}\t{tw.blen}\t{tw.s_ai}\t{tw.s_aj}\t{tw.s_bi}\t{tw.s_bj}\t{tw.s_alen}\t{tw.s_blen}"
        pOut += f"{line}\n"
        line = f"{mType}\t{pType}2\t{tm.aSeq}-{tm.aType}\t{tm.ai}\t{tm.aj}\t{tm.bSeq}-{tm.bType}\t{tm.bi}\t{tm.bj}\t{tm.peak}\t{tm.rank}\t{tm.sp_mean}\t{tm.sp_min}\t{tm.sp_max}\t{tm.mfe}\t{tm.RNA}\t{tm.structure}\t{tm.aOrt}\t{tm.bOrt}\t{tm.a_mfe}\t{tm.a_structure}\t{tm.b_mfe}\t{tm.b_structure}\t{tm.free_structure}\t{tm.free_mfe}"
        if opt["var_rep"]: line += f"\t{tm.alen}\t{tm.blen}\t{tm.s_ai}\t{tm.s_aj}\t{tm.s_bi}\t{tm.s_bj}\t{tm.s_alen}\t{tm.s_blen}"
        pOut += f"{line}\n"
    return pOut

def compareData(tlist, name, **opt):
    ## compare data sets
    ## get mutant names
    MT_names = getMutants(tlist)
    #print(MT_names)
    ## remove intra interactions
    if not opt["var_tra"]: tlist = removeIntra(tlist)
    ## separate mutant and WT interactions
    WT_list, MT_list = separateInteractions(tlist)
    #print(len(tlist),len(WT_list),len(MT_list))
    ## all mut x WT (ignore other mut)
    fin_dict = dict()
    for (mSeg, mType) in MT_names:
        WT_comp = [tx for tx in WT_list if tx.aSeq == mSeg or tx.bSeq == mSeg]
        MT_comp = [tx for tx in MT_list if tx.aType == mType or tx.bType == mType]
        #print([f"{tx.aSeq}-{tx.aType}#{tx.bSeq}-{tx.bType}" for tx in WT_comp])
        #print("#####################")
        #print([f"{tx.aSeq}-{tx.aType}#{tx.bSeq}-{tx.bType}" for tx in MT_comp])
        df_both, eq_both = list(), list()
        tm_both, tw_both = set(), set()
        ## search for overlapping interactions
        for ew,tw in enumerate(WT_comp):
            for em,tm in enumerate(MT_comp):
                if tw.aSeq == tm.aSeq and tw.bSeq == tm.bSeq and comp([tw,tm], 0):
                    if tw.peak == tm.peak: eq_both.append((tw,tm))
                    else: df_both.append((tw,tm))
                    tm_both.add(em)
                    tw_both.add(ew)
        ## add non-overlapping interactions
        tw_only = [tw for ew,tw in enumerate(WT_comp) if ew not in tw_both]
        tm_only = [tm for em,tm in enumerate(MT_comp) if em not in tm_both]
        #print(mType,len(tw_only),len(tm_only),len(df_both),len(eq_both))
        fin_dict[(mSeg,mType)] = (tw_only, tm_only, df_both, eq_both)
    ## create test lists
    tw_all, eq_all = dict(), dict()
    for (mSeg,mType),(tw_only, tm_only, df_both, eq_both) in fin_dict.items():
        for tw in tw_only:
            tw_tup = (tw.aSeq,tw.ai,tw.aj,tw.bSeq,tw.bi,tw.bj,tw.peak)
            tw_cnt = tw_all.get(tw_tup, 0)
            tw_all[tw_tup] = tw_cnt + 1
        for (tw,tm) in eq_both:
            tw_tup = (tw.aSeq,tw.ai,tw.aj,tw.bSeq,tw.bi,tw.bj,tw.peak)
            tw_cnt = eq_all.get(tw_tup, 0)
            eq_all[tw_tup] = tw_cnt + 1
    tw_test = [tw_tup for tw_tup,tw_cnt in tw_all.items() if tw_cnt >= len(MT_names)]
    eq_test = [tw_tup for tw_tup,tw_cnt in eq_all.items() if tw_cnt >= len(MT_names)]
    ## compare test lists
    out_comp = os.path.join(opt["var_pfx"], f"{name}.tsv")
    with open(out_comp, "w") as outcomp:
        header = f"mutant\tcomparison\taSeq\tai\taj\tbSeq\tbi\tbj\tpeak\trank\tsp_mean\tsp_min\tsp_max\tmfe\tRNA\tstructure\taOrt\tbOrt\ta_mfe\ta_structure\tb_mfe\tb_structure\tfree_structure\tfree_mfe"
        if opt["var_rep"]: header = f"{header}\talen\tblen\ts_ai\ts_aj\ts_bi\ts_bj\ts_alen\ts_blen"
        outcomp.write(f"{header}\n")
        for (mSeg,mType),(tw_only, tm_only, df_both, eq_both) in fin_dict.items():
            tw_only = [tw for tw in tw_only if (tw.aSeq,tw.ai,tw.aj,tw.bSeq,tw.bi,tw.bj,tw.peak) not in tw_test]
            eq_both = [(tw,tm) for (tw,tm) in eq_both if (tw.aSeq,tw.ai,tw.aj,tw.bSeq,tw.bi,tw.bj,tw.peak) not in eq_test]
            #print(mType,len(tw_only),len(tm_only),len(df_both),len(eq_both))
            outcomp.write(plotList1(tm_only, "only_mut", mType))
            outcomp.write(plotList1(tw_only, "only_wt", mType))
            outcomp.write(plotList2(df_both, "both_df", mType))
            outcomp.write(plotList2(eq_both, "both_eq", mType))
        ## get overlapping
        for (mSeg,mType),(tw_only, tm_only, df_both, eq_both) in fin_dict.items():
            tw_only = [tw for tw in tw_only if (tw.aSeq,tw.ai,tw.aj,tw.bSeq,tw.bi,tw.bj,tw.peak) in tw_test]
            eq_both = [tw for (tw,tm) in eq_both if (tw.aSeq,tw.ai,tw.aj,tw.bSeq,tw.bi,tw.bj,tw.peak) in eq_test]
            #print("WT",len(tw_only),len(eq_both))
            outcomp.write(plotList1(tw_only, "only_wt", "WT"))
            outcomp.write(plotList1(eq_both, "both_eq", "WT"))
            break

def makeDir(**opt):
    ## create directory
    dir_name, dir_base = opt["var_pfx"], opt["var_pfx"]
    if not opt["var_ovr"]:
        i = 1
        while os.path.isdir(dir_name):
            dir_name = f"{dir_base}_{i}"
            i += 1
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    opt["var_pfx"] = dir_name
    return opt

def getTime(time_s=0, name=""):
    if time_s and name:
        time_e = time.time()-time_s
        time_e = time.strftime("%H:%M:%S", time.gmtime(time_e))
        time_c = time.strftime('%x %X')
        print(f"Status: {name} finished at {time_c} in {time_e}")
    return time.time()




################################################################################
## parser
################################################################################

if __name__ == "__main__":

    ############################################################################
    ## get time and save call
    sscript = sys.argv[0]
    start_time = time.time()
    current_time = time.strftime('%x %X')
    scall = " ".join(sys.argv[1:])
    with open(f"{sscript}.log", "a") as calllog:
        calllog.write(f"Start : {current_time}\n")
        calllog.write(f"Script: {sscript}\n")
        calllog.write(f"Call  : {scall}\n")
    print(f"Call: {scall}")
    print(f"Status: Started at {current_time}")
    ############################################################################
    ## transform string into int, float, bool if possible
    def trans(s):
        if isinstance(s, str):
            try: return int(s)
            except ValueError:
                try: return float(s)
                except ValueError:
                    if s in ["True", "False"]: return s == "True"
                    else: return s
        else: return s
    ############################################################################
    ## save documentation
    rx_text = re.compile(r"\n^(.+?)\n((?:.+\n)+)",re.MULTILINE)
    rx_oneline = re.compile(r"\n+")
    rx_options = re.compile(r"\((.+?)\:(.+?)\)")
    help_dict, type_dict, text_dict, mand_list = {}, {}, {}, []
    for match in rx_text.finditer(script_usage):
        argument = match.groups()[0].strip()
        text = " ".join(rx_oneline.sub("",match.groups()[1].strip()).split())
        argopts = {"action":"store", "help":None, "default":None, "choices":None}
        for option in rx_options.finditer(text):
            key = option.group(1).strip()
            var = option.group(2).strip()
            if var == "False": argopts["action"] = "store_true"
            if var == "True": argopts["action"] = "store_false"
            if key == "choices": var = [vs.strip() for vs in var.split(",")]
            if key == "default": var = trans(var)
            argopts[key] = var
        if argopts["default"]: add_default = f" (default: {str(argopts['default'])})"
        else: add_default = ""
        argopts["help"] = rx_options.sub("",text).strip()+add_default
        argnames = argument.split(",")
        if len(argnames) > 1:
            if argopts["default"] == None:
                mand_list.append(f"var_{argnames[1][1:]}")
            type_dict[f"var_{argnames[1][1:]}"] = argopts["default"]
            argopts["argshort"] = argnames[1]
            help_dict[argnames[0]] = argopts
        else:
            text_dict[argnames[0]] = argopts["help"]
    ############################################################################
    ## get arguments
    if text_dict["dependencies"]:
        desc = f"{text_dict['description']} (dependencies: {text_dict['dependencies']})"
    else:
        desc = text_dict['description']
    p = ap.ArgumentParser(prog=sscript, prefix_chars="-", usage=text_dict["usage"],
                          description=desc, epilog=text_dict["reference"])
    p.add_argument("-v", "--version", action="version", version=text_dict["version"])
    for argname,argopts in help_dict.items():
        argshort = argopts["argshort"]
        if argopts["choices"]:
            p.add_argument(argshort, argname,            dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"],   choices=argopts["choices"])
        else:
            p.add_argument(argopts["argshort"], argname, dest=f"var_{argshort[1:]}",\
                           action=argopts["action"],     help=argopts["help"],\
                           default=argopts["default"])
    p._optionals.title = "arguments"
    opt = vars(p.parse_args())
    ############################################################################
    ## validate arguments
    if None in [opt[mand] for mand in mand_list]:
        print("Error: Mandatory arguments missing!")
        print(f"Usage: {text_dict['usage']} use -h or --help for more information.")
        sys.exit()
    for key,var in opt.items():
        if key not in mand_list:
            arg_req, arg_in = type_dict[key], trans(var)
            if type(arg_req) == type(arg_in):
                opt[key] = arg_in
            else:
                print(f"Error: Argument {key} is not of type {type(arg_req)}!")
                sys.exit()
    ############################################################################
    ## call main function
    try:
        #saved = main(**opt)
        saved = main(opt)
    except KeyboardInterrupt:
        print("Error: Interrupted by user!")
        sys.exit()
    except SystemExit:
        print("Error: System exit!")
        sys.exit()
    except Exception:
        print("Error: Script exception!")
        traceback.print_exc(file=sys.stderr)
        sys.exit()
    ############################################################################
    ## finish
    started_time = current_time
    elapsed_time = time.time()-start_time
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    current_time = time.strftime('%x %X')
    if saved:
        with open(f"{sscript}.log", "a") as calllog,\
             open(os.path.join(saved,f"call.log"), "a") as dirlog:
            calllog.write(f"Save  : {os.path.abspath(saved)}\n")
            calllog.write(f"Finish: {current_time} in {elapsed_time}\n")
            ## dirlog
            dirlog.write(f"Start : {started_time}\n")
            dirlog.write(f"Script: {sscript}\n")
            dirlog.write(f"Call  : {scall}\n")
            dirlog.write(f"Save  : {os.path.abspath(saved)}\n")
            dirlog.write(f"Finish: {current_time} in {elapsed_time}\n")
    print(f"Status: Saved at {saved}")
    print(f"Status: Finished at {current_time} in {elapsed_time}")
    sys.exit(0)
