#!/usr/bin/env python3
# script: disruptiveMutations.py
# author: Daniel Desiro'
script_usage="""
usage
    disruptiveMutations.py -f <in_fasta> -p <out_prefix> [options]

version
    disruptiveMutations.py 0.0.1 (alpha)

dependencies
    python v3.9.7, ViennaRNA v2.5.0

description
    Reduces the interaction strength between a particular area of interest in 
    the first sequence in relation to all other sequences. Ideally, the tool 
    interrupts any potential interaction between this area and all other 
    sequences. Individual predictions are performed with the RNAcofold python 
    site-package of the ViennaRNA Package 2.5.0. Example call: python 
    disruptiveMutations.py -pfx example -fsa example.fa -pss 32 -pse 96 -thr 4 
    -ovr 

################################################################

--prefix,-pfx
    output directory and prefix for result files

--fasta,-fsa
    fasta with all sequences; the first sequence has to be the short sequence 
    snipped which has to be deoptimized with all following sequences

--positionStart,-pss
    define the starting position of the sequence part which should be 
    deoptimized (default: 0)

--positionEnd,-pse
    define the ending position of the sequence part which should be deoptimized;
    a 0 will be interpreted as the end of the sequence (default: 0)

--seed,-sed
    set the seed for random (default: 0)

--reverse,-rev
    creates reverse of each strain if set (default: False)

--complement,-cmp
    creates complements of each strain if set (default: False)

--threads,-thr
    number of threads to use for RNAcofold (default: 1)

--slice,-slc
    window size for standard window search or minimum length for a seed to be
    accepted if the seed option is specified (default: 20)

--candidateIncrease,-cdi
    increases the energy of the top candidate if there is no possible mutation
    candidate (default: 1.0)

--candidateStop,-cds
    stopping energy for deoptimization (default: -10.0)

--candidateIteration,-cdt
    tries to increase the cdi value at high iteration numbers (default: 1000)

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--dangles,-dng
    use dangling ends for foldings (default: 2) (choices: 0,1,2,3)

--noLP,-nlp
    disable lonely pairs for RNAcofold (default: False)

################################################################

reference
    D. Desirò, A. Borodavka, and M. Marz.
    "DisruptiveMutations: disrupting functional long-range RNA-RNA interactions in RNA viruses."
    In Preparation, 2022.
    https://github.com/desiro/disruptiveMutations
"""

import argparse as ap
import sys
import os
import re
import traceback
import time
import pickle
try:
    from RNA import fold_compound, cvar, CONSTRAINT_DB, CONSTRAINT_DB_DEFAULT, bp_distance
except:
    print("Error: The mandatory library ViennaRNA 2.4 is missing, please read the README.md for more information!")
    exit()
import random
from operator import itemgetter
from multiprocessing import Pool
from math import ceil, floor



################################################################################
## main
################################################################################

def main(opt):
    time_s = getTime()
    ############################################################################
    ## create output folder
    print(f"Status: Create output directory ...")
    opt = makeDir(**opt)
    time_s = getTime(time_s, f"Create output directory")
    ########################################################################
    ## read fasta file
    print(f"Status: Read fasta file ...")
    data_dict, deopt, opt = readFasta(**opt)
    time_s = getTime(time_s, f"Read fasta file")
    ############################################################################
    ## create sequence snippets
    print(f"Status: Create sequence snippets ...")
    snip_set = createSnippets(data_dict, **opt)
    time_s = getTime(time_s, f"Create sequence snippets")
    ############################################################################
    ## create snippet sets
    print(f"Status: Create snippet sets ...")
    deoptimizeSnippet(snip_set, deopt, **opt)
    time_s = getTime(time_s, f"Create snippet sets")
    return opt["var_pfx"]




################################################################################
## functions
################################################################################

def getTime(time_s=0, name=""):
    if time_s and name:
        time_e = time.time()-time_s
        time_e = time.strftime("%H:%M:%S", time.gmtime(time_e))
        time_c = time.strftime('%x %X')
        print(f"Status: {name} finished at {time_c} in {time_e}")
    return time.time()

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

def readFasta(**opt):
    ## read fasta file
    data_dict, deopt = dict(), ("", "")
    first, RNA = True, ""
    with open(opt["var_fsa"], "r") as infa:
        for line in infa:
            line = line.strip()
            if re.match(r">", line):
                if RNA:
                    if first:
                        if not opt["var_pse"]: opt["var_pse"] = len(RNA)
                        opt["var_pss"] = max(1,opt["var_pss"]-floor(opt["var_slc"]/2))
                        opt["var_pse"] = min(opt["var_pse"]+ceil(opt["var_slc"]/2),len(RNA))
                        deopt = (name, revComp(RNA[opt["var_pss"]-1:opt["var_pse"]], **opt))
                        first = False
                    else: data_dict[name] = revComp(RNA, **opt)
                name, RNA = line[1:], ""
            else:
                RNA += line
        data_dict[name] = revComp(RNA, **opt)
    createFasta(data_dict, deopt, **opt)
    return data_dict, deopt, opt

def revComp(RNA, **opt):
    ## complement dictionary, or transform DNA to RNA
    RNA = RNA.upper()
    D2Rc = {"A":"U","T":"A","U":"A","C":"G","G":"C","R":"Y","Y":"R","M":"K",\
            "K":"M","S":"W","W":"S","B":"V","V":"B","D":"H","H":"D","N":"N"}
    if opt["var_cmp"]: RNA = "".join(D2Rc[i] for i in RNA)
    else:              RNA = RNA.replace("T","U")
    if opt["var_rev"]: RNA = RNA[::-1]
    return RNA

def createFasta(data_dict, deopt, **opt):
    ## create reverse complement fasta file
    ftype = ""
    if opt["var_rev"] or opt["var_cmp"]:
        if opt["var_rev"]: ftype += "_rev"
        if opt["var_rev"]: ftype += "_cmp"
        with open(os.path.join(opt["var_pfx"], f"fasta{ftype}.fa"), "w") as outfa:
            outfa.write(f">{deopt[0]}\n{deopt[1]}\n")
            for name,RNA in sorted(data_dict.items()):
                outfa.write(f">{name}\n{RNA}\n")

def createSnippets(data_dict, **opt):
    ## create sequence snippets
    l = opt["var_slc"]
    snip_list = list()
    for name,RNA in data_dict.items():
        r = len(RNA)
        snip_list.extend([RNA[i:i+l] for i in range(0,r-l+1)])
    return list(set(snip_list))

def deoptimizeSnippet(snip_set, deopt, **opt):
    ## deoptimize snippets
    pname = os.path.basename(os.path.abspath(opt["var_pfx"]))
    dRNA = deopt[1]
    r, l, n = len(dRNA), opt["var_slc"], len(snip_set)
    random.seed(opt["var_sed"])
    ## create all combinations between deoptimization candidate and other sequences
    pool_list = [(k, k+l, dRNA[k:k+l], snip, 0, opt) for k in range(0,r-l+1) for snip in snip_set]
    ## calculate structure and energy between all combinations
    with Pool(processes=opt["var_thr"]) as p:
        p_list = p.starmap(deoptimizeSnippetMulti, pool_list)
    ## deoptimize combinations, start with strongest structure
    p_list.sort(key=itemgetter(0,1))
    it_dict, cp_dict = dict(), dict()
    with open(os.path.join(opt["var_pfx"], f"{pname}_deopt.log"), "w") as outlog:
        print("Status: original mutation start end mfe mfe_mut mfe_avg iteration")
        outlog.write(f"original\tmutation\tstart\tend\tmfe\tmfe_mut\tmfe_avg\titeration\n")
        currentMin = -100.0
        deopt_orig = ""
        while currentMin < min(opt["var_cds"],0.0):
            idx = p_list.index(min(p_list, key=itemgetter(2)))
            i, j, mfe, snip, block = p_list[idx]
            it_dict[(i,j)] = it_dict.get((i,j),0) + 1
            cp_set = cp_dict.get((i,j),set())
            ## deoptimize
            min_mfe, avg_mfe, deopt_mut = deoptimize(i, j, dRNA, snip_set, **opt) # 6 seconds
            if deopt_mut == dRNA[i:j] or deopt_mut in cp_set:
                p_list[idx] = (i, j, mfe+max(opt["var_cdi"],0.01)+(it_dict[(i,j)]//opt["var_cdt"]), snip, block)
            else:
                deopt_orig = f"{dRNA}"
                dRNA = f"{dRNA[:i]}{deopt_mut}{dRNA[i+l:]}"
                c = (idx // n) * n
                x = max(0,c-(n*(ceil(l/2)-1)))
                y = min(c+(n*(floor(l/2)+1)),len(p_list))
                p_list = reoptimize(p_list, x, y, dRNA, **opt) # 3 seconds
            cp_set.add(deopt_mut)
            cp_dict[(i,j)] = cp_set
            currentMin = min_mfe
            print(f"Status: {deopt_orig[i:j]} {deopt_mut} {i} {j} {mfe:.1f} {min_mfe:.1f} {avg_mfe:.1f} {it_dict[(i,j)]}                 ", end="\r")
            outlog.write(f"{deopt_orig[i:j]}\t{deopt_mut}\t{i}\t{j}\t{mfe:.1f}\t{min_mfe:.1f}\t{avg_mfe:.1f}\t{it_dict[(i,j)]}\n")
    for (i, j, mfe, snip, block) in p_list:
        print(i, j, mfe, snip, block)
    with open(os.path.join(opt["var_pfx"], f"{pname}_deopt.fa"), "w") as outfa:
        pss, pse = opt["var_pss"], opt["var_pse"]
        RNA = ""
        for a,b in zip(deopt[1],dRNA):
            if a == b: RNA += b.upper()
            else: RNA += b.lower()
        outfa.write(f">deopt mut:{pss}-{pse} {deopt[0]}\n{RNA}\n")

def deoptimizeSnippetMulti(i, j, dsnip, snip, block, opt):
    ## create matrix and extract candidates
    RNA = f"{dsnip}&{snip}"
    constraint = f"{'<'*(j-i)}{'>'*len(snip)}"
    mfe, pattern = doCofold(RNA, constraint, **opt)
    return (i, j, float(mfe), snip, block)

def deoptimize(i, j, dRNA, snip_set, **opt):
    ## deoptimize snip
    l, t = opt["var_slc"], floor(opt["var_slc"]/2)
    mut_list = list()
    for r in ["A","C","G","U"]:
        dsnip = f"{dRNA[i:i+t]}{r}{dRNA[i+t+1:j]}"
        pool_list = [(i, j, dsnip, snip, 0, opt) for snip in snip_set]
        with Pool(processes=opt["var_thr"]) as p:
            s_list = p.starmap(deoptimizeSnippetMulti, pool_list)
        min_mfe = min(s_list, key=itemgetter(2))[2]
        avg_mfe = sum([s[2] for s in s_list])/len(s_list)
        mut_list.append((min_mfe, avg_mfe, dsnip))
    return max(mut_list, key=itemgetter(0,1))

def reoptimize(p_list, x, y, dRNA, **opt):
    ## recalculate energies
    pool_list = [(i, j, dRNA[i:j], snip, block, opt) for i,j,mfe,snip,block in p_list[x:y]]
    with Pool(processes=opt["var_thr"]) as p:
        s_list = p.starmap(deoptimizeSnippetMulti, pool_list)
    s_list.sort(key=itemgetter(0,1))
    for idx,dx in zip(range(x,y),range(0,len(s_list))):
        p_list[idx] = s_list[dx]
    return p_list

def doCofold(RNA, constraint, **opt):
    ## do Cofold
    cvar.dangles = opt["var_dng"]
    cvar.noLP = int(opt["var_nlp"])
    fc = fold_compound(RNA)
    fc.constraints_add(constraint, CONSTRAINT_DB | CONSTRAINT_DB_DEFAULT)
    pattern, mfe = fc.mfe_dimer()
    return mfe, pattern




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
