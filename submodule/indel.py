#!/usr/bin/env python
# -*- coding: utf-8 -*-
__version__ = "1.0.0"

__author__ = "\
            Keito Watano <watano.k10.yachielab@gmail.com>,\
            Naoki Konno <naoki@bs.s.u-tokyo.ac.jp>, \
            Nozomu Yachie <nzmyachie@gmail.com>"

__date__ = "2020/3/29"

import random
import numpy as np
from Bio import Phylo
import code
import sys
import argparse
import pickle
import os
import time
import re
import math

LOGO='''
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMHB""7<~_..........__<?7TWMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMB"^_.............................._MMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMMMMMB"~.....~............~.................(MMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMMMMMMB=..............~.~-((((((((-(---~..~....~.(MMMMMMMMMMM
MMMMMMMMMMMMMMMMMMMMM#"~..........~.((JgNMMMMMMMMMMMMMMMMMMMMNgJ.~..dMMMMMMMMMMM
MMMMMMMMMMMMMMMMMM#"_........~_((gMMMMMMMMMMMMMMMMMMMMD"TWMMMMMMMMNJMMMMMMMMMMMM
MMMMMMMMMMMMMMMM5~.......~((gMMMMMMMMMMMMMMMMMMMM"WMMMr......_?TMMMMMMMMMMMMMMMM
MMMMMMMMMMMMM#^......~((gMMMMMMMMMMMMMMMMMMMMMMM_  dMMr...........?WMMMMMMMMMMMM
MMMMMMMMMMM8~.....~(+MMMMMMMMMMMMMM^ ,MMMMMMMMMM\  dMMr.~.~..........?HMMMMMMMMM
MMMMMMMMM9_.....(+MMMMMMMMMMMMMMMM#   MMMMMMMMMM{  MMMMNNg&-.-.~.......(WMMMMMMM
MMMMMMMB!...._(MMMMMMMMMMMMMMMMMMMF  .MM"    7MM{ `MMMMMMMMMMMMNg,~.~....?MMMMMM
MMMMMM%....~(MMMMMMMMMMMMMMMMMMMMMF  ,MF  .,  ,M`  MMMMMMMMMMMMMMMMa,....._TMMMM
MMMM#!....(MMM!  .MMMMMM"""MMMM@!    ,M:  d#  .M.  MMMMMMMMMMMMMMMMMMN,~....?MMM
MMM#~..~-jMMMN,...M%        (M%  .,  ,M.      JM|  d""MMMMMMMMMMMMMMMMMN,....?MM
MM#_...(MMMMMM@=?TM!   .([  ,#   Mb  .M.  ..(MHMb     dMMMMMMMMMMMMMMMMMN,~...dM
MM!...(MMMMMMMb  .M    MM]  ,#   Y^   Mp  `     dNgggMMMMMMMMMMMMMMMMMMMMMp..._M
MF...(MMMMMMMMb  .N   .MM]  ,N        MMMNJ...JMMMMMMMMMM#...(MMMMMMMMMMMMM2...J
#~..(MMMMMMMMMF  -N  `JMM]  ,MNg..Nm..MMMMMMMMMMMMMMMMMMMN...(MMMMB""WMMMMMN-..(
F..(MMMMMMMMMMF  -N, .MMMN..dMMMMMMMMMMMMMMMMMMMM"HMMMMMMN_..(MM@~...._MMMMM]...
{.~JMMMMMMMMMMF  -MMMMMMMMMMMMMMMMMMMMMMMMMMMMM#~..(MMMMMM_..(MF.._d#._MMMMMN~..
.._MMMMMMMMMMMN,.dMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMNJ(MMM#=....~(M!..(M'.(MMMMM#...
..(MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM3...d#.._dMD...(_..(#..~._(dMMMMMMF..(
..(MMMMMMMMMMMMMMMMMMMMM#7WMMMMMMMMMMMMM@....(M#...(#...(M_..(@...MMM"7WMMMM>..(
;..MMMMMMMMMMMMMMMMMMMM#_..(M_.~.M#""TMF...-JMMN...JF..(MM_..JN-..T9~..(MMMF.._d
]..JMMMMMMMMMMMMMMMMMMMMMaJMM<...@....JN-.~WMMMN.~.Jr..(M#...JMN,....(+MMM#...(M
N..~HMMMMMMMMMMMMMMMMMM@...JM{..._(-.._MNe_._TMN_..J]..(=._.~JMMMMMMMMMMM#~..(MM
M[..(MMMMMMMMMMMMMMMMMM#...JM{.~.jM]...dMMNx..(M_..dN_.._(F..JMMMMMMMMMM@_.~(MMM
MN,..(MMMMMMMMMMMMMMMMM#.~.JMr..(MMb~..dMM#~..(M_..dMN&(MMh-(JMMMMMMMMM9...(dMMM
MMN_.._WMMMMMMMMMMMMMMM#..~dM]..JMMN...JMD...(MMNagMMMMMMMMMMMMMMMMMM#!...(MMMMM
MMMN-...?MMMMMMMMMMMMMM#...dMF.~(MMM-..J]...(MMMMMMMMMMMMMMMMMMMMMM#^..~.(MMMMMM
MMMMN,..._TMMMMMMMMMMMM#~..dMb..(MMML.~dMaJMMMMMMMMMMMMMMMMMMMMMM#^....(dMMMMMMM
MMMMMMp_..._TMMMMMMMMMMN_..dMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM9~...~(dMMMMMMMMM
MMMMMMMN,.....?WMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM5!.....(dMMMMMMMMMMM
MMMMMMMMMN,.~....?TMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM9=......(JMMMMMMMMMMMMMM
MMMMMMMMMMMNa-.......(7"MMMMMMMMMMMMMMMMMMMMMMMMMMW9=~.......(JMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMNJ,...~.......~?7"T""WHHH"B"""7<_..........((gMMMMMMMMMMMMMMMMMMMM
MMMMMMMMMMMMMMMMMMNJ.-........................~...._((+MMMMMMMMMMMMMMMMMMMMMMMMM
NgygkHWmgMMMMMMMMMMMMMMNagJ---_..........~_-(((JgNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

######     ######     #######     #####     #     #    #     #    #######
#     #    #     #    #          #     #    #     #    ##   ##    #
#     #    #     #    #          #          #     #    # # # #    #
######     ######     #####       #####     #     #    #  #  #    #####
#          #   #      #                #    #     #    #     #    #
#          #    #     #          #     #    #     #    #     #    #
#          #     #    #######     #####      #####     #     #    #######

prePRESUME (PRototypE of PRESUME)
feat. indel-model
github: https://github.com/yachielab/PRESUME
'''

class Lineage(Phylo.BaseTree.Clade):
    def __init__(self, branch_length=1.0, name=None, clades=None, confidence=None, color=None, width=None, 
    seq=None, mother_name=None, ROOT=False):
        super(Lineage, self).__init__(branch_length, name, clades, confidence, color, width)

        self.mother_name =mother_name
        self.branch_length = branch_length
        self.seq = seq if ROOT else self.replication(seq)

    def time_dependent_mutation(self,c,gamma,dM): # mutation of a site (GTR-Gamma model) from PRESUME
        base={'A':0, 'C':1, 'G':2, 'T': 3}
        matrix=P(dM, gamma, Al, U)
        return random.choices(['A','C','G','T'], k=1, weights=np.array(matrix[base[c]])[0])[0] # np.matrix[x] returns matrix, then the matrix is converted to array()

    def replication(self, seq):
        time = self.branch_length
        # indel
        in_length = insertion_length(time, len(seq), rambda_I, q, r)
        del_length = deletion_length(time, len(seq), rambda_D, q, r)
        both_are_occur = insertion_length != 0 and deletion_length != 0

        if both_are_occur:
            if np.random.choice([True, False]):
                newseq = insertion(seq, in_length)
                newseq = deletion(newseq, del_length)
            else:
                newseq = deletion(seq, del_length)
                newseq = insertion(newseq, in_length)     
        else:
            newseq = insertion(seq, in_length)
            newseq = deletion(newseq, del_length)
        return newseq

####################
#  some function
####################

def disturibution_of_indel_size(q, r):
    length = np.random.negative_binomial(r,q)
    return length

def insertion_length(time, L, rambda_I, q, r):
    u = disturibution_of_indel_size(q, r)
    prob = [
        1 - math.exp(-1 * rambda_I * time * (L + 1) ), # occur insertion
        math.exp(-1 * rambda_I * time * (L + 1) ) # doesn't occur insertion
    ]
    try:
        tf = np.random.choice([True, False],p=prob)
    except ValueError:
        print(prob)
    if tf:
        return u
    else:
        return 0

def deletion_length(time, L, rambda_D, q, r):
    u = disturibution_of_indel_size(q, r)
    prob = [
        1 - math.exp(-1 * rambda_D * time * ( u - 1 + L ) ), # occur deletion
        math.exp(-1 * rambda_D * time * ( u - 1 + L ) ) # doesn't occur deletion
    ]
    try:
        tf = np.random.choice([True, False],p=prob)
    except ValueError:
        print(prob)
        print(time, rambda_D, u, L)
    if tf:
        return u
    else:
        return 0

def insertion(seq, length):
    if length == 0:
        return seq
    index = len(seq)
    position = np.random.randint(0, index)
    insertseq = "".join(np.random.choice(["A", "G", "C", "T"], length))
    newseq = seq[:position] + insertseq + seq[position:] 
    # print(seq, "=>", newseq)
    return newseq

def deletion(seq, length):
    if length == 0:
        return seq
    L = len(seq)
    start_position = max(np.random.randint(-1 * length + 1, L - 1), 0)
    end_position = min(start_position + length - 1, L - 1)
    newseq = seq[:start_position] + seq[end_position + 1:]
    return newseq

def tabulate_names(tree):
    ### from https://biopython.org/wiki/Phylo_cookbook
    for idx, clade in enumerate(tree.find_clades()):
        clade.name = str(idx)
    return tree

def topology_shaper(tree):
    internal_nodes = list(tree.find_clades(terminal=False, order='preorder'))
    terminal_nodes = list(tree.find_clades(terminal=True, order='preorder'))
    topology = {} # <mother>:[<daughter_L>, <daughter_R>]
    branch_length_dict = {} # <mother.name> : <branch_length>
    root = tree.root
    branch_length_dict["root"] = root.branch_length if root.branch_length is not None else 1.0
    for clade in internal_nodes:
        if len(clade.clades) == 2:
            mother_name = clade.name
            daughter_L_name = clade.clades[0].name
            daughter_R_name = clade.clades[1].name
            mother_branch_length = clade.branch_length if clade.branch_length is not None else 1.0
            # define topology
            topology[mother_name] = [daughter_L_name, daughter_R_name]
            topology[daughter_L_name] = []
            topology[daughter_R_name] = []
            # define branch_length
            branch_length_dict[mother_name] = clade.branch_length

        else:
            raise ValueError('wrong clades! \n # of clades {}'.format(len(clade.clades)))
            pass

    for clade in terminal_nodes:
        if clade.name != root.name:
            branch_length_dict[clade.name] = clade.branch_length if clade.branch_length is not None else 1.0


    return topology, branch_length_dict

def translate_tree(topology_dict, branch_length_dict):
    initseq=''.join([np.random.choice(['A', 'G', 'C', 'T']) for i in range(1000)])
    init_clade = Lineage(branch_length=branch_length_dict["root"], name="0", seq= initseq, ROOT=True)
    newtree = Phylo.BaseTree.Tree(init_clade)
    stack=[init_clade]
    time_zero=time.time()
    cnt = 0
    while (stack !=[]):
        clade=stack.pop()
        node_name=clade.name
        mother_seq=clade.seq
        if len(topology_dict[node_name]) ==2:
            children = [
                Lineage(branch_length = branch_length_dict[topology_dict[node_name][0]], name=str(topology_dict[node_name][0]), seq=mother_seq),
                Lineage(branch_length = branch_length_dict[topology_dict[node_name][1]], name=str(topology_dict[node_name][1]), seq=mother_seq)
            ]
            clade.clades.extend(children)
            stack.extend(children)
            cnt += 1
    return newtree

def P(t, gamma, Al, U): # return transition matrix
    exp_rambda = np.diag(
        np.array([np.exp(Al[0]*t*gamma),
        np.exp(Al[1]*t*gamma),
        np.exp(Al[2]*t*gamma),
        np.exp(Al[3]*t*gamma)]))
    return np.dot(np.dot(U,exp_rambda),np.linalg.inv(U))

def PRESUME_INDEL(args):
    print(LOGO)
    argname = dir(args)
    print(argname)
    while argname != []:
        arg = argname.pop()
        find_dunder = re.match('^_', arg)
        if not find_dunder:
            if arg in {"tree", "ram_I", "ram_D", "prop_q", "int_r"} and args.__dict__[arg] is not None:
                name=arg
                val=args.__dict__[arg]
                print("PRESUMEnwk2fa: LOADED {1} as {0}".format(name, val))
    print("PRESUMEnwk2fa: other arguments are ignored(using default).")

    global rambda_I, rambda_D, q, r
    rambda_I = args.ram_I
    rambda_D = args.ram_D
    q = args.prop_q
    r = args.int_r

    ########################################################################################
    # substitution rate matrix (GTR-Gamma model) -> define substitution matrix function 
    ########################################################################################
    if args.gtrgamma is None or str(args.gtrgamma) == "default":
        MODEL="GTR{0.927000/2.219783/1.575175/0.861651/4.748809/1.000000}+FU{0.298/0.215/0.304/0.183}+G4{0.553549}"
    else:
        MODEL=str(args.gtrgamma)
    models_str=(MODEL).split("+")
    abcdef=(models_str[0].split("{"))[1].split("}")[0].split("/")
    piACGT=(models_str[1].split("{"))[1].split("}")[0].split("/")
    gamma_str=(models_str[2].split("{"))[1].split("}")[0].split("/")
    a1 = float(abcdef[0]) # A <-> C
    a2 = float(abcdef[1]) # A <-> G
    a3 = float(abcdef[2]) # A <-> T
    a4 = float(abcdef[3]) # C <-> G
    a5 = float(abcdef[4]) # C <-> T
    a6 = float(abcdef[5]) # G <-> T
    piA = float(piACGT[0])
    piC = float(piACGT[1])
    piG = float(piACGT[2])
    piT = float(piACGT[3])
    if (abs(piA+piC+piG+piT-1)>0.001): print("error piA+piC+piG+piT not equal to 1!"); sys.exit(1)
    R=np.matrix([[-(a1*piC+a2*piG+a3*piT),a1*piC,a2*piG,a3*piT],[a1*piA,-(a1*piA+a4*piG+a5*piT),a4*piG, a5*piT],[a2*piA,a4*piC,-(a2*piA+a4*piC+a6*piT),a6*piT],[a3*piA,a5*piC,a6*piG,-(a3*piA+a5*piC+a6*piG)]]) # substitution rate matrix
    print("PRESUMEnwk2fa:substitution rate matrix:")
    print(R)
    global Al, U
    Al, U = np.linalg.eig(R) # Al: eigen values (np.array), U: eigen vectors matrix : R = U * diag(Al) * U^(-1) 
    # model of site heterogeneity: get relative substitution rate gamma for each site
    MEAN=args.m
    shape = float(gamma_str[0]) #shape of gamma distribution
    global gamma
    gamma = np.random.gamma(shape,MEAN/shape,1000) # mean is args.m
    ####################################################################################################

    # processing file
    THIS_FILE_PATH=os.path.dirname(os.path.abspath(__file__)) # absolute path of PRESUME directory
    NWK = args.indel

    ### setup directory
    OUTDIR=args.output
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)
    os.chdir(OUTDIR)
    os.makedirs("PRESUMEout_from_tree", exist_ok=True)
    os.chdir("PRESUMEout_from_tree")

    ### read ###
    tree = Phylo.read(NWK, "newick")
    tree = tabulate_names(tree)

    ### translate Phylo.BaseTree.Tree into Lineage
    topology_dict, branch_length_dict = topology_shaper(tree)
    lineage_tree = translate_tree(topology_dict, branch_length_dict)
    if len(tree.get_terminals()) != len(lineage_tree.get_terminals()):
        raise ValueError('something went wrong!!!')

    fasta_data = [ [item.name ,item.seq] for item in lineage_tree.get_terminals()]
    with open("ANS_SEQ.fasta", "w") as writer:
        for item in fasta_data:
            writer.write(">{0}\n{1}\n".format(item[0], item[1]))
    # [ print(">{0}\n{1}".format(idx, item)) for idx, item in enumerate(seq)]
    Phylo.write(lineage_tree, "test.nwk", "newick")

if __name__=="__main__":
    ######## processing some args ########
    parser = argparse.ArgumentParser(description='sequene_generator_from_nwk')
    parser.add_argument("--ram_I", help="rate of insertion. (default=0)", type=float, default=0)
    parser.add_argument("--ram_D", help="rate of insertion. (default=0)", type=float, default=0)
    parser.add_argument("-q", help="q of negative binomal distribution as indel size distribution. (default=0.1, must be float)", type=float, default=0.1)
    parser.add_argument("-r", help="r of negative binomal distribution as indel size distribution. (default=1, must be integer.)", type=float, default=1)
    parser.add_argument("--output", help="output folder (default:current directory)",type=str, default=os.getcwd())
    args = parser.parse_args()

    # start simulation
    PRESUME_nwk2fa(args)

    # code.InteractiveConsole(globals()).interact()


    
