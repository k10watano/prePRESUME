#!/usr/bin/env python
# -*- coding: utf-8 -*-
__version__ = "1.0.0"

__author__ = "\
            Keito Watano <watano.k10.yachielab@gmail.com>,\
            Naoki Konno <naoki@bs.s.u-tokyo.ac.jp>, \
            Nozomu Yachie <nzmyachie@gmail.com>"

__date__ = "2020/3/29"
import subprocess
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
import shutil

LOGO='''
                               ....((++++++(-....
                         ..(zluggggHMHMMMMMNHmggezzO+..
                      .(<u&kMMHHHH@HMWWMHHHHHHHHHMMMNgsrO(.
                   ._(JMMHHH@H@@H@@`.MMQ@H@@H@@M"4HHHHHMNayz&,`    `    `    `
                 ..JMMH@H@H@HHH@HH] MHH@HH@HHM".b,H@@H@HHHMMNyzi.
              .-(MMHH@HH@HH@MHHH@HHa..(.HHHM".ma,.H@HH@H@@HHHMNmzn.  `    `
    `       .((MH@H@HH@HH@HbJH@H@H@HH@HH@@H@H@HHN.MHH@HHb.a,4HHHMmun.        `
       ` `.(jMH@HH@H@H@H@HHM+MH@HH@H@HH@HH@HH@H@HH@H@HPHH,?HH@H@HHMmzi
  `      (jdHH@H@HH@H@HH@H@HHH@HMHH@HH@HH@HH@H@H@HH@H@M+?".H@H@H@@HMNyw.  `
       .?dMHH@HH@H@HHH@HH@H@H/M#J#J@J#W@HH@HH@HH@H@HHHHHH@HHH@HMYYMHHMmz,    `
      .=qMH@H@HHH@H@H@H@HH@H@},NjbJCd1H8jMMH@HH@HH@H@H@H@H@HM4Hb.N,W@HMRz>
    `.ldMH@H@H@@HH@HH@HH@HH@H) ,H@HH#dhd6j#WH@HH@HH@H@H@HH@HM,WN HNHH@HM#v>
    .lqMH@HHHHH@HH@HH@HH@HH@H)  (HM=.H@HBudYWH@H@H@HHH@HH@H@H@m(d@H@HH@HMKr- `
 `  zAMH@H@@H@HH@HH@HH@HMY"WH\   _  HHHHMAdHYMH@HH@H@HH@HH@HHH@HH@HH@HH@HMsO.
   ,tMHHH@HH@H@H@H@H@H@HH@h,             ..dMMHH@HH@H@H@H@HH@HH@HHY"! ,HHHNtl
   rdMH@HHH@HHH@H@HH@HH@HMY=         ..gHMkQgg@HH@HH@H@HH@H@H@HHN.gHHH@H@@MKl.
  ,tdH@H@H@H@H@HHH@HHMY"`            TMH@#SQg+H@HH@HH@HH@HHH@H@H@HH@HHH@HHHN={
  (OMHH@H@HH@HH@H@HH+(JxzgQHH!     .+..7MYHaJdHH@H@H@HH@H@H@HH@HH@HH@@HH@HHMz<
  ldM@HH@HH@H@HH@H@HHmgW9dHHF.JN.` ,HHMMZWm&dH@HH@H@HH@HH@HH@HH@HH#=~..7MH@MR=
  ldHH@HH@HH@H@HH@H@H@ag9udWjMHHN. ,M#WeTmJWH@H@HH@HH@H@HH@H@H@H@M dH@Hb,HHMb=
  =dHH@H@H@HH@H@HH@HHH@QM5j#j#dBMh ,rWpdmJH@HH@H@HH@HH@H@HH@HH@HH@o?WMY^.H@MDl
  1zMHH@HH@HH@HH@HH@H@HHMHMuM1Hj@Jh,MxMNHH@H@HH@H@H@H@HH@H@HH@HH@HHNgJ+HH@HM0z
  (?M@HH@HH@HH@HH@H@HH@HH@HHMHNdNJHNHHH@HH@HH@HH@HHH@HH@HH@H@H@HH@HH@HH@HH@Mz>
  .>dH@H@H@H@H@H@H@H@HH@HH@HH@HHHH@HH@HH@HH@HH@HH@H@HH@H@HH@HH@H@MH@HH@H@HH#l:
   <?MH@HH@HH@HH@HH@H@HH@H@H@H@H@HH@H@H@HMMMMHWYY9YYYTTTTXQgHMHHM...7YMH@HM6v
   .>dHH@HH@HH@HH@HH@H@HH@HHHMMHY9Y71<;;;;;;;;;++++gQHHHHH@HHH@HHH#.M[.HHH#l:
    (;MH@H@H@H@H@H@HHMHYYT1;;;;;;;:;;;;+j+ggHH@HHHHHHH@H@HH@H@H@H@b?".H@HMO<
     <?MHHH@HHMH9YC;;;;;;;;:;;;+jjggMH@HH@HH@HH@@H@@H@H@H@HH@HH@HH@HH@HHM6<`
      <?MYYC<;;;;;;;;:;;<jjggMHH@HH@HH@H@H@HH@HH@HH@HH@HH@H@#^7MH@H@HH@M6<`
       <?z;;:;;;;;<jjgHH@HH@HH@HH@H@H@HHH@H@H@H@HH@H@HH@HH@Y.MN (WHH@HM=~`
        ~:1;;jjgHHH@HH@HH@HH@H@H@H@HH@H@HH@HHH@HH@HH@H@H@HHH@P.MF.H@M#!_
         -:?MHH@H@HH@HH@HH@HH@H@HH@HHH@H@HH@H@HH@H@HH@MH@H@HH@HLdH@M5`       `
           ~:THH@H@H@H@H@H@H@HHH@HH@@HH@H@HH@H@HH@H@HM,,MH@HH@HHHMD_
             ~<THHH@HHH@H@HH@H@HH@HH@HH@HH@HH@H@H=HH@HHx.WH@HH@H9~_
    `          ~:7MH@H@HHH@HH@H@H@HH@H@HH@H@Y7MM9.@HH@H@b.MHHM5_~`      `    `
                 _<?TMH@H@H@HH@H@H@HH@HH@HHH@h.T@.H@HH@HHHMB1<!
                    ~<<TBMHH@HHHH@H@HH@H@H@HH@Ha.JH@HHMH91<!`             `
                       `!<;?TBMMHH@H@HH@H@H@HHHHMMH96zv?!               `    `
                            _!<<>?zOTUUUU99UUwttv<?!`                `
    `                                 `````                   `    `
             `                                           `              `    `
######     ######     #######     #####     #     #    #     #    #######
#     #    #     #    #          #     #    #     #    ##   ##    #
#     #    #     #    #          #          #     #    # # # #    #
######     ######     #####       #####     #     #    #  #  #    #####
#          #   #      #                #    #     #    #     #    #
#          #    #     #          #     #    #     #    #     #    #
#          #     #    #######     #####      #####     #     #    #######
DEBUG MODE of PRESUME
feat. cassiopeia-model
github: https://github.com/yachielab/PRESUME
'''

# absolute path of python3 directory
PYTHON3 = (((
    subprocess.Popen('which python3', stdout=subprocess.PIPE, shell=True)
    .communicate()[0])
    .decode('utf-8'))
    .split('\n')
    )[0]

# absolute path of PRESUME directory
dir_str = os.path.dirname(os.path.abspath(__file__)).split("/")
dir_str.pop(-1)
PRESUME = "/".join(dir_str)

#   for Exception
class PresumeException(Exception):
    """Base class for exceptions in this module."""
    pass


class ExtinctionError(PresumeException):
    def __init__(self):
        self.message = "PRESUME: An Error Occured !!!\n \
                        doubling time of initial SEQ is 0 !"

    def __str__(self):
        return self.message


class NoChildError(PresumeException):
    def __init__(self, seqid):
        self.seqid = seqid

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        No alive children of id:{}
        """.format(self.seqid)
        return message


class CreateNewickError(PresumeException):
    def __init__(self, e):
        self.occured_error = e

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        Something went wrong when we create newick...
        ERROR MESSAGE:{}
        """.format(self.occured_error)
        return message


class UpperLimitExceededError(PresumeException):
    def __init__(self, u):
        self.upper_limit = u

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        Upper limit exceeded !
        upper limit:{}
        """.format(self.upper_limit)
        return message


class TopologyError(PresumeException):
    def __init__(self, tree):
        self.tree = tree

    def __str__(self):
        num_of_internal_nodes = len(self.tree.get_terminals())
        num_of_terminal_nodes = len(self.tree.get_terminals())
        message = """
        PRESUME: An Error Occured !!!
        PRESUME generated invalid topology tree !
        # of internal nodes:{0}
        # of terminal nodes:{1}
        """.format(num_of_internal_nodes, num_of_terminal_nodes)
        return message


class OutputError(PresumeException):
    def __init__(self, fa_cnt, tip_cnt):
        self.fa_cnt = fa_cnt
        self.tip_cnt = tip_cnt

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        number of SEQs in fasta and in newick is different!
        # of sequences:{0}
        # of tip labels:{1}
        """.format(self.fa_cnt, self.tip_cnt)
        return message


class SimulationFailureError(PresumeException):
    def __init__(self, take_count):
        self.take_count = take_count

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        Recursive limit exceeded !
        Take:{}
        """.format(self.take_count)
        return message
    
class ExceptionOfCodeOcenan(PresumeException):
    def __init__(self):
        self.message="""
        PRESUME: Because Code Ocean does not support Univa Grid Engine, 
        distributed computing mode(the option of --qsub) is disabled.
        If you want to use it, please install PRESUME from http:/github.com/yachielab/PRESUME.
        """

    def __str__(self):
        return self.message


#   for Simulation
class SEQ():
    # idM & mseq means id & sequence of mother SEQ, respectively.
    def __init__(self, SEQid, idM, mseq, CVM, rM, dM, tM, CV=False):
        self.id = SEQid  # for example: 0,1,2,...
        self.idM = idM
        self.CV = max(np.random.normal(CVM, CVM*alpha), 0)  # should be > 0
        if CV:
            self.r = self.growing_rate_dist(rM, rM*self.CV)
        else:
            self.r = self.growing_rate_dist(rM, self.CV)
        self.is_alive = self.r >= 0 and np.random.rand() > e

        if(self.is_alive):
            if self.r == 0:
                self.d = float("inf")
            else:
                self.d = 1/self.r

            self.t = tM + self.d  # time t of doubling of this SEQ
            self.seq = CAS_like_replication(mseq)
            # self.mutation_rate = compare_sequences(str(mseq), self.seq)

    def growing_rate_dist(self, mu, sigma):
        growing_rate = np.random.normal(mu, sigma)
        if growing_rate <= 0:
            growing_rate = -1
        return growing_rate

    # receive mother SEQ sequence, introduce mutations,
    # return daughter sequence.
    def daughterseq(self, seq):
        dseq = []
        for i in range(L):
            dseq = dseq + "A"
        return dseq

def CAS_like_replication(sequence):
    daughter = []
    for character in sequence:
        occur_mutation = np.random.choice(a=[True, False], p=[prob_of_edit, 1 - prob_of_edit])
        if character == 0 and occur_mutation:
            newcharacter = np.random.randint(1, state)
        else:
            newcharacter = character
            
        occur_dropout = np.random.choice(a=[True, False], p=[prob_of_dropout, 1 - prob_of_dropout])
        if occur_dropout:
            newcharacter = -1       
        daughter.append(newcharacter)
    return daughter

def compare_sequences(motherseq, myseq):
    if len(motherseq) != len(myseq):
        print('Sequences should be same length')
        return None
    score = []  # mutation: ["1"], replication: ["0"]
    for position in range(len(motherseq)):
        its_a_match = motherseq[position] == myseq[position]
        score += [int(0)] if its_a_match else [int(1)]
    return score


def argument_saver(args):
    argname = dir(args)
    arg_pair = []
    csv_L1 = []
    csv_L2 = []
    while argname != []:
        arg = argname.pop()
        find_dunder = re.match('^_', arg)
        if not find_dunder:
            if args.__dict__[arg] is not None and not(arg in {"load", "save"}):
                csv_L1.append(arg)
                csv_L2.append(args.__dict__[arg])
        arg_pair = [csv_L1, csv_L2]

        with open("args.csv", "wt") as fout:
            csvout = csv.writer(fout)
            csvout.writerows(arg_pair)


# when a child of SEQ of <SEQid> is dead,
# kill ancestors whose any descendant is dead.
def death(Lineage, SEQid):
    index = SEQid
    while(index != 0
            and Lineage[Lineage[index][1]] == ["dead"]
            and Lineage[Lineage[index][2]] == ["dead"]):
        idM = Lineage[index][0]
        Lineage[index] = ["dead"]
        index = idM


# if all SEQ died, create a file "all_SEQ_dead.out"
def all_dead(idANC):
    with open("all_SEQ_dead.out", 'w') as handle:
        handle.write(str(args.idANC)+"\n")


# count number & length of sequences in a specified fasta file
def count_sequence(in_fname):
    with open(in_fname, "r") as origin:
        input_itr = origin.readlines()
    return len(input_itr)


def create_newick(Lineage):
    try:
        init_clade = Phylo.BaseTree.Clade(name="0")
        tree = Phylo.BaseTree.Tree(init_clade)
        stack = [init_clade]
        list_of_dead = []
        while(stack != []):
            clade = stack.pop()
            SEQid = int(clade.name)
            if(len(Lineage[SEQid]) == 3):  # if the SEQ has children
                # by this, every internal node will have 2 children
                while(True):
                    both_children_are_alive = (
                        Lineage[Lineage[SEQid][1]] != ["dead"]
                        and Lineage[Lineage[SEQid][2]] != ["dead"])
                    both_children_are_dead = (
                        Lineage[Lineage[SEQid][1]] == ["dead"]
                        and Lineage[Lineage[SEQid][2]] == ["dead"])
                    if(both_children_are_alive):
                        children = [
                            Phylo.BaseTree.Clade(
                                name=str(Lineage[SEQid][i]),
                                ) for i in [1, 2]
                            ]
                        clade.clades.extend(children)
                        stack.extend(children)
                        break
                    elif(both_children_are_dead):
                        raise NoChildError(SEQid)
                        return
                    # if either of the children is dead, renew SEQid
                    # to be the id of alive child
                    elif(Lineage[Lineage[SEQid][1]] != ["dead"]):   
                        list_of_dead.append(Lineage[SEQid][2])
                        SEQid = Lineage[SEQid][1]
                    elif(Lineage[Lineage[SEQid][2]] != ["dead"]):
                        list_of_dead.append(Lineage[SEQid][1])
                        SEQid = Lineage[SEQid][2]
                    # if the clade of renewed SEQid is a terminal,
                    # rename the clade
                    if(len(Lineage[SEQid]) == 1):
                        clade.name = str(SEQid)
                        break

        # only in downstream tree of distributed computing,
        # the name of terminals is <upstream id>_<downstream id>
        if (idANC is not None):
            for clade in tree.get_terminals():
                clade_name_prefix = str(hex(idANC)).split("x")[1]
                clade_name_suffix = str(hex(int(clade.name))).split("x")[1]
                new_clade_name = "{0}_{1}". \
                    format(clade_name_prefix, clade_name_suffix)
                clade.name = new_clade_name

        # for error check
        if (len(tree.get_nonterminals()) + 1 != len(tree.get_terminals())):
            raise TopologyError(tree)
            return

        # file write in default
        if (idANC is None):
            Phylo.write(tree, "PRESUMEout.nwk", 'newick')

        else:
            # file write in case of downstream lineage
            Phylo.write(tree, "SEQ_"+str(idANC)+".nwk", 'newick')
        
        return len(tree.get_terminals()), tree, list_of_dead

    except Exception as e:
        raise CreateNewickError(e)


def survey_all_dead_lineages(Lineage):
    try:
        command = "cat intermediate/DOWN/*/all_SEQ_dead.out \
            > intermediate/all_dead.out 2> /dev/null; \
            rm intermediate/DOWN/*/all_SEQ_dead.out 2> /dev/null"
        subprocess.call(command, shell=True)

    except Exception as e:
        print(e)
        print("no lineages extinct")
        return

    with open("intermediate/all_dead.out", 'r') as handle:
        lines = handle.readlines()
        for line in lines:
            dead_SEQ_id = int(line.split("\n")[0])
            idM = Lineage[dead_SEQ_id][0]
            Lineage[dead_SEQ_id] = ["dead"]
            death(Lineage, idM)


def CombineTrees():
    top_tree = Phylo.read('PRESUMEout.nwk', 'newick')
    terminals = top_tree.get_terminals()
    TOP_CELL_CNT = len(terminals)
    index = 0
    for tip in terminals:
        index += 1
        if(index % 100 == 0):
            message = "tree assembly:{}percent finished". \
                format(index * 100 / (TOP_CELL_CNT))
            sys.stdout.write(message)
            sys.stdout.flush()
        bottom_tree_filepath = \
            "intermediate/DOWN/esu_{0}/PRESUMEout/SEQ_{0}.nwk".format(tip.name)
        bottom_tree = Phylo.read(bottom_tree_filepath, 'newick')
        tip.clades.extend(bottom_tree.clade.clades)
        newick_from = \
            "intermediate/DOWN/esu_{0}/PRESUMEout/SEQ_{0}.nwk".\
            format(tip.name)
        newick_to = \
            "intermediate/DOWN/esu_{0}/PRESUMEout/SEQ_{0}_attached.nwk".\
            format(tip.name)
        shutil.move(newick_from, newick_to)
    Phylo.write(top_tree, 'PRESUMEout_combined.nwk', 'newick')
    return len(top_tree.get_terminals())


def progress_bar(pbar, current_value):
    if pbar.n <= current_value:
        pbar.n = current_value
    else:
        pbar.n = 0
    pbar.refresh()


def sequence_writer(name, seq, file_name, overwrite_mode):
    data = [str(name)]
    for character in seq:
        data.append(str(character))

    if overwrite_mode:
        writer_mode = "a"
    else:
        writer_mode = "w"
    with open(file_name, writer_mode) as writer:
        line = "\t".join(data)
        writer.writelines(line)


def sequences_writer(list_of_sequence, file_name, no_header=False):
    # writedata = [
    # [name, character_1, character_2, ..., character_n] # header row
    # [name, 0, 0, ..., 0] # data row
    # ...
    # ]
    writedata = []

    # generate header row
    if not no_header:
        character_items = len(list_of_sequence[0].seq)
        header = ["name"]
        for index in range(character_items):
            header.append("r{}".format(index))
        header.append("\n")
        writedata.append(header)

    # generate data row
    for sequence in list_of_sequence:
        name = str(sequence.id)
        seq = sequence.seq
        entry = [name]
        for character in seq:
            entry.append(str(character))
        entry.append("\n")
        writedata.append(entry)

    with open(file_name, "w") as writer:
        for data in writedata:
            line = '\t'.join(data)
            writer.writelines(line)


def jobscript_writer(esu, serial_number, args, timelimit):
    PATH = (((
        subprocess.Popen('echo $PATH', stdout=subprocess.PIPE,
                            shell=True)
        .communicate()[0])
        .decode('utf-8'))
        .split('\n'))[0]

    LD_LIBRARY_PATH = (((
        subprocess.Popen('echo $LD_LIBRARY_PATH', stdout=subprocess.PIPE,
                        shell=True)
        .communicate()[0])
        .decode('utf-8'))
        .split('\n'))[0]

    # create intermediate file
    intermediate_file_path = \
        "intermediate/fasta/{}.seq".\
        format(str(esu.id))
    sequence_writer(esu.id, esu.seq, intermediate_file_path, True)

    with open("intermediate/shell/esu_"+str(serial_number)+".sh", 'w') as qf:
        qf.write("#!/bin/bash\n")
        qf.write("#$ -S /bin/bash\n")
        qf.write("#$ -cwd\n")
        qf.write("PATH={}\n".format(PATH))
        qf.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))
        qf.write("mkdir intermediate/DOWN/esu_"+str(esu.id)+"\n")
        qf.write("cd intermediate/DOWN/esu_"+str(esu.id)+"\n")
        qf.write("pwd\n")

        # divide until time point of (2 * timelimit)
        if args.CV:
            python_command = PYTHON3 + " " + PRESUME + "/PRESUME.py "\
                "--monitor " + str(2*timelimit)\
                + " -L "+str(args.L)\
                + " -f "+"../../../fasta/"+str(esu.id)+".seq"\
                + " -d "+str(esu.d)\
                + " -s "+str(args.s)\
                + " -T "+str(args.T)\
                + " -e "+str(args.e)\
                + " -u "+str(args.u)\
                + " --idANC "+str(esu.id)\
                + " --tMorigin "+str(esu.t-esu.d)\
                + " --CV"\
                + " --seed " + str(np.random.randint(0, args.r))\
                + " --cassiopeia"\
                + " --probedit " + str(args.probedit)\
                + " --probdout " + str(args.probdout)\
                +"\n"
        else:
            python_command = PYTHON3 + " " + PRESUME + "/PRESUME.py "\
                "--monitor " + str(2*timelimit)\
                + " -L "+str(args.L)\
                + " -f "+"../../../fasta/"+str(esu.id)+".seq"\
                + " -d "+str(esu.d)\
                + " -s "+str(args.s)\
                + " -T "+str(args.T)\
                + " -e "+str(args.e)\
                + " -u "+str(args.u)\
                + " --idANC "+str(esu.id)\
                + " --tMorigin "+str(esu.t-esu.d)\
                + " --seed " + str(np.random.randint(0, args.r))\
                + " --cassiopeia"\
                + " --probedit " + str(args.probedit)\
                + " --probdout " + str(args.probdout)\
                +"\n"

        qf.write(python_command)


def seq_reader(filepath):
    with open(filepath, "r") as reader:
        raw_row=reader.readline()
    sequence = raw_row.split("\t")[1:]
    return sequence


# main
def PRESUME_CAS(args):
    '''
    main function
            commandline argument: "python3 PRESUME.py <timelimit> <options...>"
    '''
    print(LOGO)
    global alpha, e, idANC, state, prob_of_edit, prob_of_dropout
    state = args.state
    prob_of_edit = args.probedit
    prob_of_dropout = args.probdout
    dorigin = args.d
    alpha = args.a
    growing_rate_origin = 1 / args.d
    sigma_origin = args.s
    T = args.T
    e = args.e
    idANC = args.idANC
    timelimit = args.monitor

    # initialize the pseudo-random number generator
    if args.seed is not None:
        seed = args.seed
    elif args.seed == "rand":
        seed = np.random.randint(0, args.r)
    elif args.seed is None:
        seed = 0
    else:
        seed = int(args.seed)
    np.random.seed(int(seed))
    
    # setup directory
    OUTDIR = args.output
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)
    os.chdir(OUTDIR)
    os.makedirs("PRESUMEout", exist_ok=True)
    os.chdir("PRESUMEout")

    # initital sequence specification
    if (args.f is not None):
        initseq = seq_reader(args.f)
    else:
        initseq = [0] * args.L
    
    # note: tqdm should be ver. 4.19 (conda install tqdm=4.19)
    if args.bar:
        from tqdm import tqdm
    C = args.n  # C : number of SEQs you want to get (default == None)
    if(C is not None):  # if C is specified
        if(args.qsub):  # for distributed computing
            C = int(C**(1/2))  # update C

    if args.bar:
        pbar = tqdm(range(C))

    # for safety: program definitely stops when the number of SEQs
    # exceeds UPPER_LIMIT
    UPPER_LIMIT = args.u
    delta_timelimit = timelimit
    inittimelimit = timelimit

    # next id
    i = 0

    # current existing SEQ
    SEQqueue = []

    # Lineage[j]=[k,l,m] means an SEQ whose id is j is a daughter of Lineage[k]
    # and is a mother of Lineage[l], Lineage[m]
    Lineage = [[]] * UPPER_LIMIT

    # DEBUG: for gathering mutation rate of each site
    # if args.debug:
    #     mut_rate_log = {}

    # First of all, there exits only 1 SEQ.
    if args.CV:
        SEQqueue.append(
            SEQ(i, -1, initseq, sigma_origin, growing_rate_origin,
                dorigin, args.tMorigin, True))
    else:
        SEQqueue.append(SEQ(i, -1, initseq, sigma_origin, growing_rate_origin,
                            dorigin, args.tMorigin, True))
    SEQqueue[0].seq = initseq
    i += 1
    c = 1  # current number of SEQs

    # SEQs propagation
    while(True):
        if (c == 0):
            all_dead(args.idANC)
            try:
                os.remove("ancestoral_sequences.fasta")
            except FileNotFoundError:
                pass
            if args.save:
                argument_saver(args)
            return 1

        # SEQs divide until t (time) of all of them exceed timelimit
        k = 0  # iterator of SEQqueue

        # while there remains some SEQs whose time is less than timelimit
        while(k < c):

            if (not SEQqueue[k].is_alive):
                # Note: Variable name "esu" means Evolutionally Sequence Unit
                # as container of SEQ object
                esu = SEQqueue.pop(k)
                Lineage[esu.id] = ["dead"]
                if(esu.id != 0):
                    death(Lineage, esu.idM)
                c -= 1
                if args.bar:
                    progress_bar(pbar, c)

            elif(SEQqueue[k].t < timelimit):
                esu = SEQqueue.pop(k)
                # save ancestoral sequences
                if args.viewANC:
                    sequence_writer(esu.id, esu.seq, "ancestoral_sequences.seq", True)
                # duplication
                if args.CV:
                    daughter = [SEQ(i, esu.id, esu.seq, esu.CV,
                                    esu.d, esu.r, esu.t, True),
                                SEQ(i+1, esu.id, esu.seq, esu.CV,
                                    esu.d, esu.r, esu.t, True)]
                else:
                    daughter = [SEQ(i, esu.id, esu.seq, esu.CV,
                                    esu.d, esu.r, esu.t),
                                SEQ(i+1, esu.id, esu.seq, esu.CV,
                                    esu.d, esu.r, esu.t)]

                SEQqueue.extend(daughter)

                # if args.debug:
                #     for sister in daughter:
                #         if sister.is_alive:
                #             mut_rate_log[sister.id]=sister.mutation_rate

                # [<mother>, <daughter1>, <daughter2> ]
                Lineage[esu.id] = [esu.idM, i, i+1]
                try:
                    Lineage[i] = [esu.id]
                    Lineage[i+1] = [esu.id]  # Terminal ESUs

                except IndexError:
                    print("UPPER LIMIT EXCEEDED!!!")
                    return 1

                i += 2
                c += 1
                if args.bar:
                    progress_bar(pbar, c)

            else:
                k += 1

            if(c > UPPER_LIMIT):  # for safety
                raise UpperLimitExceededError(UPPER_LIMIT)
                return 1

        if (c == 0):
            all_dead(args.idANC)
            try:
                os.remove("ancestoral_sequences.fasta")
            except FileNotFoundError:
                pass
            if args.save:
                argument_saver(args)
            return 1

        # if the required number of SEQs is not specified (default)
        if(C is None):
            break
        # if the required number of SEQs specified
        else:
            if(c < C):
                delta_timelimit = inittimelimit/(timelimit/inittimelimit)
                timelimit += delta_timelimit
            else:
                print("Number of generated sequences reached "+str(C))
                break

    # output initial sequence
    sequence_writer("root", initseq, "root.seq", False)

    # in case of distributed computing
    if args.qsub:
        # preparation for qsub
        os.mkdir("intermediate")
        os.mkdir("intermediate/DOWN")
        os.mkdir("intermediate/fasta")
        os.mkdir("intermediate/shell")

        # generate shell script
        itr = 0
        for esu in SEQqueue:
            itr += 1
            jobscript_writer(esu, itr, args, timelimit)
        
        del(SEQqueue)
        # submit job script to grid engine
        print("\ncreating bottom trees by qsub ...")
        submit_command = "qsub -l d_rt=1:00:00 -l s_rt=1:00:00 -sync y -t 1-{0} \
            {1}/exe_PRESUME.sh &> intermediate/qsub.out".\
            format(str(itr), PRESUME)
        subprocess.call(submit_command, shell=True)

        # finalize
        # remove extinct downstream lineages
        survey_all_dead_lineages(Lineage)
        tip_count, returned_tree, list_of_dead = create_newick(Lineage)
        # if args.debug:
        #     mut_rate_log_writer(mut_rate_log, list_of_dead)

        command = "cat PRESUME.e*.* > intermediate/err; \
                cat PRESUME.o*.* > intermediate/out;"
        subprocess.call(command, shell=True)

        # processing seq files
        if args.f is None:
            command = "cat intermediate/DOWN/*/PRESUMEout/PRESUMEout.seq \
                    > PRESUMEout.seq"
            subprocess.call(command, shell=True)  # combine fasta
        
        # adding header
        with open("PRESUMEout.seq","r") as reader:
            temp_data = reader.readlines()
        character_items = args.L
        header = ["name"]
        for index in range(character_items):
            header.append("r{}".format(index))
        header.append("\n")
        header_line="\t".join(header)
        writedata = [header_line] + temp_data
        with open("PRESUMEout.seq","w") as writer:
            temp_data = writer.writelines(writedata)

        fa_count = count_sequence("PRESUMEout.seq")
        tip_count = CombineTrees()  # Combine trees
        shutil.move("PRESUMEout.nwk", "intermediate")
        os.rename("PRESUMEout_combined.nwk", "PRESUMEout.nwk")
        if (not args.debug):
            shutil.rmtree("intermediate")

    else:
        # output sequences
        print("Generating a SEQ file...")
        for esu in SEQqueue:
            if(args.idANC is None):
                esu.id = str(esu.id)
            else:
                esu_name_prefix = str(hex(args.idANC)).split("x")[1]
                esu_name_suffix = str(hex(esu.id)).split("x")[1]
                new_esu_name = "{}_{}".\
                    format(esu_name_prefix, esu_name_suffix)
                esu.id = new_esu_name

        if args.n is None:
            sequences_writer(SEQqueue, "PRESUMEout.seq",True) # expect "Downstream" mode
        else:
            sequences_writer(SEQqueue, "PRESUMEout.seq")

        # output tree
        tip_count, returned_tree, list_of_dead = create_newick(Lineage)
        # save args
        if args.save:
            argument_saver(args)
    
    print("=====================================================")
    print("Simulation end time point:         "+str(2*timelimit))
    print("Number of generated sequences:     "+str(tip_count))
    print("Seed for random number generation: "+str(seed))
    print("=====================================================\n")
    return 0


def recursive_main(timelimit, limit, main_func, repeated=1):
    result = main_func(timelimit)
    if limit < repeated:
        raise SimulationFailureError(repeated)
    if result == 0:
        return repeated
    else:
        return recursive_main(timelimit, limit, main_func, repeated=repeated+1)


#   interface
if __name__ == "__main__":
    ######## processing some args ########
    parser = argparse.ArgumentParser(description='PRESUME.py', add_help=False)
    parser.add_argument(
        "--monitor",
        help="time limit (default=1)",
        type=float,
        default=1
        )

    parser.add_argument(
        "-n",
        help="required number of sequences: if you specified this parameter,\
            timelimit will be postponed until the number of the sequence reach\
            the specified number (default=None)",
        type=int
        )

    parser.add_argument(
        "-L",
        help="length of sequence (default=1000)",
        type=int,
        default=1000
        )

    parser.add_argument(
        "-d",
        help="doubling time of origin sequence (default=1)",
        type=float,
        default=1
        )

    parser.add_argument(
        "-s",
        help="sigma of doubling time of origin sequence (default=0)",
        type=float,
        default=0
        )

    parser.add_argument(
        "-a", help="CV of CV of doubling time of origin sequence (default=0)",
        type=float,
        default=0
        )

    parser.add_argument(
        "-T",
        help="Threashold of doubling time to be deleted (default=1000)",
        type=float,
        default=1000
        )

    parser.add_argument(
        "-e",
        help="random deletion probability (default=0)",
        type=float,
        default=0
        )

    parser.add_argument(
        "-u",
        help="upper limit of number of sequences (default=2^20)",
        type=int,
        default=np.power(2, 20)
        )

    parser.add_argument(
        "--output",
        help="output folder (default:current directory)",
        type=str,
        default=os.getcwd()
        )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="inactivate deletion of intermediate files"
        )

    parser.add_argument(
        "--bar",
        help="deactivate unstable functions",
        action="store_true",
        default=False
        )

    # for debug
    parser.add_argument(
        "--viewANC",
        help="generate fasta of ancestoral sequences",
        action="store_true",
        default=False
        )

    parser.add_argument(
        "--save",
        help="generate args.csv",
        action="store_true",
        default=False
        )

    parser.add_argument(
        "--CV",
        help="sigma use as CV(Coefficient Variance) of Normal Distribution",
        action='store_true',
        default=False
        )

    parser.add_argument(
        "-r",
        help="limit of retrying simulation (default=100000)",
        type=int,
        default=100000
        )

    parser.add_argument(
        "--seed",
        help="random seed used to initialize \
            the pseudo-random number generator",
        type=str,
        default=None
        )

    parser.add_argument(
        "--tMorigin",
        help="birth time of origin sequence",
        type=float,
        default=0
        )

    parser.add_argument(
        "--idANC",
        help="corresponging ancestral sequence (in upstream tree), \
            in case of distributed computing (default=None)",
        type=int
        )

    parser.add_argument(
        "--state",
        help="number of state",
        type=int,
        default=100
        )

    parser.add_argument(
        "--probedit",
        help="probability of edit",
        type=float,
        default=0
        )

    parser.add_argument(
        "--probdout",
        help="probability of dropout",
        type=float,
        default=0
        )

    args = parser.parse_args()

    # start simulation
    PRESUME_CAS(args)