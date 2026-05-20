import gzip
from contextlib import contextmanager
from itertools import product
import os

import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO


all_aas = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', 'D', 'E', 'R',
           'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q']


@contextmanager
def open_file(filename, mode="rt"):
    if filename.endswith(".gz"):
        try:
            handle = gzip.open(filename, mode)
            yield handle
        finally:
            handle.close()
    else:
        try:
            handle = open(filename, mode)
            yield handle
        finally:
            handle.close()


def calculate_conjoint_triad_composition(sequence):
    sequence = list(sequence)
    for i in range(0, len(sequence)):
        if(sequence[i] == 'U' or sequence[i] == 'X'):
            sequence[i] = 'C'
        if(sequence[i] == 'B'):
            sequence[i] = 'D'
        if(sequence[i] == 'Z'):
            sequence[i] = 'Q'
        if(sequence[i] == 'J'):
            sequence[i] = 'L'
    sequence = "".join(sequence)
    ctriad = np.zeros((7,7,7))
    aas = {'A':0,'G':0,'I':1,'L':1,'P':1,'V':0,'F':1,'W':3,'Y':2,'D':5,'E':5,'R':4,'H':3,'K':4,'S':2,'T':2,'C':6,'M':2,'N':3,'Q':3} #based on table 4 from S.Mishra et al.
    for i in range(0,len(sequence)-2):
        ctriad[aas.get(sequence[i]),aas.get(sequence[i+1]),aas.get(sequence[i+2])] += 1
    ctriad /= (len(sequence)-2)
    return ctriad