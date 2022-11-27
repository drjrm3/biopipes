#!/usr/bin/env python

import os.path
import linecache
import csv
import os
import numpy as np
import shutil
import sys


## {{{ http://code.activestate.com/recipes/65441/ (r2)
def containsAny(str, set):
    """Check whether 'str' contains ANY of the chars in 'set'"""
    return 1 in [c in str for c in set]

def containsAll(str, set):
    """Check whether 'str' contains ALL of the chars in 'set'"""
    return 0 not in [c in str for c in set]

## end of http://code.activestate.com/recipes/65441/ }}}


def contains(theString, theQueryValue):
  return theString.find(theQueryValue) > -1

""" Linear search, returns first index """
def find_first_index(lst, elem):
    ind = 0
    for l in lst:
        if str(elem)==str(l):
            return ind
    return -1


def str2bool(v):
  return v.lower() in ["yes", "true", "t", "1"]

def isExist(filename):
    if os.path.exists(filename) and os.path.isfile(filename):
        return True
    else:
        return False

def fileSize(filename):
    return int(os.path.getsize(filename))

def delete(filename):
    if os.path.exists(filename) and os.path.isfile(filename):
        os.unlink(filename)

def mkdirp(directory):
    """ makes directory if it does not exist
    """
    if not os.path.isdir(directory):
        os.makedirs(directory)

def get_column(path, c=0, r=1, sep='\t'):
    """ extracts column specified by column index
        assumes that first row as a header
    """
    try:
        reader = csv.reader(open(path, "r"), delimiter=sep)
        return [row[c] for row in reader] [r :]
    except IOError:
        print('list_rows: file "'+path +'" does not exist')
        return 'list_rows failed'


def read_one_int_col(filename):
    fh = open(filename, "r")
    values = []
    for line in fh:
        values.append(int(line.strip('\r\n')))
    return values

def read_one_float_col(filename):
    fh = open(filename, "r")
    values = []
    for line in fh:
        values.append(float(line.strip('\r\n')))
    return values

def read_one_str_col(filename):
    fh = open(filename, "r")
    values = []
    for line in fh:
        values.append(line.strip('\r\n'))
    return values

def get_index_of_col_or_row(lst, value):
    try:
        return lst.index(value)
    except:
        print('get_index_of_col_or_row: value not found "' + value + '"')
        return -1

def array2str(array, sep='\t'):
    strA = []
    for a in array:
        strA.append(str(a))
    return sep.join(strA)

def array2header(array, sep='\t'):
    strA = ["samples"]
    for a in array:
        strA.append('p'+str(a))
    return sep.join(strA)


def readindices(filename, sep='\t'):
    fh = open(filename, "r")
    values = []
    for line in fh:
        line = line.strip('\n')
        if len(line) > 0:
            if len(line.split(sep)) == 1:
                values.append(int(line))
            else:
                start=int(line.split(sep)[0])
                end=int(line.split(sep)[1])
                while start<=end:
                    values.append(start)
                    start=start+1

    return sorted(values)



def totxtfile(read_data, txtfile, compress=False, debug=True):
    b=np.array(read_data)
    np.savetxt(txtfile, b, fmt="%0.3G")
    if compress==True:
        os.system('gzip -f ' + str (txtfile) )
    if debug==True:
        print ("Written " + str(txtfile) )

def save2txt(read_data, txtfile, compress=False, debug=True):
    """
    Saves list of rows and columns in a text file
    """
    try:
        f=open(txtfile, 'w')
        tmp = array2str(array=read_data, sep='\n')
        f.write(tmp)
        if compress==True:
            os.system('gzip -f ' + str (txtfile) )
        if debug==True:
            print ("Written " + str(txtfile) )
    except IOError:
        print('save2txt: can not write to file "' + file)
        return 'save_list_of_str failed'
    finally:
        f.close()
