#!/usr/bin/env python

import os
import sys
import platform
import multiprocessing
from algbioi.com import csv


def printEnvironmentInfo():
    print("System:                %s" % platform.system())
    print("Node:                  %s" % platform.node())
    print("Number of processors:  %s" % multiprocessing.cpu_count())
    print("Platform:              %s" % platform.platform())
    print("Working dir:           %s" % os.getcwd())


def toUpper():
    for line in sys.stdin:
        print(str(line).upper())
    print("EOF")


def printDict(d):
    """
        @param d: dictionary
        @type d: dict
    """
    for k, v in d.iteritems():
        print("k: %s, v: %s" % (k,v))

def printMaps():
    a = ""
    printDict(a)


os.mkdir('tmp' )



def nestedCycles():
    """
        Nested cycls.

    """
    for i in range(10):
        for j in range(10):
            for k in range(10):
                print("%s %s %s" % (i, j, k))




















if __name__ == "__main__":
    printEnvironmentInfo()
    #toUpper()





    #print 'release  :', platform.release()
    #print 'platform.machine()', platform.machine()
    #print 'platform.version()', platform.version()
    #print 'platform.uname()', platform.uname()
    #print 'platform.system()', platform.system()
    #print 'platform.processor()', platform.processor()
    #print 'os.uname()', os.uname()
    #print 'platform.uname()', platform.uname()
    #print 'version  :', platform.version()
    #print 'machine  :', platform.machine()
    #print 'processor:', platform.processor()
    #print '------------'
    #print platform.architecture()
    #print platform.python_build()
    #print platform.python_compiler()
    #print platform.python_version()
    #print platform.python_implementation()
    #print platform.uname()
    #print platform.linux_distribution()
    #print '------------'
    #print os.getcwd()
    #raise AssertionError