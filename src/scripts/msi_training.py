#!/usr/bin/env python
# encoding: utf-8
"""
@author: jinlinfang
@contact: crawler_jinlinfang@163.com
@file: msi_traning.py
@time: 8/8/18 11:38 AM
"""
import numpy, pandas, re, subprocess, os, argparse


def func():
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="file", usage="\n"
                                                               "Usage:\n"
                                                               "")
    parser.add_argument("-n", "--normal", help="normal bam' path to construct baseline",
                        required=True)
    parser.add_argument("-s", "--mss", help="mss bam' path to construct baseline",
                        required=True)
    parser.add_argument("-h", "--msh", help="msh bam' path to construct baseline",
                        required=True)
    parser.add_argument("-b", "--bin", help="msidetector binary to call baseline data", required=True)
    parser.add_argument("-m", "--model", help="", required=True)
    parser.add_argument("-l", "--loci", help="", required=True)
    parser.add_argument("-r", "--classified_r", help="", required=True)

    args = vars(parser.parse_args())
    normal_root = args["normal"]
    mss_root = args["mss"]
    msh_root = args["msh"]
    binary = args["bin"]
    model = args["model"]


