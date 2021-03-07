#!/usr/bin/env python3

""" run all miniproject code from this script """
# Script Name: run_MiniProject.py
# Author: Jingkai Sun (ks3020@ic.ac.uk)

import subprocess as sp
from subprocess import os
import time

start = time.time()
filename = ["TPC_main.py", "model_analysis.R", "TPC_project.tex"]

# Run python and R scripts
print("%s running..." % filename[0])
print("-----------------------------------------------------")
print(os.system("python3 %s" % filename[0]))
print("%s run successfully! It takes %f minutes." % (filename[0], (time.time() - start)/60))

start2 = time.time()
print("%s running..." % filename[1])
print("-----------------------------------------------------")
print(os.system("Rscript %s" % filename[1]))
print("%s run successfully! It takes %f seconds." % (filename[1], (time.time() - start2)))

start3 = time.time()
print("Compiling LaTeX script...")
print("-----------------------------------------------------")
print(os.system("bash CompileLaTeX.sh %s" % filename[2]))
print("%s run successfully! It takes %f seconds." % (filename[2], (time.time() - start3)))
print("finished! This whole program takes %f minutes to run" % ((time.time() - start)/60))
