#!/bin/bash
# Author: jingkai.sun20@imperial.ac.uk
# Script: run_MiniProject.sh
# Desc: run all miniproject code from this script
# Arguments: none
# Date: Dec 2020

filename=$1
echo "Starting to run $filename"
Rscript $filename
echo "finished"