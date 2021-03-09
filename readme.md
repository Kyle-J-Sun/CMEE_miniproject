# CMEE Coursework MiniProject

__Author:__ Jingkai Sun

__Date:__ 2021.01.22

__Email:__ jingkai.sun20@imperial.ac.uk

## Assignment Description:
- The week8 and Week9 coursework for CMEE students at IC including: MiniProject
- Code Directory:
    - TPC_main.py: the main program of model fitting.
    - model_fitting.py: This includes all functions needed for model fitting
    - model_analysis.R: This is the R script that clean the fitted information after model fitting is done
    - run_MiniProject.py: to run the whole model fitting program in this script.
## Usage:
- In order to run all the scripts, Python 3 and R should be installed.
- All scripts in the code directory must be in the same directory for successful running.
- All necessary Python and R packages need to be installed (See next section)

>> **Warning:** The fitting process will be very slow if you choose extremely large number of iterations to find the best fit.
## Installation of Python packages:
- Please make sure all Python packages shown below have been installed:
    - scipy, numpy, pandas, lmfit, matplotlib
    - You can install the packages mentioned above using following commends in the Terminal
    ```bash
    pip3 install scipy
    pip3 install numpy
    pip3 install pandas
    pip3 install lmfit
    pip3 install matplotlib
    ```
- All R needed packages need to be installed:
    - "dplyr" and "ggplot2"
