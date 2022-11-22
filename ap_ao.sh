#!/bin/bash

#activate conda env if needed

source activate geminiconda3

#link XPA method in Py and DS9 env (needed for conda pyds9 to connect to DS9)

XPA_METHOD=localhost
export XPA_method

#applies all inputs after the wrapper call
./ap_ao.py $@ 


