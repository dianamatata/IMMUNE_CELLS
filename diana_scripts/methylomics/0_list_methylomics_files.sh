#!/bin/bash

# Goal: list all files from folder /home/users/a/avalosma/scratch/Blueprint/DNA_meth
FOLDER=/home/users/a/avalosma/scratch/Blueprint/DNA_meth
find $FOLDER -name "*.idat" > filelist_methylation.txt
 # find recursively from folder

