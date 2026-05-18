#!/usr/bin/env python3
#
# adaptML wrapper

import os
import pdb
import sys
import shutil
import subprocess as sub

print("\nwelcome to WrapLikelihood, the AdaptML wrapper.")

# read in the variables
inputs = sys.argv
tree_fn = None
outgroup = None
write_d = './'
color_fn = None
thresh_p = "0.95"

base_path = os.path.dirname(os.path.realpath(__file__)) + '/../'

for ind in range(1, len(inputs)):
    arg_parts = inputs[ind].split('=')
    code = arg_parts[0]
    arg = arg_parts[1]
    if code == 'tree':
        tree_fn = arg
    elif code == 'outgroup':
        outgroup = arg
    elif code == 'write_dir':
        write_d = arg
    elif code == 'color':
        color_fn = arg
    elif code == 'thresh':
        thresh_p = arg

###################
# get likelihoods #
###################

print("Obtaining empirical thresholds ...")

getliks_x = base_path + "clusters/getstats/GetLikelihoods.py"
getliks_l = []
getliks_l.append(sys.executable)
getliks_l.append(getliks_x)
getliks_l.append(write_d + "/emp_trees/")
getliks_l.append(write_d)
getliks_l.append(str(thresh_p))
proc = sub.Popen(getliks_l, stdout=sub.PIPE)
# wait until finished
stdout, stderr = proc.communicate()

###############
# run JointML #
###############

print("Running JointML ...")
jointml_x = base_path + "clusters/trunk/JointML.py"
jointml_l = []
jointml_l.append(sys.executable)
jointml_l.append(jointml_x)
jointml_l.append("tree=" + tree_fn)
jointml_l.append("outgroup=" + outgroup)
jointml_l.append("write=" + write_d)
jointml_l.append("habitats=" + write_d + "/habitat.matrix")
jointml_l.append("mu=" + write_d + "/mu.val")
jointml_l.append("color=" + color_fn)
jointml_l.append("thresh=" + write_d + "/thresh.file")

proc = sub.Popen(jointml_l, stdout=sub.PIPE)

# wait until finished
stdout, stderr = proc.communicate()
