#!/usr/bin/env python3
#
# adaptML wrapper that reads its settings from a single config file

import os
import pdb
import sys
import shutil
import subprocess as sub

if len(sys.argv) != 2:
    print("Please supply a single input file to read settings from")
    sys.exit(1)

# Read input file
file_location = sys.argv[1]
with open(file_location, 'r') as f:
    run_configs = f.read()

# Parse input file into dictionary
file_data = run_configs.split("\n")
run_config = {}
for line in file_data:
    line_parts = line.split('=')
    if len(line_parts) == 1:
        continue
    if len(line_parts) != 2:
        raise ValueError("Cannot parse line " + str(line_parts))
    key = line_parts[0].strip()
    value = line_parts[1].strip()
    run_config[key] = value


# Create variable defaults
tree_fn = str(None)
outgroup = str(None)
write_d = './'
color_fn = str(None)
thresh_p = str(0.95)

# Set Variables
for code, arg in run_config.items():
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
base_path = os.path.dirname(os.path.realpath(__file__)) + '/../'
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
