#!/usr/bin/python

import subprocess, sys

num_cities = 17
num_iterations = 5

def runtime(result):
    lines = result.split("\n");
    return float(lines[1][14:20])
    
infile = "input/dist%d" % (num_cities)
subprocess.call("python mkinput.py %d > %s" % (num_cities, infile), shell=True)

for cores in [128, 64, 32, 16, 8, 4, 2, 1]:
    best_speed = sys.maxint
    for i in range(num_iterations):
        result = subprocess.check_output("mpirun -np %d --hostfile hostfile wsp -i %s" % (cores, infile), shell=True).strip()
        speed = runtime(result)
        if speed < best_speed:
            best_speed = speed
    print "Processors %d -> %.3fs" % (cores, best_speed)