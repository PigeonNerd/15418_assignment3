#!/usr/bin/python

import subprocess
from prettytable import PrettyTable

num_cities = 12
num_iterations = 3

results = PrettyTable(["Seed", "1", "2", "4", "6", "12"])
results.align["Seed"] = "l"

def runtime(result):
    lines = result.split("\n");
    return float(lines[1][14:20])

for seed in [15418, 1000000, 3735928559, 111111111, 11]:
    # Generate input
    infile = "input/dist%d.%d" % (num_cities, seed)
    subprocess.call("python mkinput.py %d --seed=%d > %s" % (num_cities, seed, infile), shell=True)
    
    speeds = []
    serial = 1.0
    for cores in [1, 2, 4, 6, 12]:
        # Run parallel
        speed = 0.0
        for i in range(num_iterations):
            result = subprocess.check_output("mpirun -np %d wsp -i %s" % (cores, infile), shell=True).strip()
            speed += runtime(result)

        speed = speed / num_iterations
        
        if cores is 1:
            serial = speed
            speeds.append("%.3fs" % (speed))
        else:
            speeds.append("%.3fs (%.2fx)" % (speed, serial / speed))
        

    results.add_row([seed] + speeds)
    
print results