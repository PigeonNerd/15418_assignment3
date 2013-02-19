#!/usr/bin/python2.7

import sys
if sys.version_info < (2,7):
  raise Exception("""
Must use python 2.7 or later. 
On blacklight, run `module load python`.
On ghc, use the python in /usr/local/bin by prepending it to your PATH
""")

import argparse

import random
import math

def positive_int(value):
  ivalue = int(value)
  if ivalue <= 0:
    raise argparse.ArgumentTypeError("%s is not a positive int value" % value)
  return ivalue

parser = argparse.ArgumentParser(description="Generate input for wsp")
parser.add_argument("num_cities", 
    help="Number of cities to generate for",
    type=positive_int)
parser.add_argument("--seed",
    help="Random seed for generator",
    type=positive_int)
parser.add_argument("--size",
    help="Size of world (max value of x and y coordinat)",
    type=positive_int,
    default=100)
parser.add_argument("--print_coords", 
    help="Print city coordinates",
    action='store_true')

args = parser.parse_args()

if args.seed:
  random.seed(args.seed)

cities = []
for i in xrange(args.num_cities):
  cities.append((random.randrange(args.size), random.randrange(args.size)))

def distance(src, dst):
  return int(math.sqrt((src[0] - dst[0])**2 + (src[1] - dst[1])**2))

if args.print_coords:
  for i, city in enumerate(cities):
    print "City {} is at {}".format(i, city)
else:
  print args.num_cities
  for i, dst in enumerate(cities):
    if i == 0: continue
    print " ".join(str(distance(src,dst)) for src in cities[:i])
