#!/bin/python
import subprocess
import argparse
import os
os.chdir("build")

def run_murb_once(bodies, iters):
	"""
		Run ./bin/murb -i iters -n bodies --nv and parses the output
	"""
	
	sp = subprocess.Popen(["./bin/murb", "-i", str(iters), "-n", str(bodies), "--nv"],
		stdout=subprocess.PIPE)
	output = sp.stdout.read()
	stats = output.splitlines()[-1].split()
	ms = float(stats[3])
	fps = float(stats[5][1:])
	return ms, fps

def run_murb_many(bodies, iters, n):
	"""
		Run ./bin/murb n times to get an average
	"""

	tot_ms = 0
	tot_fps = 0

	for _ in range(n):
		ms, fps = run_murb_once(bodies, iters)
		tot_ms += ms
		tot_fps += fps
	
	return tot_ms / n, tot_fps / n

def gen_arg_parser():
	parser = argparse.ArgumentParser(
        prog="murb-bench")
	parser.add_argument("-i", help="number of iterations", type=int, default=1000)
	parser.add_argument("-b", help="number of bodies", type=int, default=300)
	parser.add_argument("-n", help="number of times to run the program", type=int, default=10)
	return parser

if __name__ == "__main__":
	parser = gen_arg_parser().parse_args()

	ms, fps = run_murb_many(parser.b, parser.i, parser.n)
	print("{} FPS, {} ms".format(fps, ms))


