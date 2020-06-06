from subprocess import Popen, PIPE
import datetime
import sys
import os

filename = os.path.join(os.getcwd(), "x64", "Release", "Sparse.exe")
algos = ["base", "custom", "blas"]
matrices = ["bcsstk01.mtx", "bcsstk06.mtx"]
repeats = 4

def one_res(name):
	proc = Popen(name, shell=True, stdout=PIPE, stderr=PIPE)
	proc.wait()
	res = proc.communicate()
	if proc.returncode:
		return res[1]
	return res[0]

def min_res(name):
	min_time = 1e6
	for i in range(repeats):
		time = float(one_res(name))
		if time < min_time:
			min_time = time
	return min_time

def get_times():
	run_name = datetime.datetime.now().strftime('%Y:%m:%d:%H:%M:%S').replace(":","_")
	os.makedirs(os.path.join(os.getcwd(), run_name), exist_ok=False)

	with open(os.path.join(os.getcwd(), run_name, "times.csv"), "w") as f:
		f.write(";") 
		for a in algos:
			f.write(a + ';')
		f.write("\n")
		for m in matrices:
			f.write(m + ";")
			for a in algos:
				string = os.path.join(os.getcwd(), filename + " " + m + " " + a)
				f.write(str(min_res(string)) + ";")
			f.write("\n")

if __name__ == "__main__":
	get_times()
