from subprocess import Popen, PIPE
import datetime
import sys
import os

filename = os.path.join(os.getcwd(), "x64", "Release", "Sparse.exe")
algos = ["base", "custom", "blas"]
matrices = os.listdir(os.path.join(os.getcwd(), "matrices"))
print(matrices)
repeats = 3

def one_res(name):
	proc = Popen(name, shell=True, stdout=PIPE, stderr=PIPE)
	proc.wait()
	res = proc.communicate()
	if proc.returncode:
		return res[1]
	return res[0]

def min_res(name):
	min_res = [1e6, 1e6]
	for i in range(repeats):
		res = (one_res(name))
		res = str(res)[2:-1]
		if float(res[0]) < float(min_res[0]):
			min_res = res
	return [float(r) for r in min_res.split(' ')]

def get_times():
	run_name = datetime.datetime.now().strftime('%Y:%m:%d:%H:%M:%S').replace(":", "_")
	os.makedirs(os.path.join(os.getcwd(), run_name), exist_ok=False)

	with open(os.path.join(os.getcwd(), run_name, "times.csv"), "w") as f:
		f.write(";") 
		for a in algos:
			f.write(a + ';' + "error;")
		f.write("\n")
		for m in matrices:
			f.write(m + ";")
			for a in algos:
				string = os.path.join(os.getcwd(), filename + " " + m + " " + a)
				res = min_res(string)
				# print(res)
				f.write(str(res[0]) + ";" + str(res[1]) + ";")
			f.write("\n")

if __name__ == "__main__":
	get_times()
