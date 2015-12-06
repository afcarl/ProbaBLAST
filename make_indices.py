import time 
import itertools
import numpy as np 
import cPickle as pickle 

thresholds = [0.925,0.95,0.975]
data_dir = './Data/'
prob_fname = data_dir+'probabilities.txt'
genome_fname = data_dir+'genome.txt'
fnames = [data_dir+'max_likelihood.pkl'] + [data_dir+'t='+str(t)+'.pkl' for t in thresholds]
w = 11 
keys = [word for word in map(''.join, itertools.product('ACTG', repeat=w))]

def make_indices(genome,probabilities):
	indices = {}
	# maximum likelihood works only with genome
	print "Building index for ML ... "
	start = time.time()
	ml_index = {word: [] for word in map(''.join, itertools.product('ACTG', repeat=w))}
	for i in range(len(genome)-w+1):
		seed = genome[i:i+w]
		ml_index[seed].append(i)
	indices['ml'] = ml_index
	end = time.time()
	print "... done"
	print "Running time:", end-start
	print
	pickle.dump(ml_index, open( fnames[0], "wb" ) )

	
	# thresholds work with both the genome and the probabilities
	for k,t in enumerate(thresholds):
		threshold = t**w
		print threshold
		print "Building index for threshold " + str(t) + " ... "
		start = time.time()
		thr_index = {word: [] for word in map(''.join, itertools.product('ACTG', repeat=w))}
		for i in range(len(genome)-w+1):
			seed = genome[i:i+w]
			prob = 1
			for j in range(i,i+w):
				prob *= probabilities[j]
			if prob >= threshold:
				thr_index[seed].append(i)
		end = time.time()
		print "... done"
		print "Running time:", end-start
		print
		pickle.dump(thr_index, open( fnames[k+1], "wb" ) )
		indices[str(t)] = thr_index
	return indices

def get_data(fname):
	with open(fname,'r') as f:
		lines = [line.rstrip('\n') for line in open(fname)]
	return lines 

def make_float_array(array):
	elts = array[0].split(' ')[:-1]
	elts = [float(e) for e in elts]
	return elts

if __name__ == "__main__":
	genome = get_data(genome_fname)
	probabilities = make_float_array(get_data(prob_fname))
	make_indices(genome[0],probabilities)
