import utils
import cPickle
import numpy as np
import matplotlib.pyplot as plt 

genome, probs = utils.get_genome_and_probs()

class Match(object):
	def __init__(self, g_ix, q_ix, length):
		self.genome_ix = g_ix
		self.query_ix = q_ix
		self.length = length
		self.score = 0
	def __str__(self):
		return 'g:%d\tq:%d\tlen:%d\tscore:%d' % (self.genome_ix, self.query_ix, self.length, self.score)
	def __eq__(self, other):
		if self.genome_ix == other.genome_ix and \
			self.query_ix == other.query_ix and \
			self.length == other.length:
			return True
		else:
			return False
	def __ne__(self, other):
		return not self.__eq__(self, other)
	def __hash__(self):
		return hash(str(self.genome_ix) + str(self.query_ix) + str(self.length))
	def __le__(self, other):
		return self.score <= other.score
	def __lt__(self, other):
		return self.score < other.score
	def __ge__(self, other):
		return self.score >= other.score
	def __gt__(self, other):
		return self.score > other.score

def load_index(fname):
	with open(fname, 'rb') as handle:
		index = cPickle.load(handle)
	return index

def get_seeds(query, index):
	seeds = []
	included = [0] * len(genome)
	for i in xrange(0, len(query)-11):
		if query[i:i+11] in index:
			genome_locs = index[query[i:i+11]]
			#seeds += [ for ix in genome_locs]
			for ix in genome_locs:
				if included[ix] != 1:
					seeds.append(Match(ix, i, 11))
					included[i:i+11] = [1] * 11

	return seeds

def ungapped_extension(query, match, threshold, delta=10):
	'''
	Given a query string and its seed, we try to extend the seed to the left and the right
	as long as the score stays above the given threshold.

	Parameters
	----------
	query: the complete query sequence
	q_ix: the index of where the seed starts within the query sequence
	q_length: the length of the seed
	g_ix: where the seed starts within the genome

	Returns
	-------
	seq_ix: the ix within query where the extended sequence starts
	seq_length: the extended sequence length
	new_g_ix: the new index in the genome where the sequence starts
	'''
	def get_score(query, m):
		seq = query[m.query_ix:m.query_ix + m.length]
		gen = genome[m.genome_ix:m.genome_ix + m.length]
		return utils.score(gen, probs[m.genome_ix:m.genome_ix+m.length], seq)
	score = get_score(query, match)
	max_score = score
	# X-Drop to the left.
	while max_score - score < delta:
		if match.query_ix < 1 or match.genome_ix < 1: break
		match.query_ix -= 1
		match.genome_ix -= 1
		match.length += 1
		score = get_score(query, match)
		if score > max_score: max_score = score

	max_score = score
	# X - drop to the right.
	while max_score - score < delta:
		if match.query_ix + match.length >= len(query) or match.genome_ix + match.length >= len(genome): break
		match.length += 1
		score = get_score(query, match)
		if score > max_score: max_score = score

	if score > threshold: 
		match.score = score
		return match
	else:
		return None

def gapped_extension():
	# TODO: Need to write NW algorithm.
	pass


if __name__ == '__main__':
	index = load_index('Data/t=0.95.pkl')
	# First we load the test sets.
	# TODO: Try test lengths of various sizes.
	t_100 = [ 25, 50, 75, 100]
	t_50 = [ 12, 25, 37, 50]
	t_25 = [ 6, 12, 19, 25]
	k = 0
	colours = ['r-', 'b-', 'g-', 'k-']
	hs = []
	xs = [i+5000*3 for i in xrange(len(genome))]
	ys = [0] * len(xs)
	ys[443302] = 100
	h = plt.plot(xs, ys, colours[3], label='Gold Standard')
	hs.append(h[0])
	for t in t_50[1:]:
		with open('test_50.txt', 'r') as handle:
			recs = []
			precs = []
			#for entry in handle.xreadlines():
			for entry in ['GATCTCCCGGCCACCAGTAGAGTATCATTATCCCCATTTTACAGGTGAGG,443302']:
				#print cur
				#cur += 1
				query, correct_ix = entry.strip().split(',')
				correct_ix = int(correct_ix)
				# First we need to find all the seeds of length 11.
				#print 'Gathering seeds...'
				seeds = get_seeds(query, index)
				# Now try to extend each seed.
				#print 'Ungapped Extension Phase...'
				extended_seeds = set()
				for seed in seeds:
					new_seed = ungapped_extension(query, seed, t)
					if not new_seed is None and new_seed not in extended_seeds:
						extended_seeds.add(new_seed)
				for seed in extended_seeds:
					#print seed
					# TODO: Call the gapped extension phase.
					pass
				if len(extended_seeds)>0:
					print query, correct_ix
					sort = sorted(extended_seeds)
					print
					xs = [i+5000*k for i in xrange(len(genome))]
					ys = [k+1] * len(genome)
					for s in sort:
						ys[s.genome_ix] = s.score + k*2
						print s
					h = plt.plot(xs, ys, colours[k], label='T=%d' %t)
					hs.append(h[0])
					k += 1
					print
				found = False
				p_count, p_total = 0. , 0.
				for match in extended_seeds:
					if correct_ix >= match.genome_ix and correct_ix <= match.genome_ix + match.length:
						p_count += 1.
						found = True
					p_total += 1.
				if found == True:
					recs.append(1.)
				else:
					recs.append(0.)
				if p_total == 0: p_total += 1.
				precs.append(p_count/p_total)
			print t
			print 'P: ', np.mean(precs)
			print 'R: ', np.mean(recs)
	plt.xlabel('Position of match in genome')
	plt.ylabel('Score of match')
	plt.legend(handles=hs)

	plt.show()			

	# TODO: Calculate precision/recall. Get ROC curves.
	# Precision is how often the match is actually where it was generated from.
	# Recall is how many matches we correctly identify. 




