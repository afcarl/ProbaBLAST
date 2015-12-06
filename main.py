import utils
import cPickle

class Match(object):
	def __init__(self, g_ix, q_ix, length):
		self.genome_ix = g_ix
		self.query_ix = q_ix
		self.length = length
	def __str__(self):
		return 'g:%d\tq:%d\tlen:%d' % (self.genome_ix, self.query_ix, self.length)
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

def load_index(fname):
	with open(fname, 'rb') as handle:
		index = cPickle.load(handle)
	return index

def get_seeds(query, index):
	seeds = []
	for i in xrange(0, len(query)-11):
		if query[i:i+11] in index:
			seeds += [Match(ix, i, 11) for ix in index[query[i:i+11]]]
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
	genome, probs = utils.get_genome_and_probs()
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
		return match
	else:
		return None

def gapped_extension():
	pass

if __name__ == '__main__':
	index = load_index('Data/max_likelihood.pkl')
	# First we load the test sets.
	with open('test.txt', 'r') as handle:
		for entry in handle.xreadlines():
			query, correct_ix = entry.strip().split(',')
			correct_ix = int(correct_ix)
			print correct_ix
			# First we need to find all the seeds of length 11.
			print 'Gathering seeds...'
			seeds = get_seeds(query, index)
			# Now try to extend each seed.
			print 'Ungapped Extension Phase...'
			extended_seeds = set()
			for seed in seeds:
				new_seed = ungapped_extension(query, seed, 10)
				if not new_seed is None and new_seed not in extended_seeds:
					extended_seeds.add(new_seed)
			for seed in extended_seeds:
				pass


			break




