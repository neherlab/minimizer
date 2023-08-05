from Bio import SeqIO
import numpy as np

# minimizer cutoff. The max is 1<<32 - 1, so with 28 uses roughly 1/16 of all kmers
cutoff = 1<<28

# from lh3
def invertible_hash(x):
    m = (1<<32) - 1
    x = (~x + (x << 21)) & m
    x = x ^ (x >> 24)
    x = (x + (x << 3) + (x << 8)) & m
    x = x ^ (x >> 14)
    x = (x + (x << 2) + (x << 4)) & m
    x = x ^ (x >> 28)
    x = (x + (x << 31)) & m
    return x

# turn a kmer into an integer
def get_hash(kmer):
    x = 0
    j = 0
    for i, nuc in enumerate(kmer):
        if i%3==2: continue # skip every third nucleotide to pick up conserved patterns
        if nuc not in 'ACGT':
            return cutoff+1 # break out of loop, return hash above cutoff
        else: # A=11=3, C=10=2, G=00=0, T=01=1
            if nuc in 'AC':
                x += 1<<j
            if nuc in 'AT':
                x += 1<<(j+1)
        j+=2

    return invertible_hash(x)

def get_minimizers(seq, k=17):
    minimizers = []
    # we know the rough number of minimizers, so we can pre-allocate the array if needed
    for i in range(len(seq) - k):
        kmer = seq[i:i+k]
        mhash = get_hash(kmer)
        if mhash<cutoff: # accept only hashes below cutoff --> reduces the size of the index and the number of look-ups
            minimizers.append(mhash)
    return np.unique(minimizers)

def make_index(fname):
    # collect minimizers for each reference sequence first
    minimizers_by_reference = list()
    for seq in SeqIO.parse(fname, 'fasta'):
        seq_str = str(seq.seq).upper().replace('-', '')
        minimizers = get_minimizers(seq_str)
        minimizers_by_reference.append({"minimizers": minimizers,
                                        "meta":{"length":len(seq.seq),
                                                "description":seq.description,
                                                "n_minimizers":len(minimizers)}})

    # construct an index where each minimizer maps to the references it contains via a bit set (here boolean np array)
    index = {"minimizers": {}, "references":[]}
    n_refs = len(minimizers_by_reference)
    for ri, minimizer_set in enumerate(minimizers_by_reference):
        for m in minimizer_set["minimizers"]:
            if m not in index["minimizers"]:
                index["minimizers"][m] = np.zeros(n_refs, dtype=bool)
            index["minimizers"][m][ri] = True # same as += 1<<ri

        # reference will be a list in same order as the bit set
        index["references"].append(minimizer_set['meta'])

    return index


if __name__=='__main__':
    index = make_index('data/references.fasta')
    normalization = np.array([x['length']/x['n_minimizers'] for x in index["references"]])

    overall_hits = np.zeros(len(index["references"]), dtype=np.int32)
    for seq in SeqIO.parse('data/queries.fasta', 'fasta'):
        seq_str = str(seq.seq).upper().replace('-', '')

        minimizers = get_minimizers(seq_str)
        hit_count = np.zeros(len(index["references"]), dtype=np.int32)
        for m in minimizers:
            if m in index["minimizers"]:
                hit_count += index["minimizers"][m]

        # we expect hits to be proportional to the length of the sequence and the number of minimizers per reference
        normalized_hits = normalization*hit_count/len(seq.seq)
        # require at least 30% of the maximal hits and at least 10 hits
        if np.max(normalized_hits)<0.3 or np.sum(hit_count)<10:
            print(seq.description, "no hit")
        else:
            overall_hits += normalized_hits>0.3
            ri = np.argmax(normalized_hits)
            print(f"{seq.description}\t best hit={normalized_hits[ri]:1.2f} to reference {index['references'][ri]['description']}")

    print("\nHits statistics:")
    for i, ref in enumerate(index["references"]):
        print(f"\t{ref['description']}\t{overall_hits[i]}")
    ## we could offer the user to run the analysis on these datasets in reverse order of the number of hits


    print(f"\nIndex statistics:")
    print(f"\tNumber of references: {len(index['references'])}")
    print(f"\tNumber of minimizers: {len(index['minimizers'])}")
    print(f"\tNumber of minimizers per kb: {1000*np.sum([x['n_minimizers'] for x in index['references']])/np.sum([x['length'] for x in index['references']]):1.2f}")
    for ref in index["references"]:
        print(f"\t\t{ref['description']}\t{ref['n_minimizers']}")

    # check uniformity of hash function
    import matplotlib.pyplot as plt
    plt.figure()
    for start in range(0, 1<<29, 1<<22):
        # compute hashes for 2^12 integers starting at start
        hashes = [invertible_hash(x) for x in range(start, start+(1<<12))]
        # sort and plot --> should be a straight line from 0 to 2^32-1
        plt.plot(sorted(hashes))
