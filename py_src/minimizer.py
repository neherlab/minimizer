from Bio import SeqIO
import numpy as np

cutoff = 1<<28
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

def get_hash(kmer):
    x = 0
    j = 0
    ambig_count= 0
    for i, nuc in enumerate(kmer):
        if i%3==2: continue
        if nuc not in 'ACGT':
            ambig_count += 1
            return cutoff+1
        else:
            if nuc in 'AC':
                x += 1<<j
            if nuc in 'AT':
                x += 1<<(j+1)
        j+=2

    return invertible_hash(x)

def get_minimizers(seq, k=17):
    minimizers = []
    for i in range(len(seq) - k):
        kmer = seq[i:i+k]
        mhash = get_hash(kmer)
        # if mhash in duplicate_minimizers:
        #     print(kmer, i)
        if mhash<cutoff:
            minimizers.append(mhash)
    return np.unique(minimizers)

def make_index(fname):
    minimizers_by_reference = list()
    for si, seq in enumerate(SeqIO.parse(fname, 'fasta')):
        seq_str = str(seq.seq).upper().replace('-', '')
        minimizers = get_minimizers(seq_str)
        minimizers_by_reference.append({"minimizers": minimizers,
                                        "meta":{"length":len(seq.seq),
                                                "description":seq.description,
                                                "n_minimizers":len(minimizers)}})

    index = {"minimizers": {}, "references":[]}
    n_refs = len(minimizers_by_reference)
    for ri, minimizer_set in enumerate(minimizers_by_reference):
        for m in minimizer_set["minimizers"]:
            if m not in index["minimizers"]:
                index["minimizers"][m] = np.zeros(n_refs, dtype=bool)
            index["minimizers"][m][ri] = True
        index["references"].append(minimizer_set['meta'])

    return index


if __name__=='__main__':
    index = make_index('data/references.fasta')

    # duplicate_minimizers = []
    # for m, refs in index.items():
    #     if len(refs)>1:
    #         duplicate_minimizers.append(m)

    normalization = np.array([x['length']/x['n_minimizers'] for x in index["references"]])
    for seq in SeqIO.parse('data/queries.fasta', 'fasta'):
        hit_count = np.zeros(len(index["references"]), dtype=np.int32)
        seq_str = str(seq.seq).upper().replace('-', '')

        minimizers = get_minimizers(seq_str)
        for m in minimizers:
            if m in index["minimizers"]:
                hit_count += index["minimizers"][m]

        res = normalization*hit_count/len(seq.seq)
        print(hit_count, np.round(res,2))
        if np.max(res)<0.3:
            print(seq.description, "no hit")
        else:
            if np.where(res>0.3)[0].shape[0]>1:
                for ri in np.where(res>0.3)[0]:
                    print(seq.description, res[ri], index["references"][ri])
