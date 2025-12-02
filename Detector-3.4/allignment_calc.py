from Bio.Align import PairwiseAligner

def alignment_calculator(human_seq, input_seq, max_bases=100000):
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0

    human_slice = str(human_seq[:max_bases])
    input_slice = str(input_seq[:max_bases])
    score = aligner.score(human_slice, input_slice)
    percentage = (score / len(input_slice)) * 100
    return percentage

# === process one input contig against human contigs ===
def process_contigs(input_seq, human_contigs):
    local_results = []
    for human_seq in human_contigs:
        pct = alignment_calculator(human_seq, input_seq)
        local_results.append((input_seq.id, human_seq.id, pct))
    return local_results
