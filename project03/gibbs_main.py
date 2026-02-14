# Import packages
import random
import numpy as np
import bamnostic as bs
import seqlogo

# Import modules
from seq_ops import reverse_complement
from motif_ops import build_pfm, build_pwm, score_kmer, pfm_ic


def GibbsMotifFinder(seqs, k, seed=None, max_iters=1000):
    """
    Gibbs sampler for motif convergence
    Args:
        seqs (list of str): DNA sequences
        k (int): motif length
        seed (int): random seed
        max_iters (int): number of iterations
    Returns:
        np.ndarray: final Position Frequency Matrix (4 x k)
    """

    # Set seed for reproducibility
    random.seed(seed)
    rng = np.random.default_rng(seed)

    def initialize_random_motifs(seqs, k):
        motifs = []
        for seq in seqs:
            if len(seq) < k:
                continue
            max_start = len(seq) - k
            start = rng.integers(0, max_start + 1)
            motifs.append(seq[start:start + k])
        return motifs

    def score_all_kmers(seq, k, pwm):
        """
        Return:
            scores: list of chosen scores (one per position) for winning motif
            strands: list of strand (0=fwd, 1=rev) selected by random, weighted probabilities
        """
        scores = []
        strands = []
        L = len(seq)
        if L < k:
            return scores, strands

        for start in range(0, L - k + 1):
            kmer = seq[start:start + k]
            rc = reverse_complement(kmer)

            # forward score
            try:
                score_fwd = score_kmer(kmer, pwm)
            except ValueError:
                score_fwd = -np.inf

            # reverse score
            try:
                score_rev = score_kmer(rc, pwm)
            except ValueError:
                score_rev = -np.inf

            # convert to probabilities
            raw = np.array([score_fwd, score_rev], dtype=float)
            raw -= np.max(raw)
            probs = np.exp(raw)
            probs /= probs.sum()

            # sample strand: 0 = forward, 1 = reverse
            strand_choice = rng.choice([0, 1], p=probs)

            # store the score associated with the chosen strand
            chosen_score = score_fwd if strand_choice == 0 else score_rev

            scores.append(chosen_score)
            strands.append(strand_choice)

        return scores, strands

    def scores_to_probabilities(scores):
        scores = np.array(scores, dtype=float)
        if scores.size == 0:
            return scores

        scores -= np.max(scores)
        exp_scores = np.exp(scores)
        total = np.sum(exp_scores)

        if total == 0:
            return np.ones_like(exp_scores) / len(exp_scores)

        return exp_scores / total

    # Initialize motifs
    motifs = initialize_random_motifs(seqs, k)
    if len(motifs) == 0:
        raise ValueError("No motifs initialized. Check sequence lengths and k.")

    # Gibbs Sampling iterations
    for it in range(max_iters):
        N = len(motifs)
        i = rng.integers(0, N)

        # Build PFM/PWM from all motifs except i
        motifs_except_i = [motifs[j] for j in range(N) if j != i]
        pfm = build_pfm(motifs_except_i, k)
        pwm = build_pwm(pfm)

        # Score all k-mers in seq[i]
        seq_i = seqs[i]
        scores, strands = score_all_kmers(seq_i, k, pwm)
        if not scores:
            continue

        probs = scores_to_probabilities(scores)
        positions = np.arange(len(scores))
        m = rng.choice(positions, p=probs)

        # choose strand for this position
        chosen_strand = strands[m]

        # update motif instance
        if chosen_strand == 0:
            motifs[i] = seq_i[m:m + k]
        else:
            motifs[i] = reverse_complement(seq_i[m:m + k])

    # Final PFM
    final_pfm = build_pfm(motifs, k)
    return final_pfm


def load_sequences_from_bam(bam_path, max_reads=10000):
    """
    Load up to max_reads sequences from BAM file
    """
    seqs = []
    with bs.AlignmentFile(bam_path) as bam:
        for idx, read in enumerate(bam):
            if idx >= max_reads:
                break
            if read.seq is None:
                continue
            seqs.append(read.seq)
    return seqs


# Run script
if __name__ == "__main__":
    bam_path = "SRR9090854.subsampled_5pct.bam"  # Input file
    k = 10  # Motif length
    seed = 42  # set seed

    # Extract sequences from BAM file
    print("Loading sequences from BAM file...")  # Print step indicator to terminal
    seqs = load_sequences_from_bam(bam_path, max_reads=500)  # limit max reads for efficiency
    print(f"Loaded {len(seqs)} sequences.")  # Print results indicator to terminal

    # Run model with sequences
    print("Running GibbsMotifFinder...")  # Print step indicator to terminal
    pfm = GibbsMotifFinder(seqs, k, seed=seed, max_iters=10000)  # Obtain freq matrix with set iterations
    print("Gibbs sampling complete.")  # Print results indicator to terminal

    # Print final Position Frequency Matrix to terminal
    print("\nFinal PFM (4 x k):")
    print(pfm)

    # Print final Position Weight Matrix to terminal
    # Calculate PWM from pfm by normalizing each column so sums to 1
    pwm = pfm / pfm.sum(axis=0)  # Use axis to indicate column
    print("\nFinal PWM (A,C,G,T rows):")  # Print results indicator to terminal
    print(pwm)  # Print actual table

    # Calculate most probable consensus sequence
    bases = np.array(["A", "C", "G", "T"])
    consensus = "".join(bases[np.argmax(pwm, axis=0)])
    print("\nMost probable consensus sequence:", consensus)

    # Print final PFM information content value
    ic = pfm_ic(pfm)
    # Print results to terminal
    print(f"\nFinal PFM information content: {ic:.3f}")  # Limit float to 3 decimals

    # Write final results to output text file
    with open("p53_results.txt", "w") as f:  # Open file in write mode
        # Write final Position Frequency Matrix
        f.write("Final PFM (4 x k):\n")  # Header
        f.write(str(pfm) + "\n\n")

        # Write final Position Weight Matrix
        f.write("Final PWM (A,C,G,T rows):\n")  # Header
        f.write(str(pwm) + "\n\n")

        # Write most probable motif string
        f.write("Consensus sequence:\n")  # Header
        f.write(consensus + "\n\n")

        # Write final information content
        f.write(f"Final PFM information content: {ic:.3f}\n")  # Limit IC float to 3 decimals

    print("Saved results to p53_results.txt")  # Print step indication to terminal

    # Create sequence logo and write to png output file
    print("Generating sequence logo...")  # Print step indication to terminal
    logo_pm = seqlogo.CompletePm(pfm=pfm.T)
    seqlogo.seqlogo(logo_pm, format='png', filename='p53_logo.png')
    print("Saved sequence logo to p53_logo.png")  # Print results indicator to terminal
