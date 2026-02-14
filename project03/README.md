# Introduction


This project applies a Gibbs sampling–based motif discovery pipeline 
to a subsampled p53 ChIP‑seq dataset from human K562 leukemia cells 
treated with daunorubicin, a chemotherapeutic agent that stabilizes and 
activates p53. Known as “the guardian of the genome,” p53 is a critical 
tumor suppressor gene located on human chromosome 17. The p53 protein has 
numerous identified binding sites within the human genome, and regulates 
cell‑cycle arrest, DNA repair, and apoptosis in response to cellular 
stress (Shield & Mirkes, 1998).

The BAM file SRR9090854.subsampled_5pct.bam contains approximately 5% of 
the original aligned read pairs from experiment SRX5865974 (run SRR9090854), 
in which a p53 antibody was used to immunoprecipitate p53‑bound chromatin 
fragments prior to Illumina NextSeq 500 sequencing. Because the reads originate 
from regions physically bound by p53, the underlying p53 DNA‑binding motif is 
expected to recur across many independent fragments. In this sense, the dataset 
is already biologically enriched for the motif, allowing motif discovery to be 
performed on a manageable subset of sequences.

The bamnostic package provides a pure‑Python interface for streaming aligned reads 
directly from the BAM file, enabling efficient extraction of nucleotide sequences 
without loading the entire alignment into memory. These sequences serve as input 
to the Gibbs sampler.

Rather than exhaustively scanning every possible k‑mer across millions of reads, 
the Gibbs sampler updates one sequence at a time. In each iteration, a single 
sequence is held out, its current motif prediction is temporarily removed, and a 
PFM/PWM is constructed from the remaining sequences’ motif instances. The algorithm 
then slides a window across the held‑out sequence, scoring every possible k‑mer under 
the current PWM. These scores are converted into a probability distribution, and a 
new motif position is sampled based on this distribution. Important to note is that
while a higher probability may be more likely to occur, it is not guaranteed each
time, which minimizes the biased selection of the most probable instance every time.
The updated motif instance is then reinserted, and a different sequence is randomly 
selected for the next iteration. Over many cycles, this process fine-tunes the PWM,
resulting in convergence upon a stable motif.

Because p53 motifs are relatively short (~10 bp), degenerate, and enriched in ChIP‑seq 
reads from p53‑bound chromatin, the sampler can converge on a high‑confidence motif 
using only a subset of the total reads. However, it is essential to treat the forward 
and reverse‑complement orientations symmetrically. p53 binds double‑stranded DNA, 
and its motif can appear in either orientation across different genomic loci. To 
avoid artificially biasing the PWM toward a single orientation, our implementation 
retains the full forward and reverse‑complement scores at each position and 
incorporates strand choice directly into the sampling step. Instead of selecting 
the strand with the highest score, the algorithm samples from the joint probability 
distribution over (position, strand) pairs. This strand randomization prevents 
the sampler from favoring an orientation early on and ensures that the final 
PWM reflects the true binding possibilities of p53.



# Pseudocode
Put pseudocode in this box:

```
Some pseudocode here
```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
