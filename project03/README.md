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


### Part 1: Setup and Initialization 
```
Initialize random number generators (use seed)

Use the code snippet Marcus gave us in an announcement to load the reduced dataset into a list called 'seqs'

Randomly initialize motifs: 

    Create an empty list motifs. 

    For each sequence in seqs: 

    Check the sequence length field and generate a random start position between 0 and seq_length – k 

        *Should also be generating random strand selection- FWD=0, REV=1 

    Append start position and strand as tuple to motifs list -or- store each motif as a tuple: (start_position, strand) 
```
### Part 2: Iterative Gibbs Sampling 
```
Repeat for max_iter iterations: 

    Randomly select a sequence index i to update. 

    Build a PWM from all motifs (99 other randomly selected ones) except the i-th: 
other_motifs = motifs[:i] + motifs[i+1:] 

        pfm = build_pfm(other_motifs, k) 

        pwm = build_pwm(pfm) 

    Score all possible k-mers in sequence i (consider both strands) based on the PWM derived from the other 99 random sequences: 

        Initialize empty lists scores and kmers. 

        For each possible start position in sequence i: 

            Extract k-mer of length k. 

            Compute reverse complement of k-mer. 

            Score both k-mer and reverse complement using PWM. 

            Take the maximum score for any given position amongst the forward and reverse strands’ scores. 

            Append score to scores and k-mer to kmers. (the kmer that gets added will be from the correct strand) 

    Convert scores to probabilities for sampling: 

        exp_scores = exp(scores - max(scores)) (raise e to the power of score-max_score to get all positive values, and subtract all scores from the max score to keep the probabilities from exceeding 1) 

        probs = exp_scores / sum(exp_scores) 

    Probabilistically select a new motif for sequence i: 

        Choose index new_index from 0 to len(kmers)-1 using probs as weights. 

        Update motifs[i] = kmers[new_index]. 

    Monitor convergence 

        Every N iterations (maybe 10,000): 

        Build current PFM from a random subset (1,000?) of motifs. 

        Calculate InfoContent using pfm_ic(pfm_current). 

        Print iteration number and InfoContent. 

        Check if the information content of the current convergence PWM is “close enough” to the information content of the previous convergence PWM 
```
### Part 3: Final Motif 
```
After all iterations, build final PFM from motifs. 

Return final_pfm. 

Plug the final convergence PWM into seqlogo and save the resulting logo to a png file :D 

```

# Successes
One of the biggest successes of this project was our team collaboration. Since we had already worked together in a previous Genomics course, we were able to communicate efficiently and divide responsibilities in a way that played to our strengths. The pseudocode development process was particularly strong; we spent significant time carefully mapping out the Gibbs Sampling algorithm before implementing it. This planning phase helped us think critically about both the conceptual MCMC framework and the practical considerations of runtime efficiency. Translating the theoretical content from lecture slides into a working mental model was challenging but ultimately very rewarding.

We are also proud of how we handled the scale of the K562 dataset. Processing millions of reads with a Gibbs Sampling approach could easily have resulted in unmanageable runtimes. Although we initially considered using MACS2 peak calls or a sliding window strategy to reduce the search space, we ultimately implemented a thoughtful sub-sampling strategy that allowed us to balance computational feasibility with biological relevance. Even though we were not able to incorporate every idea due to time constraints, the final pipeline reflects careful design decisions and strong conceptual understanding. Overall, we are very satisfied with how the project came together.

# Struggles
The primary technical challenge we faced was the sheer size of the alignment data. Gibbs Sampling requires repeated random selections and motif scoring, which can quickly become computationally expensive when applied to millions of reads. Much of our early effort focused on ensuring that our implementation would be efficient enough to run within a reasonable timeframe. This required us to be deliberate about algorithm design and mindful of how each implementation choice would affect performance.

Another major hurdle was fully grasping the MCMC framework and translating it into clear, actionable pseudocode. Converting theoretical descriptions of Gibbs Sampling into a structured algorithm required multiple rounds of discussion and refinement. On a practical level, we also encountered environment and tooling issues. In particular, syncing repositories in PyCharm led to repeated notebook corruption issues, which limited some team members’ ability to contribute directly to the scripting. Additionally, competing coursework and time constraints prevented us from implementing some of our more ambitious ideas, such as integrating MACS2 peak data. Despite these obstacles, we supported one another, adapted our approach when necessary, and were able to deliver a final product that reflects both persistence and growth.

# Personal Reflections
## Group Leader: Shameem
Group leader's reflection on the project

## Other members

### Stefanie: 
I really enjoyed working with Shameem and Justin for this project- they are both very strong teammates, and I definitely learned a lot from them. We met many times and for hours at a time, and I think we ultimately came out with a strong Gibbs sampling model. Despite the volume of data being too large for thealgorithm to manage, we were able to produce a sequence logo showing conservation in several key positions for the p53 motif. What really helped us was to learn more about the biology of p53, for instance, how p53 will bind on either forward or reverse strands, and that it is a tetramer, with each dimer recognizing one ~10bp half site. We also learned that a full p53 response unit consists of two half-sites with a 0-13bp spacer inbetween: RRRCWWGYYY – spacer – RRRCWWGYYY. Taking this into consideration when analyzing the sequence logo we obtained from running only 500 reads of data for 10,000 intervals, we were able to identify the C at position 4, purine enrichment at the 5' end, and pyrimidine enrichment at the 3' end, with significant ambiguity in the center, which is all characteristic of that p53 half-site motif. I thought it was rather impressive that we were able to build a model strong enough to detect these signals from only 500bp of data.

When we implemented the code snippet to read in the fna data for locating the Shyne Delgarno motif, we were pleased to see a strong SD signal come through, further supporting the strength of our model. Despite this, we still wanted to investigate the p53 data further to determine if there were higher quality reads we could preferentially select our data from. When we ran macs2 on the bam datafile, the cross correlation peak was shallow, indicating low signal-to-noise distinction. At this point, we felt the data had more issues than we had been trained to address and for the scope of this assignment, so we returned our focus to a successful Gibbs Sampling algorithm. 

### Justin: 
Working with Stefanie and Shameem for this project was really pleasant, and I think that already knowing each other from our Genomics class last semester really helped us work well together. I found that a lot of our issues with this project were really about the size of the data in the alignment, as processing millions of reads, especially with a Gibbs Sampling approach where we're making random selections and exhaustively scoring motifs, threatened to make our runtime unfeasibly long if we weren't careful. I think that our process of coming up with the pseudocode for this project was really interesting, as it came down to both figuring out how to implement a Gibbs Sampling algorithm and being mindful of how our implementation choices would affect our program's efficiency. While I really would've liked to have used macs2 peak data to narrow down the sequences we used, I'm still really happy with the sub-sampling approach we came up with for handling the K562 data, and we likely didn't have the time to implement everything with macs2 because of our other work that we needed to take care of. While I was very involved in the pseudocode process, I feel bad for not getting to help with the python scripting. Between a busy schedule and my computer being unreasonably uncooperative (syncing my fork while in Pycharm would cause the ipynb file to corrupt, forcing it to be viewed as plaintext, with seemingly no long-term solution other than uninstalling the IDE and deleting my fork completely), I did not get many chances to help out with the code. Thankfully, Stefanie was able to get a lot of the scripting done, and I massively am massively appreciative of the effort she put in. I'm grateful that my team was so sympathetic to the difficulties I was having getting my environment working and that I was still able to help out earlier on during the planning stages. I'm really happy both with how our project came together in the end and that my environment finally works without breaking the second it attempts to interact with a github repository. 

## Shameem:
I found this project to be a real challenge in a positive sense. Working through the pseudocode and collaborating on how to manage the data size and code were learning experiences I found very valuable. I also knew Stefanie and Justin from previous classes, so it was nice working with them again. The hardest part was definitely understanding the MCMC concept and taking the content from the slides and converting it into a mental map for our pseudocode, but after that initial hurdle, the rest was not so bad. We had a lot of ideas that bounced around for how to handle the data, such as macs2 or a sliding window, but our final product is something we are proud of, despite the challenges we faced. Personally, I found myself tested a lot on my conceptual understanding, and I really appreciate the help from my group mates. Although I didn’t do as much of the scripting as I wished I had, I spent a lot of time going through our given scripts and making a layout for our group, and that foundational work is something that reinforced the material for me. Overall, I am happy with how this project went and really grateful we had an extra week to work on this. 


# Generative AI Appendix
As per the syllabus, the AI we used was _______, and the prompts we used are shown below: 
Prompts go here (if applicable)

If no ai use: Generative AI was not used on this project
