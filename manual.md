# <samp>disruptiveMutations</samp> - manual

***

## Usage
```
disruptiveMutations.py -fsa <in_fasta> -pfx <out_prefix> [options]
```

## Version
```
disruptiveMutations.py 0.0.1 (alpha)
```

## Dependencies
```Python v3.9.7```, ```ViennaRNA v2.5.0```

## Description
<samp>disruptiveMutations</samp> reduces the interaction strength between a particular area of interest in the first sequence in relation to all other sequences. Ideally, the tool interrupts any potential interaction between this area and all other sequences. Individual predictions are performed with the RNAcofold python site-package of the ViennaRNA Package 2.5.0. Example call: python disruptiveMutations.py -pfx example -fsa example.fa -pss 32 -pse 96 -thr 4 -ovr 

## Options
```
--prefix,-pfx
    output directory and prefix for result files

--fasta,-fsa
    fasta with all sequences; the first sequence has to be the short sequence 
    snipped which has to be deoptimized with all following sequences

--positionStart,-pss
    define the starting position of the sequence part which should be 
    deoptimized (default: 0)

--positionEnd,-pse
    define the ending position of the sequence part which should be deoptimized;
    a 0 will be interpreted as the end of the sequence (default: 0)

--seed,-sed
    set the seed for random (default: 0)

--reverse,-rev
    creates reverse of each strain if set (default: False)

--complement,-cmp
    creates complements of each strain if set (default: False)

--threads,-thr
    number of threads to use for RNAcofold (default: 1)

--slice,-slc
    window size for standard window search or minimum length for a seed to be
    accepted if the seed option is specified (default: 20)

--candidateIncrease,-cdi
    increases the energy of the top candidate if there is no possible mutation
    candidate (default: 1.0)

--candidateStop,-cds
    stopping energy for deoptimization (default: -10.0)

--candidateIteration,-cdt
    tries to increase the cdi value at high iteration numbers (default: 1000)

--overwrite,-ovr
    overwrite data with folder named in prefix (default: False)

--dangles,-dng
    use dangling ends for foldings (default: 2) (choices: 0,1,2,3)

--noLP,-nlp
    disable lonely pairs for RNAcofold (default: False)
```

## References
```
D. Desirò, A. Borodavka, and M. Marz.
"DisruptiveMutations: Disrupting functional long-range RNA-RNA interactions in RNA viruses."
In Preparation, 2022.
https://github.com/desiro/disruptiveMutations
```