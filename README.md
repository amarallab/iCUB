# iCUB
Information theory based codon usage bias
This python script calculates iCUB of a given sequence.

**Parameters to run iCUB in terminal:**

        $ python nullseq.py [-h] [-s] [Sinput] [-o] [-t] [-b] [-e]
        

|  Argument  | Description | Other|
| :--- | :--- | :--- |
| -s  | Change to single nucleotide sequence mode  | Optional
| Sinput  | path/to/fasta/file or nucleotide sequence if using -m flag  |
| -o  | /output/file/Path  | Default: calculated_icub.csv |
| -t TransTable  | NCBI translation table  | Default: 11 |
| -b  | input sequenes do not inlcude start codons | Optional |
| -e  | input sequenes do not inlcude stop codons | Optional |


Example:

1. Calculate iCUB for a single nucleotide seqeunce:

        $ python iCUB.py -s nucloeotidesequence

    specifications on the sequence can be added using the -b -e and -t flags

        $ python iCUB.py -sbe nucloeotidesequence -t 10

    the above indicates that the sequence should be calculated using
    the translation table 10 and that the sequence does not have stop
    or start codon in it

2. Calculate iCUB for nucleotide sequences in a fasta file

        $ python iCUB.py /path/to/fasta/file

    specifications regarding the sequences can be added using the -b -e and -t flags

         $ python iCUB.py -be /path/to/fasta/file -t 10 -o /path/to/output/file

    the above indicates that all sequences should be calculated using
    the translation table 10 and that the sequences do not have stop
    or start codon in it and the file with the calcualted iCUB should
    be saved to /path/to/output/file

**Using iCUB as a package:**

        >>> import iCUB as iCUB
        >>> nucleotidesequence = "ATG...TGG"
        >>> iCUB.iCUB_Calculator(nucleotidesequence).get_iCUB()

If the sequences does not have start or stop codons or does not use translation
table 11, additonal attributes will need to be specified

Using other translation tables:

        >>> tanslationtable = 5
        >>> iCUB.iCUB_Calculator(nucleotidesequence, n=translationtable).get_iCUB()

Sequences without start codons:

        >>> nucleotidesequence = "GCG...TGG"
        >>> iCUB.iCUB_Calculator(nucleotidesequence, start=False).get_iCUB()

Sequences without stop codons:

        >>> nucleotidesequence = "ATG...AAA"
        >>> iCUB.iCUB_Calculator(nucleotidesequence, stop=False).get_iCUB()
