'''
File: iCUB.py
Author: Sophia Liu
Description: The contains a iCUB_Calculator object that
             Calculates the iCUB of a nucleotide sequence
'''

from Bio.Seq import Seq
import collections
from Bio.Data import CodonTable
import math
import argparse
import os.path
import csv

def Codons_for_AA(n):
    ''' Determines the codons used to code an amino acid given a codon table

    Parameters
    ----------
    n : int
        the codon table index according to NCBI

    Returns
    -------
    CodonforAA : dictionary
        Key : str
            Amino acids in Single letter notation (caps)
        Values : str
            Three letter codon (in caps)
    '''

    GivenCodonTable = CodonTable.unambiguous_dna_by_id[n]
    nucleotides = ['A', 'T', 'C', 'G']
    CodonforAA = {}
    for first in nucleotides:
        for second in nucleotides:
            for third in nucleotides:
                Codon = first + second + third
                if Codon not in CodonTable.unambiguous_dna_by_id[n].stop_codons:
                    if GivenCodonTable.forward_table[Codon] in CodonforAA.keys():
                            CodonforAA[GivenCodonTable.forward_table[Codon]].append(Codon)
                    else:
                        CodonforAA[GivenCodonTable.forward_table[Codon]] = [Codon]
                else:
                    pass
    return(CodonforAA)

def GC_Content_for_Codons(n):
    ''' Determines the number of G/C nucletides for all codons in the
    codon table

    Parameters
    ----------


    n : int
        The codon table index according to NCBI

    Returns
    -------
    CodonGCContent : dictionary
        Key : str
            Three letter codon (in caps)
        Values : int
            Number of G/C nucleotides in the codon
    '''

    CodonforAA = Codons_for_AA(n)
    CodonGCContent = collections.defaultdict(int)
    for AA in CodonforAA:
        for codon in CodonforAA[AA]:
            codonusage = collections.Counter(codon)
            CodonGCContent[codon] = codonusage['G'] + codonusage['C']

    return CodonGCContent

def Probability_Given_Beta(beta, n):
    ''' Gives the probability of using each codon given beta.
        The probabilites are determined using the bolzmann distribtion
        where P = exp(-beta*N)/Z
        N is the number of G/C nucleotides in the codon

    Parameters:
    ----------
    beta : float
        The constant beta used in the calulations of the probability

    n : int
        The codon table index according to NCBI

    Returns:
    -------
    probs : dictionary
        Key : str
            Codons  (caps)
        Values : float
            The probability of using the codon (P(codon|AA))
    '''

    probs = {}
    CodonforAA = Codons_for_AA(n)
    CodonGC = GC_Content_for_Codons(n)

    for AA in CodonforAA:
        Z = sum([math.exp(-beta*CodonGC[codon]) for codon in CodonforAA[AA]])
        for codon in CodonforAA[AA]:
            probs[codon]= math.exp(-beta*CodonGC[codon])/Z
    return probs

class iCUB_Calculator(object):
    '''
        Parameters
        ----------
        Sequence : str
            The nucleotide sequence of coding DNA sequence
            Must be in all caps
        n : int, optional
            The codon table index according to NCBI
            As defalut n is set to 11, the standard codon table
        start : bool (opt)
            True when start codon present in sequence
            Default: True
        stop : book (opt)
            True when stop codon present in sequence
            Default :True
        '''
    def __init__(self, Sequence, n=11, start=True, stop=True):
        ''' Initialization of the object

            Attributes:
            -----------
            n : int
                codon table index accoriding to NCBI
            sequences : str
                sequences
        '''

        if not isinstance(Sequence, str):
            raise TypeError("Should be a string")
        if len(Sequence)%3 != 0:
            raise ValueError("Sequence length should be mutiple of 3")

        self.n = n
        self.CodonforAA = Codons_for_AA(self.n)
        if start:
            if stop:
                self.sequence = Sequence[3:-3]
            else:
                self.sequence = Sequence[3:]
        else:
            if stop:
                self.sequence = Sequence[:-3]
            else:
                self.sequence = Sequence
        self.start = start
        self.stop = stop
        self.seq = Seq(self.sequence)
        self.AA = str(self.seq.translate(table=n))
        self.GCcontent = None
        self.AAcount = None
        self.AApercent = None
        self.CodonCount = None
        self.CodonProbabiltiy = None

    def get_AA_count(self):  # does not include first AA (Met)
        ''' Determines the number of times each amino acid is used in the gene

        Returns
        -------
        AAcount : dictionary
            Key : str
                Amino acids in Single letter notation (caps)
            Values : int
                Number of times the amino acid is used

        Notes
        -----
        Does not include start and stop codons
        '''

        if self.AAcount is None:

            CodonsforAA = self.CodonforAA
            TempAAUsageDict = collections.Counter(self.AA)
            AAUsageDict = {}
            for AA in CodonsforAA.keys():
                if AA not in TempAAUsageDict.keys():
                    AAUsageDict[AA] = 0
                else:
                    if AA in CodonsforAA.keys():
                        AAUsageDict[AA] = TempAAUsageDict[AA]
                    else:
                        pass
            self.AAcount = AAUsageDict
        else:
            pass

        return self.AAcount

    def get_AA_percent(self):
        '''Determines the percentage that each amino acid is used

        Returns
        -------
        AApercent : dictionary
            Key : str
                Amino acids in Single letter notation (caps)
            Values : float
                Percent usage of the amino acid in the sequence
        '''

        if self.AApercent is None:
            aacounts = self.get_AA_count()
            percentages = {}
            for aa in aacounts:
                percentages[aa] = aacounts[aa]/len(self.AA)
            self.AApercent = percentages
        else:
            pass

        return self.AApercent

    def get_codon_count(self):
        ''' Determines the number of times each codon is used

        Returns
        -------
        CodonCount : dictionary
            keys : str
                Three letter codon (in caps)
            values : float
                number of times the codon is used in the sequence
        '''

        if self.CodonCount is None:
            GivenCodonTable = CodonTable.unambiguous_dna_by_id[self.n]
            CodonList = [self.sequence[x:x+3] for x in range(0,len(self.sequence),3)]
            CodonUsageDict = collections.Counter(CodonList)
            for codon in GivenCodonTable.forward_table.keys():
                if codon not in CodonUsageDict.keys():
                    CodonUsageDict[codon] = 0
            self.CodonCount = CodonUsageDict

        return self.CodonCount

    def get_codon_usage_probability(self):
        ''' Determines the probability the codon is used when
        coding for its respective Amino acid (P(codon|AA))

        The probability is determined by:
        codon usage / total amino acid usage
        in a given sequence

        Returns
        -------
        CodonProbabiltiy : dictionary
            keys : str
                Three letter codon (in caps)
            values : float
                probability P(codon|AA)

        Notes:
        -----
        If the amino acid is not used, then equal codon usage is assumed
        '''

        if self.CodonProbabiltiy is None:
            GivenCodonTable = CodonTable.unambiguous_dna_by_id[self.n]
            CodonforAA = self.CodonforAA
            CodonUsage = self.get_codon_count()
            AAUsage = self.get_AA_count()
            UsageProbabilityforCodons = {}
            for Codon in GivenCodonTable.forward_table.keys():
                if AAUsage[GivenCodonTable.forward_table[Codon]] != 0:
                    UsageProbabilityforCodons[Codon] = CodonUsage[Codon]/AAUsage[GivenCodonTable.forward_table[Codon]]
                else:
                    UsageProbabilityforCodons[Codon] = 1 /len(CodonforAA[GivenCodonTable.forward_table[Codon]])
            self.CodonProbabiltiy = UsageProbabilityforCodons
        else:
            pass

        return self.CodonProbabiltiy

    def get_GC_Content(self):
        ''' Determines the G/C content of the gene

        Returns
        -------
        GCcontent : float
            The number of G/C nucleotides in the sequence

        Notes
        -----
        Does not include start and stop codons
        '''

        if self.GCcontent is None:
            NucleotideList = [self.sequence[x] for x in range(0, len(self.sequence))]
            NucleotideListDict = collections.Counter(NucleotideList)
            GCcontent =  (NucleotideListDict['G'] + NucleotideListDict['C'])
            self.GCcontent = GCcontent
        else:
            pass

        return self.GCcontent

    def compute_average_GC(self, b, given='beta'):
        '''Determines the average GC of the given sequence given codon
           probablities or a given beta

        Parameters:
        -----------
        b : float or dictionary
            When calulating for a given beta b is a float.
            When calulating for a codon distribution is a dictionary
                key : str
                    codon
                Values : float
                    P(codon)
        given : str (optional)
            When given is 'beta' the fuction calculates average GC for beta.
            When given is 'probability' the the fuction calculates average GC
            for a given codon usage proability

        Returns:
        --------
        sum(GCnumber) : float
            The expected number GC nucleotides for the given sequence given
            the conditions
        '''

        CodonGC = GC_Content_for_Codons(self.n)
        GCnumber = []
        length = (len(self.sequence)/3)-2

        if given == 'beta':
            Probability = Probability_Given_Beta(b, self.n)
            AAcount = self.get_AA_count()
            CodonforAA = self.CodonforAA

            for AA in CodonforAA:
                GCnumber.append(AAcount[AA]*sum([Probability[codon]*CodonGC[codon] for codon in CodonforAA[AA]]))

        elif given == 'probability':
            N = len(self.AA)
            GCnumber = [CodonGC[codon]*b[codon]*(N+20) for codon in b.keys()]


        return sum(GCnumber)

    def get_beta(self, given=None):
        '''Determines the value of beta given the GC content of the sequence

        Parameters:
        -----------
        given : float (optional)
            When given is none, the sequence beta is obtained from the sequence GC.
            Otherwise, a GC content can be speficied

        Returns:
        -------
        beta : float
            the constant in the equation P = exp(-beta*N)/Z
        '''

        if given is None:
            m = self.get_GC_Content()
        else:
            m = given

        rl = -20
        rr = 20

        avl = self.compute_average_GC(rl) - m
        avr = self.compute_average_GC(rr) - m

        if avl*avr > 0:
            if avl < 0:
                return rl
            else:
                return rr
        else:

            rm = (rl+rr)/2
            avv = self.compute_average_GC(rm) - m

            while math.fabs(avv) > 1e-7:
                if avv*avl > 0:
                    rl = rm
                else:
                    rr = rm
                rm =  0.5*(rl+rr)
                avv = self.compute_average_GC(rm) - m

        return rm

    def get_Info(self, prob=None):
        '''Determines raw shannon information of sequence
            H_aa = -P(codon|AA) * log(P(codon|AA),2)
            H_total = (sigma) P(AA) * H_aa

        Parameters:
        -----------
        prob : optional
            None :
                calculates the Nc_Info pertaining to the P(codon|AA) of the sequence
            dictionary :
                calculates the Nc_info values pertaining to the P(codon|AA) given

                Keys : str
                    codon in caps
                Values : float
                    P(codon|AA)

        Returns:
        --------
        H : float
            (sigma) P(AA) * H_aa
        '''

        CodonforAA = self.CodonforAA
        AAprob = self.get_AA_percent()

        if prob is None:
            CodonProb = self.get_codon_usage_probability()
        else:
            CodonProb = prob

        Haa = {}
        for AA in CodonforAA:
            Haa[AA] = sum([(-CodonProb[codon]*math.log(CodonProb[codon],2)) for codon in CodonforAA[AA] if CodonProb[codon] != 0])

        H = 0
        for AA in Haa:
            H += Haa[AA]*AAprob[AA]
        return H


    def get_iCUB(self):
        '''Calcuates the iCUB of sequence

        Returns:
        --------
        iCUB : float
        '''

        Nc = self.get_Info()
        Ncnull = self.get_Info(prob=Probability_Given_Beta(self.get_beta(), self.n))
        gamma = Nc/Ncnull
        iCUB = 20 + gamma*(61-20)

        return iCUB

def main():
    pass
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculates iCUB of nucleotide sequence')
    parser.add_argument('-s', action='store_true', help='change to single sequence mode')
    parser.add_argument('Sinput', help='file or sequence input')
    parser.add_argument('-o', default='calculated_icub.csv', help='output file path')
    parser.add_argument('-t', type=int, metavar='TransTable', default=11,
                        choices=[1, 2, 3, 4, 5,
                                 6, 9, 10, 11, 12,
                                 13, 14, 16, 21, 22,
                                 23, 24, 25],
                        help='NCBI translation Table')
    parser.add_argument('-b', action='store_false', help='input sequence/s does not include start codon')
    parser.add_argument('-e', action='store_false', help='input sequence/s does not include stop codon')
    args = parser.parse_args()

    if args.s:
        if not isinstance(args.Sinput, str):
            raise TypeError('Input nucleotide sequence must be string')
        else:
            print(iCUB_Calculator(args.Sinput, n=args.t, start=args.b, stop=args.e).get_iCUB())
    else:
        if os.path.isfile(args.Sinput):
            if args.Sinput.split('.')[-1] == 'fasta':
                with open(args.Sinput, 'r') as f:
                    lines = f.readlines()
                with open(args.o, 'w') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    csvwriter.writerow(['gene','iCUB'])
                    for i in range(len(lines)):
                        if lines[i][0] == '>':
                            gene = lines[i].split(',')[0][1:]
                            icub = iCUB_Calculator(lines[i+1].rstrip('\n'), n=args.t, start=args.b, stop=args.e).get_iCUB()
                            csvwriter.writerow([gene,icub])
                        else:
                            pass

            else:
                raise ValueError('Improper file type')
        else:
            raise ValueError('File does not exist')








