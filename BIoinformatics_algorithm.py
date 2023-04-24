
'''
# Checks if DNA sequence is valid. Returns True is sequence is valid, or False otherwise. 
#dna_seq = input("enter dna seq: ")
#print(dna_seq)
def validate_dna (dna_seq):
    seq= dna_seq.upper()
    valid= seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")
    if valid == len(seq):  
        return True
    else: 
        return False
# test
#print(validate_dna("atagagagatctcg"))
#print(validate_dna("ATAGAXTAGAT"))



# Calculates the frequency of each symbol in the sequence.Returns a dictionary
def frequency (seq):
    dic={}  
    for s in seq.upper():
        if s in dic : 
            dic[s] +=1
        else:
            dic[s] =1
    return dic


#print(frequency("atagataactcgcatag")) 
#print(frequency("MVVMKKSHHVLHSQSLIK"))



# calculate the frequency of the aminoacids in a sequence read from the users input
seq_aa= input("Protein seq: ", )
freq_seq= frequency(seq_aa)
list_f = sorted(freq_seq.items(), key=lambda x:x[1], reverse = True)
for (k,v) in list_f:
    print("Aminoacid:", k, ":", v)




# GC content (Returns percentage of G and C nucleotides in a DNA sequence)
def gc_content(dna_seq):
    gc_content=0
    for s in dna_seq:
        if s in "GCgc":
            gc_content +=1
    return gc_content / len(dna_seq)
    

# print(gc_content("atagagagatctcgatagagagatctcgatagagagatctcg"))


# computing the GC content of the non-overlapping sub-sequences of size k of the inputted sequence (Returns GC content of non−overlapping sub−sequences of size k The result is a list)
def gc_content_subseq(dna_seq, k=100):
    res =[]
    for i in range(0, len(dna_seq)-k+1, k):
        subseq = dna_seq[i:i+k]
        gc = gc_content(subseq)
        res.append(gc)
    return res

#print(gc_content_subseq("atagagagatctcgatagagagatagagagatctcgataagagagatctcgataagagagatctcgataagagagatctcgatactcgatagagagatctcg"))



# Function that computes the RNA corresponding to the transcription of the DNA sequence provided
def translation(dna_seq):
    assert validate_dna(dna_seq)
    return dna_seq.upper().replace("T" , "U")
# print(translation("atagagagatctcgatagagagatctcgatagagagatctcg"))



# Computes the reverse complement of the DNA sequence
def reverse_complement(dna_seq):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    comp =""
    for c in dna_seq.upper():
        if c == "A":
            comp = "T" + comp
        elif c == "T":
            comp = "A" + comp
        elif c == "G":
            comp = "C" + comp
        elif c== "C":
            comp = "G" + comp
    return comp



# Translates a codon into an aminoacid using an internal dictionary with the standard genetic code
def translation_codon(code):
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "TGT":"C", "TGC":"C",
    "GAT":"D", "GAC":"D",
    "GAA":"E", "GAG":"E",
    "TTT":"F", "TTC":"F",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "CAT":"H", "CAC":"H",
    "ATA":"I", "ATT":"I", "ATC":"I",
    "AAA":"K", "AAG":"K",
    "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "ATG":"M", "AAT":"N", "AAC":"N",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "TGG":"W",
    "TAT":"Y", "GGT":"Y",
    "TAA":"_", "TAG":"_", "TGA":"_"}
    if code in tc :
        return tc[code]
    else:
        return None
# print(translation_codon("GTC"))


#Translates a DNA sequence into an aminoacid sequence
# ?
def translate_seq (dna_seq, ini_pos = 0):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    seq_aa = ""
    for pos in range(ini_pos, len (seqm)-2,3):
        cod = seqm[pos:pos+3]
        seq_aa += translate_codon(cod)
    return seq_aa


# Provides the frequency of each codon encoding a given aminoacid, in a DNA sequence
#?
def codon_usage(dna_seq, aa):
    assert validate_dna.upper()
    dic={}
    total = 0
    for i in range(0, len(seqm) - 2, 3)
        cod = seqm[i:i + 3]
        if translate_codon(cod) == aa:
            if cod in dic :
                dic[cod] +=1
            else:
                dic[cod] = 1
                total =1
        elif total > 0:
            for k in dic:
                dic[k] /=total
        return dic    
    

#Computes the six reading frames of a DNA sequence (including the reverse complement
#?
def reading_frames (dna_seq):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    res.append(translate_seq(dna_seq ,0))
    res.append(translate_seq(dna_seq ,1))
    res.append(translate_seq(dna_seq ,2))
    rc = reverse_complement(dna_seq)
    res.append(translate_seq(rc,0))
    res.append(translate_seq(rc,1))
    res.append(translate_seq(rc,2))
    return res

        

#Computes all possible proteins in an aminoacid sequence. Returns list of possible proteins
#?
def all_proteins_rf (aa_seq):
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
            else :
                if aa == "M":
                    current_prot.append("")
                for i in range( len (current_prot)):
                    current_prot[i] += aa
    return proteins


#Computes all possible proteins for all open reading frames
#?
def all_orfs (dna_seq):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames(dna_seq)
    res = []
    for rf in rfs :
        prots = all_proteins_rf(rf)
        for p in prots :
            res.append(p)
    return res



#Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size
#?
def all_orfs_ord (dna_seq, minsize = 0):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames(dna_seq)
    res =[]
    for rf in rfs :
        prots = all_proteins_rf(rf)
        for p in prots:
            if len (p) > minsize: insert_prot_ord(p, res)
    return res
def insert_prot_ord (prot, list_prots):
    i=0
    while i < len (list_prots) and len (prot) < len (list_prots[i]):
        i += 1
    list_prots.insert(i, prot)


from sequences import ∗
seq = input ("Insert DNA sequence:")
if validate_dna(seq):
    print ("Valid sequence")
    print ("Transcription: ", transcription (seq))
    print ("Reverse complement:", reverse_complement(seq))
    print ("GC content (global):", gc_content(seq))
    print ("Direct translation:" , translate_seq(seq))
    print ("All proteins in ORFs (decreasing size): ", all_orfs_ord(seq))
else : 
    print ("DNA sequence is not valid")


#Reads a sequence from a multi−line text file
def read_seq_from_file(filename):
    fh = open(filename, "r")
    lines = fh.readlines()
    seq = ""
    for l in lines:
        seq += l.replace("\n","")
    fh.close()
    return seq


fname = input ("Insert input filename:")
seq = read_seq_from_file(fname)
if validate_dna (seq):
    print ("Valid sequence")
    print ("Transcription: ", transcription (seq))
    print ("Reverse complement:", reverse_complement(seq))
    print ("GC content (global):", gc_content(seq))
    print ("Direct translation:" , translate_seq(seq))
    orfs = all_orfs_ord(seq)
    i=1
    for orf in orfs:
        write_seq_to_file(orf, "orf−"+ str (i)+".txt")
        i += 1
else : 
    print ("DNA sequence is not valid")


# Class for biological sequences
#?
class MySeq:
    def __init__ ( self , seq, seq_type = "DNA"):
        self.seq = seq.upper()
        self.seq_type = seq_type

    def __len__( self ):
        return len ( self.seq)
    
    def __getitem__( self , n):
        return self.seq[n]
    def __getslice__( self , i, j):
        return self.seq[i:j]

    def __str__( self ):
        return self.seq
    
    def get_seq_biotype ( self ):
        return self.seq_type
    
    def show_info_seq ( self ):
        print ("Sequence: " + self.seq + " biotype: " + self.seq_type)

    def alphabet ( self ):
        if ( self.seq_type=="DNA"): 
            return "ACGT"
        elif ( self.seq_type == "RNA"):
            return "ACGU"
        elif ( self.seq_type=="PROTEIN"): 
            return "ACDEFGHIKLMNPQRSTVWY"
        else : 
            return None
    
    def validate ( self ):
        alp = self.alphabet()
        res = True
        i=0
        while i < len ( self.seq) and res:
            if self.seq[i] not in alp: 
                res = False
            else : i += 1
        return res

    def transcription ( self ):
        if ( self.seq_type == "DNA"):
            return MySeq( self.seq.replace("T","U"), "RNA")
        else :
            return None
        
    def reverse_comp ( self ):
        if ( self.seq_type != "DNA"): 
            return None
        comp = ""
        for c in self.seq:
            if (c == 'A'): 
                 comp = "T" + comp
            elif (c == "T"): 
                comp = "A" + comp 
            elif (c == "G"): 
                comp = "C" + comp
            elif (c== "C"): 
                comp = "G" + comp
        return MySeq(comp, "DNA")

    def translate ( self , iniPos= 0):
        if ( self.seq_type != "DNA"): 
            return None
        seq_aa = ""
        for pos in range (iniPos, len ( self.seq)-2,3):
            cod = self.seq[pos:pos+3]
            seq_aa += translate_codon(cod)
        return MySeq(seq_aa, "PROTEIN")
    

#s1 = MySeq("ATGTGATAAGAATAGAATGCTGAATAAATAGAATGACAT")
#s2 = MySeq("MKVVLSVQERSVVSLL", "PROTEIN")
#print (s1.validate(), s2.validate())
#print (s1)
#s3 = s1.transcription()
#s3.show_info_seq()
#s4 = s1.reverse_comp().translate()
#s4.show_info_seq()
'''

# Naive Algorithm for Fixed Pattern Finding
"""Both variants are implemented in the Python functions provided in the following code block.
The search_first_occ function uses an outer while cycle that finishes when the pattern is
found, or the last sub-sequence is tested. The function returns the position of the first occurrence of the pattern, or if the pattern does not occur returns −1.
The search_all_occurrences function uses a for cycle that tests for all possible subsequences (note that there are N − k + 1 possible positions where the pattern may occur),
and returns a list with all initial positions of the pattern’s occurrences (the list is empty if the
pattern does not occur)."""

def  search_first_occ(seq, pattern):
    found = False
    i = 0
    while i <= len(seq) - len(pattern) and not found:
        j =0
        while j < len (pattern) and pattern[j]==seq[i+j]:
            j=j+1
        if j== len (pattern): 
            found = True
        else : 
            i += 1
        if found: 
            return i
        else : 
            return -1


def search_all_occurrences(seq, pattern):
    res = []
    for i in range( len (seq)-len (pattern)+1):
        j=0
        while j < len (pattern) and pattern[j]==seq[i+j]:
            j=j+1
        if j == len (pattern):
            res.append(i)
    return res


#seqDNA = "ATAGAATAGATAATAGTC"
#print ( search_first_occ(seqDNA, "GAAT") )
#print ( search_first_occ(seqDNA, "TATA") )
#print ( search_all_occurrences(seqDNA, "AAT") )


def test_pat_search():
    seq = input ("Input sequence: ")
    pat = input ("Input pattern: ")
    print (pat, " occurs in the following positions:", )
    print ( search_all_occurrences(seq, pat) )


#test_pat_search()
"""
In the case of the bad-character rule, we create a dictionary with all possible symbols in the
alphabet as keys, and with values defining the rightmost position where the symbol occurs
in the pattern (−1 if the symbol does not occur). This allows to rapidly calculate the number
of positions to move forward according to this rule by calculating the offset: position of the
mismatch in the pattern – value for the symbol in the dictionary. Notice that this value might
be negative and, in this case, this means the rule does not help and it will be ignored in that
iteration. This process is done in the process_bcr function in the code given below.
The pre-processing for the good suffix rule is more complex and we will not explain here all
the details (the code is given below and we leave a detailed analysis to the interested reader).
The result of this process is to create a list that keeps the number of positions that may be
moved forward, depending on the position of the mismatch on the pattern (list index). Notice
112 Chapter 5
that in this process, both the situations illustrated above need to be taken into account. This
process is done in the process_gsr function in the code given below.
The implementation of the algorithm is given next as a Python class. The class allows to define an alphabet and a pattern in the constructor and does the pre-processing for both rules,
according to the pattern, in the function preprocess called by the constructor.
The function search_pattern allows to use an initialized object of this class to search over
target sequences for the given pattern. It is an adaptation of the naive algorithm given in the
previous section, but which makes use of the data structures from the pre-processing, using
the two rules to move forward the maximum number of allowed positions. In the worst case,
it advances a single position (as in the naive algorithm), but in other cases it can use one of
the rules to move forward more positions (the maximum of the values provided by each of the
rules).

"""



class BoyerMoore:
    def __init__( self , alphabet, pattern):
        self.alphabet = alphabet
        self.pattern = pattern
        self.preprocess()

    def preprocess( self ):
        self.process_bcr()
        self.process_gsr()

    def process_bcr( self ):
        self.occ = {}
        for symb in self.alphabet:
            self.occ[symb] = -1
        for j in range( len ( self.pattern)):
            c = self.pattern[j]
        self.occ[c] = j

    def process_gsr( self ):
        self.f = [0] * ( len ( self.pattern)+1)
        self.s = [0] * ( len ( self.pattern)+1)
        i = len ( self.pattern)
        j = len ( self.pattern)+1
        self.f[i] = j
        while i>0:
            while j<= len ( self.pattern) and self.pattern[i-1] != self.pattern[j-1]:
                if self.s[j] == 0: 
                    self.s[j] = j-i;
                j = self.f[j]
            i -= 1
            j -= 1
            self.f[i] = j
        j = self.f[0]
        for i in range( len ( self.pattern)):
            if self.s[i] == 0: 
                self.s[i] = j
            if i == j: 
                j = self.f[j]

    def search_pattern( self , text):
        res = []
        i=0
        while i <= len (text) - len ( self.pattern):
            j= len ( self.pattern)- 1
            while j>=0 and self.pattern[j]==text[j+i]: j -= 1
            if (j<0):
                res.append(i)
                i += self.s[0]
            else :
                c = text[j+i]
                i += max( self.s[j+1], j- self.occ[c])
        return res

def test():
    bm = BoyerMoore("ACTG", "ACCA")
    print (bm.search_pattern("ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC"))

# test()



# Deterministic Finite Automata
#The function overlap(s,t) prostxvides the maximum overlap between sequences s and t, as defined above, being implemented in Python in the following way (note the implementation is a simple, not an efficient one)
def overlap(s1, s2):
    maxov = min( len (s1), len (s2))
    for i in range(maxov,0,-1):
        if s1[-i:] == s2[:i]: return i
    return 0




class Automata:
    def __init__( self , alphabet, pattern):
        self.numstates = len (pattern) + 1
        self.alphabet = alphabet
        self.transition_table = {}
        self.build_transition_table(pattern)

    def build_transition_table( self , pattern ):
        for q in range( self.numstates):
            for a in self.alphabet:
                prefix = pattern[0:q] + a
                self.transition_table[(q,a)] = overlap(prefix, pattern)

    def print_automata( self ):
        print ("States: " , self.numstates)
        print ("Alphabet: " , self.alphabet)
        print ("Transition table:")
        for k in self.transition_table.keys():
            print(k[0], ",", k[1], " -> ", self.transition_table[k])

    def next_state( self , current, symbol):
        return self.transition_table.get((current, symbol))
    
    def apply_seq( self , seq):
        q=0
        res = [q]
        for c in seq:
            q = self.next_state(q, c)
            res.append(q)
        return res

    def occurrences_pattern( self , text):
        q=0
        res = []
        for i in range( len (text)):
            q = self.next_state(q, text[i])
            if q == self.numstates-1:
                res.append(i - self.numstates + 2)
        return res
    

def test():
        auto = Automata("ACGT", "ACA")
        auto.print_automata()
        print (auto.apply_seq("CACATGACATG"))
        print (auto.occurrences_pattern("CACATGACATG"))
    
# test()



# Finding Flexible Patterns: Regular Expressions
"""
import re
str = "TGAAGTATGAGA"
mo = re.search("TAT", str )

#print (mo.group())
#print (mo.span())
mo2 = re.search("TG.", str )
#print (mo2.group())
#print (mo2.span())
#print(re.findall("TA.", str ))
#print(re.findall("TG.", str ))
mos = re.finditer("TG.", str )
for x in mos:
    print (x.group())
    print (x.span())
"""

"""This program may be used to test different REs and their occurrence in biological sequences
(DNA, RNA, proteins).
One important limitation of the previous function to identify all occurrences of a pattern
(find _all_occurrences_re) is the fact that it does not consider instances of the patterns that
overlap. To illustrate this consider the following example of the previous program:
Input sequence:ATATGAAGAG
Input pattern (as a regular expression):AT.
Pattern found in position: 0
Pattern found in positions: [0] """

def find_pattern_re (seq, pat):
    from re import search
    mo = search(pat, seq)
    if (mo != None):
        return mo.span()[0]
    else :
        return -1
def find_all_occurrences_re (seq, pat):
    from re import finditer
    mos = finditer(pat, seq)
    res = []
    for x in mos:
        res.append(x.span()[0])
    return res

def test():
    seq = input ("Input sequence:")
    pat = input ("Input pattern (as a regular expression):")
    res = find_pattern_re(seq, pat)
    if res >= 0:
        print ("Pattern found in position: ", res)
    else : 
        print ("Pattern not found")

    all_res = find_all_occurrences_re(seq, pat)
    if len (all_res) > 0:
        print ("Pattern found in positions: ", all_res)
    else : 
        print ("Pattern not found")

# test()

def find_all_overlap(seq, pat):
    return find_all_occurrences_re(seq, "(?="+pat+")")

def test():
    seq = input ("Input sequence:")
    pat = input ("Input pattern (as a regular expression):")

    all_ov = find_all_overlap(seq, pat)
    if len (all_ov) > 0:
        print ("Pattern found in positions: ", all_ov)
    else :
        print ("Pattern not found")

#test()
"""
# %%
import re
# %%
seq = "AAATAGAGATGAAGAGAGATAGCGC"
rgx = re.compile ("GA.A")
rgx.search(seq).group()

# %%
rgx.findall(seq)
# %%
mo = rgx.finditer(seq)
# %%
for x in mo: 
    print (x.span())
"""

# translate_codon function

def translate_codon_re (cod):
    import re
    if re.search("GC.", cod): aa = "A"
    elif re.search("TG[TC]", cod): aa = "C"
    elif re.search("GA[TC]", cod): aa = "D"
    elif re.search("GA[AG]", cod): aa = "E"
    elif re.search("TT[TC]", cod): aa = "F"
    elif re.search("GG.", cod): aa = "G"
    elif re.search("CA[TC]", cod): aa = "H"
    elif re.search("AT[TCA]", cod): aa = "I"
    elif re.search("AA[AG]", cod): aa = "K"
    elif re.search("TT[AG]|CT.", cod): aa = "L"
    elif re.search("ATG", cod): aa = "M"
    elif re.search("AA[TC]", cod): aa = "N"
    elif re.search("CC.", cod): aa = "P"
    elif re.search("CA[AG]", cod): aa = "Q"
    elif re.search("CG.|AG[AG]", cod): aa = "R"
    elif re.search("TC.|AG[TC]", cod): aa = "S"
    elif re.search("AC.", cod): aa = "T"
    elif re.search("GT.", cod): aa = "V"
    elif re.search("TGG", cod): aa = "W"
    elif re.search("TA[TC]", cod): aa = "Y"
    elif re.search("TA[AG]|TGA", cod): aa = "_";
    else : aa = ""
    return aa

    """To finish this set of simple examples, let us recall the problem of finding a putative protein in
a sequence of aminoacids. A protein may be defined as a pattern in the sequence, that starts
with symbol “M” and ends with a “_” (the symbol representing the translation of the stop
codon). Note that in the other intermediate symbols “M” can occur, but “_” cannot. So, the
regular expression to identify a putative protein can be defined as: “M[ˆ_]*_”. This is used in
the next code block to define a function which identifies the largest possible protein contained
in an inputted aminoacid sequence.
    """

def largest_protein_re (seq_prot):
    import re
    mos = re.finditer("M[^_]∗_", seq_prot)
    sizem = 0
    lprot = ""
    for x in mos:
        ini = x.span()[0]
        fin = x.span()[1]
        s = fin - ini + 1
        if s > sizem:
            lprot = x.group()
            sizem = s
    return lprot


# finding motif
    """An example is the “Zinc finger RING-type signature” (PS00518) motif, 
    which is represtxsented by “C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]”. 
    This means a pattern starting with  aminoacid “C”, followed by any aminoacid, aminoacid “H”, 
    any aminoacid, an aminoacid in  the group “LIVMFY”, aminoacid “C”, two aminoacids, aminoacid 
    “C” and an aminoacid in  the group [LIVMYA]. 
    """


def find_zync_finger(seq):
    from re import search
    regexp = "C.H.[LIVMFY]C.{2}C[LIVMYA]"
    mo = search(regexp, seq)
    if (mo != None):
        return mo.span()[0]
    else :
        return -1
"""
de""f test():
    seq = "HKMMLASCKHLLCLKCIVKLG"
    print (find_zync_finger(seq))

test()"""



def find_prosite(seq, profile):
    from re import search
    regexp = profile.replace("-" , "")
    regexp = regexp.replace("x" , ".")
    regexp = regexp.replace("(","{")
    regexp = regexp.replace(")","}")
    mo = search(regexp, seq)
    if mo != None :
        return mo.span()[0]
    else:
        return -1
    
def test():
    seq = "HKMMLASCKHLLCLKCIVKLG"
    print (find_prosite(seq,"C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]"))

# test()

    """Given strings in this flexible alphabet, and being the purpose 
    to find their occurrences in DNA  sequences, the first task is to convert 
    strings written in this alphabet to regular expressions that  can be used to search over sequences
    """

def iub_to_RE (iub):
    dic = {"A":"A", "C":"C", "G":"G", "T":"T", "R":"[GA]", "Y":"[CT]"
    , "M":"[AC]", "K":"[GT]", "S":"[GC]", "W": "[AT]", "B":"[CGT]",
    "D":"[AGT]", "H":"[ACT]", "V":"[ACG]", "N":"[ACGT]"}

    site =  iub.replace("^","")
    regexp = ""

    for c in site:
        regexp += dic[c]
        
    return regexp

def test():
    print (iub_to_RE("G^AMTV"))
#test()

    """Given this function, we can now proceed to write functions to detect where a given enzyme
will cut a given DNA sequence, and also to calculate the resulting sub-sequences after the cut
(restriction map). These tasks are achieved by the functions cut_positions and cut_subsestxquences,
 respectively, provided in the code below.
    """

def cut_positions (enzyme, sequence):
    from re import finditer
    cutpos = enzyme.find("^")
    regexp = iub_to_RE(enzyme)

    matches = finditer(regexp, sequence)
    locs = [ ]
    for m in matches:
        locs.append(m.start() + cutpos)
    return locs

def cut_subsequences (locs, sequence):
    res = []
    positions = locs
    positions.insert(0,0)
    positions.append( len (sequence))
    for i in range( len (positions)-1):
        res.append(sequence[positions[i]:positions[i+1]])
    return res

def test():
    pos = cut_positions("G^ATTC", "GTAGAAGATTCTGAGATCGATTC")
    print (pos)
    print (cut_subsequences(pos, "GTAGAAGATTCTGAGATCGATTC"))
# test()


