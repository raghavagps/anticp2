##############################################################################
# AntiCP2 is developed for predicting, desigining and scanning anticancer    #
# peptides. It is developed by Prof G. P. S. Raghava's group. Please cite    #
# Agrawal et al., (2020) AntiCP 2.0: an updated model for predicting         #
# anticancer peptides. Briefings in Bioinformatics, doi: 10.1093/bib/bbaa153 #
# ############################################################################
import argparse  
import warnings
import pickle
import os
#import numpy as np

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-j", "--job",type=int, choices = [1,2,3], help="Job Type: 1:predict, 2:design and 3:scan, by default 1 (predict)")
parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1 by default 0.5")
parser.add_argument("-m","--model",type=int, choices = [1, 2], help="Model: 1: ACP/AMP, 2: ACP/non-ACP, by default 1")
parser.add_argument("-w","--winleng", type=int, choices =range(5, 30), help="Window Length: 5 to 30 (scan mode only), by default 10")
parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:Anticancer peptide, 2: All peptides, by default 1")
args = parser.parse_args()

#Model1
#DPC N5
#AAC N5 
# Function for generating all possible mutants
def seq_mutants(aa_seq):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    aa_seq = aa_seq.upper()
    mut_seq=[]
    for pos in range(0,len(aa_seq)):
        for aa in std:
            mut_seq += [aa_seq[:pos] + aa + aa_seq[pos+1:]]
    return mut_seq
# Function for generating pattern of a given length
def seq_pattern(aa_seq,win_len):
    aa_seq == aa_seq.upper
    seq_pat=[]
    for i1 in range(0, (len(aa_seq) + 1 - win_len)):
        i2 = i1 + int(win_len)
        seq_pat += [aa_seq[i1:i2]]
    return seq_pat


def aac_gen(seq,option,x,y):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    seq = seq.upper()
    aac=[]
    if option=='Normal':
        seq=seq
    elif option=='N':
        seq=seq[0:x]
    elif option=='C':
        seq=seq[-x:][::-1]
    elif option=='NC':
        seq=seq[0:x]+seq[-y:][::-1]
    for i in std:
        counter = seq.count(i) 
        aac+=[((counter*1.0)/len(seq))*100]
    return aac            
        
def dpc_gen(seq,option,x,y):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    seq=seq.upper()
    dpc=[]
    if option=='Normal':
        seq=seq
    elif option=='N':
        seq=seq[0:x]
    elif option=='C':
        seq=seq[-x:0][::-1]
    elif option=='NC':
        seq=seq[0:x]+seq[-y:][::-1]
    for j in std:
        for k in std:
            temp  = j+k
            count = seq.count(temp)
            dpc+=[((count*1.0)/(len(seq)-1))*100]
    return dpc
            
def bin_aac(seq,option,x=None,y=None):
    amino_acids=list("ACDEFGHIKLMNPQRSTVWY")
    Dict={}
    #print(x)
    for i,j in enumerate(amino_acids):
        Dict[j]=i
    seq=seq.upper()
    lis=np.asarray([])
    if option=='Normal':
        seq=seq
    elif option=='N':
        seq=seq[0:x]
    elif option=='C':
        seq=seq[len(seq)-x:][::-1]
    elif option=='NC':
        seq=seq[0:x]+seq[len(seq)-y:][::-1]
    for i in seq:
        a=np.zeros(len(amino_acids))    
        a[Dict[i]]=1
        lis=np.append(lis,a)
    return lis

def getVector(line,option1,option2,x=None,y=None):
    if option1=='aac':
        return aac_gen(line,option2,x,y)
    elif option1=='dpc':
        return dpc_gen(line,option2,x,y)
    elif option1=='bin':
        return bin_aac(line,option2,x,y)

def getXYforfeature(line,option1,option2,x=None,y=None):
    X=[]
    Y=[]
    X+=[getVector(line,option1,option2,x,y)]
    Y+=[+1]
    return X,Y

def adjusted_classes(y_scores, t):
    return [1 if y >= t else -1 for y in y_scores]

def Perform_testing(clf,name,X,Y,t):
    Y_test=Y
    Y_pred=clf.predict(X)
    Y_scores=[]
    if hasattr(clf,'decision_function'):
        Y_scores=clf.decision_function(X)
    else:
        Y_scores=clf.predict_proba(X)[:,1]
    Y_pred = adjusted_classes(Y_scores,t)
    return Y_pred,Y_scores 

print('##############################################################################')
print('# This program AntiCP2 is developed for predicting, desigining and scanning  #')
print('# anticancer peptides, developed by Prof G. P. S. Raghava group.             #')
print('# ############################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
 
# Output file 
 
if args.output == None:
    result_filename= "outfile.csv" 
else:
    result_filename = args.output
         
# Threshold 
if args.threshold == None:
        Threshold = 0.5
else:
        Threshold= float(args.threshold)
# Model
if args.model == None:
        Model = int(1)
else:
        Model = int(args.model)
# Job Type 
if args.job == None:
        Job = int(1)
else:
        Job = int(args.job)
# Window Length 
if args.winleng == None:
        Win_len = int(10)
else:
        Win_len = int(args.winleng)

# Display
if args.display == None:
        dplay = int(1)
else:
        dplay = int(args.display)


print('Summary of Parameters:')
print('Input File: ',Sequence,'; Model: ',Model,'; Threshold: ', Threshold,'; Job Type: ',Job)
print('Output File: ',result_filename,'; Window Length: ',Win_len,'; Display: ',dplay)


# Check above functions
#seq1 = "AGTGTGTSWARTGPTGRWK"

    
#------------------ Read input file ---------------------
def load_model(path):
    clf = pickle.load(open(path,'rb'))
    return clf

f=open(Sequence,"r")
len1 = f.read().count('>')
f.close()

f=open(Sequence,"r")
seqs=[]
seqid=[]
str1='';
header=0
line=0
if len1 >= 1: # read fasta file
#    print('Fasta File')
    for l in f:
        if l.startswith('>'):
            if header != 0:
                seqs += [str1]
                str1 = '';
            header = 1
            line+=1
            seqid += [l.rstrip()]
#            print('seqid=',seqid,'; seq=',seqs)
        else:
            str1 += l.rstrip()
 #           print('seq=',str1)
    seqs += [str1]
else: # read single line file    
    for l in f:
        if len(l) >= 5 and len(l)<= 50:
            seqs+=[l.strip()]
            seqid += ['>Seq_' + str(line)]
        line+=1

fout= open(result_filename,"w+")

i1 = 0
#======================= Prediction Module start from here =====================
if Job == 1:
    fout.write('# Sequence_ID, Sequence, Prediction, Score\n')
    if Model==1:
        clf=load_model('dpc_extra_model')
        for Sequence in seqs:
            header = seqid[i1]
            if len(Sequence) >= 5: 
                if len(Sequence) >= 51:
                    Sequence = Sequence[0:50]
                X,Y=getXYforfeature(Sequence,'dpc','Normal',0,0)
                Y_pred,Y_score=Perform_testing(clf,'svm',X,Y,Threshold)
                flag=""
                if Y_pred[0]==1:
                    flag='AntiCP'
                else:
                    flag='Non AntiCP'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
            i1 = i1 +1
    else: 
        clf=load_model('aac_extra_model')
        for Sequence in seqs:
            header = seqid[i1]
            if len(Sequence) >= 5: 
                if len(Sequence) >= 51:
                    Sequence = Sequence[0:50]
                X,Y=getXYforfeature(Sequence,'aac','Normal',0,0)
                Y_pred,Y_score=Perform_testing(clf,'svm',X,Y,Threshold)
                flag=""
                if Y_pred[0]==1:
                    flag='AntiCP'
                else:
                    flag='Non AntiCP'
                if dplay == 1:
                    if Y_pred[0]==1:
                        fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
                else:
                    fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
            i1 = i1 +1

#===================== Design Model Start from Here ======================
elif Job == 2:
    fout.write('# Sequence_ID, Sequence, Prediction, Score\n')
    print('==== Designing Peptides: Processing sequences please wait ...')
    if Model==1:
        i1 = 0
        for Sequence1 in seqs:
            pat_seq=[]
            header = seqid[i1]
            print('SeqID : ',header,'# under process ...')
            if len(Sequence1) >= 5: 
                fout.write("#Main Sequence (SequenceId, Sequence): %s,%s\n" % (header,Sequence1))
                fout.write('# Sequence_ID, pattern, Prediction, Score\n')
                if len(Sequence1) >= 51:
                    Sequence1 = Sequence1[0:50]
                pat_seq = seq_mutants(Sequence1)
                for Sequence in pat_seq:
                    clf=load_model('dpc_extra_model')
                    X,Y=getXYforfeature(Sequence,'dpc','Normal',0,0)
                    Y_pred,Y_score=Perform_testing(clf,'svm',X,Y,Threshold)
                    flag=""
                    if Y_pred[0]==1:
                        flag='AntiCP'
                    else:
                        flag='Non AntiCP'
                    if dplay == 1:
                        if Y_pred[0]==1:
                            fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
                    else:
                        fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
            i1 = i1 +1
    else:
        i1 = 0
        for Sequence1 in seqs:
            pat_seq=[]
            header = seqid[i1]
            print('SeqID : ',header,'# under process ...')
            if len(Sequence1) >= 5: 
                fout.write("#Main Sequence (SequenceId, Sequence): %s,%s\n" % (header,Sequence1))
                fout.write('# Sequence_ID, pattern, Prediction, Score\n')
                if len(Sequence1) >= 51:
                    Sequence1 = Sequence1[0:50]
                pat_seq = seq_mutants(Sequence1)
                for Sequence in pat_seq:
                    clf=load_model('aac_extra_model')
                    X,Y=getXYforfeature(Sequence,'aac','Normal',0,0)
                    Y_pred,Y_score=Perform_testing(clf,'svm',X,Y,Threshold)
                    flag=""
                    if Y_pred[0]==1:
                        flag='AntiCP'
                    else:
                        flag='Non AntiCP'
                    if dplay == 1:
                        if Y_pred[0]==1:
                            fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
                    else:
                        fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
            i1 = i1 +1

#=============== Scan Model start from here ==================

else:
    print('==== Scanning Peptides: Processing sequences please wait ...')
    if Model==1:
        i1 = 0
        for Sequence1 in seqs:
            pat_seq=[]
            header = seqid[i1]
            print('SeqID : ',header,'# under process ...')
            if len(Sequence1) >= Win_len: 
                fout.write("#Main Sequence (SequenceId, Sequence): %s,%s\n" % (header,Sequence1))
                fout.write('# Sequence_ID, pattern, Prediction, Score\n')
                pat_seq = seq_pattern(Sequence1,Win_len)
                for Sequence in pat_seq:
                    clf=load_model('dpc_extra_model')
                    X,Y=getXYforfeature(Sequence,'dpc','Normal',0,0)
                    Y_pred,Y_score=Perform_testing(clf,'svm',X,Y,Threshold)
                    flag=""
                    if Y_pred[0]==1:
                        flag='AntiCP'
                    else:
                        flag='Non AntiCP'
                    if dplay == 1:
                        if Y_pred[0]==1:
                            fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
                    else:
                        fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
            i1 = i1 +1
    else:
        i1 = 0
        for Sequence1 in seqs:
            pat_seq=[]
            header = seqid[i1]
            print('SeqID : ',header,'# under process ...')
            if len(Sequence1) >= Win_len: 
                fout.write("#Main Sequence (SequenceId, Sequence): %s,%s\n" % (header,Sequence1))
                fout.write('# Sequence_ID, pattern, Prediction, Score\n')
                pat_seq = seq_pattern(Sequence1,Win_len)
                for Sequence in pat_seq:
                    clf=load_model('aac_extra_model')
                    X,Y=getXYforfeature(Sequence,'aac','Normal',0,0)
                    Y_pred,Y_score=Perform_testing(clf,'svm',X,Y,Threshold)
                    flag=""
                    if Y_pred[0]==1:
                        flag='AntiCP'
                    else:
                        flag='Non AntiCP'
                    if dplay == 1:
                        if Y_pred[0]==1:
                            fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
                    else:
                        fout.write("%s,%s,%s,%f\n" % (header,Sequence,flag,Y_score[0]))
            i1 = i1 +1

fout.close()

print('\n======= Thanks for using AntiCP2. Your results are stored in file :',result_filename,' =====\n\n')
print('Please cite Agrawal et al. (2020) AntiCP 2.0: an updated model for predicting anticancer peptides. Briefings in Bioinformatics. doi: 10.1093/bib/bbaa153\n')



