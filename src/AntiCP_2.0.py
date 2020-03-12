import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import classification_report
from sklearn.metrics import roc_auc_score
from sklearn.metrics import cohen_kappa_score
import copy
from sklearn.metrics import confusion_matrix
from features_general import aac_gen,dpc_gen,bin_aac
from sklearn.linear_model import RidgeClassifier
from sklearn.ensemble import ExtraTreesClassifier 
import pickle
from sklearn.model_selection import cross_validate
import warnings
from sklearn.metrics import fbeta_score, make_scorer
import os
import sys
import joblib
from prettytable import PrettyTable
warnings.filterwarnings('ignore')

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

def load_model(path):
	#clf=joblib.load(path)
	clf = pickle.load(open(path,'rb'))
	return clf

#Model1
#DPC N5
#AAC N5 
tabular = PrettyTable()
tabular.field_names = ["Sequence", "Prediction", "Score"]


Sequence=input('Enter the Input Fasta File Name: ') 
Model=int(input('choose model 1 / 2 Experimental/Random :')) 
Threshold=float(input('Enter the Threshold: '))
result_filename=input('Enter the Output File Name: ')

seqs=[]
line=0
f=open(Sequence,"r")
for l in f:
	if line%2==1:
		seqs+=[l.strip()]
	line+=1

fout= open(result_filename,"w+")
fout.write('Sequence,Prediction,Score\n')
if Model==1:
	root1='./ACPs and non-ACPS' 
	clf=load_model('dpc_extra_model') 
	for Sequence in seqs:	
		X,Y=getXYforfeature(Sequence,'dpc','Normal',0,0) 
		Y_pred,Y_score=Perform_testing(clf,'svm',X,Y,Threshold)
		flag=""
		if Y_pred[0]==1:
			flag='AntiCP'
		else:
			flag='Non AntiCP' 
		tabular.add_row([Sequence,flag,Y_score[0]])
		fout.write(Sequence+","+flag+","+str(Y_score[0])+"\n")
else: 
    root1='./ACPs and random peptides' 
    clf=load_model('aac_extra_model') 
    for Sequence in seqs:
	    X,Y=getXYforfeature(Sequence,'aac','Normal',0,0) 
	    Y_pred,Y_score=Perform_testing(clf,'svm',X,Y,Threshold) 
	    flag=""
	    if Y_pred[0]==1:
	    	flag='AntiCP'
	    else:
	    	flag='Non AntiCP'
	    tabular.add_row([Sequence,flag,Y_score[0]])
	    fout.write(Sequence+","+flag+","+str(Y_score[0])+"\n")
print(tabular)
fout.close()
