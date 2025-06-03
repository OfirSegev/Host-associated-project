import os
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
import pandas as pd
#pd.set_option('display.max_columns', 20)
#pd.set_option('display.max_rows', 200)
pd.set_option('display.float_format', lambda x: '%.4f' % x)
#pd.set_option('display.max_colwidth', 1000)

import pickle
import nltk
import sys
#nltk.download('stopwords')
from nltk.stem import WordNetLemmatizer
from TrainLabelsClass import TrainLabels
#from SVMclass import TrainLabelsSVM

########################################################

#get label column to classify as input
LABEL = sys.argv[1]

#if sys.argv[2]=="repeat":
#	df = pd.read_csv("new_HvN_labels_with_prediction_LR.csv")
#	df = df[df["predicted_label"].isna()]
#	df = df.drop(columns=[LABEL+"_logistic_pred_-1",LABEL+"_logistic_pred_0",LABEL+"_logistic_pred_1","predicted_label"])
#else:
#	df = pd.read_csv("edited_mannual_labeling_sheet.csv")

#if not HvN, this means we already have a different csv to start with
columns = ['metadata processed', 'genome_id', 'meta_not_stemmed', 'source',
       'Host (1) Env(0) unknown (-1)', 'pathogen(1) commensal(0)',
       'animal (1) plant(0)', 'human(1) mouse (0)',
       'GI tract(0) Oral(1) Skin (2)', 'Arabidopsis(1) Poplar(0)',
       'roots(1) shoots(0)', 'tresterial(1) aquatic(0)',
       'inverterbrate(1) non(-1)', 'Arthropod (1) non(-1)', 'Insect(1) non(-1)' ,'Hemolymph(1) Other(0)']

if LABEL == "HvN":
	df = pd.read_csv("edited_mannual_labeling_sheet.csv", usecols=columns)
else:
	df = pd.read_csv("merged_prediction_results.csv")

print(df)


if (LABEL == "HvN"):
	label_column = "Host (1) Env(0) unknown (-1)"
elif (LABEL == "AvP"):
	label_column = "animal (1) plant(0)"
	df = df[df["predicted_label_HvN"].isin([1])]
elif (LABEL == "PvC"):
	label_column = "pathogen(1) commensal(0)"
	df = df[df["predicted_label_HvN"].isin([1])]
elif (LABEL == "HvM"):
	label_column = "human(1) mouse (0)"
	df = df[df["predicted_label_AvP"].isin([1])]
elif (LABEL == "GvOvS"):
	label_column = "GI tract(0) Oral(1) Skin (2)"
	df = df[df["predicted_label_AvP"].isin([1])]
elif (LABEL == "ArvPo"):
	label_column = "Arabidopsis(1) Poplar(0)"
	df = df[df["predicted_label_AvP"].isin([0])]
elif (LABEL == "RvS"):
	label_column = "roots(1) shoots(0)"
	df = df[df["predicted_label_AvP"].isin([0])]
elif (LABEL == "TvA"):
	label_column = "tresterial(1) aquatic(0)"
	df = df[~df["predicted_label_HvN"].isin([-1])]
elif (LABEL == "Inv"):
	label_column = "inverterbrate(1) non(-1)"
	df = df[df["predicted_label_AvP"].isin([1])]
elif (LABEL == "Art"):
	label_column = "Arthropod (1) non(-1)"
	df = df[df["predicted_label_Inv"].isin([1])]
elif (LABEL == "Ins"):
	label_column = "Insect(1) non(-1)"
	df = df[df["predicted_label_Art"].isin([1])]
elif (LABEL == "Hemo"):
	label_column = "Hemolymph(1) Other(0)"
	df = df[df["predicted_label_Inv"].isin([1])]
else:
	label_column = None
	print ("invalid label")

#change maybe to original
meta_column = "metadata processed"
#meta_column = "meta_not_stemmed"

########################################################

# Run logistic model

## balanced model when True, unbalanced when False
if LABEL == "HvN":
	to_balance = True
else:
	to_balance = False

model = TrainLabels(df, 
	label_column = label_column, 
	f = 1, 
	balance = to_balance, 
	meta_column = meta_column)

# train model
model.train()
# predct the labels
predicted = model.predictTest()
# calculate the performance scores
model.metrics()
# calculate the performance of the model on the rest of the data
# for that first
y_test, y_pred = model.predictRest()
# can use the methods of the class to calculate the performace and plot if expicitely giving the y_test and y_pred
model.metrics(y_test, y_pred)
# plot PRC and AUC
#model.plot(y_test, y_pred)
# predict unlabaled data from the input df
dfNew, predictions = model.predictNew()

print (dfNew)
print (predictions)

#dfNew.to_csv("new_labels_predicted.csv")


# the naming might diverge from the output I sent you (example from the output label_pred = 'logistic_pred_AvP')
label_pred = LABEL+'_logistic_pred'
# add the predictions to the original table
if LABEL in ["Inv","Art","Ins"]:
	df[label_pred+"_-1"] = model.predictAll()[:,0]
	df[label_pred+"_1"] = model.predictAll()[:,1]
else:
	df[label_pred+"_-1"] = model.predictAll()[:,0]
	df[label_pred+"_0"] = model.predictAll()[:,1]
	df[label_pred+"_1"] = model.predictAll()[:,2]
	if LABEL == "GvOvS":
		df[label_pred+"_2"] = model.predictAll()[:,3]

## combie with the manual labels if wished
#combined = LABEL+'_combined_manual_predicted'
#df[combined] = df[label_column]
#df.loc[df[label_column].isna(),combined] = df.loc[df[label_column].isna(),label_pred]

#df[combined].value_counts()

# export
#if sys.argv[2]=="repeat":
#	df.to_csv("repeated_new_"+LABEL+"_labels_predicted_LR.csv", index=False)
#else:
#	df.to_csv("new_"+LABEL+"_labels_predicted_LR.csv", index=False)

df = df.astype(str)
df.to_csv("new_"+LABEL+"_labels_predicted_LR.csv", index=False)

#df[combined].value_counts(bins=10)

#unclassified = df[(df[combined] > 0.3) & (df[combined] < 0.7)]
#unclassified.to_csv("unclassified_"+LABEL+".csv")

