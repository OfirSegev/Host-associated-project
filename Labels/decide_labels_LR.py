import pandas as pd
import sys
import numpy as np

#add new labels
#add condition if 2?

LABEL=sys.argv[1]

df=pd.read_csv("new_"+LABEL+"_labels_predicted_LR.csv")

def decide_label(df):

	if num_of_labels == 2:
		if(df[LABEL+"_logistic_pred_-1"]>=0.7):
			return -1
		elif(df[LABEL+"_logistic_pred_1"]>=0.7):
			return 1
		else:
			return None

	elif num_of_labels == 4:

		if(df[LABEL+"_logistic_pred_-1"]>=0.6):
			return -1
		elif(df[LABEL+"_logistic_pred_0"]>=0.6):
			return 0
		elif(df[LABEL+"_logistic_pred_1"]>=0.6):
			return 1
		elif(df[LABEL+"_logistic_pred_2"]>=0.6):
			return 2
		else:
			return None

	else:

		if(df[LABEL+"_logistic_pred_-1"]>=0.7):
			return -1
		elif(df[LABEL+"_logistic_pred_0"]>=0.7):
			return 0
		elif(df[LABEL+"_logistic_pred_1"]>=0.7):
			return 1
		else:
			return None



#switch prediction with manually labeled value if exists
def prediction_to_true_value():

	if (LABEL == "HvN"):
		original_col = "Host (1) Env(0) unknown (-1)"
	elif (LABEL == "AvP"):
		original_col = "animal (1) plant(0)"
	elif (LABEL == "PvC"):
		original_col = "pathogen(1) commensal(0)"
	elif (LABEL == "HvM"):
		original_col = "human(1) mouse (0)"
	elif (LABEL == "GvOvS"):
		original_col = "GI tract(0) Oral(1) Skin (2)"
	elif (LABEL == "ArvPo"):
		original_col = "Arabidopsis(1) Poplar(0)"
	elif (LABEL == "RvS"):
		original_col = "roots(1) shoots(0)"
	elif (LABEL == "TvA"):
		original_col = "tresterial(1) aquatic(0)"
	elif (LABEL == "Inv"):
		label_column = "inverterbrate(1) non(-1)"
	elif (LABEL == "Art"):
		label_column = "Arthropod (1) non(-1)"
	elif (LABEL == "Ins"):
		label_column = "Insect(1) non(-1)"
	elif (LABEL == "Hemo"):
		label_column = "Hemolymph(1) Other(0)"

	df[original_col] = df[original_col].astype('Int64')

	df["predicted_label_"+LABEL] = (df["predicted_label_"+LABEL]).where(df[original_col].isna(), df[original_col])

	df["predicted_label_"+LABEL] = df["predicted_label_"+LABEL].astype('Int64')

two_labels = ["Inv","Art","Ins"]
four_labels = ["GvOvS"]

if LABEL in two_labels:
	num_of_labels = 2
elif LABEL in four_labels:
	num_of_labels = 4
else:
	num_of_labels = 3

#decide labels
df["predicted_label_"+LABEL] = df.apply(decide_label, axis=1)

#add the column
prediction_to_true_value()

print (df)

df.to_csv("new_"+LABEL+"_labels_with_prediction_LR.csv", index=False)

print(df["predicted_label_"+LABEL].value_counts())
