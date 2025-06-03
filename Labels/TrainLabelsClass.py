import pandas as pd
from IPython.display import display
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, accuracy_score, f1_score
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn import metrics
from sklearn.metrics import precision_recall_fscore_support
from sklearn.model_selection import train_test_split
from imblearn.under_sampling import RandomUnderSampler

from matplotlib import pyplot as plt

class TrainLabels:
#'''
#    The class serves for the training of different logistic models
#    The parameters to vary are the label column, the fraction of data used (f) and
#    whether the data should be balanced. The meta_column is per default 2, need to change
#    to the actual column name or number.
#    To predict the label of new/unlabelled data, need to create a corpus with the words contained in it,
#    otherwise the logistic model will give an error.
#    So the class create a corpus with all data including unlabelled and work with indices for labaled and unlaballed data
#    This is expensive, so you can pass a corpus by initialization (corpus = your_corpus). 
#    Precition of new genomes (meta) might not work in this case since some words might be missing.
#'''    
    
    def __init__(self, df, meta_column = 2, label_column = -1, f = 0.1, 
                 balance = True, classifier = None, corpus = None):
        
        if not isinstance(meta_column, int):
            meta_column = df.columns.get_loc(meta_column)
        if not isinstance(label_column, int):
            label_column = df.columns.get_loc(label_column)

        self.meta_column = meta_column
        self.label_column = label_column
        self.df = df.reset_index(drop = True)
        
        # index of all instances
        index_all = self.df.index

        # can add precomputed corpus if e.g. runing same the model many times
        # the corpus has to contain all the words from training 
        # be carefull with indexing 
        if corpus is None:
            self.corpus = self.df.iloc[:,self.meta_column]
            self.corpus = self.create_tf(self.corpus)
        else:
            self.corpus = corpus
        
        # labels 
        y = self.df.iloc[:,self.label_column]
        # index for labeled and unlabaled data
        self.index_all_labeled = y[y.astype(float).isin([-1,0,1,2])].index
        self.index_unlabeled = y[~y.astype(float).isin([-1,0,1,2])].index
        
        # select labaled X and y for training
        self.y_all_labaled = y[self.index_all_labeled].reset_index(drop = True) # need to reset the index for sampling later 
        self.X_all_labeled = self.corpus[self.index_all_labeled]
        # first balance
        if balance == True:
            sampler = RandomUnderSampler()
            X, y = sampler.fit_resample(self.X_all_labeled, self.y_all_labaled)
        else:
            X, y = self.X_all_labeled, self.y_all_labaled
        # then sample
        # use the index of y since it is a Series and has an index (X is an np.array)
        y = y.sample(frac = f)
        # if balance==True and f<1 where is not used labaled data, that we can use to test the model later (indRest)
        self.indUsed = y.index
        self.indRest = self.index_all_labeled[~self.index_all_labeled.isin(self.indUsed)]
       
        print('Number of observations: ', len(self.indUsed))
        
        self.X = X[self.indUsed,:]
        self.y = y.astype(int).values
        
    def create_tf(self, data ):  
        # create the matrix for training
        self.tf = TfidfVectorizer(ngram_range=(1,2))
        return  self.tf.fit_transform(data).toarray()
    
    
    def train(self):
        
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y, test_size = 0.20)
        self.classifier = LogisticRegression(C = 5)
        self.classifier.fit(self.X_train, self.y_train)
        
        self.y_pred = self.predictTest()
        
    def predictNew(self):
        # predict unlabaled data from the input df
        dataNew = self.corpus[~self.index_all_labeled]
        #returns both the predicted label and the corresponding data
        return self.df.iloc[~self.index_all_labeled,:], self.classifier.predict(dataNew)
    
    def predictRest(self):
        # preduct the unused labeled data
        self.X_rest = self.corpus[self.indRest]
        self.y_rest = self.df.iloc[self.indRest,self.label_column].astype(int).values
        #returns both the predicted label and the true label
        
        return self.y_rest, self.classifier.predict(self.X_rest)
    
    def predictAll(self):
        return self.classifier.predict_proba(self.corpus)

    
    def predictTest(self):
        self.y_pred = self.classifier.predict(self.X_test)
        return self.y_pred 
    
    def predictTestProba(self):
        return self.y_test, self.classifier.predict_proba(self.X_test)
    
    def metrics(self,y_test = None, y_pred = None):
        # if no values are passed, will take self variables
        if (y_test is None) and (y_pred is None):
             
            y_test = self.y_test
            y_pred = self.y_pred 
            
        print(len(y_test))
        cm = confusion_matrix(y_test, y_pred)
        display("Confusion matrix: ", cm)
        precision, recall, fscore, _ = precision_recall_fscore_support(y_test,y_pred, average = 'weighted')

        print("accuracy: ", accuracy_score(self.y_test, self.y_pred).round(3))
        print("f1 score: ", f1_score(self.y_test, self.y_pred, average='weighted').round(3))
        print("precision: ", precision.round(3))
        print("recall: ", recall.round(3))
        return precision, recall, fscore

        
    def plot(self,y_test=None, y_pred=None):
        # if no values are passed, will take self variables
        if (y_test is None) and (y_pred is None):
            y_test = self.y_test
            y_pred = self.y_pred 
            
    
        metrics.RocCurveDisplay.from_predictions(y_test,y_pred)
        metrics.PrecisionRecallDisplay.from_predictions(y_test,y_pred)
    
          
    def feature_importance(self):
        # variables with the highest Beta (impact)?
        weights = pd.Series(self.classifier.coef_[0])
        features = pd.Series(self.tf.get_feature_names())
        df = pd.DataFrame([weights,features]).T
        return df.sort_values(0, ascending = False, key=abs)

