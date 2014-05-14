
# coding: utf-8

# #Cancer Survival Prediction
# We will use data that is part of an open challenge to predict survival in cancer patients from within the Cancer Genome Atlas (TCGA).  For details please see Yin Yin et al in Nature Biotechnology.  Also please consider improving on the model to derive a better predictor of survival.

# In[32]:

import os
import pandas as pd
import numpy as np
import synapseclient
from lifelines.estimation import AalenAdditiveFitter
import patsy
syn = synapseclient.login()


# ###Download data from Synapse and Load into Data Frames
# Data in Synapse is referenced by accession identifiers.  Here we define some of the identifiers that refer to data we will be using bellow.  It is possible to extract these identifiers by using the ```syn.query``` command that takes a sql like query (see ?)

# In[15]:

ACRONYM = 'KIRC'                #Kidney Renal Clear Carcinoma
trainLabelsId = 'syn1714093'   # Training bootstraps for KIRC
testLabelsId = 'syn1714090'    # Testing boostraps for KIRC
#dataId = 'syn1710306'          # for RPPA KIRC data
survivalDataId = 'syn1710303'
clinicalDataId = 'syn1715824'


# We will start by fetching the labels of the training set and testing set.  Fetching from synapse is done with ```syn.get``` given a synapse identifier such as 'syn123'.  The returned object contains information about the file in addition to the locatione where it is cached in the path variable that we use here to read the data into a pandas dataframe.

# In[16]:

#Download bootstrap labels
testLabels = pd.read_csv(syn.get(testLabelsId).path, sep='\t', header=None)
trainLabels = pd.read_csv(syn.get(trainLabelsId).path, sep='\t', header=None)


# We will also fetch the clinical covariates and survival data and read those into dataframes as well.

# In[17]:

survival=pd.read_csv(syn.get(survivalDataId).path, '\t', index_col=0)
clinical=pd.read_csv(syn.get(clinicalDataId).path, '\t', index_col=0, na_values=['[Not Available]'])
clinical = pd.concat([survival.ix[:, :-1], clinical.ix[:,:-1]], axis=1)
clinical.head()


# ##Use Lifelines package to model surivival curves
# We will use the Using a Aallen's additive model to estimate the survival times in the testing set using a model built from a training set using 100 fold sub-sampling of the data. 

# In[33]:

#Go through each training testing monteCarlo sampling and train/predict
predictions=[]
for i in range(trainLabels.shape[1]):
    X = patsy.dmatrix('age + grade + stage -1', clinical, return_type='dataframe')
    X['T'] = clinical['OS_OS']
    X['C'] = clinical['OS_vital_status']
    
    trainX = X.ix[trainLabels[i],:].reset_index()
    testX = X.ix[testLabels[i],:].reset_index()

    #Build model and train
    aaf = AalenAdditiveFitter(penalizer=1., fit_intercept=True)
    aaf.fit(trainX.drop(['index'], axis=1), duration_col='T', event_col='C',show_progress=False)
    #Predict on testing data
    median = aaf.predict_median(testX.drop(['T','C', 'index'], axis=1))
    median.index = testX['index']
    predictions.append(median.replace([np.inf, -np.inf, np.nan], 0))


# ###Saving Results to Synapse and ask Synapse to evaluate our predictions
# To document what we have done we will start by storing this code in Synapse as a file Entity.

# In[34]:

codeEntity = synapseclient.File('tcga_survival_analysis.py', parentId='syn1720423')
codeEntity = syn.store(codeEntity)


# We then save the predictions we made to a file and create a file Entity for it.

# In[95]:

#Save predictions to file
predictions = np.asarray(predictions).T
np.savetxt('predictions.csv', predictions, fmt='%.4g', delimiter='\t')
results = synapseclient.File('predictions.csv', name='GBM Aallens additive clinical model', parentId='syn1720419')


# We can the annotate the file Entity with more information such as the model we used, the input data we used (or any key-values we chose).  For this specific project it is required that dataType and cancer are specified otherwise Synapse will not evaluate the predictions.

# In[96]:

results.cancer = ACRONYM
results.dataType = 'clinical'
results.method = 'Aallen additive model'
results.normalization = 'None'
results.featureSelection='None'
results.clinicalUsed = 'Yes'


# Finnaly we push the file Entity to Synapse.  Here we use the used and executed parameters to specify the provenance of the entity as well.  This will be displayed as a provenance graph on the entity page.

# In[97]:

results = syn.store(results,  used=[trainLabelsId, testLabelsId, dataId, survivalDataId], executed=[codeEntity])


# Finally we tell Synapse that we would like it to "score" the predictions by submitting the entity to an evaluation.

# In[98]:

submission=syn.submit(1876290, results)


# In[ ]:



