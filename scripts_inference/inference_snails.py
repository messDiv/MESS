import pandas as pd
import numpy as np
import MESS as MESS


SIMOUT='/home/jluiselli/Documents/cours/internship_m1_MESS/MESS/SIMOUT_total.csv'
data = pd.read_csv('/home/jluiselli/Documents/cours/internship_m1_MESS/MESS/empirical_data/Galapagos_snails/snail_df.csv')
cla = MESS.inference.Classifier(data, SIMOUT, algorithm="rf", verbose=True)
res,proba = cla.predict(quick=True, verbose=True)
print(res)
print(proba)

res.to_csv('/home/jluiselli/Documents/cours/internship_m1_MESS/MESS/res_snails.txt')
proba.to_csv('/home/jluiselli/Documents/cours/internship_m1_MESS/MESS/proba_snails.txt')