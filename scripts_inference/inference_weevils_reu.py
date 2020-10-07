import pandas as pd
import numpy as np
import MESS as MESS


SIMOUT='/home/jluiselli/Documents/cours/internship_m1_MESS/Inference_simouts/SIMOUT_total.csv'
data = pd.read_csv('/home/jluiselli/Documents/cours/internship_m1_MESS/MESS/empirical_data/Mascarene_weevils/reunion.dat')
cla = MESS.inference.Classifier(data, SIMOUT, algorithm="rf", verbose=True)
res,proba = cla.predict(quick=True, verbose=True)
print(res)
print(proba)

res.to_csv('/home/jluiselli/Documents/cours/internship_m1_MESS/res_weevils_reu.txt')
proba.to_csv('/home/jluiselli/Documents/cours/internship_m1_MESS/proba_weevils_reu.txt')