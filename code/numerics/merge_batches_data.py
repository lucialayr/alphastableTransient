import pandas as pd
from ast import literal_eval
import sys

alpha = literal_eval(sys.argv[1])
k = literal_eval(sys.argv[2])

batches = [pd.read_csv(f"data/final_states_a{alpha}_k{round(k, 2)}_batch{0}.csv"),
           pd.read_csv(f"data/final_states_a{alpha}_k{round(k, 2)}_batch{1}.csv"),
           pd.read_csv(f"data/final_states_a{alpha}_k{round(k, 2)}_batch{2}.csv"),
           pd.read_csv(f"data/final_states_a{alpha}_k{round(k, 2)}_batch{3}.csv")]

results = pd.concat(batches)

results.to_csv(f"data/final_states_a{alpha}_k{round(k, 2)}.csv", index=False)