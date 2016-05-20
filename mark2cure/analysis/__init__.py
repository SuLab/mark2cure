import pandas as pd


synonyms_dict = pd.read_csv('mark2cure/analysis/data/synonym_dictionary.txt', sep='\t', names=['dirty', 'clean'], index_col='dirty').to_dict()['clean']
