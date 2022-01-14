import pandas as pd;
import sys

args = sys.argv
read = args[1]
Chr = args[2]
sumstat = args[3]
Newsumstat = args[4]
if read == "read_csv":
	df = pd.read_csv(sumstat, sep='\t')
else:
	df = pd.read_parquet(sumstat)

df=df[df['CHR']==int(Chr)];

df.to_csv(Newsumstat, index=False,  sep='\t',  compression='gzip');
print(df['N'].unique()[0])
#print(df[df['CHR']==int(Chr)]['N'].unique())
