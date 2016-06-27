import pandas as pd
import numpy as np

test_methods=open('test_methods.txt','w')
pep_methods=open('pep_methods.txt','w')

def define_single_dataframe(file1):
	print 'Reading in file: %s\n' %str(file1)
	usefulcols=['Sequence','Protein Group Accessions','Search Engine Rank','Precursor Area','PEP','XCorr','Intensity','Delta Mass [PPM]','First Scan']
	df=pd.read_csv(open(file1,'r'), sep='\t', usecols=usefulcols)

	print 'Processing file.\n'
	#sort dataframe, drop scan duplicates
	df=df.sort_values(['First Scan','PEP','XCorr','Search Engine Rank'],ascending=[1,1,0,1]).drop_duplicates(['First Scan'])
	df=df.rename(columns={'Protein Group Accessions':'PGA','Search Engine Rank':"Rank",'Precursor Area':'Area','Delta Mass [PPM]':'dM','First Scan':'Scan'})

	#fill in precursor values
	df['Area_mod']=np.where(df['Area']>0,df['Area'],np.where(df['Intensity']>0,df['Intensity'],100))

	#make new sequence peptide, replace modifications and replace I->L
	df['Sequence_mod']=df['Sequence'].str.upper().str.replace('I','L')

	#start building final df
	df_final=df[['Sequence_mod','PGA']]
	df_final=df_final.drop_duplicates(['Sequence_mod'])

	#calculate PSM
	df['PSM']=1
	df_final['PSM']=df.groupby(['Sequence_mod'])['PSM'].transform('count')

	#calculate average dMass
	df_final['dM']=df.groupby(['Sequence_mod'])['dM'].transform('sum')
	df_final['Avg_dM']=df_final.dM/df_final.PSM

	#calculate area for each peptide
	df=df.drop_duplicates(['Sequence_mod','Area_mod'])
	df_final['Area']=df.groupby(['Sequence_mod'])['Area_mod'].transform('sum')

	#apply PSM, AvgM, and PGA thresholds
	threshold=(df_final['PSM'] > 1) & (df_final['Avg_dM'] >-1.5) & (df_final['Avg_dM'] < 1.5) & (df_final['PGA'].str.contains('IGH'))

	#export dataframe to csv
	return df_final.ix[threshold,['Sequence_mod','PGA','PSM','Avg_dM','Area']]

def merge_injections(list_of_files):
	num=len(list_of_files)
	
	print 'Merging multiple injections.\n'
	#Merge into master list
	merge=pd.concat(list_of_files)
	
	#calculate average psm, area for each sequence
	merge['Avg_PSM']=merge.groupby(['Sequence_mod'])['PSM'].transform('sum')/num
	merge['Avg_Area']=merge.groupby(['Sequence_mod'])['Area'].transform('sum')/num
	
	#Merge files, drop sequence duplicates
	merge_final=merge[['Sequence_mod','PGA','Avg_PSM','Avg_Area']]
	merge_final=merge.drop_duplicates(['Sequence_mod'])
	
	return merge_final[['Sequence_mod','Avg_PSM','Avg_Area']]
	
def merge_FT_ELU(Elu,FT):
	
	print 'Merging data from flow through and elutions.\n'
	#Merge into master list and rename columns
	Elu.set_index(['Sequence_mod'])
	FT.set_index(['Sequence_mod'])
	ratios=pd.merge(Elu,FT, on='Sequence_mod', how='outer')
	ratios=ratios.rename(columns={'Avg_PSM_x':'Elu_PSM', 'Avg_Area_x':'Elu_Area','Avg_PSM_y':'FT_PSM','Avg_Area_y':'FT_Area'})
	ratios=ratios.fillna(0)
	
	#calculate ratios
	ratios['R_PSM']=np.where(ratios.FT_PSM>0,ratios.Elu_PSM/ratios.FT_PSM,100000)
	ratios['R_Area']=np.where(ratios.FT_Area>0,ratios.Elu_Area/ratios.FT_Area,100000)
	
	return ratios[['Sequence_mod','R_PSM','R_Area']]
	
def map_peptide_position(db_filename, merged):
	usefulcols=['ABSEQ.AA','VREGION.CDR1.AA','VREGION.CDR2.AA','CDR3.AA']
	df=pd.read_csv(open(db_filename,'r'), sep='\t', usecols=usefulcols)
	df=df.rename(columns={'ABSEQ.AA':'Seq','VREGION.CDR1.AA':'CDR1','VREGION.CDR2.AA':'CDR2','CDR3.AA':'CDR3'})
	#df=df.apply(lambda x: x.str.replace('I','L'))
	df['Seq']=df['Seq'].str.replace('I','L')
	df['CDR1']=df['CDR1'].str.replace('I','L')
	df['CDR2']=df['CDR2'].str.replace('I','L')
	df['CDR3']=df['CDR3'].str.replace('I','L')
	df=df.fillna(' ')
	
	#finds and returns starting coordinates for each CDR, returns -1 if not found
	def f(x): #x is a series, applies function on a row-by-row basis
		return x[0].find(x[1]),x[0].find(x[2]),x[0].find(x[3])
		
	#creates a virtual map of the sequence 
	#0:Seq, 1:CDR1, 2:CDR2, 3:CDR3, 4:CDR1,2,3 pos
	
	def f2(x):
		seq=[0]*len(x[0])
		for i in range(x[4][0],x[4][0]+len(x[1])):
			seq[i]=1
		for i in range(x[4][1],x[4][1]+len(x[2])):
			seq[i]=2
		for i in range(x[4][2],x[4][2]+len(x[3])):
			seq[i]=3
		return seq
		
	#map out the start locations for CDR1, CDR2, and CDR3 for each peptide within the full length sequence
	print 'Mapping peptides to CDR regions.\n'
	df['CDRmap_Coords']=df.apply(f,axis=1)
	df['CDRmap']=df.apply(f2,axis=1)
	peptide_mapping=df.set_index('Seq')['CDRmap'].to_dict()
	
	#x:Peptide sequence
	#peptide mapping (full seq: coords)
	def f3(x):
		for seq in peptide_mapping:
			if x in seq:
				pos=seq.find(x)
				end=pos+len(x)
				CDR=max(peptide_mapping[seq][pos:end])
				return CDR
			else:
				pass
		return -1 #sequence not found
			
	#merged[['Sequence_mod','R_PSM','R_Area']]
	print 'Updating CDR database.\n'
	merged['CDR']=merged['Sequence_mod'].apply(f3)
	
	df.to_csv(test_methods, sep='\t', index=None)
	merged.to_csv(pep_methods, sep='\t', index=None)
