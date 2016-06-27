import pandas as pd
import os
import numpy as np
import Pipeline_methods as pipe

db_filename='db_pepposition.txt'
CDR3_clonotyping='CDR3_clonotyping.txt'


#df['Sequence_mod','PGA','PSM','Avg_dM','Area']
Elu1=pipe.define_single_dataframe('MPER_cont3_A_1_psms.txt')
Elu2=pipe.define_single_dataframe('MPER_cont3_B_1_psms.txt')
Elu3=pipe.define_single_dataframe('MPER_cont3_C_1_psms.txt')
#Elu4=pipe.define_single_dataframe('MPER_cont3_D_1_psms.txt')
#Elu5=pipe.define_single_dataframe('MPER_cont3_E_1_psms.txt')

FT1=pipe.define_single_dataframe('MPER_cont3_A_2_psms.txt')
FT2=pipe.define_single_dataframe('MPER_cont3_B_2_psms.txt')
FT3=pipe.define_single_dataframe('MPER_cont3_C_2_psms.txt')
#FT4=pipe.define_single_dataframe('MPER_cont3_D_2_psms.txt')
#FT5=pipe.define_single_dataframe('MPER_cont3_E_2_psms.txt')

#Combine all elutions or all FT, calculate averages
Elu_combined=pipe.merge_injections([Elu1,Elu2,Elu3])
FT_combined=pipe.merge_injections([FT1,FT2,FT3])

#Combine Elution averages with FT averages, calculate ratios
merged=pipe.merge_FT_ELU(Elu_combined,FT_combined)

#Map peptides to CDR1, CDR2 or CDR3 regions
pipe.map_peptide_position(db_filename,merged)
