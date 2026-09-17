import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FixedLocator, FixedFormatter
import matplotlib.gridspec as gridspec

# Read specific sheets by name
countsdf = pd.read_excel("data/New_data/20190315_MV61_Spleen_T1, T2, MZ, MZp BPS.xlsx", sheet_name="20190315_MV61_Spleen_Cellcounts")
# remove last row from the dataframe
countsdf = countsdf.iloc[:-1]
print(countsdf.head())
df_ont = pd.read_csv('data/New_data/ontogeny_counts.csv')
host_counts_bu = pd.read_excel("data//New_data/20190712_MV_BuChi_Bcell number_MV61 MZ, T2, MZP-T2 bps.xlsx", sheet_name= 'host.cell.number')
donor_counts_bu = pd.read_excel("data//New_data/20190712_MV_BuChi_Bcell number_MV61 MZ, T2, MZP-T2 bps.xlsx", sheet_name= 'donor.cell.number')
# Calculate the total number of FM cells in the spleen
'''countsdf['FM_count'] = countsdf[['Lymphocytes/Single Cells/live/total B-cells/B220+/mature/FM/CD45.2+ | Count', 
                                'Lymphocytes/Single Cells/live/total B-cells/B220+/mature/FM/CD45.1+ | Count']].sum(axis=1)
countsdf['total_FM'] = countsdf['Cell Counts'] * countsdf['FM_count'] / countsdf['FSC-A+ | Count']'''
#Rename Columns name
countsdf.rename(columns={'Lymphocytes/Single Cells/live/total B-cells/B220+/mature/FM/CD45.2+ | Count': 'donor_FM',
                         'Lymphocytes/Single Cells/live/total B-cells/B220+/mature/FM/CD45.2+/Ki67+ | Count': 'Ki67+FM_donor',
                         'Lymphocytes/Single Cells/live/total B-cells/B220+/mature/FM/CD45.1+ | Count': 'host_FM',
                         'Lymphocytes/Single Cells/live/total B-cells/B220+/mature/FM/CD45.1+/Ki67+ | Count': 'Ki67+FM_host'}, inplace=True)
countsdf['FM_donor'] = countsdf['Cell Counts'] * countsdf['donor_FM'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_FM_donor'] = countsdf['Cell Counts'] * countsdf['Ki67+FM_donor'] / countsdf['FSC-A+ | Count']
countsdf['FM_host'] = countsdf['Cell Counts'] * countsdf['host_FM'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_FM_host'] = countsdf['Cell Counts'] * countsdf['Ki67+FM_host'] / countsdf['FSC-A+ | Count']
countsdf['total_FM'] = countsdf['FM_donor'] + countsdf['FM_host']
countsdf['fd_FM'] = countsdf['FM_donor'] / countsdf['total_FM']
countsdf['percent_total_ki67+_FM'] = (countsdf['Ki67+_FM_donor'] + countsdf['Ki67+_FM_host'])/countsdf['total_FM'] * 100
countsdf['frac_Ki67+_FM_donor'] = countsdf['Ki67+_FM_donor'] / countsdf['FM_donor']
countsdf['frac_Ki67+_FM_host'] = countsdf['Ki67+_FM_host'] / countsdf['FM_host']
countsdf = countsdf.drop(columns=['FM_donor', 'donor_FM', 'Ki67+FM_donor', 'Ki67+_FM_donor', 'host_FM', 'FM_host','Ki67+FM_host', 'Ki67+_FM_host'])
# print(countsdf.head())

# Calculate the total number of T2 MZP cells in the spleen
'''countsdf['T2_MZP_count'] = countsdf[['Lymphocytes/Single Cells/live/total B-cells/B220+/CD21hiCD24hi/T2-MZP/CD45.2+ | Count',
                                     'Lymphocytes/Single Cells/live/total B-cells/B220+/CD21hiCD24hi/T2-MZP/CD45.2+/Ki67+ | Count',
                                     'Lymphocytes/Single Cells/live/total B-cells/B220+/CD21hiCD24hi/T2-MZP/CD45.1+ | Count',
                                     'Lymphocytes/Single Cells/live/total B-cells/B220+/CD21hiCD24hi/T2-MZP/CD45.1+/Ki67+ | Count']].sum(axis=1)
countsdf['total_T2_MZP'] = countsdf['Cell Counts'] * countsdf['T2_MZP_count'] / countsdf['FSC-A+ | Count']'''
#Rename Columns name
countsdf.rename(columns={'Lymphocytes/Single Cells/live/total B-cells/B220+/CD21hiCD24hi/T2-MZP/CD45.2+ | Count': 'T2MZP_donor',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/CD21hiCD24hi/T2-MZP/CD45.2+/Ki67+ | Count': 'Ki67+_T2MZP_donor',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/CD21hiCD24hi/T2-MZP/CD45.1+ | Count': 'T2MZP_host',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/CD21hiCD24hi/T2-MZP/CD45.1+/Ki67+ | Count': 'Ki67+_T2MZP_host'}, inplace=True)

countsdf['T2_MZP_donor'] = countsdf['Cell Counts'] * countsdf['T2MZP_donor'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_T2_MZP_donor'] = countsdf['Cell Counts'] * countsdf['Ki67+_T2MZP_donor'] / countsdf['FSC-A+ | Count']
countsdf['T2_MZP_host'] = countsdf['Cell Counts'] * countsdf['T2MZP_host'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_T2_MZP_host'] = countsdf['Cell Counts'] * countsdf['Ki67+_T2MZP_host'] / countsdf['FSC-A+ | Count']
countsdf['total_T2_MZP'] = countsdf['T2_MZP_donor'] + countsdf['T2_MZP_host']
countsdf['fd_T2_MZP'] = countsdf['T2_MZP_donor'] / countsdf['total_T2_MZP']
countsdf['percent_total_ki67+_T2_MZP'] = (countsdf['Ki67+_T2_MZP_donor'] + countsdf['Ki67+_T2_MZP_host'])/countsdf['total_T2_MZP'] * 100
countsdf['frac_Ki67+_T2_MZP_donor'] = countsdf['Ki67+_T2_MZP_donor'] / countsdf['T2_MZP_donor']
countsdf['frac_Ki67+_T2_MZP_host'] = countsdf['Ki67+_T2_MZP_host'] / countsdf['T2_MZP_host']
countsdf = countsdf.drop(columns=['T2_MZP_donor', 'T2MZP_donor','Ki67+_T2_MZP_donor', 'Ki67+_T2MZP_donor', 'T2_MZP_host', 'T2MZP_host','Ki67+_T2MZP_host', 'Ki67+_T2_MZP_host'])

# Calculate the total number of T1 cells in the spleen

# countsdf['total_MZ'] = countsdf['Cell Counts'] * countsdf['MZ_count'] / countsdf['FSC-A+ | Count']
#Rename Columns name
countsdf.rename(columns={'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/MZ/CD45.2+ | Count': 'donor_MZ',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/MZ/CD45.2+/Ki67+ | Count': 'Ki67+MZ_donor',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/MZ/CD45.1+ | Count': 'host_MZ',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/MZ/CD45.1+/Ki67+ | Count': 'Ki67+MZ_host'}, inplace=True)

countsdf['MZ_donor'] = countsdf['Cell Counts'] * countsdf['donor_MZ'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_MZ_donor'] = countsdf['Cell Counts'] * countsdf['Ki67+MZ_donor'] / countsdf['FSC-A+ | Count']
countsdf['MZ_host'] = countsdf['Cell Counts'] * countsdf['host_MZ'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_MZ_host'] = countsdf['Cell Counts'] * countsdf['Ki67+MZ_host'] / countsdf['FSC-A+ | Count']
countsdf['total_MZ'] = countsdf['MZ_donor'] + countsdf['MZ_host']
countsdf['fd_MZ'] = countsdf['MZ_donor'] / countsdf['total_MZ']
countsdf['frac_Ki67+_MZ'] = (countsdf['Ki67+_MZ_donor'] + countsdf['Ki67+_MZ_host'])/ countsdf['total_MZ']
countsdf['percent_total_ki67+_MZ'] = (countsdf['Ki67+_MZ_donor'] + countsdf['Ki67+_MZ_host'])/countsdf['total_MZ'] * 100
countsdf['frac_Ki67+_MZ_donor'] = countsdf['Ki67+_MZ_donor'] / countsdf['MZ_donor']
countsdf['frac_Ki67+_MZ_host'] = countsdf['Ki67+_MZ_host'] / countsdf['MZ_host']
print('Average fraction of Ki67+ MZB host cells:', countsdf['frac_Ki67+_MZ_host'].mean())
countsdf['ratio_MZ_Ki67+_host_to_donor'] = countsdf['frac_Ki67+_MZ_host'] / countsdf['frac_Ki67+_MZ_donor']
countsdf['ratio_MZ_Ki67+_donor_to_host'] = countsdf['frac_Ki67+_MZ_donor'] / countsdf['frac_Ki67+_MZ_host']
countsdf = countsdf.drop(columns=['donor_MZ', 'Ki67+MZ_donor', 'Ki67+_MZ_donor', 'host_MZ', 'MZ_host','Ki67+MZ_host', 'Ki67+_MZ_host'])
# Calculate the total number of Transitional1 cells in the spleen
'''countsdf['T1_AA4.1_count'] = countsdf[['Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/CD21-/T1/AA4.1+/CD45.2+ | Count',
                                        'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/CD21-/T1/AA4.1+/CD45.2+/Ki67+ | Count',
                                        'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/CD21-/T1/AA4.1+/CD45.1+ | Count',
                                        'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/CD21-/T1/AA4.1+/CD45.1+/Ki67+ | Count']].sum(axis=1)
countsdf['total_T1_AA4.1+'] = countsdf['Cell Counts'] * countsdf['T1_AA4.1_count'] / countsdf['FSC-A+ | Count']'''
#Rename Columns name
countsdf.rename(columns={'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/CD21-/T1/AA4.1+/CD45.2+ | Count': 'T1AA4.1+_donor',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/CD21-/T1/AA4.1+/CD45.2+/Ki67+ | Count': 'Ki67+T1_AA4.1+_donor',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/CD21-/T1/AA4.1+/CD45.1+ | Count': 'T1AA4.1+_host',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/IgMhiCD23-/CD21-/T1/AA4.1+/CD45.1+/Ki67+ | Count': 'Ki67+T1_AA4.1+_host'}, inplace=True)

countsdf['T1_AA4.1+_donor'] = countsdf['Cell Counts'] * countsdf['T1AA4.1+_donor'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_T1_AA4.1+_donor'] = countsdf['Cell Counts'] * countsdf['Ki67+T1_AA4.1+_donor'] / countsdf['FSC-A+ | Count']
countsdf['T1_AA4.1+_host'] = countsdf['Cell Counts'] * countsdf['T1AA4.1+_host'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_T1_AA4.1+_host'] = countsdf['Cell Counts'] * countsdf['Ki67+T1_AA4.1+_host'] / countsdf['FSC-A+ | Count']
countsdf['total_T1_AA4.1+'] = countsdf['T1_AA4.1+_donor'] + countsdf['T1_AA4.1+_host']
countsdf['fd_T1_AA4.1+'] = countsdf['T1_AA4.1+_donor'] / countsdf['total_T1_AA4.1+']
countsdf['frac_Ki67+_T1_AA4.1+'] = (countsdf['Ki67+_T1_AA4.1+_donor'] + countsdf['Ki67+_T1_AA4.1+_host'])/ countsdf['total_T1_AA4.1+']
countsdf['percent_total_ki67+_T1_AA4.1+'] = countsdf['frac_Ki67+_T1_AA4.1+'] * 100
countsdf['frac_Ki67+_T1_AA4.1+_donor'] = countsdf['Ki67+_T1_AA4.1+_donor'] / countsdf['T1_AA4.1+_donor']
countsdf['frac_Ki67+_T1_AA4.1+_host'] = countsdf['Ki67+_T1_AA4.1+_host'] / countsdf['T1_AA4.1+_host']
countsdf['Nfd_MZ'] = countsdf['MZ_donor'] / (countsdf['total_MZ'] * countsdf['fd_T1_AA4.1+'])
countsdf['Nfd_FM'] = countsdf['fd_FM'] /  countsdf['fd_T1_AA4.1+']
countsdf['Nfd_T2_MZP'] = countsdf['fd_T2_MZP'] /  countsdf['fd_T1_AA4.1+']

# countsdf['Nfd_MZ'] = countsdf['fd_MZ']/ countsdf['fd_T1_AA4.1+']

countsdf = countsdf.drop(columns=['T1_AA4.1+_donor', 'T1AA4.1+_donor', 'Ki67+_T1_AA4.1+_donor', 'Ki67+T1_AA4.1+_donor', 'T1AA4.1+_host', 'Ki67+_T1_AA4.1+_host', 'Ki67+T1_AA4.1+_host','MZ_donor'])
# Calculate the total number of Immature T1 cells in the spleen
'''countsdf['T1_count'] = countsdf[['Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T1/CD45.2+ | Count',
                                        'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T1/CD45.2+/Ki67+ | Count',
                                        'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T1/CD45.1+ | Count',
                                        'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T1/CD45.1+/Ki67+ | Count']].sum(axis=1)'''
#Rename Columns name
countsdf.rename(columns={'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T1/CD45.2+ | Count': 'T1im_donor',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T1/CD45.2+/Ki67+ | Count': 'Ki67+T1_donor',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T1/CD45.1+ | Count': 'T1im_host',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T1/CD45.1+/Ki67+ | Count': 'Ki67+T1_host'}, inplace=True)

countsdf['T1_donor'] = countsdf['Cell Counts'] * countsdf['T1im_donor'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_T1_donor'] = countsdf['Cell Counts'] * countsdf['Ki67+T1_donor'] / countsdf['FSC-A+ | Count']
countsdf['T1_host'] = countsdf['Cell Counts'] * countsdf['T1im_host'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_T1_host'] = countsdf['Cell Counts'] * countsdf['Ki67+T1_host'] / countsdf['FSC-A+ | Count']
countsdf['total_T1'] = countsdf['T1_donor'] + countsdf['T1_host']
countsdf['fd_T1'] = countsdf['T1_donor'] / countsdf['total_T1']
countsdf['frac_Ki67+_T1'] = (countsdf['Ki67+_T1_donor'] + countsdf['Ki67+_T1_host'])/ countsdf['total_T1']
countsdf['frac_Ki67+_T1_donor'] = countsdf['Ki67+_T1_donor'] / countsdf['T1_donor']
countsdf['frac_Ki67+_T1_host'] = countsdf['Ki67+_T1_host'] / countsdf['T1_host']
countsdf = countsdf.drop(columns=['T1_donor', 'T1im_donor', 'Ki67+T1_donor', 'T1im_host', 'Ki67+_T1_donor', 'Ki67+T1_host', 'T1_host',  'Ki67+_T1_host'])
# Calculate the total number of Transitional2 cells in the spleen
'''countsdf['T2_count'] = countsdf[['Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T2/CD45.2+ | Count',
                                 'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T2/CD45.2+/Ki67+ | Count',
                                    'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T2/CD45.1+ | Count',
                                    'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T2/CD45.1+/Ki67+ | Count']].sum(axis=1)'''
# countsdf['total_T2'] = countsdf['Cell Counts'] * countsdf['T2_count'] / countsdf['FSC-A+ | Count']
#Rename Columns name
countsdf.rename(columns={'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T2/CD45.2+ | Count': 'T2im_donor',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T2/CD45.2+/Ki67+ | Count': 'Ki67+T2_donor',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T2/CD45.1+ | Count': 'T2im_host',
                            'Lymphocytes/Single Cells/live/total B-cells/B220+/immature/T2/CD45.1+/Ki67+ | Count': 'Ki67+T2_host'}, inplace=True)
countsdf['T2_donor'] = countsdf['Cell Counts'] * countsdf['T2im_donor'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_T2_donor'] = countsdf['Cell Counts'] * countsdf['Ki67+T2_donor'] / countsdf['FSC-A+ | Count']
countsdf['T2_host'] = countsdf['Cell Counts'] * countsdf['T2im_host'] / countsdf['FSC-A+ | Count']
countsdf['Ki67+_T2_host'] = countsdf['Cell Counts'] * countsdf['Ki67+T2_host'] / countsdf['FSC-A+ | Count']
countsdf['total_T2'] = countsdf['T2_donor'] + countsdf['T2_host']
countsdf['fd_T2'] = countsdf['T2_donor'] / countsdf['total_T2']
countsdf['frac_Ki67+_T2'] = (countsdf['Ki67+_T2_donor'] + countsdf['Ki67+_T2_host'])/ countsdf['total_T2']
countsdf['percent_total_ki67+_T2'] = countsdf['frac_Ki67+_T2'] * 100
countsdf['frac_Ki67+_T2_donor'] = countsdf['Ki67+_T2_donor'] / countsdf['T2_donor']
countsdf['frac_Ki67+_T2_host'] = countsdf['Ki67+_T2_host'] / countsdf['T2_host']
countsdf['Nfd_T2'] = countsdf['T2_donor'] / (countsdf['total_T2'] * countsdf['fd_T1_AA4.1+'])
countsdf = countsdf.drop(columns=['T2im_donor', 'Ki67+T2_donor', 'Ki67+_T2_donor', 'Ki67+T2_host', 'T2_donor', 'T2_host', 'T2im_host', 'Ki67+T2_host', 'FSC-A+ | Count', 'Cell Counts'])                                
# print(countsdf.columns)

#Create dataframe with merging T1_AA4.1+ and T2 counts 

df = countsdf[['Lamis ID', 'Age at BMT', 'Age at S1K', 'Days postBMT', 'total_T1_AA4.1+', 'total_T2', 'fd_T1_AA4.1+', 'fd_T2', 'frac_Ki67+_T1_AA4.1+_donor', 'frac_Ki67+_T2_donor', 'frac_Ki67+_T1_AA4.1+_host', 'frac_Ki67+_T2_host']]

# Merge total counts, fd, and fraction of Ki67+ donor and host cells for T1_AA4.1+ and T2 into one column
countsdf1_t1 = df[['Lamis ID', 'Age at BMT', 'Age at S1K', 'Days postBMT', 'total_T1_AA4.1+', 'fd_T1_AA4.1+', 'frac_Ki67+_T1_AA4.1+_donor', 'frac_Ki67+_T1_AA4.1+_host']].rename(
   columns={
      'total_T1_AA4.1+': 'total_count',
      'fd_T1_AA4.1+': 'fd',
      'frac_Ki67+_T1_AA4.1+_donor': 'frac_Ki67+_donor',
      'frac_Ki67+_T1_AA4.1+_host': 'frac_Ki67+_host'
   }
)
countsdf1_t1['subset'] = 'T1'

countsdf1_t2 = df[['Lamis ID', 'Age at BMT', 'Age at S1K', 'Days postBMT', 'total_T2', 'fd_T2', 'frac_Ki67+_T2_donor', 'frac_Ki67+_T2_host']].rename(
   columns={
      'total_T2': 'total_count',
      'fd_T2': 'fd',
      'frac_Ki67+_T2_donor': 'frac_Ki67+_donor',
      'frac_Ki67+_T2_host': 'frac_Ki67+_host'
   }
)
countsdf1_t2['subset'] = 'T2'

countsdf1 = pd.concat([countsdf1_t1, countsdf1_t2], ignore_index=True)

#Save as csv file 
countsdf1.to_csv('data/Transitional_cells.csv', index=False)



 

# Save the DataFrame to a CSV file

# countsdf.to_csv('data/new_data_MZB.csv', index=False)

# add agebin in total_mzcounts based on age at bmt
countsdf['agebin'] = pd.cut(countsdf['Age at BMT'], bins=[0, 70, 140], labels=['<10 weeks', '>10 weeks'])


# Specify the column name to filter out rows with any data in it
'''column_name = 'notes'

# Specify the columns to drop from the filtered DataFrames
columns_to_drop = ['notes', 'Label', 'ABCs', 'CD21intCD24hi', 'Transitional 1', 'Transitional1 AA4.1-']

# Filter out rows with any data in the specified column
filtered_host = host_counts_bu[host_counts_bu[column_name].isna()].drop(columns=columns_to_drop)
filtered_donor = donor_counts_bu[donor_counts_bu[column_name].isna()].drop(columns=columns_to_drop)



# total counts -- Filter the DataFrame to include only rows where the Ki67 column does not contain "+"
total_hostcounts = filtered_host[~filtered_host['Ki67'].str.contains('\+', na=False)].drop(columns=['Ki67']).reset_index(drop=True)
total_donorcounts= filtered_donor[~filtered_donor['Ki67'].str.contains('\+', na=False)].drop(columns=['Ki67']).reset_index(drop=True)
# print(total_hostcounts.head(10))
# print(total_donorcounts.head(10))
#total_donorcounts = total_donorcounts.set_index('Lamis.ID').reindex(total_hostcounts['Lamis.ID']).reset_index()

# Merge the two DataFrames on the specified columns (Host and Donor)
total_mzcounts = total_hostcounts.merge(total_donorcounts, on = ['Lamis.ID', 'days.post.bmt', 'age.at.S1K', 'age.at.bmt'],
                                         suffixes=("_host", "_donor"))  

# Calculate the total MZ and T1 cell counts
total_mzcounts['total_MZ'] = total_mzcounts[['MZ_host', 'MZ_donor']].sum(axis=1)
total_mzcounts['total_Transitional1 AA4.1+'] = total_mzcounts[['Transitional1 AA4.1+_host', 'Transitional1 AA4.1+_donor']].sum(axis=1)
total_mzcounts['total_FM'] = total_mzcounts[['FM_host', 'FM_donor']].sum(axis=1)
total_mzcounts['total_T2.MZP'] = total_mzcounts[['T2.MZP_host', 'T2.MZP_donor']].sum(axis=1)
total_mzcounts['total_T2'] = total_mzcounts[['T2_host', 'T2_donor']].sum(axis=1)

# Calculate the fraction of donor chimerism normalized to Transitional T1
total_mzcounts['fd_T1_AA4.1+'] = total_mzcounts['Transitional1 AA4.1+_donor'] / total_mzcounts['total_Transitional1 AA4.1+']
total_mzcounts['fd_FM'] = total_mzcounts['FM_donor'] / total_mzcounts['total_FM']
total_mzcounts['fd_T2.MZP'] = total_mzcounts['T2.MZP_donor'] / total_mzcounts['total_T2.MZP']
total_mzcounts['fd_T2'] = total_mzcounts['T2_donor'] / total_mzcounts['total_T2']
total_mzcounts['fd_mz'] = total_mzcounts['MZ_donor'] / total_mzcounts['total_MZ']



#Arrange age.at.S1K in ascending order
total_mzcounts = total_mzcounts.sort_values(by='days.post.bmt')
print(total_mzcounts.columns)


#Drop host MZ and T1 columns
# total_mzcounts = total_mzcounts.drop(columns=['MZ_host', 'T2.MZP_host', 'FM_host', 'T2_host', 'Transitional1 AA4.1+_host'])

# Ki67 counts -- Filter the DataFrame to include only rows where the Ki67 column contains '+'
ki67_hostcounts = filtered_host[filtered_host['Ki67'].str.contains('\+', na=False)].drop(columns=['Ki67']).reset_index(drop=True)
ki67_donorcounts = filtered_donor[filtered_donor['Ki67'].str.contains('\+', na=False)].drop(columns=['Ki67']).reset_index(drop=True)
ki67_hostcounts.rename(columns={'MZ': 'Ki67+ Mz', 'Transitional1 AA4.1+': 'Ki67+ T1', 'FM' : 'Ki67+ FM', 'T2' : 'Ki67+ T2', 'T2.MZP' : 'Ki67+ T2.MZP'}, inplace=True)
ki67_donorcounts.rename(columns={'MZ': 'Ki67+ Mz', 'Transitional1 AA4.1+': 'Ki67+ T1', 'FM' : 'Ki67+ FM', 'T2' : 'Ki67+ T2', 'T2.MZP' : 'Ki67+ T2.MZP'}, inplace=True)
# Merge the two DataFrames on the specified columns
Ki67counts = ki67_hostcounts.merge(ki67_donorcounts, on = ['Lamis.ID', 'days.post.bmt', 'age.at.S1K', 'age.at.bmt'],
                                         suffixes=("_host", "_donor")) 
Ki67counts['fraction_ki67T1_AA4.1_host'] = Ki67counts['Ki67+ T1_host'] / total_mzcounts['Transitional1 AA4.1+_host']
Ki67counts['fraction_ki67T1_AA4.1_donor'] = Ki67counts['Ki67+ T1_donor'] / total_mzcounts['Transitional1 AA4.1+_donor']
Ki67counts['fraction_ki67FM_host'] = Ki67counts['Ki67+ FM_host'] / total_mzcounts['FM_host']
Ki67counts['fraction_ki67FM_donor'] = Ki67counts['Ki67+ FM_donor'] / total_mzcounts['FM_donor']
Ki67counts['fraction_ki67T2_host'] = Ki67counts['Ki67+ T2_host'] / total_mzcounts['T2_host']
Ki67counts['fraction_ki67T2_donor'] = Ki67counts['Ki67+ T2_donor'] / total_mzcounts['T2_donor']
Ki67counts['fraction_ki67T2MZP_host'] = Ki67counts['Ki67+ T2.MZP_host'] / total_mzcounts['T2.MZP_host']
Ki67counts['fraction_ki67T2MZP_donor'] = Ki67counts['Ki67+ T2.MZP_donor'] / total_mzcounts['T2.MZP_donor']
Ki67counts['fraction_ki67MZ_host'] = Ki67counts['Ki67+ Mz_host'] / total_mzcounts['MZ_host']
Ki67counts['fraction_ki67MZ_donor'] = Ki67counts['Ki67+ Mz_donor'] / total_mzcounts['MZ_donor']
#Arrange days post bmt  in ascending order
Ki67counts = Ki67counts.sort_values(by='days.post.bmt')
print(Ki67counts.columns)

# Filter out Nfd values that seem unphysiological (>1.5)
# total_mzcounts = total_mzcounts[total_mzcounts['total_MZ'] <= 1e7]
# total_mzcounts = total_mzcounts[total_mzcounts['Nfd'] <= 1.2]'''

# Create a figure with a custom grid layout
fig = plt.figure(figsize=(14, 8))
gs = gridspec.GridSpec(2, 4)  # Define a 2x4 grid

# Top row: Two plots
ax1 = fig.add_subplot(gs[0, 0:2])  # First plot spans columns 0 and 1
ax2 = fig.add_subplot(gs[0, 2:4])  # Second plot spans columns 2 and 3

# Bottom row: One centered plot
ax3 = fig.add_subplot(gs[1, 1:3])  # Centered plot spans columns 1 and 2


# Plot 1: Precursor cell counts
sns.scatterplot(ax=ax1, data=countsdf, x='Age at S1K', y='total_MZ', hue='agebin', s=100)
ax1.set_yscale('log')
ax1.set_ylim(1e5, 1e7)
ax1.set_xscale('log')
ax1.set_xlim(40, 800)
ax1.set_xticks([75, 150, 300, 600])
ax1.set_xticklabels([75, 150, 300, 600])
ax1.set_ylabel('')
ax1.set_xlabel('Host age (days)', fontsize=16)
ax1.set_title('Variation in pool size across the lifetime', fontsize=17)
ax1.tick_params(axis='both', labelsize=15)
ax1.legend(title='Age at BMT', fontsize=12)
ax1.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot 2: Donor fraction
sns.scatterplot(ax=ax2, data=countsdf, x='Days postBMT', y='Nfd_MZ', hue='agebin', s=100, legend=False)
ax2.set_yticks(np.arange(0, 1.25, 0.25))
ax2.set_xscale('log')
ax2.set_xlim(10, 750)
ax2.set_xticks([10, 50, 150, 300, 600])
ax2.set_xticklabels([10, 50, 150, 300, 600])
ax2.set_xlabel('Days post BMT', fontsize=16)
ax2.set_title('Donor fraction in MZ B normalized to chimerism in T1', fontsize=17)
ax2.tick_params(axis='both', labelsize=15)
ax2.set_ylabel('')
ax2.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot 3: Proportion of Ki67+ MZ Donor and Host cells
sns.scatterplot(ax=ax3, data=countsdf, x='Days postBMT', y='frac_Ki67+_MZ_donor', alpha=0.85, s=100, color='brown', label='Donor')
sns.scatterplot(ax=ax3, data=countsdf, x='Days postBMT', y='frac_Ki67+_MZ_host', alpha=0.85, s=100, color='blue', label='Host')
ax3.set_yticks(np.arange(0, 1.25, 0.25))
ax3.set_xscale('log')
ax3.set_xlim(10, 750)
ax3.set_xticks([10, 50, 150, 300, 600])
ax3.set_xticklabels([10, 50, 150, 300, 600])
ax3.set_xlabel('Days post BMT', fontsize=16)
ax3.tick_params(axis='both', labelsize=14)
ax3.set_ylabel('')
ax3.set_title('Proportion of Ki67+ MZ cells', fontsize=17)
ax3.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
ax3.legend()

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig('Custom_MZ1_plots.png')
plt.close()

# Plotting FM cells
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
sns.scatterplot(ax=axes[0, 0], data=countsdf, x='Age at S1K', y='total_FM', s=100)
# sns.scatterplot(ax=axes[0, 0], data=total_mzcounts, x='age.at.S1K', y='total_FM', s=100)
axes[0, 0].set_yscale('log')  # sets the y-axis to a logarithmic scale
axes[0, 0].set_ylim(5*1e6, 1e8)
axes[0, 0].set_ylabel('')
axes[0,0].set_xscale('log')
axes[0,0].set_xlim(50, 800)
axes[0,0].set_xticks([75, 150,  300,  600])
axes[0,0].set_xticklabels([75, 150,  300,  600])
axes[0, 0].set_xlabel('Host Age (days)')
axes[0, 0].set_title('Total FM population size')
axes[0, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
axes[0,0].legend(labels=['315_MV61', '712_mv_BuChi'])


# Plot donor fraction 
sns.scatterplot(ax=axes[0, 1], data=countsdf, x='Days postBMT', y='fd_FM', s=100)
# sns.scatterplot(ax=axes[0, 1], data=total_mzcounts, x='days.post.bmt', y='fd_FM', s=100)
axes[0, 1].set_yticks(np.arange(0, 1.25, 0.25))
axes[0,1].set_xlim(10, 400)
axes[0,1].set_xscale('log')
axes[0, 1].set_xticks([30, 100, 300])
axes[0, 1].set_xticklabels([30, 100, 300])
axes[0, 1].set_xlabel('Days post BMT')    
axes[0, 1].set_title('Donor fraction in FM compartment')
axes[0, 1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot Propotion of Ki67+ MZ Donor cells
sns.scatterplot(ax=axes[1, 0], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_FM_donor', alpha= 0.85, s=100)
# sns.scatterplot(ax=axes[1, 0], data=Ki67counts, x='age.at.S1K', y= 'fraction_ki67FM_donor', alpha= 0.85, s=100)
axes[1, 0].set_yticks(np.arange(0, 1.25, 0.25))
axes[1, 0].set_ylabel('')
axes[1,0].set_xscale('log')
axes[1,0].set_xlim(50, 800)
axes[1,0].set_xticks([75, 150,  300,  600])
axes[1,0].set_xticklabels([75, 150,  300,  600])
axes[1, 0].set_xlabel('Host Age (days)')
axes[1, 0].set_title('Proportion of Ki67+ FM Donor cells')
axes[1, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot Propotion of Ki67+ MZ Donor cells
sns.scatterplot(ax=axes[1, 1], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_FM_host', alpha= 0.85, s=100)
# sns.scatterplot(ax=axes[1, 1], data=Ki67counts, x='age.at.S1K', y= 'fraction_ki67FM_host', alpha= 0.85, s=100)
axes[1,1].set_yticks(np.arange(0, 1.25, 0.25))
axes[1,1].set_ylabel('')
axes[1,1].set_xscale('log')
axes[1,1].set_xlim(50, 800)
axes[1,1].set_xticks([75, 150,  300,  600])
axes[1,1].set_xticklabels([75, 150,  300,  600])
axes[1,1].set_xlabel('Host Age (days)')
axes[1,1].set_title('Proportion of Ki67+ FM Host cells')
axes[1,1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


plt.tight_layout()
plt.savefig('Compare_FM_plots.png')
plt.close()

# Plotting T2.MZP cells
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
sns.scatterplot(ax=axes[0, 0], data=countsdf, x='Age at S1K', y='total_T2_MZP', s=100)
# sns.scatterplot(ax=axes[0, 0], data=total_mzcounts, x='age.at.S1K', y='total_T2.MZP', s=100)
axes[0, 0].set_yscale('log')  # sets the y-axis to a logarithmic scale
axes[0, 0].set_ylim(1e5, 1e7)
axes[0, 0].set_ylabel('')
axes[0,0].set_xscale('log')
axes[0,0].set_xlim(50, 800)
axes[0,0].set_xticks([75, 150,  300,  600])
axes[0,0].set_xticklabels([75, 150,  300,  600])
axes[0, 0].set_xlabel('Host Age (days)')
axes[0, 0].set_title('Total T2 MZP population')
axes[0, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
axes[0,0].legend(labels=['315_MV61', '712_mv_BuChi'])


# Plot donor fraction 
sns.scatterplot(ax=axes[0, 1], data=countsdf, x='Days postBMT', y='fd_T2_MZP', s=100)
# sns.scatterplot(ax=axes[0, 1], data=total_mzcounts, x='days.post.bmt', y='fd_T2.MZP', s=100)
axes[0, 1].set_yticks(np.arange(0, 1.25, 0.25))
axes[0,1].set_xlim(10, 300)
#axes[0,1].set_xscale('log')
axes[0, 1].set_xticks([0, 100, 300])
axes[0, 1].set_xticklabels([ 30, 100, 300])
axes[0, 1].set_xlabel('Days post BMT')    
axes[0, 1].set_title('Donor fraction in T2 MZP compartment')
axes[0, 1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot Propotion of Ki67+ MZ Donor cells
sns.scatterplot(ax=axes[1, 0], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_MZP_donor', alpha= 0.85, s=100)
# sns.scatterplot(ax=axes[1, 0], data=Ki67counts, x='age.at.S1K', y= 'fraction_ki67T2MZP_donor', alpha= 0.85, s=100)
axes[1, 0].set_yticks(np.arange(0, 1.25, 0.25))
axes[1, 0].set_ylabel('')
axes[1,0].set_xscale('log')
axes[1,0].set_xlim(50, 800)
axes[1,0].set_xticks([75, 150,  300,  600])
axes[1,0].set_xticklabels([75, 150,  300,  600])
axes[1, 0].set_xlabel('Host Age (days)')
axes[1, 0].set_title('Proportion of Ki67+ T2 MZP Donor cells')
axes[1, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot Propotion of Ki67+ MZ Donor cells
sns.scatterplot(ax=axes[1, 1], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_MZP_host', alpha= 0.85, s=100)
# sns.scatterplot(ax=axes[1, 1], data=Ki67counts, x='age.at.S1K', y= 'fraction_ki67T2MZP_host', alpha= 0.85, s=100)
axes[1,1].set_yticks(np.arange(0, 1.25, 0.25))
axes[1,1].set_ylabel('')
axes[1,1].set_xscale('log')
axes[1,1].set_xlim(50, 800)
axes[1,1].set_xticks([75, 150,  300,  600])
axes[1,1].set_xticklabels([75, 150,  300,  600])
axes[1,1].set_xlabel('Host Age (days)')
axes[1,1].set_title('Proportion of Ki67+ T2 MZP Host cells')
axes[1,1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


plt.tight_layout()
plt.savefig('Compare_T2.MZP_plots.png')
plt.close()

# Plotting T2 cells
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
mean = countsdf['total_T2'].mean()
median = countsdf['total_T2'].median()
# print('Mean of T2 precursor population:', mean)
print('Median of T2 precursor population:', median)
sns.scatterplot(ax=axes[0, 0], data=countsdf, x='Age at S1K', y='total_T2', s=100)
# sns.scatterplot(ax=axes[0, 0], data=total_mzcounts, x='age.at.S1K', y='total_T2', s=100)
# axes[0,0].axhline(y=mean, linestyle='solid', color= 'red', linewidth=1.5)
axes[0,0].axhline(y=median, linestyle='solid', linewidth=1.5)
axes[0, 0].set_yscale('log')  # sets the y-axis to a logarithmic scale
axes[0, 0].set_ylim(1e5, 1e7)
axes[0, 0].set_ylabel('')
axes[0,0].set_xscale('log')
axes[0,0].set_xlim(50, 800)
axes[0,0].set_xticks([75, 150,  300,  600])
axes[0,0].set_xticklabels([75, 150,  300,  600])
axes[0, 0].set_xlabel('Host Age (days)')
axes[0, 0].set_title('Total T2 population size')
axes[0, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
axes[0,0].legend(labels=['315_MV61', '712_mv_BuChi'])


# Plot donor fraction 
fd_value = countsdf['fd_T2']
fd_value = fd_value[fd_value > 0.25]
mean_fd= fd_value.mean()
median_fd = fd_value.median()
mean= countsdf['fd_T2'].mean()
median = countsdf['fd_T2'].median()
print('Mean of donor fraction in T2 precursor population:', mean_fd)
print('Median of precursor population:', median_fd)
sns.scatterplot(ax=axes[0, 1], data=countsdf, x='Days postBMT', y='fd_T2', s=100)
# sns.scatterplot(ax=axes[0, 1], data=total_mzcounts, x='days.post.bmt', y='fd_T2', s=100)
axes[0,1].axhline(y=mean_fd, linestyle='solid',  color= 'red', linewidth=1.5)
axes[0,1].axhline(y=median_fd, linestyle='solid', linewidth=1.5)
axes[0, 1].set_yticks(np.arange(0, 1.25, 0.25))
axes[0,1].set_xlim(0, 300)
# axes[0,1].set_xscale('log')
axes[0, 1].set_xticks([30, 100, 300])
axes[0, 1].set_xticklabels([ 30, 100, 300])
axes[0, 1].set_xlabel('Days post BMT')    
axes[0, 1].set_title('Donor fraction in T2 compartment')
axes[0, 1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot Propotion of Ki67+ T2 Donor cells
mean= countsdf['frac_Ki67+_T2_donor'].mean()
median = countsdf['frac_Ki67+_T2_donor'].median()
# print('Mean of K67+ fraction of donor in T2 precursor population:', mean)
print('Median of fraction of Ki67+ donor cells in T2 precursor population:', median)
sns.scatterplot(ax=axes[1, 0], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_donor', alpha= 0.85, s=100)
# sns.scatterplot(ax=axes[1, 0], data=Ki67counts, x='age.at.S1K', y= 'fraction_ki67T2_donor', alpha= 0.85, s=100)
# axes[1,0].axhline(y=mean, linestyle='solid',  color= 'red', linewidth=1.5)
axes[1,0].axhline(y=median, linestyle='solid', linewidth=1.5)
axes[1, 0].set_yticks(np.arange(0, 1.25, 0.25))
axes[1, 0].set_ylabel('')
axes[1,0].set_xscale('log')
axes[1,0].set_xlim(50, 800)
axes[1,0].set_xticks([75, 150,  300,  600])
axes[1,0].set_xticklabels([75, 150,  300,  600])
axes[1, 0].set_xlabel('Host Age (days)')
axes[1, 0].set_title('Proportion of Ki67+ T2 Donor cells')
axes[1, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot Propotion of Ki67+ T2 Host cells
mean= countsdf['frac_Ki67+_T2_host'].mean()
median = countsdf['frac_Ki67+_T2_host'].median()
# print('Mean of Ki67+ fraction of host in T2 precursor population:', mean)
print('Median of fraction of Ki67+ host cells in T2 precursor population:', median)
sns.scatterplot(ax=axes[1, 1], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_host', alpha= 0.85, s=100)
# sns.scatterplot(ax=axes[1, 1], data=Ki67counts, x='age.at.S1K', y= 'fraction_ki67T2_host', alpha= 0.85, s=100)
# axes[1,1].axhline(y=mean, linestyle='solid',  color= 'red', linewidth=1.5)
axes[1,1].axhline(y=median, linestyle='solid', linewidth=1.5)
axes[1,1].set_yticks(np.arange(0, 1.25, 0.25))
axes[1,1].set_ylabel('')
axes[1,1].set_xscale('log')
axes[1,1].set_xlim(50, 800)
axes[1,1].set_xticks([75, 150,  300,  600])
axes[1,1].set_xticklabels([75, 150,  300,  600])
axes[1,1].set_xlabel('Host Age (days)')
axes[1,1].set_title('Proportion of Ki67+ T2 Host cells')
axes[1,1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


plt.tight_layout()
plt.savefig('Compare_T2_plots.png')
plt.close()

# Plotting T1_AA4.1+ cells
mean= countsdf['total_T1_AA4.1+'].mean()
median = countsdf['total_T1_AA4.1+'].median()
print('Mean of T1 precursor population:', mean)
# print('Median of precursor population:', median)

fig, axes = plt.subplots(2, 2, figsize=(15, 10))
sns.scatterplot(ax=axes[0, 0], data=countsdf, x='Age at S1K', y='total_T1_AA4.1+', s=100)
# sns.scatterplot(ax=axes[0, 0], data=total_mzcounts, x='age.at.S1K', y='total_Transitional1 AA4.1+', s=100)
axes[0,0].axhline(y=mean, linestyle='solid', linewidth=1.5)
# axes[0,0].axhline(y=median, linestyle='solid', linewidth=1.5)
axes[0, 0].set_yscale('log')  # sets the y-axis to a logarithmic scale
axes[0, 0].set_ylim(1e5, 1e7)
axes[0, 0].set_ylabel('')
axes[0,0].set_xscale('log')
axes[0,0].set_xlim(50, 800)
axes[0,0].set_xticks([75, 150,  300,  600])
axes[0,0].set_xticklabels([75, 150,  300,  600])
axes[0, 0].set_xlabel('Host Age (days)')
axes[0, 0].set_title('Total T1_IgMhi population size')
axes[0, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
axes[0,0].legend(labels=['315_MV61', 'Mean'])


# Plot donor fraction 
fd_value = countsdf['fd_T1_AA4.1+']
fd_value = fd_value[fd_value > 0.25]
mean_fd= fd_value.mean()
median_fd = fd_value.median()
print('Mean of donor fraction in T1 precursor population:', mean_fd)
sns.scatterplot(ax=axes[0, 1], data=countsdf, x='Days postBMT', y='fd_T1_AA4.1+', s=100)
# sns.scatterplot(ax=axes[0, 1], data=total_mzcounts, x='days.post.bmt', y='fd_T1_AA4.1+', s=100)
axes[0, 1].axhline(y=mean_fd, linestyle='solid', linewidth=1.5)
# axes[0, 1].axhline(y=median_fd, color = 'red', linestyle='solid', linewidth=1.5)
axes[0, 1].set_yticks(np.arange(0, 1.25, 0.25))
axes[0,1].set_xlim(0, 300)
#axes[0,1].set_xscale('log')
axes[0, 1].set_xticks([30, 100, 300])
axes[0, 1].set_xticklabels([ 30, 100, 300])
axes[0, 1].set_xlabel('Days post BMT')    
axes[0, 1].set_title('Donor fraction in T1_IgMhi compartment')
axes[0, 1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot Propotion of Ki67+ Donor cells

mean = countsdf['frac_Ki67+_T1_AA4.1+_donor'].mean()
median = countsdf['frac_Ki67+_T1_AA4.1+_donor'].median()
print('Mean of Ki67+ T1 precursor population:', mean)
sns.scatterplot(ax=axes[1, 0], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_AA4.1+_donor', alpha= 0.85, s=100)
# sns.scatterplot(ax=axes[1, 0], data=Ki67counts, x='age.at.S1K', y= 'fraction_ki67T1_AA4.1_donor', alpha= 0.85, s=100)
axes[1,0].axhline(y=mean, linestyle='solid', linewidth=1.5)
# axes[1,0].axhline(y=median, linestyle='solid', linewidth=1.5)
axes[1, 0].set_yticks(np.arange(0, 1.5, 0.25))
axes[1, 0].set_ylabel('')
axes[1,0].set_xscale('log')
axes[1,0].set_xlim(50, 800)
axes[1,0].set_xticks([75, 150,  300,  600])
axes[1,0].set_xticklabels([75, 150,  300,  600])
axes[1, 0].set_xlabel('Host Age (days)')
axes[1, 0].set_title('Proportion of Ki67+ T1_IgMhi Donor cells')
axes[1, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot Propotion of Ki67+ Host cells
mean = countsdf['frac_Ki67+_T1_AA4.1+_host'].mean()
median = countsdf['frac_Ki67+_T1_AA4.1+_host'].median()
print('Median of Ki67+ T1 precursor population:', median)
sns.scatterplot(ax=axes[1, 1], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_AA4.1+_host', alpha= 0.85, s=100)
# sns.scatterplot(ax=axes[1, 1], data=Ki67counts, x='age.at.S1K', y= 'fraction_ki67T1_AA4.1_host', alpha= 0.85, s=100)
axes[1,1].axhline(y=median, linestyle='solid', linewidth=1.5)
axes[1,1].set_yticks(np.arange(0, 1.25, 0.25))
axes[1,1].set_ylabel('')
axes[1,1].set_xscale('log')
axes[1,1].set_xlim(50, 800)
axes[1,1].set_xticks([75, 150,  300,  600])
axes[1,1].set_xticklabels([75, 150,  300,  600])
axes[1,1].set_xlabel('Host Age (days)')
axes[1,1].set_title('Proportion of Ki67+ T1_IgMhi Host cells')
axes[1,1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


plt.tight_layout()
plt.savefig('Compare_T1 AA4.1+_plots.png')
plt.close()

# Plotting T1 cells
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
sns.scatterplot(ax=axes[0, 0], data=countsdf, x='Age at S1K', y='total_T1', s=100)
sns.scatterplot(ax=axes[0, 0], data=countsdf, x='Age at S1K', y='total_T1_AA4.1+', s=100)
axes[0, 0].set_yscale('log')  # sets the y-axis to a logarithmic scale
axes[0, 0].set_ylim(1e5, 1e7)
axes[0, 0].set_ylabel('')
axes[0,0].set_xscale('log')
axes[0,0].set_xlim(50, 800)
axes[0,0].set_xticks([75, 150,  300,  600])
axes[0,0].set_xticklabels([75, 150,  300,  600])
axes[0, 0].set_xlabel('Host Age (days)')
axes[0, 0].set_title('Total T1(diff subset) population size')
axes[0, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
axes[0,0].legend(labels=['T1_IgMHi', 'T1_IgDlo'])


# Plot donor fraction 
sns.scatterplot(ax=axes[0, 1], data=countsdf, x='Days postBMT', y='fd_T1', s=100)
sns.scatterplot(ax=axes[0, 1], data=countsdf, x='Days postBMT', y='fd_T1_AA4.1+', s=100)
axes[0, 1].set_yticks(np.arange(0, 1.25, 0.25))
axes[0,1].set_xlim(10, 400)
axes[0,1].set_xscale('log')
axes[0, 1].set_xticks([30, 100, 300])
axes[0, 1].set_xticklabels([ 30, 100, 300])
axes[0, 1].set_xlabel('Days post BMT')    
axes[0, 1].set_title('Donor fraction in T1 compartment')
axes[0, 1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot Propotion of Ki67+ Donor cells
sns.scatterplot(ax=axes[1, 0], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_donor', alpha= 0.85, s=100)
sns.scatterplot(ax=axes[1, 0], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_AA4.1+_donor', alpha= 0.85, s=100)
axes[1, 0].set_yticks(np.arange(0, 1.25, 0.25))
axes[1, 0].set_ylabel('')
axes[1,0].set_xscale('log')
axes[1,0].set_xlim(50, 800)
axes[1,0].set_xticks([75, 150,  300,  600])
axes[1,0].set_xticklabels([75, 150,  300,  600])
axes[1, 0].set_xlabel('Host Age (days)')
axes[1, 0].set_title('Proportion of Ki67+ T1 Donor cells')
axes[1, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot Propotion of Ki67+ Host cells
sns.scatterplot(ax=axes[1, 1], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_host', alpha= 0.85, s=100)
sns.scatterplot(ax=axes[1, 1], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_AA4.1+_host', alpha= 0.85, s=100)
axes[1,1].set_yticks(np.arange(0, 1.25, 0.25))
axes[1,1].set_ylabel('')
axes[1,1].set_xscale('log')
axes[1,1].set_xlim(50, 800)
axes[1,1].set_xticks([75, 150,  300,  600])
axes[1,1].set_xticklabels([75, 150,  300,  600])
axes[1,1].set_xlabel('Host Age (days)')
axes[1,1].set_title('Proportion of Ki67+ T1 Host cells')
axes[1,1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


plt.tight_layout()
plt.savefig('Compare_T1_plots.png')
plt.close()

#Plot T1+T2 

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
sns.scatterplot(ax=axes[0,0], data=countsdf, x='Age at S1K', y='total_T1_AA4.1+', s=75, label='T1')
sns.scatterplot(ax=axes[0,0], data=countsdf, x='Age at S1K', y='total_T2', s=75, label='T2')
axes[0,0].axhline(y=mean, linestyle='solid', linewidth=1.5)
# axes[0,0].axhline(y=median, linestyle='solid', linewidth=1.5)
axes[0,0].set_yscale('log')  # sets the y-axis to a logarithmic scale
axes[0,0].set_ylim(1e5, 1e7)
axes[0, 0].set_ylabel('')
axes[0,0].set_xscale('log')
axes[0,0].set_xlim(50, 800)
axes[0,0].set_xticks([75, 150,  300,  600])
axes[0,0].set_xticklabels([75, 150,  300,  600])
axes[0,0].set_xlabel('Host Age (days)')
axes[0,0].set_title('Total Transitional (T1 and T2) population size')
axes[0,0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
# axes[0,0].legend(labels=['315_MV61', 'Mean'])


# Plot donor fraction 
fd_value = countsdf1['fd']
fd_value = fd_value[fd_value > 0.25]
mean_fd= fd_value.mean()
median_fd = fd_value.median()
print('Mean of donor fraction in T1 precursor population:', mean_fd)
sns.scatterplot(ax=axes[0, 1], data=countsdf, x='Days postBMT', y='fd_T1_AA4.1+', s=75)
sns.scatterplot(ax=axes[0, 1], data=countsdf, x='Days postBMT', y='fd_T2', s=75)
# axes[0, 1].axhline(y=mean_fd, linestyle='solid', linewidth=1.5)
# axes[0, 1].axhline(y=median_fd, color = 'red', linestyle='solid', linewidth=1.5)
axes[0, 1].set_yticks(np.arange(0, 1.25, 0.25))
axes[0,1].set_xlim(0, 300)
#axes[0,1].set_xscale('log')
axes[0, 1].set_xticks([30, 100, 300])
axes[0, 1].set_xticklabels([ 30, 100, 300])
axes[0, 1].set_xlabel('Days post BMT')    
axes[0, 1].set_title('Donor fraction in transitional compartment')
axes[0, 1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot Propotion of Ki67+ Donor cells

# mean = countsdf['frac_Ki67+_T1_AA4.1+_donor'].mean()
# median = countsdf['frac_Ki67+_T1_AA4.1+_donor'].median()
# print('Mean of Ki67+ T1 precursor population:', mean)
sns.scatterplot(ax=axes[1, 0], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_AA4.1+_donor', alpha= 0.85, s=75)
sns.scatterplot(ax=axes[1, 0], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_donor', alpha= 0.85, s=75)
# axes[1,0].axhline(y=mean, linestyle='solid', linewidth=1.5)
# axes[1,0].axhline(y=median, linestyle='solid', linewidth=1.5)
axes[1, 0].set_yticks(np.arange(0, 1.5, 0.25))
axes[1, 0].set_ylabel('')
axes[1,0].set_xscale('log')
axes[1,0].set_xlim(50, 800)
axes[1,0].set_xticks([75, 150,  300,  600])
axes[1,0].set_xticklabels([75, 150,  300,  600])
axes[1, 0].set_xlabel('Host Age (days)')
axes[1, 0].set_title('Fraction of Ki67+ donor cells')
axes[1, 0].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot Propotion of Ki67+ Host cells
# mean = countsdf['frac_Ki67+_T1_AA4.1+_host'].mean()
# median = countsdf['frac_Ki67+_T1_AA4.1+_host'].median()
# print('Median of Ki67+ T1 precursor population:', median)
sns.scatterplot(ax=axes[1, 1], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_AA4.1+_host', alpha= 0.85, s=75)
sns.scatterplot(ax=axes[1, 1], data=countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_host', alpha= 0.85, s=100)
# axes[1,1].axhline(y=median, linestyle='solid', linewidth=1.5)
axes[1,1].set_yticks(np.arange(0, 1.25, 0.25))
axes[1,1].set_ylabel('')
axes[1,1].set_xscale('log')
axes[1,1].set_xlim(50, 800)
axes[1,1].set_xticks([75, 150,  300,  600])
axes[1,1].set_xticklabels([75, 150,  300,  600])
axes[1,1].set_xlabel('Host Age (days)')
axes[1,1].set_title('Fraction of Ki67+ host cells')
axes[1,1].grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


plt.tight_layout()
plt.savefig('Compare_T1+T2.png')
plt.close()

# Plot total MZ B cells for young + old hosts
# plot 1 graph

'''fig = plt.figure(figsize=(10, 6))
sns.scatterplot(data=countsdf, x='Age at S1K', y='total_MZ', s=100, alpha=0.85, label = 'Busulfan Chimeras')
sns.scatterplot(data=df_ont, x='age.at.S1K', y='MZ_total', alpha=0.85, s=100, color = 'red', label='Young WT mice')
plt.yscale('log')  # sets the y-axis to a logarithmic scale
plt.ylim(5*1e3, 1e7)
plt.ylabel('')
plt.xscale('log')
plt.xlim(0, 800)
plt.xticks([10, 30, 100, 300, 600], [10, 30, 100, 300, 600])
# axes.set_xticks([75, 150, 300, 600])
# axes.set_xticklabels([75, 150, 300, 600])
plt.xlabel('Host Age (days)')
plt.title('Total MZ B population size')
plt.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
# plt.legend(labels=['315_MV61', '712_mv_BuChi'])
plt.tight_layout()
plt.savefig('Total_MZ_plot.png')
plt.close()'''
