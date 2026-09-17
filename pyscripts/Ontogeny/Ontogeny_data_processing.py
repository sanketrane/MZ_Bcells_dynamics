import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read specific sheets by name
countsdf = pd.read_excel("data/Ontogeny_MZT1_AA41.xlsx", sheet_name="SPLEEN_Cell number") # counts for the host compartment MZ B cell population


countsdf = countsdf.drop(columns=['T1 AA4.1-', 'Transitional1'])
print(countsdf.head(10))


# Create two dataframes for Ki67 with + and without +
df_Ki67_pos = countsdf[countsdf['Ki67'].str.contains('\+', na=False)].reset_index(drop=True)
df_total = countsdf[~countsdf['Ki67'].str.contains('\+', na=False)].reset_index(drop=True)


#Drop Ki67 column
df_Ki67_pos = df_Ki67_pos.drop(columns=['Ki67'])
df_total= df_total.drop(columns=['Ki67'])


# Merge the two DataFrames on the specified columns
ont_df = df_Ki67_pos.merge(df_total, on = ['age.at.S1K', 'mouse'], suffixes=("_Ki67+", "_total"))
print(ont_df.head(10))


# Calculate fraction of Ki67+ Mz cells
ont_df['frac_MZ_Ki67+'] = (ont_df['MZ_Ki67+'] / ont_df['MZ_total'])
print(ont_df['frac_MZ_Ki67+'])

# Calculate the fraction of T1 cells that are Ki67+
ont_df['frac_T1_AA4.1+_Ki67+'] = (ont_df['T1 AA4.1+_Ki67+'] / ont_df['T1 AA4.1+_total'])

#Drop values with Ki67 value less than 0.1

ont_df = ont_df[ont_df['frac_MZ_Ki67+'] > 0.01]
ont_df = ont_df[ont_df['frac_T1_AA4.1+_Ki67+'] > 0.01]


# Calculate mean & median of fraction of Ki67+ T1 cells
mean_T1_Ki67 = ont_df['frac_T1_AA4.1+_Ki67+'].mean()
median_T1_Ki67 = ont_df['frac_T1_AA4.1+_Ki67+'].median()
# Save the DataFrame to a CSV file
ont_df.to_csv("ontogeny_counts.csv", index=False)

# Fit a model to T1 cell counts

# Plot the total MZ and T1 cell counts as subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 10))

# Plot total T1 cell counts
sns.scatterplot(ax=ax1, data=ont_df, x='age.at.S1K', y='T1 AA4.1+_total', s=100, label='Young WT mice', legend=False)
ax1.set_yscale('log')  # sets the y-axis to a logarithmic scale
ax1.set_ylim(1e5, 1e8)
ax1.set_xscale('log')
ax1.set_xlim(10, 200)
ax1.set_xticks([10, 30, 100, 200])
ax1.set_xticklabels([10, 30, 100, 200])
# ax1.set_ylabel('Total MZ B cell numbers', fontsize=14)
ax1.set_xlabel('Mouse Age (days)')
ax1.set_ylabel('')
ax1.set_title('T1 AA4.1 B cells', fontweight='bold')
ax1.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot total MZ cell counts
sns.scatterplot(ax=ax2, data=ont_df, x='age.at.S1K', y='MZ_total', s=100, label='Young WT mice', legend=False)
ax2.set_yscale('log')  # sets the y-axis to a logarithmic scale
ax2.set_ylim(5*1e3, 1e7)
ax2.set_xscale('log')
ax2.set_xlim(10, 200)
ax2.set_xticks([10, 30, 100, 200])
ax2.set_xticklabels([10, 30, 100, 200])
ax2.set_xlabel('Mouse Age (days)')
ax2.set_ylabel('')
ax2.set_title('MZ B cells', fontweight='bold')
ax2.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot percentage of Ki67+ MZ cells
sns.scatterplot(ax=ax4, data=ont_df, x='age.at.S1K', y='frac_MZ_Ki67+',  s=100, label='Young WT mice', legend=False)
ax4.set_ylim(0, 1)
# ax3.set_yscale('logit')
ax4.set_xscale('log')
ax4.set_xlim(10, 200)
ax4.set_xticks([10, 30, 100, 200])
ax4.set_xticklabels([10, 30, 100, 200])
ax4.set_xlabel('Mouse Age (days)')
ax4.set_ylabel('')
ax4.set_title('Fraction of Ki67+ MZ B cells', fontweight='bold')
ax4.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot fraction of Ki67+ T1 cells
sns.scatterplot(ax=ax3, data=ont_df, x='age.at.S1K', y='frac_T1_AA4.1+_Ki67+', s =100, label='Young WT mice', legend=False)
ax3.set_ylim(0, 1)
# ax4.axhline(mean_T1_Ki67, color='r', linestyle='solid', label='Mean')
# ax4.axhline(median_T1_Ki67, color='b', linestyle='dashed', label='Median')
# ax4.set_yscale('log')
ax3.set_xscale('log')
ax3.set_xlim(10, 200)
ax3.set_xticks([10, 30, 100, 200])
ax3.set_xticklabels([10, 30, 100, 200])
ax3.set_xlabel('Mouse Age (days)')
ax3.set_ylabel('')
ax3.set_title('Fraction of Ki67+ T1 AA4.1 B cells', fontweight='bold')
ax3.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)



plt.tight_layout()
plt.savefig('ontogeny_plot.png')










