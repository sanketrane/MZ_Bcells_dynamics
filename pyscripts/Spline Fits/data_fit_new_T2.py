import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
# from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.preprocessing import MinMaxScaler, StandardScaler

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# Read specific sheets by name
Total_countsdf = pd.read_csv("/home/apoorva/Desktop/MZB-new-analysis/data/New_data/new_data_MZB.csv") # counts for the host compartment MZ B cell population

# Ignore fraction of Ki67 donor and host < 0.7
Total_countsdf1 = Total_countsdf[Total_countsdf['frac_Ki67+_T2_donor'] >= 0.7]

Total_countsdf3 = Total_countsdf[Total_countsdf['frac_Ki67+_T2_host'] >= 0.7]

# Remove an obervation with days postBmt 117 and 82 with Nfd <0.25
Total_countsdf2 = Total_countsdf[Total_countsdf['Days postBMT'] != 117]
Total_countsdf2 = Total_countsdf[~((Total_countsdf['Days postBMT'] == 82) & (Total_countsdf['fd_T2'] < 0.25))]


#Curve fitting function
def exponential_model(x, a, b):
    return a*(np.exp(b* (x-72)))
def linear_model(x, a, b):
    return a*(x-40) + b
def lin_model(x):
    return np.exp(0.002*(x-40) + 17.19)
def gaussian_model(x, a, b):
    return  np.exp(-b*(x - 72)) + a
def gau_function(x):
    return  np.exp(-0.5*(x-69)) + 0.17
def exp_func(x, a, b):
    return   a*(1- np.exp(-b * (x-10)))
def nfd_func(x):
    return 0.83*(1 - np.exp(-0.08*(x-13))) 

x_data = np.array(Total_countsdf1['Age at S1K'])
x_data1 = np.array(Total_countsdf3['Age at S1K'])
x_bmt = np.array(Total_countsdf1['Age at BMT'])
mean_age_at_BMT = x_bmt.mean()
median_age_at_BMT = Total_countsdf1['Age at BMT'].median()
mode_age_at_BMT = Total_countsdf1['Age at BMT'].mode()

x_fd = np.array(Total_countsdf['Days postBMT'])
# y_precursor = np.array(Total_countsdf1['total_counts'])
y_fd = np.array(Total_countsdf['fd_T2'])
y_donor= np.array(Total_countsdf1['frac_Ki67+_T2_donor'])
y_host= np.array(Total_countsdf3['frac_Ki67+_T2_host'])

# Initial guess for the parameters
initial_guess_p = [1, 0.05]
#initial_guess_d= [ 0.02, 0.25]

# Scale the input data
scaler = MinMaxScaler()

 # Perform curve fitting

popt, pcov = curve_fit(linear_model, x_data, y_donor, p0=initial_guess_p, maxfev=10000)
popt1, pcov1 = curve_fit(linear_model, x_data1, y_host, p0=initial_guess_p, maxfev=10000)
# popt3, pcov3 = curve_fit(exponential_model, x_data, y_precursor_log, p0=initial_guess_p, maxfev=10000)
# popt1, pcov1 = curve_fit(gaussian_model, x_data, y_donor, maxfev=10000)
popt2, pcov2 = curve_fit(exp_func, x_fd, y_fd, maxfev=10000)
# # popt, pcov = curve_fit(exponential_model, x_data_scaled, y_host, p0=initial_guess, maxfev=10000)

print("Fitted parameters for FD(T2):", popt2)
# print("Covariance matrix for FD:", pcov2)   
# Print the fitted parameters
# print("Fitted parameters for Precursor population:", popt)
print("Fitted parameters for Ki67 Donor Fraction:", popt)
print("Fitted parameters for Ki67 Host Fraction:", popt1)
# print("Fitted parameters for FD:", popt2)

#unsacle the data for log plot
#x_fit_scaled = scaler.transform(x_data.reshape(-1, 1)).flatten()
x_dense = np.linspace(59, 750, 1000)
x_bmt_dense = np.linspace(41, 101, 1000)
x_densefd = np.linspace(13, 690, 1000)
# x_fit_scaled = scaler.transform(x_data.reshape(-1, 1)).flatten()
# y_fit_log = exponential_model(x_data, *popt3)
# y_fit1_log = linear_model(x_data, *popt)
# y_fit1_scaled = np.exp(y_fit1_log)

# Inverse transform the scaled y values to get them back to the original scale
#y_fitp = scaler.inverse_transform(y_fit_scaled.reshape(-1, 1)).flatten()

# Create a 3x2 grid of subplots
fig, axes = plt.subplots(2, 2, figsize=(15, 10))


# Plot precursor cell counts
mean= Total_countsdf['total_T2'].mean()
median = Total_countsdf['total_T2'].median()
# variance_T2 = Total_countsdf['total_T2'].var()
# variance_T1 = Total_countsdf['total_T1_AA4.1+'].var()
# std_dev_t1 = Total_countsdf['total_T1_AA4.1+'].std()
# std_dev_t2 = Total_countsdf['total_T2'].std()

# print('Ratio of stds of T2 and T1 precursor populations:', std_dev_t2/std_dev_t1)
# print('Variance of T2 precursor population:', variance_T2)
# print('Variance of T1 precursor population:', variance_T1)
print('Mean of T2 precursor population:', mean)
print('Median of T2 precursor population:', median)

fig, axes = plt.subplots(1, 3, figsize=(23, 5))
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]
sns.scatterplot(ax=axes[0], data=Total_countsdf, x='Age at S1K', y='total_T2', s=75, alpha= 0.85, linewidth=0.5,color='black')
# sns.scatterplot(ax=axes[0, 0], data=total_mzcounts, x='age.at.S1K', y='total_Transitional1 AA4.1+', s=100)
ax1.axhline(y=mean, linestyle='solid', linewidth=1.5, color='black')
ax1.axhline(y=median, color='red', linestyle='dashed')
# axes[0,0].axhline(y=median, linestyle='solid', linewidth=1.5)
ax1.set_yscale('log')  # sets the y-axis to a logarithmic scale
ax1.set_ylim(1e5, 1e7)
ax1.set_ylabel('')
ax1.set_xscale('log')
ax1.set_xlim(50, 800)
ax1.set_xticks([75, 150,  300,  600])
ax1.set_xticklabels([75, 150,  300,  600])
ax1.set_xlabel('Host Age (days)', fontsize=16)
ax1.set_title('Total T2 population size')
ax1.tick_params(axis='both') 
ax1.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot donor fraction FO compartment
# fd_mean= Total_countsdf['fd_T2'].mean()
# fd_median = Total_countsdf['fd_T2'].median()
# print('Mean of donor fraction in T2 precursor population:', fd_mean)
# print('Median of donor fraction in T2 precursor population:', fd_median)


sns.scatterplot(ax=axes[1], data=Total_countsdf, x='Age at S1K', y='fd_T2', alpha=0.85, s=75, linewidth=0.5,color='black')
y_fitd= exp_func(x_densefd, *popt2)
y_fitd1 = nfd_func(x_densefd)
x_dense1 = x_densefd + mean_age_at_BMT
# ax2.plot(x_dense1, y_fitd1, color='black')
ax2.set_yticks(np.arange(0, 1.25, 0.25)[[0, 2, 4]])
ax2.set_xlim(50, 800)
ax2.set_xscale('log')
ax2.set_xticks([75, 150, 300, 600])
ax2.set_xticklabels([75, 150, 300, 600])
# ax2.set_xlabel('Days post BMT')  
ax2.set_ylabel('')  
ax2.set_title('Donor fraction')
ax2.tick_params(axis='both')
ax2.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot Propotion of Ki67+ FO cells
# mean_d= Total_countsdf['frac_Ki67+_T2_donor'].mean()
# median_d = Total_countsdf['frac_Ki67+_T2_donor'].median()
# # print('Mean of Ki67+ donor T2 precursor population:', mean_d)
# print('Median of Ki67+ donor T2 precursor population:', median_d)
# mean_h = Total_countsdf['frac_Ki67+_T2_host'].mean()
# print('Mean of Ki67+ host T2 precursor population:', mean_h)
# median_h = Total_countsdf['frac_Ki67+_T2_host'].median()
# print('Median of Ki67+ host T2 precursor population:', median_h)

sns.scatterplot(ax=axes[2], data=Total_countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_host', alpha= 0.85, s=75,  linewidth=0.5,color='#ab2239')
y_fitd= linear_model(x_dense, *popt)
y_fith= linear_model(x_dense, *popt1)
sns.scatterplot(ax=axes[2], data=Total_countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_donor',alpha= 0.85, s=75,  linewidth=0.5,color='#3b9bb3')

# ax3.axhline(y=median_d, color='red', linestyle='solid', label= 'Fitted donor curve', linewidth=0.85)
ax3.plot(x_dense, y_fitd, linestyle='solid', label= 'Fitted donor curve', linewidth=0.85, color='#3b9bb3')

# ax3.axhline(y=0.825, color='blue', linestyle='solid', label= 'Fitted host curve', linewidth=0.85)
ax3.plot(x_dense, y_fith, color='#ab2239', linestyle='solid', label= 'Fitted host curve', linewidth=0.85)
ax3.set_yticks(np.arange(0, 1.25, 0.25))
ax3.set_ylim(0.25, 1.05)
ax3.set_ylabel('')
ax3.set_xscale('log')
# ax3.set_yscale('log')
ax3.set_xlim(50, 750)
ax3.set_xticks([75, 150,  300, 600])
ax3.set_xticklabels([75, 150,  300, 600])
ax3.set_xlabel('Host Age (days)', fontsize=16)
ax3.set_title('Proportion of Ki67+ T2 cells', fontsize=17)
ax3.tick_params(axis='both')
ax3.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
ax3.legend(labels=['Host', 'Donor'])
plt.tight_layout()
plt.savefig('Precursor_T2_new_plots.pdf', dpi=300)

                                       


