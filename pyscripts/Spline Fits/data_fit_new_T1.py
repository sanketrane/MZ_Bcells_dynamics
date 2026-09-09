import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.preprocessing import MinMaxScaler, StandardScaler

# PROJECT_ROOT = '/Users/apoorvasingh/Desktop/MZB-new-analysis'
# os.chdir(PROJECT_ROOT)
# print('Working directory:', os.getcwd())

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Read specific sheets by name
Total_countsdf = pd.read_csv("/home/apoorva/Desktop/MZB-new-analysis/data/New_data/new_data_MZB.csv") # counts for the host compartment MZ B cell population
# Ignore fraction of Ki67 donor and host < 0.7
Total_countsdf1 = Total_countsdf[Total_countsdf['frac_Ki67+_T2_donor'] >= 0.7]
# Total_countsdf1 = Total_countsdf1[Total_countsdf1['Age at S1K'] <= 320]
Total_countsdf2 = Total_countsdf1[Total_countsdf['frac_Ki67+_T2_host'] >= 0.7]
# Total_countsdf2 = Total_countsdf2[Total_countsdf2['Age at S1K'] <= 320]
# Remove an obervation with days postBmt 117 and 82 with Nfd <0.25
# Total_countsdf1 = Total_countsdf[Total_countsdf['Days postBMT'] != 117]
# Total_countsdf1 = Total_countsdf[~((Total_countsdf['Days postBMT'] == 82) & (Total_countsdf['fd_T1_AA4.1+'] < 0.25))]
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
    return  np.exp(-0.12*(x-69)) + 0.17
def exp_func(x, a, b):
    return   a*(1- np.exp(-b * (x-10)))
def exp_func1(x):
    return 0.85*(1- np.exp(-0.1 * (x-10)))
def nfd_func(x):
    return 0.925*(1 - np.exp(-0.155*(x)))


x_data = np.array(Total_countsdf1['Age at S1K'])
x_data1 = np.array(Total_countsdf2['Age at S1K'])
x_fd = np.array(Total_countsdf['Days postBMT'])

y_fd = np.array(Total_countsdf['fd_T1_AA4.1+'])
y_fdt2 = np.array(Total_countsdf['fd_T2'])

y_donor= Total_countsdf1['frac_Ki67+_T1_AA4.1+_donor']
y_host= Total_countsdf2['frac_Ki67+_T1_AA4.1+_host']

#logit transform for donor fraction
y_fd_logit = np.log(y_fd / (1 - y_fd))

# Initial guess for the parameters
initial_guess_p = [1, 0.05, 1]
#initial_guess_d= [ 0.02, 0.25]

# Scale the input data
scaler = MinMaxScaler()

# Perform curve fitting
popt, pcov = curve_fit(linear_model, x_data, y_donor, maxfev=10000)
popt1, pcov1 = curve_fit(linear_model, x_data1, y_host, maxfev=10000)
popt2, pcov2 = curve_fit(exp_func, x_fd, y_fd, maxfev=10000)
popt3, pcov3 = curve_fit(exp_func, x_fd, y_fdt2, maxfev=10000)

# Print the fitted parameters
print("Fitted parameters for Ki67 Donor Fraction:", popt)
print("Fitted parameters for Ki67 Host Fraction:", popt1)
print("Fitted parameters for Fd(T1):", popt2)
# print("Covariance matrix for FD(T1):", pcov2)   

print("Fitted parameters for Fd(T2):", popt3)
# print("Covariance matrix for FD(T2):", pcov3)   

#unsacle the data for log plot
#x_fit_scaled = scaler.transform(x_data.reshape(-1, 1)).flatten()
x_dense = np.linspace(59, 731, 1000)
x_densefd = np.linspace(10, 690, 1000)

# Plot precursor cell counts
mean= Total_countsdf['total_T1_AA4.1+'].mean()
median = Total_countsdf['total_T1_AA4.1+'].median()
variance = Total_countsdf['total_T1_AA4.1+'].var()
# print('Variance of T1 precursor population:', variance)
print('Mean of T1 precursor population:', mean)
print('Median of precursor population:', median)

#plot 4*4 grid 

fig, axes = plt.subplots(1, 3, figsize=(18, 4))
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]
# Plot total MZ cell counts
sns.scatterplot(ax=ax1, data=Total_countsdf, x='Age at S1K', y='total_T1_AA4.1+', alpha= 0.85,s=75, linewidth=0.5,color='black')
ax1.axhline(y=mean, linestyle='solid', color='black')
ax1.set_yscale('log')  # sets the y-axis to a logarithmic scale
ax1.set_ylim(1e5, 1e7)
ax1.set_ylabel('')
ax1.set_xscale('log')
ax1.set_xlim(50, 800)
ax1.set_xticks([75, 150,  300,  600])
ax1.set_xticklabels([75, 150,  300,  600])
ax1.set_xlabel('Host Age (days)')
ax1.set_title('Total T1 population size')
ax1.tick_params(axis='both') 
ax1.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)

# Plot donor fraction FO compartment

sns.scatterplot(ax=ax2, data=Total_countsdf, x='Days postBMT', y='fd_T1_AA4.1+', alpha=0.85, s=75,  linewidth=0.5,color='black')
# sns.scatterplot(ax=axes[1], data=Total_countsdf, x='Days postBMT', y='fd_T2', color= 'red',alpha=0.85, s=100)
# y_fitd= exp_func(x_densefd, *popt2)
# y_fitd= exp_func(x_densefd, *popt2)
# y_fitt2= exp_func(x_densefd, *popt3)
# y_fitd1 = nfd_func(x_densefd)
# y_fitd2 = nfd_func1(x_densefd)
y_fit_T1 = exp_func1(x_densefd)
# y_fit_T1 = exp_func1(x_densefd, *popt2)


ax2.plot(x_densefd, y_fit_T1, color='black')
# axes[1].plot(x_densefd, y_fitd, label='T1')
# axes[1].plot(x_densefd, y_fitt2, linestyle='dashed', color='blue', label='T2')
# axes[1].plot(x_densefd, y_fitt2, color='red',label='Fitted curve')
ax2.set_yticks(np.arange(0, 1.25, 0.25)[[0, 2, 4]])
ax2.set_xlim(9, 700)
ax2.set_xscale('log')
ax2.set_xticks([10, 30, 100, 300, 600])
ax2.set_xticklabels([10, 30, 100, 300, 600])
ax2.set_xlabel('Days post BMT')
ax2.set_ylabel('')
ax2.set_title('Donor fraction')
ax2.tick_params(axis='both')
ax2.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)


# Plot Propotion of Ki67+ FO cells
mean_d= Total_countsdf['frac_Ki67+_T1_AA4.1+_donor'].mean()
median_d = Total_countsdf['frac_Ki67+_T1_AA4.1+_donor'].median()
print('Mean of Ki67+ donor T1 precursor population:', mean_d)
mean_h = Total_countsdf['frac_Ki67+_T1_AA4.1+_host'].mean()
median_h = Total_countsdf['frac_Ki67+_T1_AA4.1+_host'].median()
print('Median of Ki67+ host T1 precursor population:', median_h)

y_fitd= linear_model(x_dense, *popt)
y_fith= linear_model(x_dense, *popt1)

y_pred_host = linear_model(x_dense, a=-0.0002, b=0.98)

sns.scatterplot(ax=ax3, data=Total_countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_AA4.1+_host', alpha= 0.85, s=75,  linewidth=0.5,color='#ab2239')
sns.scatterplot(ax=ax3, data=Total_countsdf, x='Age at S1K', y= 'frac_Ki67+_T1_AA4.1+_donor', alpha= 0.85, s=75,  linewidth=0.5,color='#3b9bb3')
# axes[2].axhline(y=median_h, color='blue', linestyle='solid', label= 'Fitted host curve', linewidth=0.85)
# axes[2].axhline(y=mean_d, linestyle='solid', label= 'Fitted donor curve', linewidth=0.85, color='black')
ax3.plot(x_dense, y_fitd, linestyle='solid', label= 'Fitted donor curve', linewidth=0.85, color='#3b9bb3')
ax3.plot(x_dense, y_fith, linestyle='solid', label= 'Fitted host curve', linewidth=0.85, color='#ab2239')
# ax3.plot(x_dense, y_pred_host, linestyle='dashed', label= 'Fitted host curve', linewidth=0.85, color='#ab2239')
ax3.set_yticks(np.arange(0, 1.25, 0.25))
ax3.set_ylabel('')
ax3.set_ylim(0.25, 1.05)
ax3.set_xscale('log')
#axes[2].set_yscale('log')
ax3.set_xlim(50, 800)
ax3.set_xticks([75, 150,  300,  600])
ax3.set_xticklabels([75, 150,  300,  600])
ax3.set_xlabel('Host Age (days)')
ax3.set_title('Proportion of Ki67+ cells')
ax3.tick_params(axis='both')
ax3.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
ax3.legend(labels=['Host', 'Donor'])

plt.tight_layout()
plt.savefig('Precursor_T1_new_plots.pdf', dpi=300)


# Plot residual for donor fraction
# y_fitfd= exp_func(x_fd, *popt2)
# residuals_d = Total_countsdf['fd_T1_AA4.1+'] - y_fitfd
# plt.figure(figsize=(8, 5))
# plt.scatter(Total_countsdf['Days postBMT'], residuals_d, color='purple', s=100)
# plt.axhline(0, color='black', linestyle='--', linewidth=1)
# plt.xlabel('Days post BMT', fontsize=16)
# plt.ylabel('Residuals', fontsize=16)
# plt.title('Residuals of Donor Fraction Fit', fontsize=17)
# plt.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.savefig('Residuals_FD_T1.png')
#plt.show()


