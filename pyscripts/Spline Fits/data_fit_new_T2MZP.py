import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.preprocessing import MinMaxScaler, StandardScaler

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Read specific sheets by name
Total_countsdf = pd.read_csv("/Users/apoorvasingh/Desktop/MZB-new-analysis/data/new_data_MZB.csv") # counts for the host compartment MZ B cell population

#Curve fitting function
def exponential_model(x, a, b):
    return np.exp(b* (x-40)+a)
def linear_model(x, a):
    return a*(x-10) 
def line_func(x):
    return np.exp(0.0035*(x-10)) - 0.9865
def logistic_model(x):
    # return 0.7 / (1 + np.exp(- 0.1* (x-10)))
    return 0.73 / (1 + np.exp(- 0.08* (x-59)))

    # return 0.9 / (1 + (60/(x-13)))

def lin_model(x):
    return np.exp(0.00045*(x-40) + 16.98)
def gaussian_model(x, a, b):
    return  np.exp(-b*(x - 40)) + a
# def gau_function(x):
#     return  np.exp(-0.0465*(x-40)) + 0.095
def gau_function(x):
    return  (0.4*np.exp(-0.048*(x-59))) +0.095
def exp_func(x, a, b):
    return   a*(1- np.exp(-b * (x-10)))
def nfd_func(x):
    return 0.8*(1 - np.exp(-0.0003*(x-10)**2)) 

# drop the last two rows to get a better fit
Total_countsdf1 = Total_countsdf.drop(Total_countsdf.tail(2).index)
print(Total_countsdf1.head())
x_data = np.array(Total_countsdf['Age at S1K'])
x_data1 = np.array(Total_countsdf1['Age at S1K'])
x_fd = np.array(Total_countsdf['Days postBMT'])
y_precursor = np.array(Total_countsdf['total_T2_MZP'])
y_precursor1 = np.array(Total_countsdf1['total_T2_MZP'])
y_fd = np.array(Total_countsdf['fd_FM'])
y_donor= np.array(Total_countsdf['frac_Ki67+_T2_MZP_donor'])
y_host= np.array(Total_countsdf['frac_Ki67+_T2_MZP_host'])

# Initial guess for the parameters
initial_guess_p = [0.5, 0.5]
#initial_guess_d= [ 0.02, 0.25]

# Scale the input data
scaler = MinMaxScaler()
# scaler1 = StandardScaler()
x_fd_scaled = scaler.fit_transform(x_fd.reshape(-1, 1)).flatten()
# x_data_scaled = scaler.fit_transform(x_data.reshape(-1, 1)).flatten()
# x_data_scaled1 = scaler1.fit_transform(x_data.reshape(-1, 1)).flatten()

#scale presursor data
# y_precursor_scaled = scaler.fit_transform(y_precursor.reshape(-1, 1)).flatten()
y_precursor_log = np.log(y_precursor)

# # Perform curve fitting
# popt, pcov = curve_fit(linear_model, x_data, y_precursor_log, p0=initial_guess_p, maxfev=10000)
# popt3, pcov3 = curve_fit(exponential_model, x_data, y_precursor_log, p0=initial_guess_p, maxfev=10000)
popt1, pcov1 = curve_fit(gaussian_model, x_data, y_donor, maxfev=10000)
# popt3, pcov3 = curve_fit(logistic_model, x_data, y_donor, maxfev=10000)

popt2, pcov2 = curve_fit(exp_func, x_fd_scaled, y_fd, maxfev=10000)
# # popt, pcov = curve_fit(exponential_model, x_data_scaled, y_host, p0=initial_guess, maxfev=10000)

# Print the fitted parameters
# print("Fitted parameters for Precursor population:", popt)
print("Fitted parameters for Donor Fraction:", popt1)
# print("Fitted parameters for Donor Fraction:", popt3)
print("Fitted parameters for FD:", popt2)

#unsacle the data for log plot
#x_fit_scaled = scaler.transform(x_data.reshape(-1, 1)).flatten()
x_dense = np.linspace(59, 750, 1000)
x_densefd = np.linspace(10, 650, 1000)

# x_fit_scaled = scaler.transform(x_data.reshape(-1, 1)).flatten()
# y_fit_log = exponential_model(x_data, *popt3)
# y_fit1_log = linear_model(x_data, *popt)
y_fit_fd = exp_func(x_densefd, *popt2)
# y_fit1_scaled = np.exp(y_fit1_log)

# Inverse transform the scaled y values to get them back to the original scale
#y_fitp = scaler.inverse_transform(y_fit_scaled.reshape(-1, 1)).flatten()

# Create a 3x2 grid of subplots
fig, axes = plt.subplots(2, 2, figsize=(15, 10))


# Plot precursor cell counts
mean= Total_countsdf['total_T2_MZP'].mean()
median = Total_countsdf['total_T2_MZP'].median()
print('Mean of T2 MZP precursor population:', mean)
print('Median of T2 MZP precursor population:', median)

fig, axes = plt.subplots(1, 3, figsize=(23, 5))
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]
sns.scatterplot(ax=axes[0], data=Total_countsdf, x='Age at S1K', y='total_T2_MZP', s=75, alpha=0.85, color='black')
ax1.axhline(y=mean, linestyle='solid', linewidth=1.5, color='black')
ax1.set_yscale('log')  # sets the y-axis to a logarithmic scale
ax1.set_ylim(1e5, 1e7)
ax1.set_ylabel('')
ax1.set_xscale('log')
ax1.set_xlim(50, 800)
ax1.set_xticks([75, 150,  300, 600])
ax1.set_xticklabels([75, 150,  300, 600])
ax1.set_xlabel('Host Age (days)')
ax1.set_title('Total T2 MZP population size')
ax1.tick_params(axis='both') 
ax1.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)



# Plot donor fraction FO compartment
sns.scatterplot(ax=axes[1], data=Total_countsdf, x='Days postBMT', y='fd_T2_MZP', alpha=0.85, s=75, color='black')
# y_fitd= exp_func(x_densefd, *popt2)
# y_fitfd1 = linear_model(x_fd_scaled, *popt2)
# y_fitfd1 = logistic_model(x_densefd, *popt3)
# y_fitfd1 = logistic_model(x_densefd)
y_fitfd2 = line_func(x_densefd)
y_fitd1 = nfd_func(x_densefd)
ax2.plot(x_densefd, y_fitd1, color='black', label='Fitted curve')
# ax2.plot(x_densefd, y_fitfd1, color='black', label='Fitted curve')

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
mean_h = Total_countsdf['frac_Ki67+_T2_MZP_host'].mean()
median_h = Total_countsdf['frac_Ki67+_T2_MZP_host'].median()
print('Mean of Ki67+ host T2 MZP precursor population:', mean_h)
print('Median of Ki67+ host T2 MZP precursor population:', median_h)
y_fitd = gaussian_model(x_dense,*popt1)
y_fitd1 = gau_function(x_dense)
# y_fitd1 = gau_function(x_data)

sns.scatterplot(ax=axes[2], data=Total_countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_MZP_host', alpha= 0.85, s=75, color='#ab2239')
sns.scatterplot(ax=axes[2], data=Total_countsdf, x='Age at S1K', y= 'frac_Ki67+_T2_MZP_donor', alpha= 0.85, s=75, color='#3b9bb3')
ax3.plot(x_dense, y_fitd1, color='#3b9bb3',label='Fitted donor curve', linewidth=1.5)
ax3.axhline(y=median_h, color='#ab2239', label= 'Fitted host curve', linewidth=1.5)
ax3.set_yticks(np.arange(0, 1.25, 0.25))
ax3.set_ylabel('')
ax3.set_xscale('log')
ax3.set_xlim(50, 800)
ax3.set_xticks([75, 150,  300,  600])
ax3.set_xticklabels([75, 150,  300,  600])
ax3.set_xlabel('Host Age (days)')
ax3.set_title('Proportion of Ki67+ cells')
ax3.tick_params(axis='both', labelsize=15) 
ax3.grid(True, which='both', ls="solid", linewidth=0.5, alpha=0.3)
ax3.legend(labels=['Host', 'Donor'])



plt.tight_layout()
plt.savefig('Precursor_T2MZP_new_plots.pdf', dpi=300)
#plt.show()'''
                                       

                                       


