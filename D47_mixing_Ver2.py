import math
import matplotlib.pyplot as plt

# Script to forward model D47 mixing using the non-linerity calculation of Defliece and Lohmann (2015).
# Runs in Python version 3. There are plotting issues in Python version 2...
# Requires matplotlib package to make plots.

# USAGE:
# Assign isotopic compositions for end member compositions. Default d13C is in PDB and d18O is in SMOW.
# 		If you want to input d18O in PDB than uncomment conversions below endmember inputs.
# Sample data should be in a .csv file with the following columns:
#		Col 1: Sample names
#		Col 2: d13C in VPDB
#		Col 3: d18O in VSMOW (uncomment conversions for VPDB)
#		Col 4: D47 values from polly ARF file
# Constants for ETF (lines 65,66), EGL (line 67), and reference gas (lines 72, 73) are from an old ARF
#		output file, you may want to verify and update for your data.
# Script outputs a 4-panel plot that helps visualize how the end members relate to your data and overlays
#		the mixing model onto your data. Sample data shown as black dots, end members shown as red dots
#		with linear mixing line, mixing model shown as blue line with orange linear trend for comparison.
#		Panel 1: d13C vs d18O cross plot
#		Panel 2: d13C vs D47
#		Panel 3: d18O vs D47
#		Panel 4: mixing ratio vs D47

# Things to consider:
# - The mixing ratio for samples is currently determined using d13C relative to assigned end members.
#		To calculate mixture using d18O, change 'd13C_mixRatios[:]' on line 175 to 'd18O_mixRatios[:]'.
# - 17O abundances are calculated using the equations in Defliece and Lohmann (2015), which follow Gonfiantini
#		et al. (1993). In other words, updated values from Brand et al. (2010) are not used.
# - If you are getting modeled mixes that look inverted relative to your data, try swapping the end member values.
# - Acid correction is set to 0, corresponding to the ARF file output. This can be changed on line 68.

# NOTE: output d13C values in VPDB, d18O values in VSMOW.

# load sample data ()
samples = np.loadtxt("exampleData.csv",delimiter=',',dtype=np.str, usecols=[0], skiprows=1)
data = np.loadtxt("exampleData.csv",delimiter=',', usecols=[1,2,3], skiprows=1)

# End-member sample compositions (to input d18O as PDB, uncomment conversions below):

# End member 1:
d13C_e1 = 2
d18O_e1 = 25
D47_e1 = 0.65

# End member 2:
d13C_e2 = -4.5
d18O_e2 = 17
D47_e2 = 0.58

# Optional VSMOW to VPDB conversions (uncomment to use):
#d18O_e1 = d18O_e1*1.03091 + 30.91
#d18O_e2 = d18O_e2*1.03091 + 30.91


# Mass ratio of end member 1:
mix = np.linspace(0,1,11)

#calc mixing ratios for samples
d13C_mixRatios=(data[:,0]-d13C_e2)/(d13C_e1-d13C_e2)
d18O_mixRatios=(data[:,1]-d18O_e2)/(d18O_e1-d18O_e2)

#Constants from polly ARF document:
ETF_int = 0.92402
ETF_slope = 1.0237
EGL_slope = -0.0005259
ac = 0


# Polly ref. gas composition, ratios, abundances, and stochastic distributions. Note VPDB to VSMOW conversion for d18O_rg
d13C_rg = -10.3
d18O_rg = -6.44*1.03091 + 30.91
R13_rg = (d13C_rg/1000 + 1) * 0.0112372
R18_rg = (d18O_rg/1000 + 1) * 0.0020052
R17_rg = (R18_rg/0.0020052)**0.5164 * 0.0003799
C12_rg = 1/(1+R13_rg)
C13_rg = R13_rg/(1+R13_rg)
O16_rg = 1/(1+R17_rg+R18_rg)
O17_rg = R17_rg/(1+R17_rg+R18_rg)
O18_rg = R18_rg/(1+R17_rg+R18_rg)
R45_rg_stoc = (C12_rg*O16_rg*O16_rg + 2*C12_rg*O16_rg*O17_rg)/(C12_rg*O16_rg*O16_rg)
R46_rg_stoc = (2*C12_rg*O16_rg*O18_rg + C12_rg*O17_rg*O17_rg + 2*C13_rg*O16_rg*O17_rg)/(C12_rg*O16_rg*O16_rg)
R47_rg_stoc = (2*C13_rg*O16_rg*O18_rg + C13_rg*O17_rg*O17_rg + 2*C12_rg*O17_rg*O18_rg)/(C12_rg*O16_rg*O16_rg)

# Get isotope ratios, abundances and stochastic distributions for endmember 1:
R13_e1 = (d13C_e1/1000 + 1) * 0.0112372
R18_e1 = (d18O_e1/1000 + 1) * 0.0020052
R17_e1 = (R18_e1/0.0020052)**0.5164 * 0.0003799
C12_e1 = 1/(1+R13_e1)
C13_e1 = R13_e1/(1+R13_e1)
O16_e1 = 1/(1+R17_e1+R18_e1)
O17_e1 = R17_e1/(1+R17_e1+R18_e1)
O18_e1 = R18_e1/(1+R17_e1+R18_e1)
R45_e1_stoc = (C12_e1*O16_e1*O16_e1 + 2*C12_e1*O16_e1*O17_e1)/(C12_e1*O16_e1*O16_e1)
R46_e1_stoc = (2*C12_e1*O16_e1*O18_e1 + C12_e1*O17_e1*O17_e1 + 2*C13_e1*O16_e1*O17_e1)/(C12_e1*O16_e1*O16_e1)
R47_e1_stoc = (2*C13_e1*O16_e1*O18_e1 + C13_e1*O17_e1*O17_e1 + 2*C12_e1*O17_e1*O18_e1)/(C12_e1*O16_e1*O16_e1)
d45_e1 = (R45_e1_stoc/R45_rg_stoc - 1)*1000
d46_e1 = (R46_e1_stoc/R46_rg_stoc - 1)*1000
D47_e1_RF = D47_e1 - ac
D47_e1_vsWG = (D47_e1_RF - ETF_int)/ETF_slope
d47_e1 = ((D47_e1_vsWG + 1000)*R47_e1_stoc - 1000*R47_rg_stoc)/(R47_rg_stoc - EGL_slope*R47_e1_stoc)

# Repeat for endmember 2:
R13_e2 = (d13C_e2/1000 + 1) * 0.0112372
R18_e2 = (d18O_e2/1000 + 1) * 0.0020052
R17_e2 = (R18_e2/0.0020052)**0.5164 * 0.0003799
C12_e2 = 1/(1+R13_e2)
C13_e2 = R13_e2/(1+R13_e2)
O16_e2 = 1/(1+R17_e2+R18_e2)
O17_e2 = R17_e2/(1+R17_e2+R18_e2)
O18_e2 = R18_e2/(1+R17_e2+R18_e2)
R45_e2_stoc = (C12_e2*O16_e2*O16_e2 + 2*C12_e2*O16_e2*O17_e2)/(C12_e2*O16_e2*O16_e2)
R46_e2_stoc = (2*C12_e2*O16_e2*O18_e2 + C12_e2*O17_e2*O17_e2 + 2*C13_e2*O16_e2*O17_e2)/(C12_e2*O16_e2*O16_e2)
R47_e2_stoc = (2*C13_e2*O16_e2*O18_e2 + C13_e2*O17_e2*O17_e2 + 2*C12_e2*O17_e2*O18_e2)/(C12_e2*O16_e2*O16_e2)
d45_e2 = (R45_e2_stoc/R45_rg_stoc - 1)*1000
d46_e2 = (R46_e2_stoc/R46_rg_stoc - 1)*1000
D47_e2_RF = D47_e2 - ac
D47_e2_vsWG = (D47_e2_RF - ETF_int)/ETF_slope
d47_e2 = ((D47_e2_vsWG + 1000)*R47_e2_stoc - 1000*R47_rg_stoc)/(R47_rg_stoc - EGL_slope*R47_e2_stoc)

# Calculate mixed compositions between end members:

d13C_mix = mix*d13C_e1 + (1-mix)*d13C_e2
d18O_mix = mix*d18O_e1 + (1-mix)*d18O_e2
d45_mix = mix*d45_e1 + (1-mix)*d45_e2
d46_mix = mix*d46_e1 + (1-mix)*d46_e2
d47_mix = mix*d47_e1 + (1-mix)*d47_e2
R13_mix = (d13C_mix/1000 + 1) * 0.0112372
R18_mix = (d18O_mix/1000 + 1) * 0.0020052
R17_mix = (R18_mix/0.0020052)**0.5164 * 0.0003799
R45_mix = (d45_mix/1000 + 1)*R45_rg_stoc
R46_mix = (d46_mix/1000 + 1)*R46_rg_stoc
R47_mix = (d47_mix/1000 + 1)*R47_rg_stoc
C12_mix = 1/(1+R13_mix)
C13_mix = R13_mix/(1+R13_mix)
O16_mix = 1/(1+R17_mix+R18_mix)
O17_mix = R17_mix/(1+R17_mix+R18_mix)
O18_mix = R18_mix/(1+R17_mix+R18_mix)
R45_mix_stoc = (C12_mix*O16_mix*O16_mix + 2*C12_mix*O16_mix*O17_mix)/(C12_mix*O16_mix*O16_mix)
R46_mix_stoc = (2*C12_mix*O16_mix*O18_mix + C12_mix*O17_mix*O17_mix + 2*C13_mix*O16_mix*O17_mix)/(C12_mix*O16_mix*O16_mix)
R47_mix_stoc = (2*C13_mix*O16_mix*O18_mix + C13_mix*O17_mix*O17_mix + 2*C12_mix*O17_mix*O18_mix)/(C12_mix*O16_mix*O16_mix)
D47_mix_vsWG = ((R47_mix/R47_mix_stoc - 1) - (R46_mix/R46_mix_stoc - 1) - (R45_mix/R45_mix_stoc -1))*1000
D47_mix_vsWGo = D47_mix_vsWG - d47_mix*EGL_slope
D47_mix = D47_mix_vsWGo*ETF_slope + ETF_int


#plot data

plt.subplot(2,2,1)
plt.title("bulk isotopes")
plt.plot(data[:,0],data[:,1], 'k.')
plt.plot([d13C_e1, d13C_e2],[d18O_e1, d18O_e2],'r-o')
plt.xlabel("d13C VPDB")
plt.ylabel("d18O VSMOW")

plt.subplot(2,2,2)
plt.title("D47 vs d13C")
plt.plot(data[:,0],data[:,2], 'k.')
plt.plot([d13C_e1, d13C_e2],[D47_e1, D47_e2],'r-o')
plt.xlabel("d13C VPDB")
plt.ylabel("D47")

plt.subplot(2,2,3)
plt.title("D47 vs d18O")
plt.plot(data[:,1],data[:,2], 'k.')
plt.plot([d18O_e1, d18O_e2],[D47_e1, D47_e2],'r-o')
plt.xlabel("d18O VSMOW")
plt.ylabel("D47")

plt.subplot(2,2,4)
plt.title("Mixing")
plt.plot(mix, D47_mix)
plt.plot([0,1],[D47_e2,D47_e1])
plt.plot(d13C_mixRatios[:],data[:,2], 'k.')
plt.xlabel("mix ratio")
plt.ylabel("D47")

plt.subplots_adjust(hspace=0.5,wspace=0.25)

plt.show()





