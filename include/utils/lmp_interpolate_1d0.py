import numpy as np
import scipy.interpolate as spinterp
from matplotlib import pyplot as plt

coeff = [-6.65303514e-09, 2.95224385e-04, -6.19561357e-01]
C = 28.5 + np.log10(3600);

time_scale = 3600
time = np.asarray([65,50.6,44,45,56,10])
T = 840
P = lambda t: T*(C+np.log10(t))
pval = P(time)

sfun = np.poly1d(coeff)
stress = 10**sfun(pval)

print(stress[::-1], time[::-1]*3600, )

plt.plot(pval, stress)
plt.show()
