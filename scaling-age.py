import numpy as np 
from uncertainties import ufloat

# Enter some data for an example red giant whose age we want to estimate 
# Use ufloat(0,0) for any measurement that is not available
nu_max   = ufloat( 140.87,   0.83 ) # muHz, second argument is uncertainty 
Delta_nu = ufloat(  12.715,  0.021) # muHz
Teff     = ufloat(5054,     95    ) # K
Fe_H     = ufloat(  -0.6,    0.1  ) # dex

# Calibrated exponents from Table 2
# P       =       [   alpha,    beta,  gamma,  delta]
P_age = np.array([[- 9.760 , 13.08  , -6.931, 0.4894],  # 1
                  [- 7.778 , 10.77  , -11.05,      0],  # 2
                  [-12.19  , 15.86  ,      0,  1.027],  # 3
                  [       0,  1.396 , -22.32, -1.046],  # 4
                  [  1.084 ,       0, -23.28, -1.165],  # 5
                  [- 8.837 , 11.73  ,      0,      0],  # 6
                  [       0,  0.9727, -14.64,      0],  # 7
                  [  0.6424,       0, -13.82,      0]]) # 8
sigma_sys = np.array([0.25, 0.32, 0.34, 0.82, 0.92, 0.86, 1.2, 1.3])

# Apply the scaling relation 
def scaling_age(nu_max, Delta_nu, Teff, Fe_H, 
        nu_max_ref=104.5, Delta_nu_ref=9.25, Teff_ref=4790, age_ref=4.3):
    # Determine which row of the table to use by checking which entries are 0
    star = np.array([nu_max!=0, Delta_nu!=0, Teff!=0, Fe_H!=0])
    found = False
    for idx, exponents in enumerate(P_age):
        found = np.array_equal(np.nonzero(exponents)[0], np.nonzero(star)[0])
        if found: 
            break 
    
    if not found: # No applicable scaling relation 
        return np.nan 
    
    # Equation 1, plus the systematic error of the corresponding relation 
    alpha, beta, gamma, delta = exponents
    return ((nu_max   /   nu_max_ref) ** alpha  * 
            (Delta_nu / Delta_nu_ref) ** beta   * 
            (Teff     /     Teff_ref) ** gamma  * age_ref * 
            (np.e**Fe_H             ) ** delta) + ufloat(0, sigma_sys[idx])

age = scaling_age(nu_max, Delta_nu, Teff, Fe_H)
print('Age:', '{:.2u}'.format(age), '[Gyr]')