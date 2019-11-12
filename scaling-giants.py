import numpy as np 
from uncertainties import ufloat
from uncertainties.unumpy import isnan 
from warnings import warn 

# Enter some data for an example red giant whose age we want to estimate 
# Use ufloat(0,0) for any measurement that is not available
nu_max   = ufloat( 140.87,   0.83 ) # muHz, second argument is uncertainty 
Delta_nu = ufloat(  12.715,  0.021) # muHz
Teff     = ufloat(5054,     95    ) # K
Fe_H     = ufloat(  -0.6,    0.1  ) # dex

refs = [4.3, 1.33, 6.6] # age [Gyr], mass [solar units], radius [solar units]

# Calibrated exponents from Table 2
# P = [   alpha,    beta,  gamma,  delta]
P_age_full = np.array([
      [- 9.760 , 13.08  , -6.931, 0.4894],  #  1
      [- 7.778 , 10.77  , -11.05,      0],  #  2
      [-12.19  , 15.86  ,      0,  1.027],  #  3
      [       0,  1.396 , -22.32, -1.046],  #  4
      [  1.084 ,       0, -23.28, -1.165],  #  5
      [- 8.837 , 11.73  ,      0,      0],  #  6
      [       0,  0.9727, -14.64,      0],  #  7
      [  0.6424,       0, -13.82,      0],  #  8
      [       0,       0,      0,      0],  #  9
      [       0,       0,      0,      0]]) # 10
sigma_sys_age = np.array([0.25, 0.32, 0.34, 0.82, 0.92, 
                          0.86, 1.2,  1.3,  0,    0])

# Calibrated exponents from Table 3
# P = [  alpha,    beta, gamma,   delta]
P_mass_full = np.array([
      [ 2.901 , -3.876 , 1.621,       0],  #  1
      [ 2.901 , -3.876 , 1.621,       0],  #  2
      [ 3.546 , -4.619 ,     0, -0.1457],  #  3
      [      0, -0.3845, 5.740,  0.4290],  #  4
      [-0.2976,       0, 5.935,  0.4594],  #  5
      [  3.056, -4.015 ,     0,       0],  #  6
      [      0,       0,     0,       0],  #  7
      [      0,       0,     0,       0],  #  8
      [      0,       0,     0,       0],  #  9
      [      0,       0,     0,       0]]) # 10
sigma_sys_M = np.array([0.023, 0.023, 0.023, 0.10, 0.11,  
                        0.046, 0.15,  0,     0,    0])

# Calibrated exponents from Table 4
# P = [  alpha,    beta,  gamma,  delta]
P_radius_full = np.array([
      [ 0.9570, -1.955 , 0.6288,      0],  #  1
      [ 0.9570, -1.955 , 0.6288,      0],  #  2
      [ 1.008 , -1.999 ,      0,      0],  #  3
      [      0, -0.8048, 2.062 , 0.1378],  #  4
      [-0.6593,       0, 2.953 , 0.2283],  #  5
      [ 1.008 , -1.999 ,      0,      0],  #  6
      [      0, -0.7362, 0.8088,      0],  #  7
      [-0.5591,       0, 0.7857,      0],  #  8
      [      0, -0.7038,      0,      0],  #  9
      [-0.5353,       0,      0,      0]]) # 10
sigma_sys_R = np.array([0.037, 0.037, 0.075, 0.16, 0.25, 
                        0.075, 0.24,  0.38,  0.24, 0.36])

# now stack all these tables together 
P         = np.array([P_age_full,    P_mass_full, P_radius_full])
sigma_sys = np.array([sigma_sys_age, sigma_sys_M, sigma_sys_R])

# Apply the scaling relation. 
# Returns age, mass, and radius in Gyr, solar masses, and solar radii. 
def scaling_giants(nu_max, Delta_nu, Teff, Fe_H, 
        nu_max_ref=104.5, Delta_nu_ref=9.25, Teff_ref=4790, 
        check_bounds=True, warn_bounds=True, warn_combo=True, star_name=''):
    
    result = [np.nan, np.nan, np.nan] # nannannannannan batman! 
    
    # check that the data are in bounds 
    if check_bounds: 
        if (nu_max   != 0 and nu_max   <   27.9  or nu_max   >   255.6  or
            Delta_nu != 0 and Delta_nu <    3.73 or Delta_nu >    17.90 or 
            Teff     != 0 and Teff     < 4520    or Teff     > 5120     or 
            Fe_H     != 0 and Fe_H     <   -1.55 or Fe_H     >     0.50): 
            if warn_bounds:
                warn("Input data out of range of training data. " + star_name)
            return result 
    
    # Determine which row of the table to use by checking which entries are 0
    star = np.array([nu_max!=0, Delta_nu!=0, Teff!=0, Fe_H!=0])
    found = False
    for combo_idx in range(len(P_age_full)):
        exponents = P_age_full[combo_idx,]
        if not np.any(exponents): # relation 9 or 10
            exponents = P_radius_full[combo_idx,]
        found = np.array_equal(np.nonzero(exponents)[0], np.nonzero(star)[0])
        if found: 
            break 
    
    if not found: # No applicable scaling relation 
        if warn_combo:
            warn("No applicable relation for input combination. " + star_name)
        return result 
    
    # Equation 1, plus the systematic error of the corresponding relation 
    for var_idx, exponents in enumerate(P):
        if not np.any(exponents[combo_idx,]):
            continue 
        alpha, beta, gamma, delta = exponents[combo_idx,]
        sigma = sigma_sys[var_idx, combo_idx]
        result[var_idx] = (
            (nu_max   /   nu_max_ref) ** alpha  * 
            (Delta_nu / Delta_nu_ref) ** beta   * 
            (Teff     /     Teff_ref) ** gamma  * refs[var_idx] * 
            (np.e**Fe_H             ) ** delta) + ufloat(0, sigma)
    
    return result 

age, mass, radius = scaling_giants(nu_max, Delta_nu, Teff, Fe_H)

if not isnan(age):
    print('Age:', '{:.2u}'.format(age), '[Gyr]')

if not isnan(mass):
    print('Mass:', '{:.2u}'.format(mass), '[solar masses]')

if not isnan(radius):
    print('Radius:', '{:.2u}'.format(radius), '[solar radii]')