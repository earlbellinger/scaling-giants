## A seismic scaling relation for stellar age II: The red giant branch 

[Bellinger, E. P.](https://earlbellinger.com), "A seismic scaling relation for stellar age II: The red giant branch", *MNRAS Letters* submitted. 

The file `scaling_age.py` is identical to that in Appendix A of that paper. 

```
$ python3 scaling_age.py
Age: 7.7+/-1.2 [Gyr]
```

The file `scaling_giants.py` additionally provides calibrated scaling mass and radius relations as well as the scaling age relations. 

```
$ python3 scaling_giants.py \
    --nu_max 140.87 0.83 \
    --Delta_nu 12.715 0.021 \
    --Teff 5054 95 \
    --Fe_H -0.6 0.1 
Age: 7.7+/-1.2 [Gyr]
Mass: 1.005+/-0.042 [solar masses]
Radius: 4.877+/-0.075 [solar radii]
```

Some additional functionality is built into the scaling_giants routine. It optionally checks that the star's values are within the ranges of training data (default: `True`), and optionally warns the user when no applicable scaling relation is found (default: `True`). It can be loaded as a library with `from scaling_giants import scaling_giants`. Also, as is evident, it provides a command line interface: 

```
$ python3 scaling_giants.py -h
usage: scaling_giants.py [-h] [-n NU_MAX NU_MAX] [-d DELTA_NU DELTA_NU]
                         [-t TEFF TEFF] [-f FE_H FE_H]
                         [-c SUPPRESS_CHECK_BOUNDS] [-wb SUPPRESS_WARN_BOUNDS]
                         [-wc SUPPRESS_WARN_COMBO] [-s STAR_NAME]

Input data: value and uncertainty.

optional arguments:
  -h, --help            show this help message and exit
  -n NU_MAX NU_MAX, --nu_max NU_MAX NU_MAX
                        frequency at maximum power in microHertz
  -d DELTA_NU DELTA_NU, --Delta_nu DELTA_NU DELTA_NU
                        large frequency separation in microHertz
  -t TEFF TEFF, --Teff TEFF TEFF
                        effective temperature in Kelvin
  -f FE_H FE_H, --Fe_H FE_H FE_H
                        metallicity [Fe/H]

Additional options:
  -c SUPPRESS_CHECK_BOUNDS, --suppress_check_bounds SUPPRESS_CHECK_BOUNDS
                        don't enforce that inputs are within training data
                        bounds (not recommended)
  -wb SUPPRESS_WARN_BOUNDS, --suppress_warn_bounds SUPPRESS_WARN_BOUNDS
                        don't raise warning when star is rejected for being
                        out of bounds
  -wc SUPPRESS_WARN_COMBO, --suppress_warn_combo SUPPRESS_WARN_COMBO
                        don't raise warning when no applicable scaling
                        variable combo is found
  -s STAR_NAME, --star_name STAR_NAME
                        name of star
```
