## A seismic scaling relation for stellar age II: The red giant branch 

[Bellinger, E. P.](https://earlbellinger.com), "A seismic scaling relation for stellar age II: The red giant branch", *MNRAS Letters* submitted. 

The file `scaling-age.py` is identical to that in Appendix A of that paper. 

```
$ python3 scaling-age.py
Age: 7.7+/-1.2 [Gyr]
```

The file `scaling-giants.py` additionally provides calibrated scaling mass and radius relations as well as the scaling age relations. 

```
$ python3 scaling-giants.py 
Age: 7.7+/-1.2 [Gyr]
Mass: 1.005+/-0.042 [solar masses]
Radius: 4.877+/-0.075 [solar radii]
```

Some additional functionality is built into the scaling-giants routine. It optionally checks that the star's values are within the ranges of training data (default: `True`), and optionally warns the user when no applicable scaling relation is found (default: `True`). 