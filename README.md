# SNoIn
A collection of Python modules for developing programs of **S**hell **No**nlocal **In**teraction with holey substrates (pronounced "snowing"). This library relies on [ShNaPr](https://github.com/david-kamensky/ShNAPr), [tIGAr](https://github.com/david-kamensky/tIGAr), and their associated dependencies, chiefly [FEniCS](https://fenicsproject.org/).

The library was written to support the following paper on JMPS:
```
@article{Shell-HoleySub,
title = {Atomistically-informed continuum modeling and isogeometric analysis of 2D materials over holey substrates},
journal = {Journal of the Mechanics and Physics of Solids},
volume = {170},
pages = {105100},
year = {2023},
issn = {0022-5096},
doi = {https://doi.org/10.1016/j.jmps.2022.105100},
url = {https://www.sciencedirect.com/science/article/pii/S0022509622002770},
author = {Moon-ki Choi and Marco Pasetto and Zhaoxiang Shen and Ellad B. Tadmor and David Kamensky},
}

```

# Installation 
It is recommended to use this library in conjunction with the advanced form compiler TSFC. For its installation, see Kamensky's [answer](https://fenicsproject.discourse.group/t/quadrature-representation/2025/2).
```
pip3 install git+https://github.com/blechta/tsfc.git@2018.1.0
pip3 install git+https://github.com/blechta/COFFEE.git@2018.1.0
pip3 install git+https://github.com/blechta/FInAT.git@2018.1.0
pip3 install singledispatch networkx pulp
```

# Remarks
* The materials cosidered in the problem are MoS2 and Si3N4. However, users can introduce different materials as they desire. 
* For new substrates, users should follow the format in [SNoIn/Substrate.py](https://github.com/Zhaoxiang-Shen/SNoIn/SNoIn/Substrate.py). Nevertheless, one should keep in mind the limits of the given tengent plane approximation, see Section 3.2.2 in the paper.
