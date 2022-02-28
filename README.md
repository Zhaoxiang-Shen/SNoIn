# SNoIn
A collection of Python modules for developing programs of **S**hell **No**nlocal **In**teraction with holey substrates (pronounced "snowing"). This library relies on [ShNaPr](https://github.com/david-kamensky/ShNAPr), [tIGAr](https://github.com/david-kamensky/tIGAr), and their associated dependencies, chiefly [FEniCS](https://fenicsproject.org/).

The library was written to support the following paper:
```
@article{Shell-HoleySub,
title = "Continuum modeling and isogeometric analysis of 2D materials over holey substrates",
journal = "TBD",
author = "MK,ZS,MC,ET,DK"
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
* The parallel running option `mpirun` fails when the number of element is larger than 100X100 in the circular case (i.e. [demo/cosSub.py](https://github.com/Zhaoxiang-Shen/SNoIn/demo/cosSub.py)), while the library performs stably in serial.  
* The materials cosidered in the problem are MoS2 and Si3N4. However, users can introduce different materials as they desire. 
* For new substrates, users should follow the format in [SNoIn/Substrate.py](https://github.com/Zhaoxiang-Shen/SNoIn/SNoIn/Substrate.py). Nevertheless, one should keep in mind the limits of the given tengent plane approximation, i.e. TODO. See the paper for details.
