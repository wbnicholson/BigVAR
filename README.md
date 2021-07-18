BigVAR
======



Tools for modeling sparse high-dimensional multivariate time series

# R Package

For a demonstration of the package's capabilities, see the recently updated [BigVAR Tutorial](http://www.wbnicholson.com/BigVAR.html), the [Shiny App](http:/bigvar.ddns.net:3838/BigVAR/), or the slightly out of date user guide available [on Arxiv](https://arxiv.org/abs/1702.07094).

Note: This package utilizes C++11, so it requires a compiler with C++11 support (which should include most modern compilers) and a later version of R (version 3.1 is the oldest that I can confirm works).

To install the development version of BigVAR, after installing the devtools package, run the following commands

```R
library(devtools)

install_github("wbnicholson/BigVAR/BigVAR")
```

The stable version is available on [cran](https://cran.r-project.org/package=BigVAR).


If you experience any bugs or have feature requests, contact me at wbn8@cornell.edu.

# Python Package

A minimalist Python implementation (partially inspired by this [abandoned effort](https://github.com/josh-alley/BigVARPython)) has been released.  Currently, it only has the capability to fit VARs Basic or Elastic Net penalty structures.  Feel free to suggest other functionality or submit pull requests.

## Installation

In order to install the Python implementation, clone the repository, navigate to the python directory and run 

```bash
pip install -e .
```

### Usage

An example script is below

```python

import numpy as np
from BigVAR.BigVARSupportFunctions import MultVARSim, CreateCoefMat
from BigVAR.BigVARClass import BigVAR,rolling_validate

# example coefficient matrix
B1=np.array([[.4,-.02,.01],[-.02,.3,.02],[.01,.04,0.3]])
B2=np.array([[.2,0,0],[0,.3,0],[0,0,0.13]])
B=np.concatenate((B1,B2),axis=1)
B=np.concatenate((B,np.zeros((k,2*k))),axis=1)
k=3;p=4
A=CreateCoefMat(B,p,k)
Y=MultVARSim(A,p,k,0.01*np.identity(3),500)
VARX={}

# construct BigVAR object:
# Arguments:
# Y T x k multivariate time series
# p: lag order
# penalty structure (only Basic and BasicEN supported)
# granularity (depth of grid and number of gridpoints)
# T1: Start of rolling validation
# T2: End of rolling validation
# alpha: elastic net alpha candidate
# VARX: VARX specifications as dict with keys k (number of endogenous series), s (lag order of exogenous series)

mod=BigVAR(Y,p,"Basic",[50,10],50,80,alpha=0.4,VARX=VARX)

res=rolling_validate(mod)

# coefficient matrix
res.B

# out of sample MSFE
res.oos_msfe

#optimal lambda
res.opt_lambda
```



