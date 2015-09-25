BigVAR
======

Dimension Reduction Methods for Multivariate Time Series


For a demonstration of the package's capabilities, see the shiny app: http://ec2-107-23-247-210.compute-1.amazonaws.com:3838 or the user guide available [on my website](http://www.wbnicholson.com/BigVAR.pdf).

Note: This package utilizes C++11, so it requires a compiler with C++11 support (which should include most modern compilers) and a later version of R (version 3.1 is the oldest that I can confirm works).

To install BigVAR, after installing the devtools package, run the following commands

```R
if(!require("zoo"))
{install.packages("zoo")}

library(devtools)

install_github("wbnicholson/BigVAR/BigVAR")
```

The package source is available [here](http://www.wbnicholson.com/BigVAR_1.0.tar.gz).

