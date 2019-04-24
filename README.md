BigVAR
======

Tools for modeling sparse high-dimensional multivariate time series

For a demonstration of the package's capabilities, see the recent updated [BigVAR Tutorial](http://www.wbnicholson.com/BigVAR.html) or the slightly out of date user guide available [on Arxiv](https://arxiv.org/abs/1702.07094).  The shiny app is available [here](http://ec2-54-226-7-230.compute-1.amazonaws.com:3838/BigVAR/).

Note: This package utilizes C++11, so it requires a compiler with C++11 support (which should include most modern compilers) and a later version of R (version 3.1 is the oldest that I can confirm works).

To install the development version of BigVAR, after installing the devtools package, run the following commands

```R
library(devtools)

install_github("wbnicholson/BigVAR/BigVAR")
```

The stable version is available on [cran](https://cran.r-project.org/package=BigVAR).


If you experience any bugs or have feature requests, contact me at wbn8@cornell.edu.

