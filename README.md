# SpikeInference <img src="spike_inference_hex.png" align="right" width="150px"/>
Quantifying uncertainty for spikes estimated from calcium imaging data

`SpikeInference` is an `R` package for quantifying uncertainty (i.e., obtaining valid p-values and confidence intervals) for spikes estimated from calcium imaging data via an <img src="https://render.githubusercontent.com/render/math?math=\ell_0"> penalized algorithm described in Jewell and Witten (2018) and Jewell et al. (2019).

[![Build Status](https://travis-ci.com/yiqunchen/SpikeInference.svg?token=quzzuXpzUN1XM57uHTXX&branch=main)](https://travis-ci.com/yiqunchen/SpikeInference)

## Installation

To download the SpikeInference package, use the code below.
``` r
# install.packages("devtools")
devtools::install_github("yiqunchen/SpikeInference")
library(SpikeInference)
```

Note that `SpikeInference` imports the package `Rcpp`. If the installation process fails on MAC OS due to issues related to `R` and `cpp` compiler tools, this [post](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/) might provide some useful information.

## Use
See [here]( https://yiqunchen.github.io/SpikeInference/) for tutorials and examples.  Please [file an issue](https://github.com/yiqunchen/SpikeInference/issues) if you have a request for a tutorial that is not currently included.

## Citation

If you use `SpikeInference` for your analysis, please cite our manuscript:

Yiqun Chen, Sean Jewell, and Daniela Witten. (2021). Uncertainty quantification for
spikes estimated from calcium imaging data. (insert arxiv link)

## Bug Reports / Change Requests

If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/yiqunchen/SpikeInference/issues).

## Additional References
Jewell, S. and Witten, D. (2018). Exact spike train inference via l0 optimization. Ann. Appl. Stat., 12(4):2457–2482.

Jewell, S. W., Hocking, T. D., Fearnhead, P., and Witten, D. M. (2019b). Fast nonconvex deconvolution of calcium imaging data. Biostatistics.

![](https://github.com/yiqunchen/SpikeInference/blob/main/man/figures/combined_plot_exp_7_cell_29_paper_example.png)
