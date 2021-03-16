# SpikeInference <img src="spike_inference_hex.png" align="right" width="150px"/>

`SpikeInference` is an `R` package for quantifying uncertainty (i.e., obtaining valid p-values and confidence intervals) for spikes estimated from calcium imaging data via an <img src="https://render.githubusercontent.com/render/math?math=\ell_0"> penalized algorithm described in Jewell and Witten (2018) and Jewell et al. (2020).

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

Chen YT, Jewell SW, Witten DM. (2021) [Quantifying uncertainty in spikes estimated from calcium imaging data](https://arxiv.org/abs/2103.07818v1
). arXiv:2103.0781 [statME].

## Bug Reports / Change Requests

If you encounter a bug or would like to make a change request, please file it as an issue [here](https://github.com/yiqunchen/SpikeInference/issues).

## Additional References
Jewell, S. and Witten, D. (2018). Exact spike train inference via l0 optimization. Ann. Appl. Stat., 12(4):2457–2482.

Jewell SW, Hocking TD, Fearnhead P, Witten DM. Fast nonconvex deconvolution of calcium imaging data. Biostatistics. 2020;21(4):709-726. doi:10.1093/biostatistics/kxy083

![](https://github.com/yiqunchen/SpikeInference/blob/main/man/figures/combined_plot_exp_7_cell_29_paper_example.png)
