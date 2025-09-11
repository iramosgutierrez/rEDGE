
<img src="man/figures/EDGEcalc_logo.png" width="140px" align="right"/>

# EDGEcalc

Package to calculate EDGE scores

Documentation of website under current development…

## The EDGE index

## EDGE1

Based on Isaac et al., 2007.

Following formula:

log(1+ED)+GE∗log(2)

## EDGE2

Based on Gumbs et al., 2023

2 functions: `calculate_EDGE2` calculates 1 iteration
`calculate_EDGE2_multiple` allows to calculate several (`n.iter`)
replicates, even in parallel using multiple CPU cores.
