# API

## Solar position
The following functions implement the Solar Position Algorithm, as described in []. It follows the NREL's [C implementation](https://midcdmz.nrel.gov/spa/), with the only difference of using Julia's [`datetime2julian`](https://docs.julialang.org/en/v1/stdlib/Dates/#Dates.julian2datetime) for the calculation of the Julian day. 

## References
