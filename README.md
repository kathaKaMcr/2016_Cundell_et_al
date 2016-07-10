# 2016_Cundell_et_al
Custom code for a large-scale analysis phospho-proteome-dynamics at mitotic exit.
Dependencies:
- **[XPPjl](https://github.com/novakgroupoxford/XPPjl)**, a modelling platform written in julia
- **[PPSIM](https://github.com/novakgroupoxford/PPSIM)**
- **Julia-Packages:**
  - DataFrames
  - JSON
  - DataStructures
  - Optim
  - PyPlot

...make sure these are available to you.


## Useage
Used in conjunction with [PPSIM](https://github.com/novakgroupoxford/PPSIM) and [XPPjl](https://github.com/novakgroupoxford/XPPjl).

PPSIM provides a general interface for specifying the structure of the data of interest, a dynamical model to fit to the data, and for building data-processing pipelines.

It requires the presence of the following `*.yml`-files:

- `specifications.yml` the main file for setting up the data processing, including locaton of the data, regular expressions for data wrangling and a

- `model.yml`, a file specifying the mathematical model, model states and data bindings

- `stream.yml`, a placeholder file for specifying the pipleine (not yet working)

Look at each in detail to see what they contain.
PPSIM generates an object that enables easy access to attributes and setings used in the construction of pipelines.
The specifications-file allows for external resources, such as custom code to be included in the processing.
At the moment, the pipeline-feature is deprecated, `script.jl` describes the pipeline in pure julia code.
To process the data, simply modify the `specifications.yml`-file to point to the right data-source, specify an appropriate path for the folder to hold the processed data, and run `scrip.jl`
