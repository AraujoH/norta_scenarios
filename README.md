# Scenario Generation with NORTA
## Author: Hugo Santarém de Araújo

## Overview

This project implements scenario generation for power system operations using the NORTA (NORmal To Anything) method with Gaussian copulas. It generates spatially and temporally correlated stochastic scenarios for electricity load, solar power, and wind power based on ERCOT historical data and probabilistic forecasts.

The code supports both day-ahead and intraday market scenario generation. For intraday markets, scenarios are conditioned on realized values from earlier hours, enabling analysis of multi-stage stochastic optimization problems under uncertainty.

**Key Features:**
- Gaussian copula-based scenario generation preserving spatial and temporal correlations
- Landing probability estimation from probabilistic forecasts
- Intraday market scenario generation with proper conditioning on past realizations
- Support for ERCOT load, solar, and wind data
- Configurable scenario dimensions (length, number of scenarios, sheets, iterations)

## Changelog

### 2026-02-11
**Bug Fixes:**
- Fixed critical bug in IDM scenario conditioning where scenarios were incorrectly conditioned on their own iteration instead of the specified base iteration. Added `iteration_index` parameter to `generate_probability_IDM_scenarios_cube!` to properly control which base iteration is used for conditioning past values.

**Performance Improvements:**
- Consolidated output files from 180+ intermediate files to 18 consolidated files by implementing 5D/4D array storage
- Made save path parameters optional in `generate_probability_IDM_scenarios_cube!` and `convert_land_prob_cube_to_data` to eliminate unnecessary intermediate file generation
- Pre-allocated large arrays in main script to hold all iterations and hours in memory before serialization

**Code Quality:**
- Improved parameter naming clarity: renamed `number_of_iterations` to `number_of_iterations_IDM` in IDM generation function
- Removed confusing "normal" terminology from parameter names (e.g., `normal_scenario_iteration` → `iteration_index`)
- Updated file naming conventions to use `iter_$(iteration)` instead of `normal_$(iter)`
- Added `iteration` parameter to `write_scenarios_to_file` for better file tracking

### 2025-12-22
- The code now generates scenarios as required by the "Intraday Market under price consistency" project. Instead of generating independent scenarios and rearranging them to have the past realizations fixed, the code now generates scenarios that fix and condition new realizations on past ones.


This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> norta_scenarios

It is authored by Hugo Santarem de Araujo.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "norta_scenarios"
```
which auto-activate the project and enable local path handling from DrWatson.
