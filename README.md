# Scenario Generation with NORTA
## Author: Hugo Santarém de Araújo


**Dec. 22, 2025**
- The code now generates scenarios as required by the "Intraday Market under 
price consistency" project. Instead of generating independent scenarios and 
rearrange them to have the past realizations fixed, the code now generates 
scenarios that fix and condition new realizations on past ones. 

--- To be edited
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
