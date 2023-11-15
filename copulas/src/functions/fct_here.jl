#= Function HERE
========================================================================
Take a path 

--- Input:
        path: strings to directory and/or file. The root directory is 
              the directory of the directory /src
--- Output:
        string with full desired path name.
=======================================================================#

function here(path...)
    root_dir = "C:\\Users\\ks885\\Documents\\aa_research\\Modeling\\norta_scenarios\\copulas"
    joinpath(root_dir, path...)
end

