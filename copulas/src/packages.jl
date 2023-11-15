# Load packages
const required_packages = [
    "AWSS3",
    "CSV",
    "DataFrames",
    "Dates",
    "Distributions",
    "HDF5",
    "JuliaFormatter",
    "LinearAlgebra",
    "LinearSolve",
    "Random",
    "Statistics",
    "StatsBase",
    "Tables",
    "TSFrames",
    "TimeZones"

]


# for pkg in required_packages
#     if haskey(Pkg.installed(), Symbol(pkg))
#         using $(Meta.parse(pkg))
#         println("Package $pkg was already installed and has been loaded.")
#     else
#         Pkg.add(pkg)
#         using $(Meta.parse(pkg))
#         println("Package $pkg was not installed but has been installed and loaded.")
#     end
# end

# Install package if necessary. Load all packages.
for pkg in required_packages
    if pkg ∉ keys(Pkg.installed())
        Pkg.add(pkg)
    end
    using (Meta.parse(pkg))
end
