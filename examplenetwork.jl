# Import needed modules
include("FunctionalMotifs.jl")
include("ReactivityResults.jl")

using Random,
LinearAlgebra,
CSV,
DataFrames,
.FunctionalMotifs,
.ReactivityResults

# Network setup -----------------------------------------
# Read the example interaction matrix 
df = CSV.read("ex_intmat.csv", DataFrame; header = false)
interactions = Matrix{Int64}(df)   

# Initializing ------------------------------------------
# Number of iterations (determines the number of different parametrizations or communities)
n = 10^4;

# Find and analyze motifs -------------------------------
# Note: Version right now checks for all motifs. Future updates will include motifs to be picked individually
results = FunctionalMotifs.func_react_motifs(n, interactions)

# Create histograms -------------------------------------
ReactivityResults.hist_probdens2x2([
    (data = results.props[:pp], title = "Predator-Prey"),
    (data = results.props[:ac], title = "Apparent Competition"),
    (data = results.props[:ec], title = "Exploitative Competition"),
    (data = results.props[:ch3], title = "Tri-Trophic Chain"),
])