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
# Number of species
N = 15;

# Find and analyze motifs -------------------------------
# Note: Version right now checks for all motifs. Future updates will include motifs to be picked individually
results = FunctionalMotifs.func_react_motifs(n, interactions);

# Repeat process for topologies made with the niche model
results_niche = FunctionalMotifs.react_motifs_niche(n, N);

# Clean data --------------------------------------------
# Necessary step because the niche model sometimes produces topologies where we don't find all motifs. 
# Empty values have to be removed for the plot
cl_prop_pp = filter(!=(-999), results_niche.props[:pp])
cl_prop_ac = filter(!=(-999), results_niche.props[:ac])
cl_prop_ec = filter(!=(-999), results_niche.props[:ec])
cl_prop_ch3 = filter(!=(-999), results_niche.props[:ch3])

# Create histograms -------------------------------------
ReactivityResults.hist_probdens2x2_comp([
    (data1 = cl_prop_pp, data2 = results.props[:pp], title = "Predator-Prey"),
    (data1 = cl_prop_ac, data2 = results.props[:ac], title = "Apparent Competition"),
    (data1 = cl_prop_ec, data2 = results.props[:ec], title = "Exploitative Competition"),
    (data1 = cl_prop_ch3, data2 = results.props[:ch3], title = "Tri-Trophic Chain")
],
maxval = 6.5,
outpath = "comparison.png",
leg_labels = ("Different Topologies", "Example Network"))