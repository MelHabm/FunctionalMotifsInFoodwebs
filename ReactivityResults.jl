module ReactivityResults

using Plots,
Plots.PlotMeasures,
StatsPlots

export hist_probdens, hist_probdens2x2, hist_probdens_comp, hist_probdens2x2_comp, get_maxmotif

"""
Author:       Melanie Habermann
Date:         04.10.2024
Last changes: 29.07.25
"""

""" ------------------------------------------
hist_probdens()

Creates a single histogram that plots the probability density against the 
proportion of system reactivity.

Input:
data      : data set that should be used for the histogram

Optional: 
titlename : a string with the title for the plot (default is no title)
n_bins    : integer that determines the number of bins in the histogram (default is 30)
maxval    : max. value to be used for the y-axis (only use if you want to set that yourself)
hist_col  : determines the colour of the histogram (default is orange)
"""
function hist_probdens(data; 
                        titlename::String = "", 
                        n_bins::Int64 = 30, 
                        maxval::Number = 10.0, 
                        hist_col = :orange)

    # Set a range for bins
    b_range = range(0, 1.0, length = n_bins)

    histogram(data, 
        bins = b_range, 
        normalize = :pdf, 
        xlabel = "Proportion of System Reactivity", 
        ylabel = "Probability Density",
        title = titlename, 
        legend = false, 
        ylim = (0.0, maxval), 
        titlefont = font(18),
        ytickfont = font(20),
        guidefont = font(20),
        alpha = 0.5, 
        color = hist_col,
        framestyle = :box,
        xtickfontcolor = :white,
        grid = false)

end # end function

"""
hist_probdens_comp()

Creates a single plot with two histograms that plots the probability density against the 
proportion of system reactivity and compares the results of two sets of data.

Input:
datasets  : Tuple including data and label that should be displayed in the legend

Optional: 
titlename : a string with the title for the plot (default is no title)
n_bins    : sets the number of bins in the histogram (default is 30)
maxval    : max. value to be used for the y-axis (only use if you want to set that yourself)

Example:
hist_probdens_comp([
    (data = prob_pp_ex, label = "Example Network"), 
    (data = prob_pp_dt, label = "Different Topologies")
])
"""
function hist_probdens_comp(datasets; 
                            titlename::String = "", 
                            n_bins::Int64 = 30, 
                            maxval::Number = 10.0)
    
    # safety first
    if length(datasets) != 2
        error("hist_probdens_comp() expects exactly two datasets.")
    end

    # Set range for bins
    b_range = range(0, 1.0, length = n_bins)

    # Set default colours
    colors = [:orange, :blue]

     # Start with empty plot
    plt = plot() 

    # Fill with the histograms that should be displayed
    for (k, item) in enumerate(datasets)
        # Use default color if not provided
        hcol = get(item, :color, colors[k])

        histogram!(plt, item.data;
            bins = b_range,
            normalize = :pdf,
            alpha = 0.5,
            color = hcol,
            label = item.label)
    end

        # Add plot settings
    plot!(plt;
        title = titlename,
        xlabel = "Proportion of System Reactivity",
        ylabel = "Probability Density",
        titlefont = font(18),
        ytickfont = font(20),
        guidefont = font(20),
        legend = true,
        ylim = (0.0, maxval),
        framestyle = :box,
        xtickfontcolor = :white,
        grid = false)

end # end function

"""
hist_probdens2x2()

Creates a plot with 2x2 subplots full of histograms that plot the probability density against the 
proportion of system reactivity.

Input:
datasets : motif datasets to plot (proportion of system reactivity)

Optional:
n_bins   : Int, number of bins in the histogram (default 30)
maxval   : Int or Float, max value for the y-axis (default 10.0)
outpath  : String, path or name of the plot for save file (default hist_grid.png)
"""
function hist_probdens2x2(datasets;
                            n_bins::Int64 = 30,
                            maxval::Number = 10.0,
                            outpath::String = "hist_grid.png")

    # safety first
    if length(datasets) != 4
        error("hist_probdens_grid() expects exactly four datasets.")
    end

    # Bin range for all plots
    b_range = range(0, 1.0, length = n_bins)

    # Create 2x2 layout
    plt = plot(layout = (2, 2), size = (1080, 720),
               left_margin = 10mm, right_margin = 10mm,
               top_margin = 1mm, bottom_margin = 3mm)

    # Plot each subplot
    for (k, item) in enumerate(datasets)
        histogram!(plt, item.data;
            bins = b_range,
            normalize = :pdf,
            xlabel = (k > 2) ? "Proportion of System Reactivity" : "",
            ylabel = (k % 2 == 1) ? "Probability Density" : "",
            title = item.title,
            titlefont = font(18),
            legend = false,
            ylim = (0.0, maxval),
            subplot = k,
            alpha = 0.5,
            color = get(item, :color, :blue),
            framestyle = :box,
            xtickfont = (k > 2) ? font(20) : nothing,
            ytickfont = (k % 2 == 1) ? font(20) : nothing,
            xtickfontcolor = (k <= 2) ? :white : :black,
            ytickfontcolor = (k % 2 == 0) ? :white : :black,
            guidefont = font(20),
            grid = true,
            bottom_margin = (k > 2) ? 3mm : 1mm)
    end

    # Save to file
    savefig(plt, outpath)

    # Show plot
    display(plt)

end # end function

""" ------------------------------------------
hist_probdens2x2_comp()

Creates a plot with 2x2 subplots full of histograms that plot the probability density against the 
proportion of system reactivity and compares therein two sets of data.

Input: 
panels  : NamedTupels with (data1, data2, label) that should be used in the plot

Optional:
n_bins  : number of bins in the histogram (default 30)
maxval  : max. value for the y-axis (default 10)
colors  : Tuple with the colours that should be used for the comparison of the histograms - should be chosen colourblind friendly (default orange and blue)
outpath : string with the path or name as which the plot should be saved (default comparison_2x2.png)
"""
function hist_probdens2x2_comp(panels;
                                n_bins::Int64 = 30,
                                maxval::Number = 10.0,
                                colors::Tuple = (:blue, :orange),
                                outpath::String = "comparison_2x2.png",
                                leg_labels = ("Set 1", "Set 2"))

    # Safety first
    if length(panels) != 4
        error("Function expects exactly four panels (NamedTuples with data1, data2, and title).")
    end

    # Set bin range
    b_range = range(0, 1.0, length = n_bins)

    # Create layout
    plt = plot(layout = (2, 2), size = (1080, 720),
               left_margin = 10mm, right_margin = 10mm,
               top_margin = 1mm, bottom_margin = 3mm)

    for (i, panel) in enumerate(panels)
        # Extract data
        d1 = panel.data1
        d2 = panel.data2
        ttl = panel.title

        # Plot both histograms into subplot i
        histogram!(plt, d1;
            bins = b_range,
            normalize = :pdf,
            alpha = 0.5,
            color = colors[1],
            label = leg_labels[1],
            subplot = i)

        histogram!(plt, d2;
            bins = b_range,
            normalize = :pdf,
            alpha = 0.5,
            color = colors[2],
            label = leg_labels[2],
            subplot = i)

        # Styling each subplot
        plot!(plt;
            subplot = i,
            title = ttl,
            xlabel = (i > 2 ? "Proportion of System Reactivity" : ""),
            ylabel = (i == 1 || i == 3 ? "Probability Density" : ""),
            titlefont = font(18),
            legend = false,
            ylim = (0.0, maxval),
            guidefont = font(20),
            ytickfont = (i == 2 || i == 4 ? font(1) : font(20)),
            ytickfontcolor = (i == 2 || i == 4 ? :white : :black),
            xtickfont = (i > 2 ? font(20) : font(1)),
            xtickfontcolor = (i > 2 ? :black : :white),
            framestyle = :box,
            grid = true)
    end

    # Save and display
    savefig(plt, outpath)
    display(plt)

end # end function

"""
get_maxmotif()

Use for motifs of the same size!

The function combines the results for different motifs of the same size (i.e., motifs with the same number of nodes) and 
determines which of these motifs has the highest reactivity.

Input:
vectors... : vectors with the reactivity values of the motifs to analyze

Output:
vals       : reactivity of the most reactive motifs
idxs       : indices of the corresponding motifs
"""
function get_maxmotif(vectors::AbstractVector...)
    # Combine the vectors column-wise into a matrix.
    # Each row now has the reactivity of the most reactive instance of each motif for a single parametrization
    A = hcat(vectors...)

    # For each row, find the max and the index. The index gives information about the kind of motif that is the one with
    # the highest reactivity of the bunch
    results = map(eachrow(A)) do row
        findmax(row)
    end
    vals = getindex.(results, 1)
    idxs = getindex.(results, 2)

    return vals, idxs
end

"""
count_mot_appearance()

Counts how often a motif appears as the most reactive one.

Input:
idxs : Vector{Int} with the indices of the motifs to count
n    : number of motifs
"""
function count_mot_appearance(idxs::Vector{Int}, n_motifs::Int)
    return [sum(idxs .== i) for i in 1:n_motifs]
end

"""
box_motifs()

Creates a boxplot that compares the most reactive motif of a certain size (for example two and three node motifs) to the 
reactivity of the entire system.

Input:
vectors : vectors to plot

Optional:
ticklab : ticklabels in a Vector{String} (default: strings in the form "i nodes" where i is 1:number of vectors)
outpath : path or name to save the plot as (default: maxz_boxplot.png)

Example:
box_motifs(vals1, vals2, vals3, ticklabels = ["3-node", "4-node", "System"])
"""
function box_motifs(vectors...; 
                    ticklab::Vector{String} = ["$i Nodes" for i in 1:length(vectors)],
                    outpath::String = "maxz_boxplot.png")

    # Determine how many datasets will be compared in the boxplot
    n = length(vectors)

    # Cretae plot values
    plot_vals = [vectors...]

    if length(ticklab) != n
        error("Length of 'ticklab' must match number of input vectors.")
    end

    plt = boxplot(plot_vals,
        ylabel = "Reactivity",
        legend = false,
        xticks = (1:n, ticklab),
        xlabel = "")

    savefig(plt, outpath)
    display(plt)
    
end

"""
count_species(indices)

Counts how often a species appears in the vector indices.

Input: indices = Vector{Any}, vector that includes the indices of the most reactive motifs

Output: num_spec = Dict{Int64, Int64}, frequencies of each species
"""
function count_species(indices::Vector{Any})

    # Flatten the vector
    flat_v = vcat(indices);

    # Step 1: Count the frequency of each species
    num_spec = Dict{Int64, Int64}()
    for element in flat_v
        if haskey(num_spec, element)
            num_spec[element] += 1
        else
            num_spec[element] = 1
        end
    end

    return num_spec
        
end # end function

""" ------------------------------------------
bar_species_frequencies()

Creates a bar plot with the number of appearances 
"""
function bar_species_frequencies(num_spec::Dict{Int64, Int64})
    # Prepare data for plotting
    unique_elements = collect(keys(num_spec))
    frequencies = collect(values(num_spec))
    
    # Create a bar plot or histogram
    bar(unique_elements, frequencies, 
        xlabel = "Node", 
        ylabel = "Frequency", 
        title = "Element Frequency Distribution")

end # end function


end # end module
