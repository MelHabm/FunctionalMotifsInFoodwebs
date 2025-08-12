# FunctionalMotifsInFoodwebs
Code is written using Julia v1.11.1. Make sure the correct version is installed on your computer. 

Run examplenetwork.jl to read the interaction matrix from ex_intmat.csv and calculate the most reactive motifs and create a histogram for the example that shows the distribution of the proportion of system reactivity that is explained by the most reactive instance of different motifs. 

The whole project is divided into the following modules:

MotifFinder.jl:
    Use interaction matrices to find a number of different motifs in food webs that have been associated with foodweb stability in previous studies and combinations of those motifs. 
    The motifs included are:
    
    2 species:
      predator-prey
      mutualism (in apparent competition motifs)
    
    3 species:
      apparent competition (two prey species that share the same predator)
      exploitative competition (two predator species that share the same prey)
      omnivory (one species eats both species in a predator-prey pair)
      tri-trophic foodchains
    
    4 species:
      four species foodchains
      bifan (combination of two exploitative competition motifs with the same predators)
      diamond (combination of apparent competition and exploitative competition motifs)
      combined tri-trophic foodchains (two tri-trophic chains that share the same top predator)

FunctionalMotifs.jl:
    Indentifies the most reactive instances of the motifs in question and calculates for each the reactivity value z in each community, the indices of this motif inds, and the proportion of system reactivity explained by this motif. Different motifs can be called using the corresponding keys (see function func_react_motifs() for more explanation of the keys in question).

Reactivity.jl:
    Calculates the reactivity of a symmetric Jacobian and motifs of interest.

GeneralizedModel.jl: 
    Generate Jacobian matrices of stable systems based on a Generalized Foodweb Model.

ReactivityResults.jl: 
    Module for plotting the results. Includes all functions used to plot histograms that compare the probability density of the proportion of system reactivity that is explained by a single motif.

ex_intmat.csv: 
    The interaction matrix of the example network used in the paper "Functional Motifs in Foodwebs and Networks". (i,j)th element = 1, if i eats j.
