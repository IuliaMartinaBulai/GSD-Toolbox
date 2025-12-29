# GSD-Toolbox
Graph Spectral Dimension Toolbox

(c) Developed by Iulia Martina Bulai, Isidoros Iakovidis, Tim Steger

Description:
This toolbox contains MATLAB code for estimating the spectral dimension of graphs using stochastic eigenvalue counting methods.
The approach is based on Chebyshev polynomial approximations and Hutchinson trace estimators.
The code supports both synthetic and real-world graph datasets.

The numerical experiments implemented here correspond to those presented in the associated paper. 
Iulia Martina Bulai, Isidoros Iakovidis, Tim Steger "Efficient Calculation of Graph Spectral Dimension"

Datasets:
The toolbox supports six datasets:
-ring (synthetic, generated in code)
-path (synthetic, generated in code)
-swiss_roll (synthetic, generated in code)
-bunny (real point cloud dataset, generated in code)
-minnesota (road network graph)
-brain (large-scale brain graph)

Scripts:
-graph_spectral_dimension_main.m
Main script for spectral dimension estimation.
Computes eigenvalues when feasible, performs stochastic counting, fits the scaling law, and produces figures.

-brain_demo.m
Stochastic-only analysis for large brain graph.
Avoids full eigenvalue computation and saves numerical results and figures.

-brain_data_manipulation.m
Post-processing script for brain experiments.
Reloads saved results and regenerates figures without recomputation.

-filters_plot.m
Plots the filter functions used.

-select_filter.m
Defines the spectral filters used in the counting procedure.

-get_dataset.m
Loads or generates datasets depending on the dataset name.

-Brain_results 
Folder containing the results of the brain

Utils
Folder containing necessary functions from sgwt package

This toolbox is released under the GNU General Public License (GPL), version 3 or later.

The software is free to use, modify, and redistribute under the terms of the GPL.
It is distributed without any warranty, including implied warranties of merchantability or fitness for a particular purpose.

A copy of the GNU General Public License should be included with this toolbox.
If not, see http://www.gnu.org/licenses/.
