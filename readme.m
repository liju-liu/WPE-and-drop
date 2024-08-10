% README: Code Usage Instructions
% This guide provides an overview of the scripts used to analyze BOLD signal data 
% and compute the complexity properties within brain regions.
% 
% Step 1: Extract BOLD Signal Time Series
% The first step involves extracting the BOLD signal time series from a brain atlas template. 
% This is achieved using the extract_Timematrix.m script.
% 
% Script: extract_Timematrix.m
% Input: Brain atlas template and BOLD signal files.
% Output: Time matrix containing the average BOLD signal time series for each brain region (ROI).
% 
% Step 2: Calculate Weighted Permutation Entropy (WPE)
% Once the time series is extracted, the next step is to calculate the Weighted Permutation Entropy (WPE)
% for each region. This is done using the WPE.m script.
% Script: WPE.m
% Input: BOLD signal time series extracted from the previous step.
% Output: WPE matrix, representing the complexity of BOLD signals over time.

% Step 3: Calculate Complexity Drop Propagation Properties
% The final step is to analyze the complexity drop propagation and compute the node probability spatial topology for the four key time points (start, half, peak, and end) of each cascade. This involves running the script that calculates these properties based on the WPE matrices.
% % Script: complexity drop_propagation.m
% Input: WPE matrix from Step 2.
% Output:
% Complexity drop propagation properties.
% Node probability spatial topology for the start, half, peak, and end time points at group level.
