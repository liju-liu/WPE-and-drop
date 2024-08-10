% % This script calculates the Weighted Permutation Entropy (WPE) matrix 
% for BOLD signal time series using the sliding window method.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clearvars
close all
clc

%specify windowing parameters for time-resolved complexity analysis
window_size         =   20                 ;
slide_samples       =   1                  ;

% specify parameters for computation of weighted permutation entropy
m                   =   3                  ;
tau                 =   1                   ;
ties                =   'sequence'          ;

% the BOLD timeseries to analyze 
path_data=('E:\MRI_246\HC\NoZscore\Timematrix');
%WPE Savepath
Savepath=('E:\MRI_246\HC\NoZscore\20_1');


mkdir([Savepath,'\individual_WPE_coarsening_3\']);

temp=dir(path_data);
temp=temp(3:end); 

for j=1:length(temp)
    
    Mask_Timematrix=importdata(fullfile(path_data,temp(j).name));
    
    % get number of regions and number of windows
    n_regions           =   size(Mask_Timematrix,2)                                          ;
    n_windows           =   length(1:slide_samples:(size(Mask_Timematrix,1)-window_size))    ;
    
    % preallocate array of WPE timeseries
    WPE_mat             =   nan(n_windows, n_regions);
    
    % loop over regions (can run as parfor but left sequential here)
    for i = 1:n_regions
    
        % BOLD timeseries of region i
        x                   =   Mask_Timematrix(:,i)                             ;

        % run complexity estimation on BOLD signal of region i
        WPE_timeseries      =   timeresolved_permEntropy(window_size, ...
                                                         slide_samples, ...
                                                         x, m, tau, ties)   ;

        % collect normalized WPE values over windows
        WPE_mat(:,i)        =   [WPE_timeseries.WPE_norm]' ;
        
        
                                     
    % keyboard
    end
      %Save each individual's WPE matrix (dimension: number of windows ¡Á number of brain regions).
      save([Savepath,'\individual_WPE\',temp(j).name],'WPE_mat');
      fprintf([temp(j).name, '\n']);
end
     disp('Congratulation!!!!!!');




    
 

