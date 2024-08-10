function [PE_timeseries] = timeresolved_permEntropy(window_size, slide_samples, x, m, tau, ties)

% This function computes permutation entropy measures in a time- resolved 
% fashion. Specifically, for a timeseries 'x', windows of size
% 'window_size' are considered to compute the entropy measures over that
% current window.
% Windows are consecutively slid by 'slide_samples' sample points.
%
% This functions calls permEntropy.m. For a description of the parameters
% for permutation entropy calculation (i.e., m, tau, and ties), see the
% function header therein.
%
% Input:
%           window_size - (scalar) length of timeseries window
%           slide_samples - (scalar) number of sample points the window is
%                            shifted
%           x - (vector) the signal timeseries
%           m - (scalar) the motif length
%           tau - (scalar) the time delay for entropy calculation
%           ties - (string) method to deal with the theoretical possibility
%                   of rank ties (essentially non-existent in real-valued
%                   signals); can be 'sequence' or 'GaussianNoise' (see 
%                   permEntropy.m for details)
%
% Output:
%
%           PE_timeseries - (structure) contains the timeseries of the
%           entropy measures (unweighted as well as weighted permutation
%           entropy in both raw and normalized form) and some descriptive
%           parameters of the computation such as number of windows.
%
%
% Example:
%
%   x               = rand(1000,1)      ;
%   m               = 3                 ;
%   tau             = 1                 ;
%   ties            = 'sequence'        ;
%   window_size     = 80                ;
%   slide_samples   = 10                ;
% 
% 
%   [PE_timeseries] = timeresolved_permEntropy(window_size, slide_samples, x, m, tau, ties)
%
% SKrohn, August 2019
% -------------------------------------------------------------------------

sample_length   =   length(x)                                       ; % length of total timeseries

idx_vec         =   1:slide_samples:sample_length-window_size       ; % vector of indicies where each sliding window starts

nr_windows      =   length(idx_vec)                                 ; % number of sliding windows


% initialize output by tracking computation parameters
PE_timeseries.sample_length     =       sample_length                   ;
PE_timeseries.window_size       =       window_size                     ;
PE_timeseries.slide_samples     =       slide_samples                   ;
PE_timeseries.perc_overlap      =       1- slide_samples/window_size    ;
PE_timeseries.nr_windows        =       nr_windows                      ;

% set last entry of permutation entropy measures for implicit
% preallocation
PE_timeseries(nr_windows).PE    =       []                              ;


% loop over of windows of the timeseries
for i = 1:nr_windows
    
    current_window                          =   x(idx_vec(i): idx_vec(i) + window_size-1)       ; % get the current window of the timeseries
    
    [PE_i, PE_norm_i, WPE_i, WPE_norm_i]    =   permEntropy(current_window, m, tau, ties)       ; % compute the permutation entropy measures over the current window
    
    % assign permutation entropy measures over current window to output
    PE_timeseries(i).PE             =   PE_i            ;
    PE_timeseries(i).PE_norm        =   PE_norm_i       ;
    PE_timeseries(i).WPE            =   WPE_i           ;
    PE_timeseries(i).WPE_norm       =   WPE_norm_i      ;
    
end

end
