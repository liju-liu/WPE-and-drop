function [PE, PE_norm, WPE, WPE_norm]     =   permEntropy(x, m, tau, ties)

% This function computes permutation entropy measures for an input
% timeseries x.
%
% Input:
%           x           (vector) - the timeseries
%           m           (scalar) - the motif length
%           tau         (scalar) - the time delay to partition x
%           ties        (string) - argument to specify which method to
%                                  apply to account for tied ranks in the
%                                  ordinal pattern (e.g., [1.2 3.1 1.2]
%                                  would result in a pattern of [1 2 1];
%                                  can be 'sequence' (i.e., the first
%                                  occurence will get the lower rank) or
%                                  'GaussianNoise' (i.e., minimal Gaussian
%                                  noise is added to resolve the ties);
%                                  although a theoretical possibility,
%                                  early studies showed this to be
%                                  extremely rare in real-valued signals
%
% Recommended input values are tau = 1 and m in [3, ..., 7], cf.
% Bandt, C., & Pompe, B. (2002). Permutation entropy: a natural complexity
% measure for time series. Physical review letters, 88(17), 174102.
%
% Output:
%
%           PE          (scalar) - permutation entropy in bits
%           PE_norm     (scalar) - permutation entropy in bits normalized
%                                  to [0,1]
%           WPE         (scalar) - weighted permutation entropy after
%                                  Fadlallah et al. (2013) using the
%                                  variance in each partition to factor in
%                                  amplitude information
%           WPE_norm    (scalar) - weighted permutation entropy in bits
%                                  normalized to [0,1]
%
% Reference:
%
% Fadlallah, B., Chen, B., Keil, A., & Príncipe, J. (2013).
% Weighted-permutation entropy: A complexity measure for time series
% incorporating amplitude information. Physical Review E, 87(2), 022911.
%
%
% Example:
%   x       =   [3.5 2.1 5.3 1.2 4.2 5.7 6.5]   ;   % real-valued signal vector
%   m       =   3                               ;   % motif length
%   tau     =   1                               ;   % time lag
%   ties    =   'sequence'                      ;   % account for tied ranks
% 
% [PE, PE_norm, WPE, WPE_norm]     =   permEntropy(x, m, tau, ties)
%
%
% SKrohn, August 2019
% -------------------------------------------------------------------------

% input checks
if nargin<4
    error('Missing input arguments.')
end


possible_perms          =   perms(1:m)                  ; % all possible permutations of rank patterns given a motif length of m

ncols_partition         =   length(x) - (m-1)*tau       ; % number of ways to split up the time series into chunks of length m, moving on tau steps in the data; this is equivalent to index 'j' in Fadlallah et al. (2013), see paragraph above eq. (1) there

partition_mat           =   nan(m, ncols_partition)     ; % preallocate a matrix where each column corresponds to a chunk of time series x of length m

rank_mat                =   nan(m, ncols_partition)     ; % preallocate the matrix of oridinal patterns

pattern_freg            =   [possible_perms zeros(size(possible_perms,1),1)]      ; % count the frequency of the occurring patterns (column 4)

weighted_pattern_freq   =   [possible_perms zeros(size(possible_perms,1),1)]      ; % count the weighted frequency of the occurring patterns (column 4)

pattern_mat             =   nan(2, ncols_partition)     ; % a matrix to collect the means (row 1) and weights (row 2) for every partition snippet of the time series in partition_mat


% split up time series into chunks of length m, progressing tau steps at a
% time
for i = 1:ncols_partition
    
    partition_mat(:,i)              =   x(i:(i+m-1))                                            ; % partition the time series
    
    [~,idx_tmp, sequence_rank]      =   unique(partition_mat(:,i))                              ; % find the ordinal pattern of the i'th chunk
    
    
    % this is the case when the ranks are tied
    if size(sequence_rank,1)~=size(unique(sequence_rank),1)    
        
        [count, ~, idxcount] =  histcounts(sequence_rank,numel(idx_tmp))    ; % count unique values; cf. https://www.mathworks.com/matlabcentral/answers/175086-finding-non-unique-values-in-an-array
        idx_non_unique       =  count(idxcount)>1                           ; % index non-unique entries in partitioned snippet
        
        % which method is used for accommodating the ties?
        switch ties
                        
            % resolve ties by adding minimal Gaussian noise to the ties
            % data
            case 'GaussianNoise'
                
                partition_mat(idx_non_unique,i) =   partition_mat(idx_non_unique,i) + rand(length(partition_mat(idx_non_unique,i)),1)*10^(-10)      ; % add minimal Gaussian noise to the non-unique entries             
                [~,~, sequence_rank]            =   unique(partition_mat(:,i))                                                                      ; % repeat the ranking with resolved ties
               
            % resolve ties by order of occurrence in the snippet
            case 'sequence'
                              
                ties_sorted          =   sort(partition_mat(:,i))                    ; % make sure the relative ranks across ties are accurate (e.g., [1.2 3 1.2] must be [1 3 2] not [1 2 3], i.e. global maxima must remain global)
                [~, seq_rank_tmp]    =   ismember(partition_mat(:,i), ties_sorted)   ; % cf. https://www.mathworks.com/matlabcentral/answers/33296-ranking-ordering-values-with-repeats
                
                
                seq_rank_tmp(idx_non_unique)    =   nan                             ; % set non-unique ranks to NaN
                
                possible_ranks                  =   1: size(seq_rank_tmp,1)                     ; % possible ranks in the motif
                missing_ranks                   =   setdiff(possible_ranks, seq_rank_tmp)       ; % those ranks not assigned (the former ties)
                
                seq_rank_tmp(isnan(seq_rank_tmp))   =   missing_ranks                           ; % assign ranks in consecutive order (i.e., the first tied value will get the lower rank in the pattern)
                
                sequence_rank                   =   seq_rank_tmp                                ; % assign rank sequence with resolved ties back to rank sequence variable            
                          
        end        
    end
    
    
    rank_mat(:,i)               =   sequence_rank                                           ; % keep in rank sequence in rank matrix
    
    [~, pattern_idx]            =   ismember(rank_mat(:,i)', possible_perms, 'rows')        ; % check which of the possible patterns the i'th partition snippet is equivalent to
    
    pattern_freg(pattern_idx,4) =   pattern_freg(pattern_idx,4) +    1                      ; % count the occurrences of each pattern
    
    pattern_mat(1,i)            =   mean(partition_mat(:,i))                                ; % mean of the i'th timeseries partition
    pattern_mat(2,i)            =   (1/m)*sum((partition_mat(:,i)-pattern_mat(1,i)).^2)     ; % weighting factor for each pattern based on the snippet variance (eq. 7 in Fadlallah et al., 2013)
    
    weighting_factor_i          =   pattern_mat(2,i)                                        ; % reassign weighting factor for clarity
    
    weighted_pattern_freq(pattern_idx,4)  =   weighted_pattern_freq(pattern_idx,4) + 1*weighting_factor_i           ; % count the occurrences of each pattern weighted by the weighting factor of the i'th snippet
    
    % % note: the inbuilt variance function leads to the same WPE output (taking var(partition_mat(:,i)) returns the normalized variance, which is identical to eq. 7 in Fadlallah et al., 2013)
    % variance_i                            =   var(partition_mat(:,i))                                     ; % variance of the i'th partition
    % weighted_pattern_freq(pattern_idx,4)  =   weighted_pattern_freq(pattern_idx,4) + 1*variance_i         ; % compute the variance-weighted occurrence of the i'th partition
    
end


% compute the permutation entropy on x 
prob_vec    =   pattern_freg(:,4)/ncols_partition                               ;   % the (unweighted) probability vector for each pattern occurrence

PE          =    -sum(prob_vec(prob_vec>0).*log2(prob_vec(prob_vec>0)))         ;   % the permutation entropy is the negative sum over all pattern probabilities (p_i) * the logarithm of that probability; zero entries discarded (term is zero by definition)
PE_norm     =    PE / log2(size(possible_perms,1))                              ;   % normalize permutation entropy to the interval [0 1]; the denominator is equivalent to log2(factorial(m))


% compute the weighted permutation entropy on x
weighted_rank_patterns      =   weighted_pattern_freq(weighted_pattern_freq(:,4)>0,4)       ; % weighted rank patterns with zeros discarded (term is zero in entropy calculation by definition)
weighted_prob_vec           =   weighted_rank_patterns./sum(weighted_rank_patterns)         ; % weighted probability vector; sums to unity

WPE         =     -sum(weighted_prob_vec.*log2(weighted_prob_vec))              ;   % the weighted permutation entropy is the Shannon entropy over the weighted probability vector
WPE_norm    =     WPE / log2(size(possible_perms,1))                            ;   % normalize permutation entropy to the interval [0 1]; the denominator is equivalent to log2(factorial(m))


%  keyboard

end