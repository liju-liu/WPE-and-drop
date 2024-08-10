function [ Timematrix ] = TimesExtract_Faster(Datapath ,Maskpath, varargin )

% Input variables
% Datapath: Location of the subject's BOLD signal file
% Maskpath: Location of the template file for extracting time series

% Output variable
% TimeMatrix: The average BOLD time series within the template for each subject, extracted based on the template

DefaultOpts=struct('datatype',0);
Args=parseInOpts(DefaultOpts,varargin);
if Args.datatype == 1
    Data = Datapath;
    Mask = Maskpath;
else
    Matlab_version=version('-release');
    Matlab_version=Matlab_version(1:end-1);
    if str2num(Matlab_version)>2017
        Readfun='niftiread';
    else
         Readfun='y_Read';
    end
    eval(['Mask=',Readfun,'(Maskpath);']);
    eval(['Data = ',Readfun,'(Datapath);']);
end
Timematrix = zeros(size(Data,4),max(Mask(:)));
[N1,N2,N3,N4] = size(Data);
[M1,M2,M3] = size(Mask);
if N1 ~= M1 || N2 ~= M2 || N3 ~= M3
    error('The dimension of mask is not equal with dimension of data \n');
end
Reshdata = reshape(Data,N1*N2*N3,N4);
for loop=1:max(Mask(:))
    Regiondata = Reshdata(Mask==loop,:);
    Timematrix(:,loop) = nanmean(Regiondata);
end
end

