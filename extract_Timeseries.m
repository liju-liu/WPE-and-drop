 %Extract average BOLD time series for each ROI. 

clc 
clear all

%the path of Brainnetome atlas
maskpath = '...\BN_Atlas_246_3mm.nii';


%Output file save path.
Savepath='the save path';

mkdir([Savepath,'\Timematrix\']);


%Path to save the preprocessed data.
path_data=('...\FunImgARWSDCF');
temp=dir(path_data);
temp=temp(3:end);   

for i=1:length(temp)
    %Extract average BOLD time series for each ROI.   
    Mask_Timematrix = TimesExtract_Faster(fullfile(path_data,temp(i).name,'Filtered_4DVolume.nii'),maskpath );
    %Save the extracted ROI time series.
    save([Savepath,'\Timematrix\',temp(i).name,'.mat'],'Mask_Timematrix');
    fprintf([temp(i).name, '\n']);
     
end
        
     disp('Congratulation!!!!!!');