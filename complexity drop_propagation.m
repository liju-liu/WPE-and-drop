% This script processes BOLD signal data to analyze and visualize complexity drops in brain regions.
% 
% Key Steps:
% 1. Load and concatenate WPE matrices across subjects to compute a global drop threshold.
% 2. Create a binary Drop_Affinity_matrix to identify regions with complexity drops based on the threshold.
% 3. Analyze the Drop_Affinity_matrix to identify and record cascades, including their start,half, peak, and end windows.
% 4. Calculate node probability vectors for the start, half, peak, and end layers of each cascade.
% 5. Save individual and group-level results, including duration and node probability vectors.
% 6. Convert and save node probability vectors as NIfTI (.nii) files for visualization.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
clc
clear all;

data_path_wpematrix=('the save path of WPE_mat');
Savepath_Drop_Affinity_matrix=('Set a save path');
temp1=dir(data_path_wpematrix);
temp1=temp1(3:end); 

for q=1:length(temp1)
    WPE_mat=importdata(fullfile(data_path_wpematrix,temp1(q).name));
    if(q==1)
        All_WPE_mat=WPE_mat; 
    else
        %Concatenate the WPE matrices longitudinally to calculate the drop threshold (1%).
        All_WPE_mat=[All_WPE_mat;WPE_mat];                       
    end
end

%calculate the drop threshold (1%).
drop_thresh = prctile(All_WPE_mat(:),1); 


%------------------------------------------------------------------------------%
%For each time window, determine the brain regions that experience a complexity drop,
%set the corresponding values to 1, and otherwise to 0. This results in a binary matrix
% with dimensions: number of windows ¡Á number of brain regions, named Drop_Affinity_matrix.
temp2=dir(data_path_wpematrix);
temp2=temp2(3:end); 
for j=1:length(temp2)        
        Drop_Affinity_matrix=importdata(fullfile(data_path_wpematrix,temp2(j).name));
        Drop_Affinity_matrix(Drop_Affinity_matrix<=drop_thresh)=1;
        Drop_Affinity_matrix(Drop_Affinity_matrix<1)=0;
        save([Savepath_Drop_Affinity_matrix,'\',temp2(j).name],'Drop_Affinity_matrix');
end
%------------------------------------------------------------------------------%

%% 
%------------------------------------------------------------------------------%
data_path=Savepath_Drop_Affinity_matrix;%Load the Drop_Affinity_matrix
Savepath=('Set a save path');
mask246_nii='the path of Brainnetome atlas';%the path of Brainnetome atlas
% load('BNA_labels.mat');
% BNA_Network_labels=table2array(BNA_labels(:,3));


mkdir([Savepath,'\Propagation\']);
mkdir([Savepath,'\transfer_steps\']);
mkdir([Savepath,'\cascades_indiv\']);
mkdir([Savepath,'\layer_regions\']);
mkdir([Savepath,'\origion_layer_indiv']);
mkdir([Savepath,'\peak_layer_indiv']);
mkdir([Savepath,'\end_layer_indiv']);
mkdir([Savepath,'\half_layer_indiv']);
mkdir([Savepath,'\cascades_indiv']);


temp=dir(data_path);
temp=temp(3:end);   

for i=1:length(temp)
   name=strsplit(temp(i).name,'.');
   data=importdata(fullfile(data_path,temp(i).name));
   drop_regions=sum(data,2);%Calculate the number of brain regions that experience a complexity drop in each time window.
   [n_windows,n_regions]=size(data);
   cascades=0;
   peak_layer=[];
   start_layer=[];
   end_layer=[];
   j=2;
   
   %record the four key time windows(start  layer, half layer, peak layer, and end layer) of each cascade.
        while( j<=n_windows-2)

            if(drop_regions(j)>=10 && drop_regions(j)>drop_regions(j+1) &&drop_regions(j)>drop_regions(j-1))
                              
                cascades=cascades+1;
                peak_layer(cascades)=j;
                m=j-1;
                n=j+1;               
                if (m==1)
                    start_layer(cascades)=m;
                else
                    while(m>=2 && drop_regions(m)>drop_regions(m-1)&&drop_regions(m-1)~=0)
                        m=m-1;
                    end
                    start_layer(cascades)=m;
                end
                while(n<=n_windows-1 && drop_regions(n)>drop_regions(n+1)&&drop_regions(n+1)~=0)
                    n=n+1;
                end
                end_layer(cascades)=n;
                
                j=n;
                
            end
            j=j+1;
        end
%-----------------------------------------------------------------------------%
%Calculate the node probability values for the cascade's start  layer, 
% half layer, peak layer, and end layer at the group level.
start_layer_probability_indiv=zeros(1,n_regions);
peak_layer_probability_indiv=zeros(1,n_regions);
end_layer_probability_indiv=zeros(1,n_regions);
half_layer_probability_indiv=zeros(1,n_regions);    
    if(i==1)
        start_layer_probability_group=zeros(1,n_regions);
        peak_layer_probability_group=zeros(1,n_regions);
        end_layer_probability_group=zeros(1,n_regions);
        half_layer_probability_group=zeros(1,n_regions);
        sum_indiv=0;        
%         origion_layer_regions=zeros(length(temp),1);
%         peak_layer_regions=zeros(length(temp),1);
%         end_layer_regions=zeros(length(temp),1);   
    end

    half_layer=zeros(1,cascades);
    cascades_half=0;
    for k=1:cascades
        start_layer_probability_indiv=start_layer_probability_indiv+data(start_layer(k),:);
        peak_layer_probability_indiv=peak_layer_probability_indiv+data(peak_layer(k),:);
        end_layer_probability_indiv=end_layer_probability_indiv+data(end_layer(k),:);
        Duration_start_peak(i,k)=peak_layer(k)-start_layer(k);%the duration for the start-peak 
        Duration_peak_end(i,k)=end_layer(k)-peak_layer(k);%the duration for the peak-end 
        if Duration_start_peak(i,k)>1
            half_layer(1,k)=fix(Duration_start_peak(i,k)/2)+start_layer(k);
            half_layer_probability_indiv=half_layer_probability_indiv+data(half_layer(1,k),:);
            cascades_half=cascades_half+1;
        end
           
%         origion_layer_regions(i,1)= origion_layer_regions(i,1)+drop_regions(start_layer(k),1);
%         peak_layer_regions(i,1)= peak_layer_regions(i,1)+drop_regions(peak_layer(k),1);
%         end_layer_regions(i,1)= end_layer_regions(i,1)+drop_regions(end_layer(k),1);
        
    end  
    
    if(cascades~=0)
        start_layer_probability_indiv=start_layer_probability_indiv/cascades;
        peak_layer_probability_indiv=peak_layer_probability_indiv/cascades;
        end_layer_probability_indiv=end_layer_probability_indiv/cascades;
        if (cascades_half~=0)
        half_layer_probability_indiv=half_layer_probability_indiv/cascades_half;
        end
        sum_indiv=sum_indiv+1;  
        save([Savepath,'\origion_layer_indiv\',temp(i).name],'start_layer_probability_indiv');
        save([Savepath,'\half_layer_indiv\',temp(i).name],'half_layer_probability_indiv');
        save([Savepath,'\peak_layer_indiv\',temp(i).name],'peak_layer_probability_indiv');
        save([Savepath,'\end_layer_indiv\',temp(i).name],'end_layer_probability_indiv');
     end       
    
    start_layer_probability_group=start_layer_probability_group+start_layer_probability_indiv;
    half_layer_probability_group=half_layer_probability_group+half_layer_probability_indiv;
    peak_layer_probability_group=peak_layer_probability_group+peak_layer_probability_indiv;
    end_layer_probability_group=end_layer_probability_group+end_layer_probability_indiv;
    
    %Record the number of cascade for each individual.
    cascades_indiv(i)=cascades;
    
end

%Calculate the duration (in terms of the number of time windows) for the start-peak 
%and peak-end phases of each cascade. The start-end duration is the sum of the two durations.
Duration_start_peak(Duration_start_peak==0)=nan;
Duration_peak_end(Duration_peak_end==0)=nan;
transfer_steps_up_mean=nanmean(Duration_start_peak,2);
transfer_steps_down_mean=nanmean(Duration_peak_end,2);

save([Savepath,'\transfer_steps\','transfer_steps_up_mean'],'transfer_steps_up_mean');
save([Savepath,'\transfer_steps\','transfer_steps_down_mean'],'transfer_steps_down_mean');
save([Savepath,'\cascades_indiv\','cascades_indiv'],'cascades_indiv');

%The node probability vectors corresponding to the four time windows (start, half, peak, and end)
Average_probability_start=start_layer_probability_group/sum_indiv;
Average_probability_peak=peak_layer_probability_group/sum_indiv;
Average_probability_end=end_layer_probability_group/sum_indiv;
Average_probability_half=half_layer_probability_group/sum_indiv;

%%
%Convert the node probability vectors corresponding to the four time windows into NIfTI (.nii) files.
[mask1,info1] = y_Read(mask246_nii);
img1=(Average_probability_start);
for n=1:n_regions
    mask1(find(mask1==n))=img1(n);
end
info1.dt=[16,0];
y_Write(mask1,info1,fullfile(Savepath,'\Propagation\','Average_cascade_origion.nii'));

[mask2,info2] = y_Read(mask246_nii);
img2=(Average_probability_peak);
for n=1:n_regions
    mask2(find(mask2==n))=img2(n);
end
info2.dt=[16,0];
y_Write(mask2,info2,fullfile(Savepath,'\Propagation\','Average_cascade_peak.nii'));

[mask3,info3] = y_Read(mask246_nii);
img3=(Average_probability_end);
for n=1:n_regions
    mask3(find(mask3==n))=img3(n);
end
info3.dt=[16,0];
y_Write(mask3,info3,fullfile(Savepath,'\Propagation\','Average_cascade_end.nii'));

[mask4,info4] = y_Read(mask246_nii);
img4=(Average_probability_half);
for n=1:n_regions
    mask4(find(mask4==n))=img4(n);
end
info4.dt=[16,0];
y_Write(mask4,info4,fullfile(Savepath,'\Propagation\','Average_cascade_half.nii'));

disp('Congratulation!!!!!!');

