function [scout_to_import, scouts, CM_loc, Time_ext, Space_ext, ...
     Speed, Regions_ons] = evolving_scouts_creation_anatomy(ROIs, source_map, threshold, win)

% This function is part of the study
% Spike Propagation Mapping Reveals Effective Connectivity 
% and Predicts Surgical Outcome in Epilepsy
% Brain, 2023, Matarrese M.A.G. et al.

% Author: Matarrese M.A.G.
% contact: mag.matarrese@gmail.com
 
% Function to create a propagation in source space
% INPUTS: 
% ROIs is a struct (exported from Brainstorm) containing the AAL3
% scouts created on the specific patient anatomy. See the data uploaded for
% further details. This file allows identifying the anatolical gyrus from
% from where the propagation begins
% source_map is the event of interest in the source domain to analyze. This
% file can be obtained directly from Brainstorm. It is crucial that this
% struct contains the field "ImageGridAmp", i.e. the source map computed on
% Brainstorm has to be a Full Result matrix.
% Threshold is the activation cut-off to consider one source as active. In
% our study we set this parameter to 0.7
% win is the duration to average the activation values whitin this
% interval. In our study we set this to 5 ms.

% OUTPUTS:
% scout_to_import is a scruct containing the scouts of the propagation in
% the source space that can be imported in Brainstorm for visualization.
% The scouts are colorcoded from red (onset) to blue (last scout in a
% propagation). 
% scout is a struct containing the coordinates in the 3D MRI space for each
% VRA, the corresponding center of mass, the boundary and convex hull
% definition, the timing of occurrence in the propagation, the area, the
% volume and the anatomical localization
% Time_ext, Space_ext, Speed are respectively the Time duration of the
% propagatio, the Spatial dispacement and the speed of the propagation.
% ROI_num, 
% Regions_ons is a vector cointaining the index(es) of the anatomical ROI(s) 
% where the sources of the onset are located

% To export from brainstorm the uncostrained source map (it should be a
% Full results with ImagiGridAmp matrix) and call it "source_map"
ROIs = table({ROIs.Vertices}.', {ROIs.Label}.', 'VariableNames', {'Vertices', 'Label'});
ROIs.ind=[1:size(ROIs,1)]';

ROIs_expand = dataset;
for r=1:size(ROIs,1)
    temp = dataset;
    temp.Vertices = [ROIs.Vertices{r,1}]';
    temp.ind = ROIs.ind(r,1)*ones(size( temp.Vertices));
    Labels = cell(size(temp.Vertices));
    Labels(:) = {ROIs.Label{r,1}};
    temp.Labels = Labels;
    ROIs_expand = [ROIs_expand; temp];
end


dip_amp=source_map.ImageGridAmp;
volume_loc=source_map.GridLoc;
scout_to_import=struct;

x_dip=dip_amp(1:3:end,:);
y_dip=dip_amp(2:3:end,:);
z_dip=dip_amp(3:3:end,:);


MRI_amp=sqrt(x_dip.^2+y_dip.^2+z_dip.^2);
MRI_amp_th=MRI_amp./max(max(MRI_amp));

time_real=source_map.Time;
sample0 = find(source_map.Time==0);
samF=round(mean(diff(time_real).^-1));
sample_vec=1:size(MRI_amp_th,2);

% In order to get a minimum number of scout, the amplitude of 5ms (10
% sample for 2 kHz windows) of source map is averaged to get one single scout

win = round(win*10^-3*(samF+2));
% win = round(win*10^-3*(samF+2));
jump=1:win:size(MRI_amp_th,2);
jump_sec=time_real(jump);


MRI_amp_th_group=[];
CM_loc=[];
for kk=1:length(jump)
    if jump(kk)+win<size(MRI_amp_th,2)
        temp_MRI_amp=mean(MRI_amp_th(:,jump(kk):jump(kk)+win)')';
        
    elseif jump(kk)==size(MRI_amp_th,2)
        
    else
        temp_MRI_amp=mean(MRI_amp_th(:,jump(kk):size(MRI_amp_th,2))')';
    end
    vertices=find(temp_MRI_amp>threshold);
    scouts(kk).Vertices=vertices';
    
    scouts(kk).coord_vertices=volume_loc(vertices,:);
    MRI_amp_th_group=[MRI_amp_th_group temp_MRI_amp];
    
    if isempty(vertices)
        warning('No scout creation for this time windows')
    else
        scout_to_import(kk).Vertices=vertices';
        
        if length(vertices)==1
            seed_loc=volume_loc(vertices,:);
        else
            seed_loc=mean(volume_loc(vertices,:));
        end
        
        CM_loc=[CM_loc; seed_loc jump_sec(kk)];
        scouts(kk).CM_loc=seed_loc;
        scouts(kk).Time=jump_sec(kk);
        [~, seed]=minimum_distance(seed_loc,volume_loc);
        scout_to_import(kk).Seed=seed;
        scout_to_import(kk).Color=[0.8,0.4,0];
        label=strcat(num2str(round(jump_sec(kk)*1e+3)),'_',num2str(threshold));
        scout_to_import(kk).Label=label;
        scout_to_import(kk).Function='Mean';
        scout_to_import(kk).Region='RO';
        scout_to_import(kk).Handles=[]';
    end
    
end


%%%%%% To delete initial empty row
aaa=[];
for ll=1:length(scout_to_import)
    aaa=[aaa; isempty(scout_to_import(ll).Vertices)];
end
scout_to_import(aaa==1)=[];

aaa=[];
for ll=1:length(scouts)
    aaa=[aaa; isempty(scouts(ll).coord_vertices)];
end
scouts(aaa==1)=[];


% Volume of scouts
for jj=1:length(scouts)
    if length(scouts(jj).coord_vertices)>3
        patchFaces = boundary(scouts(jj).coord_vertices, 0.7);
        try
            patchFaces1 = convhulln(scouts(jj).coord_vertices);
        catch
            patchFaces1=patchFaces;
        end
        [totalVolume,totalArea] = stlVolumeNormals(scouts(jj).coord_vertices',patchFaces1'); %[m3 m2]
        scouts(jj).Faces_boundary=patchFaces;
        scouts(jj).Faces_convexhull=patchFaces1;
        scouts(jj).totalVolume = totalVolume*1e+6;
        scouts(jj).totalArea = totalArea*1e+4;
    end
end

%%%%%% Define color map for scouts and identify ROI of overlap
colors=flipud(jet(size(scout_to_import,2)));
for c=1:size(scout_to_import,2)
    scout_to_import(c).Color=colors(c,:);
    Vert_to_test = scouts(c).Vertices;
    Regions = ROIs_expand.ind(ismember(ROIs_expand.Vertices,Vert_to_test));
    if c==1
        Regions_ons = Regions;
    end
    [GC,GR,GP] = groupcounts(Regions);
    % Percentage = GC*100/size(Vert_to_test,2);
    Labels = unique(ROIs_expand.Labels(ismember(ROIs_expand.Vertices,Vert_to_test)));
    scouts(c).Vertices_Region = Regions;
    scouts(c).ROI_ind = GR;
    scouts(c).ROI_Labels = Labels;
    scouts(c).Occurrences = GC;
    scouts(c).ROI_percentage = GP;
end

Time_ext=(scouts(end).Time-scouts(1).Time)*1000;
for i=2:1:size(CM_loc,1)
    dist(i-1)=sqrt((CM_loc(i,1)-CM_loc(i-1,1))^2+(CM_loc(i,2)-CM_loc(i-1,2))^2+(CM_loc(i,3)-CM_loc(i-1,3))^2);
end
Space_ext=sum(dist)*100; % space [cm]
Speed=Space_ext/Time_ext; % cm/ms

scouts([scouts.Time].'<=0)=[];

end
