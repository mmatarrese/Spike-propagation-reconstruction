% This code is part of the study
% Spike Propagation Mapping Reveals Effective Connectivity 
% and Predicts Surgical Outcome in Epilepsy
% Brain, 2023, Matarrese M.A.G. et al.

% Author: Matarrese M.A.G.
% contact: mag.matarrese@gmail.com

% This MAIN reproduces part of Figure 4 of the above research article. This
% first part is necessary to load and define the correct files to process
% in order to get a propagation in source domain. These files can be
% replaced with others to analyze a new patient. 
path=cd;
cd('Data')
files = dir('*.mat');
for i=1:length(files)
    load(files(i).name);
end
cd(path)
source_map.ImageGridAmp = [ImageGridAmp1 ImageGridAmp2 ImageGridAmp3...
    ImageGridAmp4 ImageGridAmp5 ImageGridAmp6...
    ImageGridAmp7 ImageGridAmp8 ImageGridAmp9];
source_map.GridLoc=GridLoc;
source_map.GoodChannel=GoodChannel;
source_map.Time=Time;
clearvars -except ROIs resection volume source_map
%%
win = 5; % average source map activity within this window
tol = 10; % spatial tolerance to include a point in resection 5mm
threshold=0.7;

[scout_to_import, scouts, CM_loc, Time_ext, ...
    Space_ext, Speed, Regions_ons] = ...
    evolving_scouts_creation_anatomy(ROIs, source_map, threshold, win);

ons_percent=10;
[Vol_ons, Vol_espr, Vol_spr, Vol_RES, coord_ESI_ons, coord_ESI_spr_early, coord_ESI_spr, scout_ESI_zones] = ...
    ESI_three_propagating_zones(source_map, scouts, scout_to_import, resection, volume, ons_percent);


%% Plot propagation in MRI space:
% Press one button on the keyboard to get the next scout
Coord_Res=volume.GridLoc(resection.Vertices,:);
Coord_Res=unique(Coord_Res,'rows'); 
Faces_boundary = boundary(Coord_Res, 0.7);
open 'Cortex.fig'
hold on
trisurf(Faces_boundary,Coord_Res(:,1),Coord_Res(:,2),Coord_Res(:,3),'Facecolor',[0.77 0.33 0.88],'FaceAlpha',0.4,'EdgeColor', 'none','FaceLighting','gouraud')

color=flipud(jet(size(scouts,2)));
for pp=1:size(scouts,2)
    
    scouts_one = scouts(pp);
    Coord = scouts_one.coord_vertices;
    if size(Coord,1)~=1
        Coord = rmoutliers(scouts_one.coord_vertices);
    end
    
    if isempty(scouts_one.Faces_boundary)
        hold on
        scatter_scout=scatter3(Coord(:,1),Coord(:,2),Coord(:,3),30,'filled','MarkerFaceColor',color(pp,:));
    else
        Faces_boundary = boundary(Coord,0.7);
        hold on
        surf_scout=trisurf(Faces_boundary,Coord(:,1),Coord(:,2),Coord(:,3),'Facecolor',color(pp,:),'FaceAlpha',0.7,'EdgeColor', 'none','FaceLighting','gouraud');
    end

    pause
    set(surf_scout,'Facecolor','none')
    if isempty(scouts_one.Faces_boundary)
        set(scatter_scout,'MarkerFaceColor','none')
    end
    
end

%% Plot ESI zones in MRI space
open 'Cortex.fig'
hold on
Coord_Res=volume.GridLoc(resection.Vertices,:);
Coord_Res=unique(Coord_Res,'rows'); 
Faces_boundary = boundary(Coord_Res, 0.7);
trisurf(Faces_boundary,Coord_Res(:,1),Coord_Res(:,2),Coord_Res(:,3),'Facecolor',[0.77 0.33 0.88],'FaceAlpha',0.4,'EdgeColor', 'none','FaceLighting','gouraud')

A=coord_ESI_ons;
Faces_ons = boundary(A, 0.4);
trisurf(Faces_ons,A(:,1),A(:,2),A(:,3),'Facecolor',[1 0 0],'FaceAlpha',0.7,'EdgeColor', 'none','FaceLighting','gouraud')

A=coord_ESI_spr_early;
Faces_espr = boundary(A, 0.4);
trisurf(Faces_espr,A(:,1),A(:,2),A(:,3),'Facecolor',[1 0.4 0.2],'FaceAlpha',0.7,'EdgeColor', 'none','FaceLighting','gouraud')

A=coord_ESI_spr;
Faces_lspr = boundary(A, 0);
trisurf(Faces_lspr,A(:,1),A(:,2),A(:,3),'Facecolor',[0 0 1],'FaceAlpha',0.7,'EdgeColor', 'none','FaceLighting','gouraud')



