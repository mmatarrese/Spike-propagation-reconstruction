load('Data1.mat')
load('Data2.mat')
win = 5; % average source map activity within this window
tol = 10; % spatial tolerance to include a point in resection 5mm
threshold=0.7;

[scout_to_import, scouts, CM_loc, VRA, Time_ext, ...
    Space_ext, Speed, Propagation_scouts, ROI_num,...
    Regions_ons] = evolving_scouts_creation_anatomy(ROIs, source_map, threshold, win);

ons_percent=10;

[Vol_ons, Vol_spr_early, Vol_RES, coord_ESI_ons, coord_ESI_spr_early, scout_ESI_zones] = ...
    ESI_propagating_zones(source_map, scouts, scout_to_import, resection, volume, ons_percent, 2*ons_percent);

[~, Vol_spr, ~, ~, coord_ESI_spr, scout_ESI_zones_late] = ...
    ESI_propagating_zones(source_map, scouts, ...
    scout_to_import, resection, volume, ons_percent, 100-ons_percent);

scout_ESI_zones(2).Label = 'ESI_early_Spread_zone';
scout_ESI_zones(3)=scout_ESI_zones_late(2);
scout_ESI_zones(2).Color=[1 0.5 0];


%%
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

%%
open 'Cortex.fig'
hold on
Coord_Res=volume.GridLoc(resection.Vertices,:);
Coord_Res=unique(Coord_Res,'rows'); 
Faces_boundary = boundary(Coord_Res, 0.7);
trisurf(Faces_boundary,Coord_Res(:,1),Coord_Res(:,2),Coord_Res(:,3),'Facecolor',[0.77 0.33 0.88],'FaceAlpha',0.4,'EdgeColor', 'none','FaceLighting','gouraud')

A=rmoutliers(coord_ESI_ons);
Faces_ons = boundary(A, 0.7);
trisurf(Faces_ons,A(:,1),A(:,2),A(:,3),'Facecolor',[1 0 0],'FaceAlpha',0.7,'EdgeColor', 'none','FaceLighting','gouraud')

A=rmoutliers(coord_ESI_spr_early);
Faces_espr = boundary(A, 0.7);
trisurf(Faces_espr,A(:,1),A(:,2),A(:,3),'Facecolor',[1 0.4 0.2],'FaceAlpha',0.7,'EdgeColor', 'none','FaceLighting','gouraud')

A=rmoutliers(coord_ESI_spr);
Faces_lspr = boundary(A, 0.7);
trisurf(Faces_lspr,A(:,1),A(:,2),A(:,3),'Facecolor',[0 0 1],'FaceAlpha',0.7,'EdgeColor', 'none','FaceLighting','gouraud')



