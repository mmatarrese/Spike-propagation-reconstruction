function [Vol_ons, Vol_espr, Vol_spr, Vol_RES, coord_ESI_ons, coord_ESI_espr, coord_ESI_spr, scout_ESI_zones] = ...
    ESI_three_propagating_zones(source_map, scouts, scout_to_import, resection, volume, ons_percent)

% This function is part of the study
% Spike Propagation Mapping Reveals Effective Connectivity 
% and Predicts Surgical Outcome in Epilepsy
% Brain, 2023, Matarrese M.A.G. et al.

% Author: Matarrese M.A.G.
% contact: mag.matarrese@gmail.com


% Function to get the coordinates of onset, early spread, spead and 
% resection zones for one patient.

% INPUTS: 
% source_map is the event of interest in the source domain to analyze. This
% file can be obtained directly from Brainstorm. It is crucial that this
% struct contains the field "ImageGridAmp", i.e. the source map computed on
% Brainstorm has to be a Full Result matrix.
% scout_to_import is a scruct containing the scouts of the propagation in
% the source space that can be imported in Brainstorm for visualization.
% scout is a struct containing the coordinates in the 3D MRI space for each
% VRA, the corresponding center of mass, the boundary and convex hull
% definition, the timing of occurrence in the propagation, the area, the
% volume and the anatomical localization
% resection and volume are two files from brainstorm containing the 
% resected volume definition and the Grid of points (volume) on which the
% resection is built. 
% ons_percent is the percentage of propagation for defining the onset. In
% our study we defined onset as the first 10% of the total duration of the
% propagation, the following 10% was the early spread zone and all the rest
% was considered as late spread area.

% OUTPUTS: Vol_ons, Vol_espr, Vol_spr, Vol_RES are the volumes in cm3 of
% onset, early-spread, late-spread and resection respectively
% coord_ESI_ons, coord_ESI_espr, coord_ESI_spr, scout_ESI_zones are the
% coordinated in the MRI space of onset, early-spread and late-spread

% Definition of ESI Onset-zone and ESI spread-zone from scouts
time_ESI = [scouts.Time].';
time_ons=(ons_percent*(scouts(end).Time-scouts(1).Time))/100;
[~, ind_ons]=min(abs(time_ESI-time_ons));

%time_spr = 2*time_ons;
time_espr=(2*ons_percent*(scouts(end).Time-scouts(1).Time))/100;
[~, ind_espr]=min(abs(time_ESI-time_espr));

coord_ESI_ons = rmoutliers(unique(cell2mat({scouts(1:ind_ons).coord_vertices}.'),'row'));
coord_ESI_espr = unique(cell2mat({scouts(ind_ons:ind_espr).coord_vertices}.'),'row');
coord_ESI_spr = unique(cell2mat({scouts(ind_espr:end).coord_vertices}.'),'row');

% to define the spread zone as the remaining 85% in coord_ESI_spr ind_spr
% must be changed with end
Patch_ons = delaunayTriangulation(coord_ESI_ons);
Patch_espr = delaunayTriangulation(coord_ESI_espr);
Patch_spr = delaunayTriangulation(coord_ESI_spr);

if size(Patch_ons) ~= 0
    Patch_ons_convhull = convexHull(Patch_ons);
    [Vol_ons] = stlVolumeArea(coord_ESI_ons',Patch_ons_convhull'); %[cm3 cm2]
else
    Patch_ons_convhull = [];
    Vol_ons=NaN;
end

if size(Patch_espr) ~= 0
    Patch_espr_convhull = convexHull(Patch_espr);
    [Vol_espr] = stlVolumeArea(coord_ESI_espr',Patch_espr_convhull'); %[cm3 cm2]
    
else
    Patch_espr_convhull = [];
    Vol_espr=NaN;
end

if size(Patch_spr) ~= 0
    Patch_spr_convhull = convexHull(Patch_spr);
    [Vol_spr] = stlVolumeArea(coord_ESI_spr',Patch_spr_convhull'); %[cm3 cm2]
    
else
    Patch_spr_convhull = [];
    Vol_spr=NaN;
end

% To load resection and volume whitin it is defined to get all the vertices
% 3D coordinates

Vertices=resection.Vertices;
Vol_Grid=volume.GridLoc;
Coord_Res=Vol_Grid(Vertices,:);
Coord_Res=unique(Coord_Res,'rows'); %% added to modify res volume resampling
Patch_RES = delaunayTriangulation(Coord_Res);
Patch_RES_convhull = convexHull(Patch_RES);

[Vol_RES] = stlVolumeArea(Coord_Res',Patch_RES_convhull'); %[cm3 cm2]

% Struct to import in brainstorm with ESI Onset and Spread Zone
volume_loc=source_map.GridLoc;
seed_loc_ons=mean(coord_ESI_ons);
Vertices_ESI_ons = unique(cell2mat({scout_to_import(1:ind_ons).Vertices}));
[~, seed_ons]=minimum_distance(seed_loc_ons,volume_loc);
scout_ESI_zones(1).Vertices=Vertices_ESI_ons;
scout_ESI_zones(1).Seed=seed_ons;
scout_ESI_zones(1).Color=[1 0 0];
scout_ESI_zones(1).Label='ESI_Onset_zone';
scout_ESI_zones(1).Function='Mean';
scout_ESI_zones(1).Region='RO';
scout_ESI_zones(1).Handles=[];

seed_loc_espr=mean(coord_ESI_espr);
Vertices_ESI_espr = unique(cell2mat({scout_to_import(ind_ons:ind_espr).Vertices}));
[~, seed_espr]=minimum_distance(seed_loc_espr,volume_loc);
scout_ESI_zones(2).Vertices=Vertices_ESI_espr;
scout_ESI_zones(2).Seed=seed_espr;
scout_ESI_zones(2).Color=[0 0.2 0.8];
scout_ESI_zones(2).Label='ESI_early_Spread_zone';
scout_ESI_zones(2).Function='Mean';
scout_ESI_zones(2).Region='RO';
scout_ESI_zones(2).Handles=[];

seed_loc_spr=mean(coord_ESI_spr);
Vertices_ESI_spr = unique(cell2mat({scout_to_import(ind_espr:end).Vertices}));
[~, seed_spr]=minimum_distance(seed_loc_spr,volume_loc);
scout_ESI_zones(3).Vertices=Vertices_ESI_spr;
scout_ESI_zones(3).Seed=seed_spr;
scout_ESI_zones(3).Color=[0 0.2 0.8];
scout_ESI_zones(3).Label='ESI_Spread_zone';
scout_ESI_zones(3).Function='Mean';
scout_ESI_zones(3).Region='RO';
scout_ESI_zones(3).Handles=[];

end