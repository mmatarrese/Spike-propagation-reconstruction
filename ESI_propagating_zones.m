function [Vol_ons, Vol_spr, Vol_RES, coord_ESI_ons, coord_ESI_spr, scout_ESI_zones] = ...
    ESI_propagating_zones(source_map, scouts, scout_to_import, resection, volume, ons_percent, spr_percent)


% Function to get the volume of onset, spead and resection zones and the
% coordinates of all the points and of internal points.

% Definition of ESI Onset-zone and ESI spread-zone from scouts
time_ESI = [scouts.Time].';
% ons_percent=15;
time_ons=(ons_percent*(scouts(end).Time-scouts(1).Time))/100;
[~, ind_ons]=min(abs(time_ESI-time_ons));

%time_spr = 2*time_ons;
time_spr=(spr_percent*(scouts(end).Time-scouts(1).Time))/100;
[~, ind_spr]=min(abs(time_ESI-time_spr));

coord_ESI_ons = unique(cell2mat({scouts(1:ind_ons).coord_vertices}.'),'row');
coord_ESI_spr = unique(cell2mat({scouts(ind_ons:ind_spr).coord_vertices}.'),'row');
% to define the spread zone as the remaining 85% in coord_ESI_spr ind_spr
% must be changed with end
Patch_ons = delaunayTriangulation(coord_ESI_ons);
Patch_spr = delaunayTriangulation(coord_ESI_spr);

if size(Patch_ons) ~= 0
    Patch_ons_convhull = convexHull(Patch_ons);
    [Vol_ons] = stlVolumeArea(coord_ESI_ons',Patch_ons_convhull'); %[cm3 cm2]
else
    Patch_ons_convhull = [];
    Vol_ons=NaN;
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

seed_loc_spr=mean(coord_ESI_spr);
Vertices_ESI_spr = unique(cell2mat({scout_to_import(ind_ons:ind_spr).Vertices}));
[~, seed_spr]=minimum_distance(seed_loc_spr,volume_loc);
scout_ESI_zones(2).Vertices=Vertices_ESI_spr;
scout_ESI_zones(2).Seed=seed_spr;
scout_ESI_zones(2).Color=[0 0.2 0.8];
scout_ESI_zones(2).Label='ESI_Spread_zone';
scout_ESI_zones(2).Function='Mean';
scout_ESI_zones(2).Region='RO';
scout_ESI_zones(2).Handles=[];


end