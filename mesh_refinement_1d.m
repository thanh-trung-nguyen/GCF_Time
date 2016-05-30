function newmesh = mesh_refinement_1d(oldmesh,indicationfun, threshold)

% refine a mesh based on an indication function 
% refine mesh at grid points such that the absolute value of the indication
% function is larger than threshold*maximum of the indication function
% we first find the grid points with large value of the indication
% function, then refine the mesh on the left and right of these points

indicationfun = abs(indicationfun);
NrOld = length(oldmesh);

refine_grid_points = 0*oldmesh;
refine_grid_points(indicationfun > threshold*max(indicationfun)) = 1; % find all grid points to refine the mesh

Nr_refine_points = sum(refine_grid_points);

NrNew = NrOld + 2*Nr_refine_points; %maximum number of grid points in the new mesh
newmesh = zeros(1,NrNew);

idx = 1; NrNewUpdate = 1; 
while idx <= NrOld - 1
    if (refine_grid_points(idx) == 0)
        newmesh(NrNewUpdate) = oldmesh(idx);
        NrNewUpdate = NrNewUpdate+ 1;
    else
        newmesh(NrNewUpdate) = (oldmesh(idx-1) + oldmesh(idx))/2;
        newmesh(NrNewUpdate+1) = oldmesh(idx);
        newmesh(NrNewUpdate+2) = (oldmesh(idx) + oldmesh(idx+1))/2;
        NrNewUpdate = NrNewUpdate + 3;
        if refine_grid_points(idx+1) == 1
            NrNewUpdate = NrNewUpdate - 1;
        end
    end
    idx = idx + 1;
end
idx = NrOld;
    if (refine_grid_points(idx) == 0)
        newmesh(NrNewUpdate) = oldmesh(idx);
    else
        newmesh(NrNewUpdate) = (oldmesh(idx-1) + oldmesh(idx))/2;
        newmesh(NrNewUpdate+1) = oldmesh(idx);
        NrNewUpdate = NrNewUpdate + 1;
    end

newmesh = newmesh(1:NrNewUpdate);


