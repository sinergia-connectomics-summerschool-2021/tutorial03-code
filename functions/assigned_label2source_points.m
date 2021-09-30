function [label_vector, ROIs_elec]= ...
    assigned_label2source_points(sourcemodel, nifti_atlas, tbl_idx2label)

% get all voxels which are not outside the brain, in white matter or
% cerebellum
I = find(nifti_atlas.anatomy+1); % add 1 to include non-labeled part of the MRI
% convert the indices into voxel coordinates
[x,y,z] = ind2sub(nifti_atlas.dim,I);
% convert voxel coordinates into coord in mm in real space
coord_mm = nifti_atlas.transform *[x,y,z,ones(length(x),1)]';
% find the 26 closest voxel around each source points
[idx] = knnsearch(coord_mm(1:3,:)',sourcemodel.pos(sourcemodel.inside,:),'K',26);
% get the most present tissue around the source point
idx_assgnTissue = mode(nifti_atlas.anatomy(I(idx)),2);
% assigned the label of the most present tissue around the source point
label_vector = repmat({'Outside'},length(sourcemodel.inside),1);
bln_zeros = idx_assgnTissue == 0;
assigned_label_vector = repmat({'Outside'},sum(sourcemodel.inside),1);
assigned_label_vector(~bln_zeros) = ...
    tbl_idx2label.abbreviation(idx_assgnTissue(~bln_zeros));
label_vector(sourcemodel.inside) = assigned_label_vector;
% get the centroids of the ROIs
centroid_coord = zeros(height(tbl_idx2label),3);
for id_abrv = 1:height(tbl_idx2label)
    bln = nifti_atlas.anatomy(:) == id_abrv;
    centroid_coord(id_abrv,:) = mean(coord_mm(1:3,bln),2)';
end

centroid_coord = cellfun(@(str)  mean(sourcemodel.pos(strcmp(str,label_vector),:),1),...
    tbl_idx2label.abbreviation, 'uni',0);


centroid_coord = cell2mat(centroid_coord);
ROIs_elec = struct('elecpos', centroid_coord,...
    'chanpos',centroid_coord, 'label',{tbl_idx2label.abbreviation});




% % for debugging
% pos = sourcemodel.pos;
% inside = sourcemodel.inside;
% bln = strcmp(label_vector, 'L Hipp');
% figure; scatter3(pos(inside,1),pos(inside,2),pos(inside,3))
% hold on;
% scatter3(pos(bln,1),pos(bln,2),pos(bln,3), 'filled')
% 
% bln = strcmp(label_vector, 'L STG1');
% scatter3(pos(bln,1),pos(bln,2),pos(bln,3), 'filled')
% bln = strcmp(label_vector, 'L STG2');
% scatter3(pos(bln,1),pos(bln,2),pos(bln,3), 'filled')
% axis vis3d