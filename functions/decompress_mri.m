function [mri_filename, temp_folder] = decompress_mri(mri_filename)

if ispc
    guid = System.Guid.NewGuid();
    temp_folder = char(guid.ToString());
else
    temp_folder = datestr(now,'dd-mmm-yyyy_HH-MM-SS');
end
mri_filename = gunzip(mri_filename, temp_folder);
mri_filename = mri_filename{1};