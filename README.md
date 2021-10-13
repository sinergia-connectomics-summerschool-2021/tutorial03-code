# tutorial03-code
In this tutorial, we will give an overview of the preprocessing steps to project the scalp EEG data into the source space and resume them at the level of regions of interest (ROIs) using Fieldtrip. To do so we will you the first subject of the VEPCON: https://openneuro.org/datasets/ds003505/versions/1.0.2 described in: https://doi.org/10.1101/2021.03.16.435599

run the scripts in the following order:
* eeg_preprocessing.m
* mri_preprocessing.m
* source_modelling.m
* ROI_time_courses.m

Useful Fieldtrip Tutorials:
* https://www.fieldtriptoolbox.org/example/use_independent_component_analysis_ica_to_remove_eog_artifacts/
* https://www.fieldtriptoolbox.org/example/use_independent_component_analysis_ica_to_remove_ecg_artifacts/
* https://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/
* https://www.fieldtriptoolbox.org/tutorial/sourcemodel/
* https://www.fieldtriptoolbox.org/workshop/oslo2019/forward_modeling/
* https://www.fieldtriptoolbox.org/faq/how_can_i_fine-tune_my_bem_volume_conduction_model/#converting_segmentation_to_segmentation