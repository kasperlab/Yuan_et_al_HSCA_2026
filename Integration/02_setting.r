py_run_string("merge_by = ['DonorID', 'SampleID', 'LibraryPlatform', 'ExperimentID']")
py_run_string("factors = ['StudyID', 'DonorID_short', 'Age', 'Sex', 'Ethnicity1', 'SubjectType',
'SampleID', 'SampleType', 'SampledSiteCondition', 'Tissue', 'LibraryPlatform', 'StrandSequence', 'SequencingPlatform', 'ReferenceGenome', 'SampleStatus',
'AnatomicalRegionLevel1', 'AnatomicalRegionLevel2', 'AnatomicalRegionLevel3',
'Ethnicity2', 'EthnicityDetail', 'BiologicalUnit', 'SampleCultured', 'DonorID', 'ExperimentID']")
py_run_string("factors_only_binary = ['DonorID_short', 'SampleID', 'DonorID']")
py_run_string("annotation = ['OriginalAnnotation', 'SeparateAnnotation_1', 'SeparateAnnotation_2', 'IntegratedAnnotation_2_v35', 'IntegratedAnnotation_2_v65', 'IntegratedAnnotation_3_v65']")
############
# Signatures
############
markers_list = py$marker_dict_JID
