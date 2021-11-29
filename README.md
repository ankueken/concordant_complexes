# concordant_complexes
 
Requirements: 

- MATLAB 2020b
- CobraToolbox v3 (Heirendt et. al., Creation and analysis of biochemical constraint-based models: the COBRA Toolbox v3.0, Nature Protocols, volume 14, pages 639–702, 2019 doi.org/10.1038/s41596-018-0098-2.)
- SFAnalysis (https://www.nature.com/articles/s41467-019-08746-5, https://github.com/adbroido/SFAnalysis)


Usage:

Find concordant complexes under two scenarios:

(i) Considering irreversibility constraints
- Run concordant_pairs_irrev.m

(ii) Considering irreversibility and growth optimality constraints
- Run concordant_pairs_objective.m

1. models are preprocessed (removal of blocked reactions, split reversible reactions into two irreversible reactions)
2. balanced complex detection
3. concordant complex detection
4. find concordance modules


Files for analysis:

Number of concordant complex pairs and complexes in concordance relation 
- Results_TableS1.m


Degree, effective degree, module size distributions
- Results_write_degseq_TableS3.m
- Results_effective_deg_TableS3.m
- Results_TableS4_metabolite_degree.m
- Results_overview_power_law_TableS3.m

- To analyze the fit to power law we use the code of Broido et al. (SFAnalysis (https://www.nature.com/articles/s41467-019-08746-5, https://github.com/adbroido/SFAnalysis)) - -code not copied to this repository


Modules and metabolic subsystems
- Results TableS5_S6_classification.m
- Results_module_pathways_*.m
- Results_hypergeom_test_modules_*.m


Reducibility index
- Results_reducibility_index.m


Random network variants
- folder 'randomization\' includes random network variants
- folder 'Results_concordant_randomization_irreversibility_considered\' 
  includes Results of calculating concordant complexes and modules running
  function 'concordant_randomized_networks' for 58 random network variants
- Results_statistics_random_variants.m 
