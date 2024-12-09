# FARM [http://www.farm-model.org/](http://www.farm-model.org/)

**Multi-grid Eulerian model** for dispersion, transformation, and deposition of airborne pollutants in gas and aerosol phases.

### Suited for:
- **Case studies**, episodes reconstruction, and studies on pollutants formation and accumulation.
- **Regional air quality assessments**.
- **Analysis of scenarios** and the effects of emissions abatement policies.
- **Pollutants forecast at multiple scales**, fed by prognostic meteorological models.
# REPOSITORY STRUCTURE

.. code-block:: text

    ├── model_mod.f90
    ├── python_code
    │   ├── bash_templates
    │   │   ├── new_run_submit.template.bsh
    │   │   ├── run_filter.template.bsh
    │   │   ├── run_submit_ens.template.bsh
    │   │   ├── run_submit.template.bsh
    │   │   ├── submit_ens.template.bsh
    │   │   ├── submit_filter.template.bsh
    │   │   └── submit.template.bsh
    │   ├── data_assimilation_cameo_path_manager.py
    │   ├── logs_orchestrator
    │   │   └── farm_to_dart_full_logs
    │   ├── orchestrator_farm_dart.py
    │   ├── orchestrator.py
    │   ├── orchestrator_utils.py
    │   ├── perturbation
    │   │   ├── bash_emissions
    │   │   │   ├── perturb_fields_emission.bsh
    │   │   │   └── run_pertub.bsh
    │   │   ├── delete_missing_value.py
    │   │   ├── perturb_emi_multi.py
    │   │   ├── perturb_emi_opt.py
    │   │   ├── perturb_emi.py
    │   │   ├── perturb_icbc_opt.py
    │   │   ├── replace_perturb_into_original_emission.py
    │   │   └── submit_replace.sh
    │   ├── post_processing
    │   │   ├── area_km.py
    │   │   ├── ellipsoid.py
    │   │   ├── handle_time_conventions.py
    │   │   ├── infograph.py
    │   │   ├── obs_avg.py
    │   │   ├── pdf_fit.py
    │   │   ├── plot_avg_maps.py
    │   │   ├── plot_maps_check.py
    │   │   ├── plot_maps.py
    │   │   ├── plot_obs_space_new.py
    │   │   ├── plot_obs_space.py
    │   │   ├── plot_obs_space_sub.py
    │   │   ├── plot_rms_maps.py
    │   │   ├── plot_total_st.py
    │   │   ├── post_processing_maps.py
    │   │   ├── process_assim_nc.py
    │   │   ├── query_station.py
    │   │   ├── split_times.py
    │   │   ├── submit.py
    │   │   └── test.py
    │   ├── __pycache__
    │   ├── run_orchestrator.bsh
    │   └── test.py
    ├── sandbox
    │   ├── test_24h
    │   │   └── shift.py
    │   └── test_time_split
    │       └── split_times.py
    └── work
        ├── handle_time_conventions.py
        ├── input.nml
        ├── input_template.nml
        ├── input_template_noise.nml
        ├── pdf_fit.py
        ├── plot_avg_maps.py
        ├── plot_maps_check.py
        ├── plot_maps.py
        ├── plot_obs_space.py
        ├── plot_obs_space_sub.py
        ├── plot_total_st.py
        ├── process_assim_nc.py
        ├── query_station.py
        ├── quickbuild.sh
        ├── split_times.py
        └── submit.py

