# MIMESI Orchestrator (Operational Guide)

<div class="mimesi-hero-logo">
  <img alt="MIMESI logo" src="assets/Logo_mimesi.svg" width="220">
</div>

This repository runs a time-stepped workflow driven by a single YAML configuration file.

You typically use it to:

1. Prepare a run directory (paths, inputs, templates).
2. Point the orchestrator to a YAML config.
3. Monitor logs until completion.

## What is included here

Only the **CHIMERE 2017** pipeline is enabled in this repository. If your configuration asks for other pipelines or model types, the program will stop at startup with an error. For the implementation of other models please contact
alessandro.dausilio@suez.com

## What you need to know to run it

- You run `main.py` from a shell: `python -u main.py -c path/to/config.yaml`
- All run locations (data folders, templates, listing files) are provided through the YAML config
- Logs are written automatically (see [Quickstart](quickstart.md) and [Troubleshooting](troubleshooting.md))

## Where to go next

- [Quickstart](quickstart.md): minimal step-by-step run
- [Installation](installation.md): environment prerequisites
- [Example Config (ITA7)](config_mimesi_ITA7.md): explanation of the recommended starting config
- [Troubleshooting](troubleshooting.md): common startup/runtime failures
