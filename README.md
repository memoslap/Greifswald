# Greifswald
This repository is related to the script we use for neuroimaging analysis in MeMoSLAP Greifswald

# Simulation_scripts_P9
Preparation | Install SimNIBS and Python with all related packages.
00b-charm.sh is a script that runs the SIMNIBS brain reconstruction. To run go to terminal and write ./00b-charm.sh
02a-optimized-sim.py is a script that runs the simulations based on projects specific region of interests and settings found in Memoslap folder. To run it, go to the terminal and write simnibs_python 02a-optimized-sim.py.
02b-optimized_collect_vals.py is a script that collects values from all simulations and creates one pickle and excel file for all subjects and all projects. To run it, go to the terminal and write python 02b-optimized-collect-vals.py
02c-optimized_fig.py and 02d-optimized_violin are a scripts that create the figures for the simulations. To run it, go to the terminal and write python 02c-optimized_fig.py or python 02d-optimized_violin.py
02e-optimized_report.py is a script that creates screenshots of all simulations, e-field magnitude with ROI and electrode positioning. To run it go to the terminal and write simnibs_python 02e-optimized_report.py
