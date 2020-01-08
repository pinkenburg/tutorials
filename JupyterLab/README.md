This is an example [sPHENIX](https://github.com/sPHENIX-Collaboration/) [ROOT Jupyter Notebook](https://root.cern.ch/notebooks/HowTos/HowTo_ROOT-Notebooks.html) to run the sPHENIX anaysis on the [BNL SDCC Jupyter Lab](https://jupyter.sdcc.bnl.gov/). Please use lauch the HTC Jupyter hub for this tutorial: 

![SDCC HTC Jupyter Hub](https://jupyter.sdcc.bnl.gov/images/jupyterhub-logo-card-htc.png)

Jupyter Notebook provide an analysis interface that can integrate documentation, code and output plot in a single file. It supports all sPHENIX environment. It is suitable for final analysis step, while large scale analysis and reconstruction should go through the compiled `c++` analysis module.  

Two examples are provided: 
* run the analysis in ROOT C++ macros: [example_tpc_testbeam2019_root.ipynb](./example_tpc_testbeam2019_root.ipynb)
* `pyROOT` macro perform similar task if you would like to use `python` for analysis:  [example_tpc_testbeam2019_python.ipynb](./example_tpc_testbeam2019_python.ipynb)
