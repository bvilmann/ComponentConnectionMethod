# Component Connection Method
The Component Connection Method (CCM) composes a single small-signal model from predefined and interconnected dynamical systems in state space using Python and Excel. The method is similar to the `connect()` function known from [MATLAB](https://se.mathworks.com/help/control/ref/dynamicsystem.connect.html). 

The motivation of the implementation is to extend options for handling interconnected dynamical systems in Python-based analysis frameworks. Using Excel was a practical need for a tabular editing view. This dependency can easily be changed within the code.

# Future work
We are currently working on incorporating the CCM functionality for nonlinear dynamical models together with a load flow solver for better small-signal sensitivity analysis features in Python-based investigations.

# License
The repository is distributed under the MIT License.

# Research and development purposes
If you are using the code for research, please cite the paper, [Composing Small-Signal Models via the Component Connection Method with Python and Excel](https://ieeexplore-ieee-org.proxy.findit.cvt.dk/document/10892445):

|Citation||
|----|----|
|Text|B. Vilmann, O. Saborío-Romano and D. Müller, "Composing Small-Signal Models via the Component Connection Method with Python and Excel," 2024 59th International Universities Power Engineering Conference (UPEC), Cardiff, United Kingdom, 2024, pp. 1-6, doi: 10.1109/UPEC61344.2024.10892445. keywords: {Power engineering;Codes;Frequency-domain analysis;Electricity;Power system stability;Stability analysis;Generators;Power generation;component connection method;state-space;small-signal model;applied engineering;open-source;grid code requirements},|
|BIBTEX|@INPROCEEDINGS{10892445, author={Vilmann, Benjamin and Saborío-Romano, Oscar and Müller, Daniel},   booktitle={2024 59th International Universities Power Engineering Conference (UPEC)},    title={Composing Small-Signal Models via the Component Connection Method with Python and Excel},    year={2024},   volume={},   number={},    pages={1-6},   keywords={Power engineering;Codes;Frequency-domain analysis;Electricity;Power system stability;Stability analysis;Generators;Power generation;component connection method;state-space;small-signal model;applied engineering;open-source;grid code requirements}, doi={10.1109/UPEC61344.2024.10892445}}|
