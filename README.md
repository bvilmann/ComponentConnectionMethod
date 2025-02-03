# Component Connection Method
The Component Connection Method (CCM) composes a single small-signal model from predefined and interconnected dynamical systems in state space using Python and Excel. The method is similar to the `connect()` function known from [MATLAB](https://se.mathworks.com/help/control/ref/dynamicsystem.connect.html). 

The motivation of the implementation is to extend options of handling interconnected dynamical systems in Python-based analysis frameworks, as there is an increasing demand for small-signal models (SSMs) compatible with Python from transmission system operators [1]. The choice of using Excel was a practical need for a tabular editing view. This can easily be changed within the code

# Usage
## Research and development purposes
If you are using the code for research, you must cite the paper [2].

### BIBTEX

# Future work
- We are currently working on incorporating the CCM functionality with load flows for better small-signal sensitivity analysis features in Python-based investigations.

# License
The repository is distributed under the MIT License.

# References
[1] Energinet (2023), "Small-signal model requirement of DC-connected power park module", internal document, document number 23/04536-18, 

[2] Vilmann, B (2024), "Composing Small-Signal Models via the Component Connection Method with Python and Excel", IEEE, in publishing
