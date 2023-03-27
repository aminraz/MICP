
# Contents



The Python scripts presented in this repository can simulate
microbially induced calcium carbonate precipitation (MICP) experiments
reported in the work of [Martinez et al (2013)](https://ascelibrary.org/doi/abs/10.1061/%28ASCE%29GT.1943-5606.0000787). Each python file
simulates one experiment that is mentioned in the docstring at the top of
each file. For the full description of each experiment refer to the
work of [Martinez et al. (2013)](https://ascelibrary.org/doi/abs/10.1061/%28ASCE%29GT.1943-5606.0000787) and [Martinez et al. (2014)](https://linkinghub.elsevier.com/retrieve/pii/S0266352X14000214).

A geochemical solver called [Phreeqc](https://www.usgs.gov/software/phreeqc-version-3) is used for geochemical
calculations. [Phreeqpy](https://www.phreeqpy.com/) is a python interface for Phreeqc that is used
here. Phreeqc handles various geochemical calculations using different
keyword blocks. A keyword block is a code segment in the input file
processed by Phreeqc. To develop a model for conducting MICP in a
column, we utilized the following keyword blocks for specific
purposes:

-   `Solution`: for defining initial and inlet solutions.
-   `Rate`: to define the rate of ureolysis and precipitation.
-   `Kinetics`: to add the kinetics using rate equations to the solutions.
-   `Transport`: to make a column with specific number of cells and initiate flow.
-   `Selected_Output`: to collect output data.

The full description of each of the aforementioned keyword blocks is
provided in the [Phreeqc manual](https://pubs.usgs.gov/tm/06/a43/). 

