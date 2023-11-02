# env_global_epist

Code to generate the panels in the figures of the paper:

Juan Diaz-Colunga, Alvaro Sanchez & C. Brandon Ogbunugafor (2022). Environmental modulation of global epistasis in a drug resistance fitness landscape.

The ./data directory contains data files obtained from their original source [[1](https://doi.org/10.1111/evo.14428)].

The ./scripts directory contains the R scripts to produce the plots in the main paper and supplementary material. Plots are saved under the ./plots directory. Files description:

**eFunctions.R**   Contains auxiliary functions for data preprocessing (described in [[2](https://doi.org/10.1101/2022.06.21.496987)]).

**analyze_data.R**   Produces the panels for all main and supplementary figures in the paper.

[1] C. Brandon Ogbunugafor (2022). The mutation effect reaction norm (mu-rn) highlights environmentally dependent mutation effects and epistatic interactions. *Evolution*. DOI: [10.1111/evo.14428](https://doi.org/10.1111/evo.14428)

[2] Juan Diaz-Colunga, Abigail Skwara, Jean C. C. Vila, Djordje Bajic, Alvaro Sanchez (2022). Emergent ecosystem functions follow simple quantitative rules. *biorXiv*. DOI: [10.1101/2022.06.21.496987](https://doi.org/10.1101/2022.06.21.496987)
