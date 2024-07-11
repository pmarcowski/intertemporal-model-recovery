# Model recovery analysis for intertemporal choice
This repository contains the data and analysis scripts for model recovery in intertemporal choice.

Intertemporal choice refers to decisions involving tradeoffs between costs and benefits occurring at different points in time, such as choosing between a smaller reward available immediately and a larger reward available after a delay.

The primary goal of this project is to validate our model selection procedure and assess the reliability of the best-fitting intertemporal choice model. This is achieveed through a model recovery analysis, which involves:

- Model comparison to identify the best-fitting model for explaining empirical choices.
- Simulation of experimental procedures using agents based on the best-fitting model.
- Model recovery analysis to determine if the best-fitting model and its true parameters can be accurately recovered from simulated data.

This approach allows us to:

- Verify the validity of our model selection procedure
- Assess the identifiability of the best-fitting model
- Evaluate the accuracy of parameter estimates
- Understand potential limitations in inferring decision processes from choice data

# Installation
No installation is required beyond having R and the necessary libraries. To run the script, simply clone this repository or download the script files and execute them within your R environment.

Clone this repository:
```
git clone https://github.com/username/intertemporal-choice-recovery.git
```

Navigate to the project directory:
```
cd path/to/intertemporal-choice-recovery
```

Refer to the header of `recover_model.R` script for a complete list of necessary libraries.

# Usage
Set the working directory to the project folder in your R environment:
```
setwd("path/to/intertemporal-choice-recovery")
```

Run the main analysis script `recover_model.R`.

# Project Structure
The project's key components are organized as follows:

- Model functions: `R/models/`
  Contains implementations of the intertemporal choice models.
- Helper functions: `R/tools/`
  Includes utility functions for data preprocessing, model fitting, and result visualization.
- Data: `data/`
  Stores the empirical data used in the analysis.
- Output: `output/`
  Contains generated results and visualizations.
- Checkpoints: `output/checkpoints/`
  Stores intermediate results for efficient analysis resumption.

## Models

Six intertemporal choice models are evaluated, implemented in `R/models/`:

Discounting models (with power choice rule):
- Exponential discounting model (EXPO)
- Hyperboloid model (HYPER2)
- Dual-Exponential model (DEXPO)

Heuristic models (with exponential choice rule):
- Intertemporal Choice Heuristic (ITCH)
- Difference-Ratio-Interest-Finance-Time model (DRIFT)
- Trade-off model (TRADE)

Model definitions are based on:
Wulff, D. U., & van den Bos, W. (2018). Modeling Choices in Delay Discounting. *Psychological Science, 29(11)*. [https://doi.org/10.1177/0956797616664342](https://doi.org/10.1177/0956797616664342)

## Data

The empirical data used for initial model comparison is stored in `data/` and is from:
Marzilli Ericson, K. M., White, J. M., Laibson, D., & Cohen, J. D. (2015). Money Earlier or Later? Simple Heuristics Explain Intertemporal Choices Better Than Delay Discounting Does. *Psychological Science, 26(6)*, 826–833. [https://doi.org/10.1177/0956797615572232](https://doi.org/10.1177/0956797615572232)

The data was collected in an experiment that presents participants with a series of binary choices between smaller-sooner and larger-later monetary rewards, varying the amounts and time delays to assess individual intertemporal preferences.

## Checkpoints

The analysis uses checkpoints to save intermediate results in `output/checkpoints/`:

1. Cross-validation results: "checkpoints/cross_validation.rds"
2. Agent parameters: "checkpoints/agent_parameters.rds"
3. Simulations: "checkpoints/simulations.rds"
4. Recovered parameters: "checkpoints/recovered_parameters.rds"

These checkpoints allow for efficient resumption of the analysis if interrupted. To re-run specific sections of the analysis, simply delete the corresponding checkpoint files in `output/checkpoints/`.

# Analysis
1. Model comparison using Monte Carlo Cross-Validation (MCCV) to identify the best-fitting model
2. Parameter estimation for the best-fitting model
3. Simulation of choices using agents based on the best-fitting model and estimated parameters
4. Model recovery analysis on simulated data:
- Recover the best-fitting model
- Estimate parameters and compare to true values
5. Bland-Altman analysis for assessing agreement between true and recovered parameters
6. Statistical inference using permutation tests
7. Visualization of results

Please note that the modeling procedure, particularly cross-validation, require relatively large computational power and may take some time.

# Results
The analysis generates several outputs:

- Visualization of model errors (`output/model_error.png`)
- Parameter recovery plots, including Bland-Altman plots (`output/parameter_recovery.png`)
- Combined figure for model evaluation and recovery (`output/figure.png`)

These results provide insights into the reliability of the model selection procedure, the identifiability of the best-fitting intertemporal choice model, and the accuracy of parameter recovery.

# License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/pmarcowski/intertemporal-model-recovery/blob/main/LICENSE) file for details.

# Contact
For any questions or feedback, please contact the author directly.
