# Analysis code

## Folders

- **Estimates** - R scripts with the selected effect estimates from extracted data; one file per cancer. Analysis is based on final data in these files.
-- Stage III-IV - R scripts with the effect estimates for the sensitivity analysis using stage III/IV as a late stage definition across all cancers.
-- Timepoints - R scripts with the effect estimates for looking at the impact of timing.
- **Bootstrap** - R scripts with the bootstrap functions for the calculation of the 95% confidence intervals of the summary measures.

## Functions 

- Allcancer_functions.R - Functions used in estimate selection and analysis.

## R Markdown files - Analyses and Results

1. Primary Analysis
- Primary_surrogate1.Rmd
- Primary_surrogate3.Rmd

2. Subgroup Analyses
- Trials with surrogate measured earlier than mortality:
   - Subgroup_EarlyFU_surrogate1.Rmd
   - Subgroup_EarlyFU_surrogate3.Rmd

3. Test type:
   - Subgroup_Test_surrogate1.Rmd
   - Subgroup_Test_surrogate3.Rmd

4. Stage III-IV as late stage definition:
   - Sensitivity_StageDefinition_surrogate1.Rmd
   - Sensitivity_StageDefinition_surrogate3.Rmd

5. Impact of timing:
   - Timepoints.Rmd

6. Exploratory Analyses
   - Exploratory - FE.Rmd

# License

GNU GPL v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Copyright 2025, Nefeli Kouppa, Adam Brentnall

