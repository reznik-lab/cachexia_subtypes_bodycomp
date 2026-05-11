# Anatomical dynamics define cancer cachexia subtypes and identify systemic inflammation as a marker of lethal wasting

This repository contains the core analysis code for cachectic subtypes discovery from body composition data. 
To identify cachectic weight loss episodes from longitudinal BMI data, see https://github.com/reznik-lab/cachexia.

# Data requirements 
Input files are not included in the repository 

## Body composition file
A table of body composition with two time points per patient including:
- A unique patient identifier (e.g. `MRN`)
- A `[SCAN_TYPE]` columns which designates either 'Start-Cachexia' or 'Post-Cachexia'
- Columns corresponding to changes in body composition  (including SAT, VAT, SKM, pancreas, liver, spleen, kidney, BMD and all associated densities)

# Running code

## Preprocessing 
To process body composition data that is outside the physiological range, and remove large outliers use the 
`./processing/clean_bodycomp_data_func.R` or run 
`./processing/process_bodycomp.R`

## Calculating body composition changes
Use `Supplementary Table 1`, which is already pre-processed body composition values
This function calculates % change for all areas and volumes, and raw changes for all densities
`./processing/create_deltas_file.R` includes additional filtering to remove outliers
Run the `./processing/calculate_deltas_func.R` from `./processing/create_deltas_file.R` to quantify changes in body composition between Start-Cachexia and Post-Cachexia

## Subtype discovery 
Use `Supplementary Table 2`, or the outputs from `./processing/create_deltas_file.R`. 
Clustering using deltas file is available at `./processing/cluster_deltas.R`