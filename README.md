# Anatomical dynamics define cancer cachexia subtypes and identify systemic inflammation as a marker of lethal wasting

This repository contains the core analysis code for cachectic subtypes discovery from body composition data. 
To identify cachectic weight loss episodes, see https://github.com/reznik-lab/cachexia.

# Data requirements 
Input files are not included in the repository 

## Body composition file
A table of body composition deltas per patient including:
- A unique patient identifier
- A `[SCAN_TYPE]` columns which designates either 'Start-Cachexia' or 'Post-Cachexia'
- Columns corresponding to changes in body composition  (including SAT, VAT, SKM, pancreas, liver, spleen, kidney, BMD and all associated densities)

# Running code

## Preprocessing 
To process body composition data that is outside the physiological range, and remove large outliers use the 
`./processing/clean_bodycomp_data_func.R` or run 
`./processing/process_bodycomp.R`

## Calculating body composition changes
Run the `./processing/calculate_deltas_func.R` from `./processing/create_deltas_file.R` to quantify changes in body composition between Start-Cachexia and 'Post-Cachexia

## Subtype discovery 
Clustering using deltas file is available at `./processing/cluster_deltas.R`