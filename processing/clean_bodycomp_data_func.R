# function to qc OSCAR outputs
# physiological ranges available in supplementary table

clean_bodycomp_data <- function(bodycomp, l3 = TRUE) {
  
  # make sure columns are numeric 
  to_numeric <- function(df, cols) {
    for (col in cols) {
      df[[col]] <- as.numeric(df[[col]])
    }
    return(df)
  }
  
  # bodycomp columns 
  num_cols <- c(
    "MuscleValues.L1TotalMuscleMeanHU", "MuscleValues.L3TotalMuscleMeanHU",
    "MuscleValues.L1TotalMuscleMedianHU", "MuscleValues.L3TotalMuscleMedianHU",
    "MuscleValues.L1MuscleOnlyMedianHU", "MuscleValues.L3MuscleOnlyMedianHU",
    "MuscleValues.L1MuscleOnlyMeanHU", "MuscleValues.L3MuscleOnlyMeanHU",
    "CalciumScoring.AbdominalAgatston", "CalciumScoring.TotalAorticAgatston",
    "LiverValues.LiverMedianHU", "LiverValues.LiverVolume",
    "SpleenValues.SpleenMedianHU", "SpleenValues.SpleenVolume",
    "L3FatValues.L3VATSATRatio", "L1FatValues.L1VATSATRatio",
    "KidneyValues.KidneyMedianHU", "KidneyValues.KidneyVolume",
    "PancreasValues.PancreasMedianHU", "PancreasValues.PancreasVolume",
    "L3FatValues.L3VATMedian", "MuscleValues.L3IMATMeanHU", "MuscleValues.L3IMATMedianHU",
    "L1FatValues.L1VATMedian", "MuscleValues.L1IMATMeanHU", "MuscleValues.L1IMATMedianHU",
    "T10FatValues.T10VATMedian", "MuscleValues.T10IMATMeanHU", "MuscleValues.T10IMATMedianHU",
    "T12FatValues.T12VATMedian", "MuscleValues.T12IMATMeanHU", "MuscleValues.T12IMATMedianHU",
    "MuscleValues.T10TotalMuscleMeanHU", "MuscleValues.T12TotalMuscleMeanHU",
    "MuscleValues.T10TotalMuscleMedianHU", "MuscleValues.T12TotalMuscleMedianHU",
    "MuscleValues.T10MuscleOnlyMedianHU", "MuscleValues.T12MuscleOnlyMedianHU",
    "MuscleValues.T10MuscleOnlyMeanHU", "MuscleValues.T12MuscleOnlyMeanHU"
  )
  bodycomp <- to_numeric(bodycomp, num_cols)

  bodycomp <- bodycomp %>%
    mutate(across(contains("BMD"), ~ as.numeric(as.character(.)))) %>%
    mutate(across(contains("BMD"), ~ ifelse(. >= -50 & . <= 1200, ., NA_real_)))
  
  qc_rules <- list(
    "MuscleValues.L3MuscleOnlyMeanHU"   = c(-50, 200),
    "MuscleValues.L3MuscleOnlyMedianHU" = c(-50, 200),
    "MuscleValues.L3MuscleMeanHU"       = c(-50, 200),
    "MuscleValues.L3MuscleMedianHU"     = c(-50, 200),
    "MuscleValues.L3TotalMuscleMeanHU"  = c(-50, 200),
    "MuscleValues.L3TotalMuscleMedianHU"= c(-50, 200),
    "MuscleValues.L3MuscleOnlyArea"     = c(25, 500),
    "MuscleValues.L3TotalMuscleArea"    = c(25, 500),
    "L3FatValues.L3TATArea"             = c(0.1, 1500),
    "L3FatValues.L3VATArea"             = c(0.1, 1200),
    "L3FatValues.L3SATArea"             = c(0.1, 1000),
    "L3FatValues.L3SATMedian"           = c(-120, -30),
    "L3FatValues.L3VATMedian"           = c(-120, -60),
    
    "MuscleValues.T10MuscleOnlyMeanHU"   = c(-50, 200),
    "MuscleValues.T10MuscleOnlyMedianHU" = c(-50, 200),
    "MuscleValues.T10MuscleMeanHU"       = c(-50, 200),
    "MuscleValues.T10MuscleMedianHU"     = c(-50, 200),
    "MuscleValues.T10TotalMuscleMeanHU"  = c(-50, 200),
    "MuscleValues.T10TotalMuscleMedianHU"= c(-50, 200),
    "MuscleValues.T10MuscleOnlyArea"     = c(25, 500),
    "MuscleValues.T10TotalMuscleArea"    = c(25, 500),
    "T10FatValues.T10TATArea"            = c(0.1, 1500),
    "T10FatValues.T10VATArea"            = c(0.1, 1200),
    "T10FatValues.T10SATArea"            = c(0.1, 1000),
    "T10FatValues.T10SATMedian"          = c(-120, -30),
    "T10FatValues.T10VATMedian"          = c(-120, -60),
    
    "MuscleValues.T12MuscleOnlyMeanHU"   = c(-50, 200),
    "MuscleValues.T12MuscleOnlyMedianHU" = c(-50, 200),
    "MuscleValues.T12MuscleMeanHU"       = c(-50, 200),
    "MuscleValues.T12MuscleMedianHU"     = c(-50, 200),
    "MuscleValues.T12TotalMuscleMeanHU"  = c(-50, 200),
    "MuscleValues.T12TotalMuscleMedianHU"= c(-50, 200),
    "MuscleValues.T12MuscleOnlyArea"     = c(25, 500),
    "MuscleValues.T12TotalMuscleArea"    = c(25, 500),
    "T12FatValues.T12TATArea"            = c(0.1, 1500),
    "T12FatValues.T12VATArea"            = c(0.1, 1200),
    "T12FatValues.T12SATArea"            = c(0.1, 1000),
    "T12FatValues.T12SATMedian"          = c(-120, -30),
    "T12FatValues.T12VATMedian"          = c(-120, -60),
    
    "MuscleValues.L1MuscleOnlyMeanHU"    = c(-50, 200),
    "MuscleValues.L1MuscleOnlyMedianHU"  = c(-50, 200),
    "MuscleValues.L1TotalMuscleMeanHU"   = c(-50, 200),
    "MuscleValues.L1TotalMuscleMedianHU" = c(-50, 200),
    "MuscleValues.L1MuscleOnlyArea"      = c(25, 500),
    "MuscleValues.L1TotalMuscleArea"     = c(25, 500),
    "L1FatValues.L1TATArea"              = c(0.1, 1500),
    "L1FatValues.L1VATArea"              = c(0.1, 1200),
    "L1FatValues.L1SATArea"              = c(0.1, 1000),
    "L1FatValues.L1SATMedian"            = c(-120, -30),
    "L1FatValues.L1VATMedian"            = c(-120, -60),
    
    "CalciumScoring.AbdominalAgatston"   = c(0, 40000),
    "CalciumScoring.TotalAorticAgatston" = c(0, 40000),
    
    "LiverValues.LiverMedianHU"          = c(-50, 180),
    "LiverValues.LiverVolume"            = c(100, 5000),
    
    "SpleenValues.SpleenMedianHU"        = c(10, 500),
    "SpleenValues.SpleenVolume"          = c(50, 6000),
    
    "L3FatValues.L3VATSATRatio"          = c(0.02, 10),
    "L1FatValues.L1VATSATRatio"          = c(0.02, 10),
    
    "KidneyValues.KidneyMedianHU"        = c(5, 300),
    "KidneyValues.KidneyVolume"          = c(50, 750),
    
    "PancreasValues.PancreasMedianHU"    = c(-10, 260),
    "PancreasValues.PancreasVolume"      = c(20, 140),
    
    "BMDL3Values.BMDL3HighSensitivityHU" = c(-50, 1200),
    "BMDL3Values.BMDL3StandardHU"        = c(-50, 1200),
    
    "BMDL1Values.BMDL1HighSensitivityHU" = c(-50, 1200),
    "BMDL1Values.BMDL1StandardHU"        = c(-50, 1200)
  )
  
  for (col in names(qc_rules)) {
    if (col %in% names(bodycomp)) {
      range <- qc_rules[[col]]
      bodycomp[[col]] <- ifelse(bodycomp[[col]] < range[1] | bodycomp[[col]] > range[2], NA, bodycomp[[col]])
    }
  }
  return(bodycomp)
}