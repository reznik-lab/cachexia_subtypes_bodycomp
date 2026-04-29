# AUTHOR: Sonia Boscenco
# Function to calculate body composition changes across cachexia
calculate_body_composition_deltas <- function(df) {
  df %>%
    group_by(MRN) %>%
    reframe(
      
      # densities
      delta_MuscleDensity = MuscleValues.L3TotalMuscleMedianHU[SCAN_TYPE == "Post-Cachexia"] - MuscleValues.L3TotalMuscleMedianHU[SCAN_TYPE == "Start-Cachexia"],
      delta_SATDensity    = L3FatValues.L3SATMedian[SCAN_TYPE == "Post-Cachexia"] - L3FatValues.L3SATMedian[SCAN_TYPE == "Start-Cachexia"],
      delta_VATDensity    = L3FatValues.L3VATMedian[SCAN_TYPE == "Post-Cachexia"] - L3FatValues.L3VATMedian[SCAN_TYPE == "Start-Cachexia"],
      delta_KidneyDensity = KidneyValues.KidneyMedianHU[SCAN_TYPE == "Post-Cachexia"] - KidneyValues.KidneyMedianHU[SCAN_TYPE == "Start-Cachexia"],
      delta_PancreasDensity = PancreasValues.PancreasMedianHU[SCAN_TYPE == "Post-Cachexia"] - PancreasValues.PancreasMedianHU[SCAN_TYPE == "Start-Cachexia"],
      delta_SpleenDensity   = SpleenValues.SpleenMedianHU[SCAN_TYPE == "Post-Cachexia"] - SpleenValues.SpleenMedianHU[SCAN_TYPE == "Start-Cachexia"],
      delta_LiverDensity    = LiverValues.LiverMedianHU[SCAN_TYPE == "Post-Cachexia"] - LiverValues.LiverMedianHU[SCAN_TYPE == "Start-Cachexia"],
      
      # volumes
      delta_LiverArea    = (LiverValues.LiverVolume[SCAN_TYPE == "Post-Cachexia"] - LiverValues.LiverVolume[SCAN_TYPE == "Start-Cachexia"]) / LiverValues.LiverVolume[SCAN_TYPE == "Start-Cachexia"] * 100,
      delta_KidneyVolume = (KidneyValues.KidneyVolume[SCAN_TYPE == "Post-Cachexia"] - KidneyValues.KidneyVolume[SCAN_TYPE == "Start-Cachexia"]) / KidneyValues.KidneyVolume[SCAN_TYPE == "Start-Cachexia"] * 100,
      delta_PancreasVolume = (PancreasValues.PancreasVolume[SCAN_TYPE == "Post-Cachexia"] - PancreasValues.PancreasVolume[SCAN_TYPE == "Start-Cachexia"]) / PancreasValues.PancreasVolume[SCAN_TYPE == "Start-Cachexia"] * 100,
      delta_SpleenVolume   = (SpleenValues.SpleenVolume[SCAN_TYPE == "Post-Cachexia"] - SpleenValues.SpleenVolume[SCAN_TYPE == "Start-Cachexia"]) / SpleenValues.SpleenVolume[SCAN_TYPE == "Start-Cachexia"] * 100,
      
      # log ratio
      log_LiverVolume    = log(LiverValues.LiverVolume[SCAN_TYPE == "Post-Cachexia"]) - log(LiverValues.LiverVolume[SCAN_TYPE == "Start-Cachexia"]),
      log_KidneyVolume = log(KidneyValues.KidneyVolume[SCAN_TYPE == "Post-Cachexia"]) - log(KidneyValues.KidneyVolume[SCAN_TYPE == "Start-Cachexia"]),
      log_PancreasVolume = log(PancreasValues.PancreasVolume[SCAN_TYPE == "Post-Cachexia"]) - log(PancreasValues.PancreasVolume[SCAN_TYPE == "Start-Cachexia"]),
      log_SpleenVolume   = log(SpleenValues.SpleenVolume[SCAN_TYPE == "Post-Cachexia"]) - log(SpleenValues.SpleenVolume[SCAN_TYPE == "Start-Cachexia"]),
      log_MuscleArea = log(MuscleValues.L3TotalMuscleArea[SCAN_TYPE == "Post-Cachexia"]) - log(MuscleValues.L3TotalMuscleArea[SCAN_TYPE == "Start-Cachexia"]),
      log_SAT        = log(L3FatValues.L3SATArea[SCAN_TYPE == "Post-Cachexia"]) - log(L3FatValues.L3SATArea[SCAN_TYPE == "Start-Cachexia"]),
      log_VAT        = log(L3FatValues.L3VATArea[SCAN_TYPE == "Post-Cachexia"]) - log(L3FatValues.L3VATArea[SCAN_TYPE == "Start-Cachexia"]),
      log_IMAT       = log(MuscleValues.L3IMATArea[SCAN_TYPE == "Post-Cachexia"]) - log(MuscleValues.L3IMATArea[SCAN_TYPE == "Start-Cachexia"]),
      
      # other
      delta_BMDLStandard  = BMDL3Values.BMDL3StandardHU[SCAN_TYPE == "Post-Cachexia"] - BMDL3Values.BMDL3StandardHU[SCAN_TYPE == "Start-Cachexia"],
      delta_TotalAortic   = CalciumScoring.TotalAorticAgatston[SCAN_TYPE == "Post-Cachexia"] - CalciumScoring.TotalAorticAgatston[SCAN_TYPE == "Start-Cachexia"],
      
      # areas
      delta_MuscleArea = (MuscleValues.L3TotalMuscleArea[SCAN_TYPE == "Post-Cachexia"] - MuscleValues.L3TotalMuscleArea[SCAN_TYPE == "Start-Cachexia"]) / MuscleValues.L3TotalMuscleArea[SCAN_TYPE == "Start-Cachexia"] * 100,
      delta_SAT        = (L3FatValues.L3SATArea[SCAN_TYPE == "Post-Cachexia"] - L3FatValues.L3SATArea[SCAN_TYPE == "Start-Cachexia"]) / L3FatValues.L3SATArea[SCAN_TYPE == "Start-Cachexia"] * 100,
      delta_VAT        = (L3FatValues.L3VATArea[SCAN_TYPE == "Post-Cachexia"] - L3FatValues.L3VATArea[SCAN_TYPE == "Start-Cachexia"]) / L3FatValues.L3VATArea[SCAN_TYPE == "Start-Cachexia"] * 100,
      delta_IMAT       = (MuscleValues.L3IMATArea[SCAN_TYPE == "Post-Cachexia"] - MuscleValues.L3IMATArea[SCAN_TYPE == "Start-Cachexia"]) / (MuscleValues.L3IMATArea[SCAN_TYPE == "Start-Cachexia"] + 1) * 100,
      
      # time
      scans_days = as.Date(SCAN_DATE[SCAN_TYPE == "Post-Cachexia"]) - as.Date(SCAN_DATE[SCAN_TYPE == "Start-Cachexia"])
    )
}