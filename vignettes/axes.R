ax <- spec_axes(
  outcome    = c(cdc_any = "y_cdc_any", cdc_severe = "y_cdc_severe", efp3 = "y_efp3", efp2 = "y_efp2"),
  covariates = list(crude = character(0),
                    demo  = c("Age", "Gender", "Race"),
                    socio = c("Age", "Gender", "Race", "Educational_level", "PIR", "Marriage_ststus",
                              "Smoking_status", "Drinking_status", "Flossing"),
                    full  = c("Age", "Gender", "Race", "Educational_level", "PIR", "Marriage_ststus",
                              "Smoking_status", "Drinking_status", "Flossing", "BMI", "Diabetes",
                              "Hypertension", "Cardiovascular_disease", "Hyperlipidemia", "Arthritis", "HbA1c")),
  coding     = c("per_sd", "log2"),
  sample     = list(all = TRUE, age40 = quote(Age >= 40)),
  weighting  = "design")
