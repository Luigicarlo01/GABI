####################################################################################################
# GIS-FA MAIZE
# ASREML FA(3) + PLS
# LOEO for all environments
# ncomp = 1:10
#
# Evaluation:
#   1. BLUE observed in TargetEnv
#   2. BLUP from FULL FA
#
# Outputs:
#   - Metrics per environment
#   - Mean and SD across environments
#   - Pooled metrics
#   - Individual predictions
####################################################################################################

rm(list = ls())

options(echo = TRUE)


####################################################################################################
# 0. PACKAGES
####################################################################################################

library(dplyr)
library(asreml)
library(pls)


####################################################################################################
# 0.1 ASREML OPTIONS
####################################################################################################

asreml.options(
  maxit = 1000,
  workspace = "8gb",
  pworkspace = "16gb",
  ai.penalty = 50,
  ai.loadings = 10,
  step.size = 0.05,
  update.step.size = 0.05
)


####################################################################################################
# 0.2 FA OUTPUT FUNCTION
####################################################################################################

source(
  "https://raw.githubusercontent.com/saulo-chaves/May_b_useful/refs/heads/main/fa_outs.R"
)


####################################################################################################
# 0.3 OUTPUT DIRECTORY
####################################################################################################

dir.create(
  "resMod",
  showWarnings = FALSE,
  recursive = TRUE
)


####################################################################################################
# 1. INPUT DATA
####################################################################################################


# --------------------------------------------------------------------------------------------------
# Environmental covariates
# --------------------------------------------------------------------------------------------------

tmp_covamb = tempfile(
  fileext = ".RData"
)

download.file(
  "https://raw.githubusercontent.com/Luigicarlo01/GABI/main/covamb.RData",
  destfile = tmp_covamb,
  mode = "wb",
  quiet = TRUE
)

load(
  tmp_covamb
)

unlink(
  tmp_covamb
)


# --------------------------------------------------------------------------------------------------
# Phenotypic data
# --------------------------------------------------------------------------------------------------

data = read.csv2(
  "https://raw.githubusercontent.com/Luigicarlo01/GABI/main/data_EmbI.csv",
  sep = ","
)


# --------------------------------------------------------------------------------------------------
# Genomic relationship matrix
# --------------------------------------------------------------------------------------------------

tmp_GK = tempfile(
  fileext = ".RData"
)

download.file(
  "https://raw.githubusercontent.com/Luigicarlo01/GABI/main/GK.RData",
  destfile = tmp_GK,
  mode = "wb",
  quiet = TRUE
)

load(
  tmp_GK
)

unlink(
  tmp_GK
)


# --------------------------------------------------------------------------------------------------
# BLUEs
# --------------------------------------------------------------------------------------------------

Gblues = read.table(
  "https://raw.githubusercontent.com/Luigicarlo01/GABI/main/Gblues.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)


# --------------------------------------------------------------------------------------------------
# Objects used in analysis
# --------------------------------------------------------------------------------------------------

covamb = matcovamb

EnvMarkers = as.data.frame(
  covamb
)

bluesYT = as.data.frame(
  Gblues
)

relMat = as.matrix(
  GK
)


colnames(bluesYT) = c(
  "Unit",
  "GenID",
  "YT"
)


####################################################################################################
# 2. ENVIRONMENTS
####################################################################################################

envSet = c(
  "AL_E01", "AL_E02",
  "BE_E03", "BE_E04", "BE_E05",
  "BI_E06", "BI_E07", "BI_E08",
  "BV_E09",
  "CM_E10", "CM_E11", "CM_E12", "CM_E13",
  "CQ_E14",
  "GO_E15", "GO_E16",
  "IP_E17",
  "LO_E18", "LO_E19",
  "MA_E20",
  "MR_E21", "MR_E22", "MR_E23", "MR_E24",
  "PA_E25", "PA_E26", "PA_E27",
  "PT_E28", "PT_E29",
  "PL_E30", "PL_E31", "PL_E32",
  "PG_E33",
  "RB_E34",
  "SR_E35", "SR_E36", "SR_E37", "SR_E38",
  "SP_E39",
  "SL_E40", "SL_E41", "SL_E42", "SL_E43",
  "SI_E44",
  "SO_E45",
  "TE_E46",
  "VI_E47"
)


env_map = setNames(
  envSet,
  sub(
    ".*_",
    "",
    envSet
  )
)


####################################################################################################
# 3. PHENOTYPIC DATA
####################################################################################################

myY = data[, c(
  10,
  11,
  15,
  16,
  17,
  19
)]


colnames(myY)[5] =
  "GenID"


myY$EnvID = env_map[
  sub(
    ".*_",
    "",
    myY$env
  )
]


myY = myY %>%
  na.omit() %>%
  filter(
    EnvID %in% envSet
  ) %>%
  droplevels()


####################################################################################################
# 4. BLUE DATA
####################################################################################################

bluesYT$EnvID = env_map[
  bluesYT$Unit
]


bluesYTdf = bluesYT %>%
  na.omit() %>%
  filter(
    EnvID %in% envSet
  ) %>%
  mutate(
    
    GenID =
      as.character(
        GenID
      ),
    
    EnvID =
      as.character(
        EnvID
      )
  )


####################################################################################################
# 5. GENOMIC RELATIONSHIP MATRIX
####################################################################################################

relMat = relMat[
  rownames(relMat) %in%
    bluesYTdf$GenID,
  
  colnames(relMat) %in%
    bluesYTdf$GenID,
  
  drop = FALSE
]


relMat = relMat[
  order(
    rownames(relMat)
  ),
  
  order(
    colnames(relMat)
  ),
  
  drop = FALSE
]


myY = myY[
  myY$GenID %in%
    colnames(relMat),
  ,
  drop = FALSE
] %>%
  droplevels()


bluesYTdf = bluesYTdf[
  bluesYTdf$GenID %in%
    colnames(relMat),
  ,
  drop = FALSE
]


####################################################################################################
# 6. ENVIRONMENTAL COVARIATES
####################################################################################################

EnvMarkers = as.data.frame(
  covamb
)


env_original = rownames(
  EnvMarkers
)


# --------------------------------------------------------------------------------------------------
# Convert environment names to E01, E02, ..., E47
#
# Supported examples:
#
# E01
# EMB1_01
# AL_E01
# --------------------------------------------------------------------------------------------------

env_code = ifelse(
  
  grepl(
    "^E[0-9]{2}$",
    env_original
  ),
  
  env_original,
  
  ifelse(
    
    grepl(
      "_E[0-9]{2}$",
      env_original
    ),
    
    sub(
      ".*_",
      "",
      env_original
    ),
    
    paste0(
      "E",
      sub(
        ".*_",
        "",
        env_original
      )
    )
  )
)


# --------------------------------------------------------------------------------------------------
# Convert:
#
# E01 -> AL_E01
# E02 -> AL_E02
# ...
# E47 -> VI_E47
# --------------------------------------------------------------------------------------------------

new_names = unname(
  env_map[
    env_code
  ]
)


rownames(EnvMarkers) =
  new_names


# --------------------------------------------------------------------------------------------------
# Remove unmapped environments
# --------------------------------------------------------------------------------------------------

EnvMarkers = EnvMarkers[
  !is.na(
    rownames(
      EnvMarkers
    )
  ),
  ,
  drop = FALSE
]


# --------------------------------------------------------------------------------------------------
# Keep only environments used in analysis
# --------------------------------------------------------------------------------------------------

EnvMarkers = EnvMarkers[
  rownames(EnvMarkers) %in%
    envSet,
  ,
  drop = FALSE
]


# --------------------------------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------------------------------

stopifnot(
  nrow(EnvMarkers) ==
    length(envSet)
)


stopifnot(
  all(
    envSet %in%
      rownames(EnvMarkers)
  )
)


# --------------------------------------------------------------------------------------------------
# Exact order
# --------------------------------------------------------------------------------------------------

EnvMarkers = EnvMarkers[
  envSet,
  ,
  drop = FALSE
]


stopifnot(
  identical(
    rownames(EnvMarkers),
    envSet
  )
)


####################################################################################################
# 7. FULL FA(3)
####################################################################################################

myY$EnvID = factor(
  myY$EnvID,
  levels = envSet
)


full_data = droplevels(
  myY
)


full_data$GenID = factor(
  full_data$GenID,
  levels = colnames(
    relMat
  )
)


assign(
  "full_data_asreml",
  full_data,
  envir = .GlobalEnv
)


assign(
  "relMat_asreml",
  relMat,
  envir = .GlobalEnv
)


FA_full = asreml(
  
  fixed =
    yield ~
    rep:EnvID +
    EnvID,
  
  random =
    ~ fa(
      EnvID,
      3
    ):vm(
      GenID,
      relMat_asreml
    ) +
    diag(
      EnvID
    ):block,
  
  residual =
    ~ dsum(
      ~ id(units) | EnvID
    ),
  
  data =
    full_data_asreml,
  
  maxiter =
    1000
)


####################################################################################################
# 7.1 UPDATE FULL FA
####################################################################################################

FA_full = update(
  FA_full
)


####################################################################################################
# 7.2 FULL FA OUTPUTS
####################################################################################################

fa_full = fa.outs(
  FA_full,
  name.env = "EnvID",
  name.gen = "GenID"
)


####################################################################################################
# 8. FULL-FA BLUP REFERENCE
####################################################################################################

bl = fa_full$blups


nm = names(
  bl
)


env_col = nm[
  grepl(
    "^env$|^EnvID$|environment",
    nm,
    ignore.case = TRUE
  )
][1]


gen_col = nm[
  grepl(
    "^gen$|^GenID$|genotype",
    nm,
    ignore.case = TRUE
  )
][1]


mar_col = nm[
  grepl(
    "marginal",
    nm,
    ignore.case = TRUE
  )
][1]


BLUP_full = bl %>%
  transmute(
    
    EnvID =
      as.character(
        .data[[env_col]]
      ),
    
    GenID =
      as.character(
        .data[[gen_col]]
      ),
    
    BLUP_ref =
      as.numeric(
        .data[[mar_col]]
      )
  )


####################################################################################################
# 9. LOEO
####################################################################################################

results_list = list()

prediction_list = list()


counter_results = 1

counter_predictions = 1


####################################################################################################
# ONLY FA1, FA2 AND FA3
####################################################################################################

fa_names = c(
  "fa1",
  "fa2",
  "fa3"
)


####################################################################################################
# LOEO LOOP
####################################################################################################

for(TargetEnv in envSet){
  
  
  ##################################################################################################
  # 9.1 TRAINING PHENOTYPES
  ##################################################################################################
  
  train_data = myY %>%
    filter(
      EnvID != TargetEnv
    ) %>%
    droplevels()
  
  
  train_data$GenID = factor(
    train_data$GenID,
    levels = colnames(
      relMat
    )
  )
  
  
  assign(
    "train_data_asreml",
    train_data,
    envir = .GlobalEnv
  )
  
  
  ##################################################################################################
  # 9.2 FA(3) WITHOUT TargetEnv
  ##################################################################################################
  
  FA_train = asreml(
    
    fixed =
      yield ~
      rep:EnvID +
      EnvID,
    
    random =
      ~ fa(
        EnvID,
        3
      ):vm(
        GenID,
        relMat_asreml
      ) +
      diag(
        EnvID
      ):block,
    
    residual =
      ~ dsum(
        ~ id(units) | EnvID
      ),
    
    data =
      train_data_asreml,
    
    maxiter =
      1000
  )
  
  
  ##################################################################################################
  # 9.3 UPDATE FA
  ##################################################################################################
  
  FA_train = update(
    FA_train
  )
  
  
  ##################################################################################################
  # 9.4 FA OUTPUTS
  ##################################################################################################
  
  fa_res = fa.outs(
    FA_train,
    name.env = "EnvID",
    name.gen = "GenID"
  )
  
  
  GenScores = as.matrix(
    fa_res$rot.scores
  )
  
  
  EnvLoadings = as.matrix(
    fa_res$rot.loads
  )
  
  
  ##################################################################################################
  # 9.5 ENVIRONMENTAL TRAINING AND TEST SETS
  ##################################################################################################
  
  X_train = EnvMarkers[
    rownames(
      EnvMarkers
    ) != TargetEnv,
    ,
    drop = FALSE
  ]
  
  
  X_test = EnvMarkers[
    rownames(
      EnvMarkers
    ) == TargetEnv,
    ,
    drop = FALSE
  ]
  
  
  ##################################################################################################
  # 9.6 ALIGN ENVIRONMENTAL COVARIATES AND FA LOADINGS
  ##################################################################################################
  
  common_env = intersect(
    rownames(
      EnvLoadings
    ),
    rownames(
      X_train
    )
  )
  
  
  EnvLoadings = EnvLoadings[
    common_env,
    ,
    drop = FALSE
  ]
  
  
  X_train = X_train[
    common_env,
    ,
    drop = FALSE
  ]
  
  
  ##################################################################################################
  # 9.7 REMOVE ENVIRONMENTAL VARIABLES WITHOUT VARIATION
  ##################################################################################################
  
  train_mean = sapply(
    X_train,
    mean,
    na.rm = TRUE
  )
  
  
  train_sd = sapply(
    X_train,
    sd,
    na.rm = TRUE
  )
  
  
  keep_variables =
    is.finite(
      train_mean
    ) &
    is.finite(
      train_sd
    ) &
    train_sd > 0
  
  
  X_train = X_train[
    ,
    keep_variables,
    drop = FALSE
  ]
  
  
  X_test = X_test[
    ,
    keep_variables,
    drop = FALSE
  ]
  
  
  train_mean = train_mean[
    keep_variables
  ]
  
  
  train_sd = train_sd[
    keep_variables
  ]
  
  
  ##################################################################################################
  # 9.8 MISSING VALUES
  ##################################################################################################
  
  for(
    j in seq_len(
      ncol(
        X_train
      )
    )
  ){
    
    X_train[
      is.na(
        X_train[, j]
      ),
      j
    ] =
      train_mean[j]
    
    
    X_test[
      is.na(
        X_test[, j]
      ),
      j
    ] =
      train_mean[j]
  }
  
  
  ##################################################################################################
  # 9.9 STANDARDIZATION BASED ONLY ON TRAINING ENVIRONMENTS
  ##################################################################################################
  
  X_train = sweep(
    as.matrix(
      X_train
    ),
    2,
    train_mean,
    "-"
  )
  
  
  X_train = sweep(
    X_train,
    2,
    train_sd,
    "/"
  )
  
  
  X_test = sweep(
    as.matrix(
      X_test
    ),
    2,
    train_mean,
    "-"
  )
  
  
  X_test = sweep(
    X_test,
    2,
    train_sd,
    "/"
  )
  
  
  ##################################################################################################
  # 9.10 KEEP FA1-FA3
  ##################################################################################################
  
  EnvLoadings = EnvLoadings[
    ,
    fa_names,
    drop = FALSE
  ]
  
  
  GenScores = GenScores[
    ,
    fa_names,
    drop = FALSE
  ]
  
  
  ##################################################################################################
  # 9.11 NUMBER OF PLS COMPONENTS
  ##################################################################################################
  
  max_ncomp = min(
    10,
    nrow(
      X_train
    ) - 1,
    ncol(
      X_train
    )
  )
  
  
  ##################################################################################################
  # 9.12 PLS DATA
  ##################################################################################################
  
  data_pls = data.frame(
    
    lambda =
      I(
        EnvLoadings
      ),
    
    CovAmb =
      I(
        X_train
      )
  )
  
  
  datTest = data.frame(
    
    CovAmb =
      I(
        X_test
      )
  )
  
  
  ##################################################################################################
  # 9.13 BLUE REFERENCE
  ##################################################################################################
  
  BLUE_target = bluesYTdf %>%
    filter(
      EnvID == TargetEnv
    ) %>%
    transmute(
      
      GenID =
        as.character(
          GenID
        ),
      
      BLUE_ref =
        as.numeric(
          YT
        )
    )
  
  
  ##################################################################################################
  # 9.14 FULL-FA BLUP REFERENCE
  ##################################################################################################
  
  BLUP_target = BLUP_full %>%
    filter(
      EnvID == TargetEnv
    ) %>%
    select(
      GenID,
      BLUP_ref
    )
  
  
  ##################################################################################################
  # 9.15 PLS GRID
  ##################################################################################################
  
  for(nc in 1:max_ncomp){
    
    
    ################################################################################################
    # PLS
    ################################################################################################
    
    pls_fit = pls::plsr(
      
      lambda ~ CovAmb,
      
      data =
        data_pls,
      
      ncomp =
        nc
    )
    
    
    ################################################################################################
    # Predict TargetEnv loadings
    ################################################################################################
    
    pred_raw = predict(
      pls_fit,
      newdata = datTest,
      ncomp = nc
    )
    
    
    pred_pls = matrix(
      drop(
        pred_raw
      ),
      nrow = 1
    )
    
    
    colnames(pred_pls) =
      fa_names
    
    
    rownames(pred_pls) =
      TargetEnv
    
    
    ################################################################################################
    # Reconstruct genotype values
    ################################################################################################
    
    EBV =
      pred_pls[
        1,
        fa_names,
        drop = FALSE
      ] %*%
      t(
        GenScores[
          ,
          fa_names,
          drop = FALSE
        ]
      )
    
    
    EBV = as.numeric(
      EBV
    )
    
    
    predictions = data.frame(
      
      Env =
        TargetEnv,
      
      ncomp =
        nc,
      
      GenID =
        as.character(
          rownames(
            GenScores
          )
        ),
      
      EBV =
        EBV
    )
    
    
    ################################################################################################
    # Join BLUE and BLUP references
    ################################################################################################
    
    ModelPred = predictions %>%
      left_join(
        BLUE_target,
        by = "GenID"
      ) %>%
      left_join(
        BLUP_target,
        by = "GenID"
      )
    
    
    ################################################################################################
    # BLUE METRICS
    ################################################################################################
    
    N_BLUE = sum(
      complete.cases(
        ModelPred$BLUE_ref,
        ModelPred$EBV
      )
    )
    
    
    Cor_BLUE = cor(
      ModelPred$BLUE_ref,
      ModelPred$EBV,
      use = "pairwise.complete.obs"
    )
    
    
    RMSE_BLUE = sqrt(
      mean(
        (
          ModelPred$BLUE_ref -
            ModelPred$EBV
        )^2,
        na.rm = TRUE
      )
    )
    
    
    ################################################################################################
    # BLUP METRICS
    ################################################################################################
    
    N_BLUP = sum(
      complete.cases(
        ModelPred$BLUP_ref,
        ModelPred$EBV
      )
    )
    
    
    Cor_BLUP = cor(
      ModelPred$BLUP_ref,
      ModelPred$EBV,
      use = "pairwise.complete.obs"
    )
    
    
    RMSE_BLUP = sqrt(
      mean(
        (
          ModelPred$BLUP_ref -
            ModelPred$EBV
        )^2,
        na.rm = TRUE
      )
    )
    
    
    ################################################################################################
    # STORE RESULTS
    ################################################################################################
    
    results_list[[counter_results]] = data.frame(
      
      Env =
        TargetEnv,
      
      Factor =
        "FA1-FA3",
      
      ncomp =
        nc,
      
      N_BLUE =
        N_BLUE,
      
      Cor_BLUE =
        Cor_BLUE,
      
      RMSE_BLUE =
        RMSE_BLUE,
      
      N_BLUP =
        N_BLUP,
      
      Cor_BLUP =
        Cor_BLUP,
      
      RMSE_BLUP =
        RMSE_BLUP
    )
    
    
    prediction_list[[counter_predictions]] =
      ModelPred
    
    
    counter_results =
      counter_results + 1
    
    
    counter_predictions =
      counter_predictions + 1
  }
  
  
  ##################################################################################################
  # CLEAN MEMORY
  ##################################################################################################
  
  gc()
}


####################################################################################################
# 10. COMBINE RESULTS
####################################################################################################

results_all = bind_rows(
  results_list
)


predictions_all = bind_rows(
  prediction_list
)


####################################################################################################
# 11. RESULTS PER ENVIRONMENT
####################################################################################################

write.table(
  results_all,
  file =
    "resMod/GISFA_FA3_PLS_per_environment.tsv",
  sep =
    "\t",
  quote =
    FALSE,
  row.names =
    FALSE
)


####################################################################################################
# 12. INDIVIDUAL PREDICTIONS
####################################################################################################

write.table(
  predictions_all,
  file =
    "resMod/GISFA_FA3_PLS_predictions.tsv",
  sep =
    "\t",
  quote =
    FALSE,
  row.names =
    FALSE
)


####################################################################################################
# 13. MEAN AND SD ACROSS ENVIRONMENTS
####################################################################################################

summary_environment = results_all %>%
  group_by(
    ncomp
  ) %>%
  summarise(
    
    N_Envs =
      n_distinct(
        Env
      ),
    
    Mean_Cor_BLUE =
      mean(
        Cor_BLUE,
        na.rm = TRUE
      ),
    
    SD_Cor_BLUE =
      sd(
        Cor_BLUE,
        na.rm = TRUE
      ),
    
    Mean_RMSE_BLUE =
      mean(
        RMSE_BLUE,
        na.rm = TRUE
      ),
    
    SD_RMSE_BLUE =
      sd(
        RMSE_BLUE,
        na.rm = TRUE
      ),
    
    Mean_Cor_BLUP =
      mean(
        Cor_BLUP,
        na.rm = TRUE
      ),
    
    SD_Cor_BLUP =
      sd(
        Cor_BLUP,
        na.rm = TRUE
      ),
    
    Mean_RMSE_BLUP =
      mean(
        RMSE_BLUP,
        na.rm = TRUE
      ),
    
    SD_RMSE_BLUP =
      sd(
        RMSE_BLUP,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


write.table(
  summary_environment,
  file =
    "resMod/GISFA_FA3_PLS_mean_SD.tsv",
  sep =
    "\t",
  quote =
    FALSE,
  row.names =
    FALSE
)


####################################################################################################
# 14. POOLED METRICS
####################################################################################################

summary_pooled = predictions_all %>%
  group_by(
    ncomp
  ) %>%
  summarise(
    
    N_BLUE =
      sum(
        complete.cases(
          BLUE_ref,
          EBV
        )
      ),
    
    Cor_BLUE_Pooled =
      cor(
        BLUE_ref,
        EBV,
        use =
          "pairwise.complete.obs"
      ),
    
    RMSE_BLUE_Pooled =
      sqrt(
        mean(
          (
            BLUE_ref -
              EBV
          )^2,
          na.rm = TRUE
        )
      ),
    
    N_BLUP =
      sum(
        complete.cases(
          BLUP_ref,
          EBV
        )
      ),
    
    Cor_BLUP_Pooled =
      cor(
        BLUP_ref,
        EBV,
        use =
          "pairwise.complete.obs"
      ),
    
    RMSE_BLUP_Pooled =
      sqrt(
        mean(
          (
            BLUP_ref -
              EBV
          )^2,
          na.rm = TRUE
        )
      ),
    
    .groups =
      "drop"
  )


write.table(
  summary_pooled,
  file =
    "resMod/GISFA_FA3_PLS_pooled.tsv",
  sep =
    "\t",
  quote =
    FALSE,
  row.names =
    FALSE
)


####################################################################################################
# END
####################################################################################################
