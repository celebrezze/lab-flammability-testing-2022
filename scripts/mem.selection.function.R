
# MINOR FUNCTIONS

# Mallow's CP
# note: opted not to really use this in this model selection, still useful to include in case we want to consider it in future model selections
mallows.cp <- function(model, k, n) {
  return(AIC(model) + 2*k / (n - k - 1))
}

### Model Performance Function (AIC, BIC, Nakagawa's R2 values)
mem.performance <- function(model){
  aic <- AIC(model)
  bic <- BIC(model)
  cond.r2 <- r2_nakagawa(model)[[1]]
  marg.r2 <- r2_nakagawa(model)[[2]]
  print(paste('AIC =', aic, 'BIC =', bic, 'Cond. Rsq =', cond.r2, 'Marg. Rsq =', marg.r2))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRIMARY MIXED EFFECTS MODEL SELECTION FUNCTION (JUST FIXED EFFECTS, NO INTERACTIONS)
# NOTE: y.var and rem.str must be a character vector, predictors must be a dataframe
mem.selection <- function(y.var, predictors, df, rem.str = '(1|species/plant_id)'){
  # FIRST: List models
  # empty list
  mod.list <- list()
  call.vec <- c()
  
  # setting up a nested for loop to list all possible models including those predictors
  for(i in 2:ncol(predictors)){
    call <- colnames(predictors) %>%
      combinations(n = ncol(predictors), r = i, repeats.allowed = F) %>% 
      apply(1, paste0, collapse = ' + ') # all possible combinations of models from 2 - # of predictors we're interested in (note: cannot include only 1 predictor, as this breaks multicollinearity for loop below; also, we're probably not interested in a model with only one explanatory variable)
    for(j in (1+length(call.vec)):(length(call)+length(call.vec))){ # adding linear mixed effects model to mod.list
      mod.list[[j]] <- lmer(as.formula(paste(y.var, '~', call[j-length(call.vec)], '+', rem.str, sep = '')), REML = F, data = df)
    }
    call.vec <- append(call.vec, call) # to index where to put model into model list
  }
  
  # CREATE DATAFRAME
  model.df <- data.frame(model.id = c(1:length(mod.list)),
                         call = c(rep(NA, length(mod.list))),
                         multicollinearity = c(rep(NA, length(mod.list))),
                         AIC = c(rep(NA, length(mod.list))),
                         BIC = c(rep(NA, length(mod.list))),
                         CP = c(rep(NA, length(mod.list))),
                         Cond_R2 = c(rep(NA, length(mod.list))),
                         Marg_R2 = c(rep(NA, length(mod.list))),
                         formula = c(rep(NA, length(mod.list))))
  
  # SECOND: FLAG ANY MODELS THAT BREAK MULTICOLINEARITY ASSUMPTIONS
  for(i in 1:length(mod.list)){
    vif.df <- multicollinearity(mod.list[[i]])
    ifelse(max(vif.df$VIF) > 5, model.df$multicollinearity[i] <- 'yes', 
           model.df$multicollinearity[i] <- 'no')
  }
  
  for(i in 1:length(mod.list)){
    if(model.df$multicollinearity[i] == 'yes'){
      mod.list[[i]] <- NA
    }}
  
  # THIRD: CALCULATE AIC, BIC, MALLOWS CP
  for(i in 1:length(mod.list)) {
    if(model.df$multicollinearity[i] == 'yes')
    {model.df[i, 3:6] <- NA}
    else{
      model.df$call[i] <- paste(colnames(mod.list[[i]]@frame), collapse = ', ')
      model.df$AIC[i] <- AIC(mod.list[[i]])
      model.df$BIC[i] <- BIC(mod.list[[i]])
      model.df$Cond_R2[i] <- r2_nakagawa(mod.list[[i]])[[1]]
      model.df$Marg_R2[i] <- r2_nakagawa(mod.list[[i]])[[2]]
      model.df$CP[i] <- mallows.cp(mod.list[[i]], k = length(mod.list[[i]]@beta - 1), n = nrow(df))
      model.df$formula[i] <- as.character(mod.list[[i]]@call)[2]
    }
  }
  output <- list(model.df, mod.list)
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODEL SELECTION TABLE

mem.selection.table <- function(df, mod.list, output){
  kable.df <- df %>% 
    na.omit() %>% 
    filter(AIC < (min(AIC)+2) | BIC < (min(BIC)+2)) %>%
    arrange(AIC) %>% 
    select(model.id, call, AIC, BIC)
  
  for(i in 1:nrow(kable.df)){
    kable.df$cond.r2[i] <- r2_nakagawa(mod.list[[kable.df$model.id[i]]])[[1]]
    kable.df$marg.r2[i] <- r2_nakagawa(mod.list[[kable.df$model.id[i]]])[[2]]
  }
  kable.df %>% 
    select(-model.id) %>% 
    kable(format = 'html', escape = F, col.names = c('Model Structure', 'AIC', 'BIC', 'Conditional Rsq', 'Marginal Rsq')) %>% 
    kable_styling(bootstrap_options = c('hover', 'bordered', 'condensed'), fixed_thead = T) %>% 
    save_kable(here('figures', 'MEM', output))
  return(kable.df)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# UNIVARIATE MIXED EFFECTS MODEL SELECTION
# predictors must be character vector
mem.univariate.selection <- function(y.var, predictors, df, rem.str = '(1|species/plant_id)'){
  # FIRST: List models
  # Create an empty list to store model results
  mod.list <- list()
  
  # Loop through the interaction combinations
  for (i in 1:length(predictors)) {
    # Create the formula for the model (e.g., outcome ~ predictor1 * predictor2)
    formula <- as.formula(paste(y.var, "~", predictors[i], "+", rem.str))
    
    # Fit the mixed-effects model
    model <- lmer(formula, REML = F, data = mem.df)
    
    # Store the model result in the list
    mod.list[[i]] <- model
  }
  
  # CREATE DATAFRAME
  model.df <- data.frame(model.id = c(1:length(mod.list)),
                         call = c(rep(NA, length(mod.list))),
                         multicollinearity = c(rep(NA, length(mod.list))),
                         AIC = c(rep(NA, length(mod.list))),
                         BIC = c(rep(NA, length(mod.list))),
                         CP = c(rep(NA, length(mod.list))),
                         Cond_R2 = c(rep(NA, length(mod.list))),
                         Marg_R2 = c(rep(NA, length(mod.list))),
                         formula = c(rep(NA, length(mod.list))))
  
  # SECOND: CALCULATE AIC, BIC, MALLOWS CP
  for(i in 1:length(mod.list)) {
    model.df$call[i] <- paste(colnames(mod.list[[i]]@frame), collapse = ', ')
    model.df$AIC[i] <- AIC(mod.list[[i]])
    model.df$BIC[i] <- BIC(mod.list[[i]])
    model.df$Cond_R2[i] <- r2_nakagawa(mod.list[[i]])[[1]]
    model.df$Marg_R2[i] <- r2_nakagawa(mod.list[[i]])[[2]]
    model.df$CP[i] <- mallows.cp(mod.list[[i]], k = length(mod.list[[i]]@beta - 1), n = nrow(df))
    model.df$formula[i] <- as.character(mod.list[[i]]@call)[2]
  }
  output <- list(model.df, mod.list)
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FIRST ORDER INTERACTIONS MIXED EFFECTS MODEL SELECTION
# NOTE: y.var, predictors, and rem.str must be a character vector, df must be a dataframe
mem.firstorder.int.selection <- function(y.var, predictors, df, rem.str = '(1|species/plant_id)'){
  # FIRST: List models
  # Generate all possible combinations of two predictors for interaction
  interaction_combinations <- combn(predictors, 2, simplify = FALSE)
  
  # Create an empty list to store model results
  mod.list <- list()
  
  # Loop through the interaction combinations
  for (combo in interaction_combinations) {
    # Create the formula for the model (e.g., outcome ~ predictor1 * predictor2)
    formula <- as.formula(paste(y.var, "~", paste(combo, collapse = "*"), "+", rem.str))
    
    # Fit the mixed-effects model
    model <- lmer(formula, REML = F, data = mem.df)
    
    # Store the model result in the list
    mod.list[[paste(combo, collapse = "_")]] <- model
  }
  
  # CREATE DATAFRAME
  model.df <- data.frame(model.id = c(1:length(mod.list)),
                         call = c(rep(NA, length(mod.list))),
                         multicollinearity = c(rep(NA, length(mod.list))),
                         AIC = c(rep(NA, length(mod.list))),
                         BIC = c(rep(NA, length(mod.list))),
                         CP = c(rep(NA, length(mod.list))),
                         Cond_R2 = c(rep(NA, length(mod.list))),
                         Marg_R2 = c(rep(NA, length(mod.list))),
                         formula = c(rep(NA, length(mod.list))))
  
  # SECOND: FLAG ANY MODELS THAT BREAK MULTICOLINEARITY ASSUMPTIONS
  for(i in 1:length(mod.list)){
    vif.df <- multicollinearity(mod.list[[i]])
    ifelse(max(vif.df$VIF) > 5, model.df$multicollinearity[i] <- 'yes', 
           model.df$multicollinearity[i] <- 'no')
  }
  
  for(i in 1:length(mod.list)){
    if(model.df$multicollinearity[i] == 'yes'){
      mod.list[[i]] <- NA
    }}
  
  # THIRD: CALCULATE AIC, BIC, MALLOWS CP
  for(i in 1:length(mod.list)) {
    if(model.df$multicollinearity[i] == 'yes')
    {model.df[i, 3:6] <- NA}
    else{
      model.df$call[i] <- paste(colnames(mod.list[[i]]@frame), collapse = ', ')
      model.df$AIC[i] <- AIC(mod.list[[i]])
      model.df$BIC[i] <- BIC(mod.list[[i]])
      model.df$Cond_R2[i] <- r2_nakagawa(mod.list[[i]])[[1]]
      model.df$Marg_R2[i] <- r2_nakagawa(mod.list[[i]])[[2]]
      model.df$CP[i] <- mallows.cp(mod.list[[i]], k = length(mod.list[[i]]@beta - 1), n = nrow(df))
      model.df$formula[i] <- as.character(mod.list[[i]]@call)[2]
    }
  }
  output <- list(model.df, mod.list)
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FIRST ORDER INTERACTIONS AND OTHER FIXED EFFECTS MODEL SELECTION
# NOTE: y.var, interaction, and rem.str must be character vectors, predictors and df must be dataframes
mem.int.fixed.selection <- function(y.var, interaction, predictors, df, rem.str = '(1|species/plant_id)'){
  
  # FIRST: List models
  
  # Create an empty list to store model results
  mod.list <- list()
  call.vec <- c()
  
  # Remove terms from interaction term from predictors list so they are not duplicated in the function call; interaction terms are automatically included as single predictors as well as first-order interactions in lmer, so there is no need to have instances with lmer(y ~ a*b + a + b)
  predictors <- predictors %>% 
    select(-interaction)
  
  # Loop through the interaction combinations
  for(i in 1:ncol(predictors)){
    call <- colnames(predictors) %>%
      combinations(n = ncol(predictors), r = i, repeats.allowed = F) %>% 
      apply(1, paste0, collapse = ' + ') # all possible combinations of models from 2 - # of predictors we're interested in (note: cannot include only 1 predictor, as this breaks multicollinearity for loop below; also, we're probably not interested in a model with only one explanatory variable)
    for(j in (1+length(call.vec)):(length(call)+length(call.vec))){ # adding linear mixed effects model to mod.list
      mod.list[[j]] <- lmer(as.formula(paste(y.var, '~',  paste(interaction, collapse = "*"), '+', call[j-length(call.vec)], '+', rem.str, sep = '')), REML = F, data = df)
    }
    call.vec <- append(call.vec, call) # to index where to put model into model list
  }
  
  # CREATE DATAFRAME
  model.df <- data.frame(model.id = c(1:length(mod.list)),
                         call = c(rep(NA, length(mod.list))),
                         multicollinearity = c(rep(NA, length(mod.list))),
                         AIC = c(rep(NA, length(mod.list))),
                         BIC = c(rep(NA, length(mod.list))),
                         CP = c(rep(NA, length(mod.list))),
                         Cond_R2 = c(rep(NA, length(mod.list))),
                         Marg_R2 = c(rep(NA, length(mod.list))),
                         formula = c(rep(NA, length(mod.list))))
  
  # SECOND: FLAG ANY MODELS THAT BREAK MULTICOLINEARITY ASSUMPTIONS
  for(i in 1:length(mod.list)){
    vif.df <- multicollinearity(mod.list[[i]])
    ifelse(max(vif.df$VIF) > 5, model.df$multicollinearity[i] <- 'yes', 
           model.df$multicollinearity[i] <- 'no')
  }
  
  for(i in 1:length(mod.list)){
    if(model.df$multicollinearity[i] == 'yes'){
      mod.list[[i]] <- NA
    }}
  
  # THIRD: CALCULATE AIC, BIC, MALLOWS CP
  for(i in 1:length(mod.list)) {
    if(model.df$multicollinearity[i] == 'yes')
    {model.df[i, 3:6] <- NA}
    else{
      model.df$call[i] <- paste(colnames(mod.list[[i]]@frame), collapse = ', ')
      model.df$AIC[i] <- AIC(mod.list[[i]])
      model.df$BIC[i] <- BIC(mod.list[[i]])
      model.df$Cond_R2[i] <- r2_nakagawa(mod.list[[i]])[[1]]
      model.df$Marg_R2[i] <- r2_nakagawa(mod.list[[i]])[[2]]
      model.df$CP[i] <- mallows.cp(mod.list[[i]], k = length(mod.list[[i]]@beta - 1), n = nrow(df))
      model.df$formula[i] <- as.character(mod.list[[i]]@call)[2]
    }
  }
  output <- list(model.df, mod.list)
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GENERALIZED LINEAR MIXED EFFECTS MODEL SELECTION

glmm.selection <- function(y.var, predictors, df, mod.family, zero.inflation, rem.str = '(1|plant_id)'){
  # FIRST: List models
  # empty list
  mod.list <- list()
  call.vec <- c()
  
  # setting up a nested for loop to list all possible models including those predictors
  for(i in 2:ncol(predictors)){
    call <- colnames(predictors) %>%
      combinations(n = ncol(predictors), r = i, repeats.allowed = F) %>% 
      apply(1, paste0, collapse = ' + ') # all possible combinations of models from 2 - # of predictors we're interested in (note: cannot include only 1 predictor, as this breaks multicollinearity for loop below; also, we're probably not interested in a model with only one explanatory variable)
    for(j in (1+length(call.vec)):(length(call)+length(call.vec))){ # adding linear mixed effects model to mod.list
      mod.list[[j]] <- glmmTMB(as.formula(paste(y.var, '~', call[j-length(call.vec)], '+', rem.str, sep = '')), data = df, family = mod.family, ziformula = as.formula(paste0('~', zero.inflation, sep = '')))
    }
    call.vec <- append(call.vec, call) # to index where to put model into model list
  }
  
  # CREATE DATAFRAME
  model.df <- data.frame(model.id = c(1:length(mod.list)),
                         call = c(rep(NA, length(mod.list))),
                         multicollinearity = c(rep(NA, length(mod.list))),
                         AIC = c(rep(NA, length(mod.list))),
                         BIC = c(rep(NA, length(mod.list))),
                         R2_marg = c(rep(NA, length(mod.list))),
                         R2_cond = c(rep(NA, length(mod.list))),
                         formula = c(rep(NA, length(mod.list))))
  
  # SECOND: CALCULATE AIC, BIC, MALLOWS CP, R2
  for(i in 1:length(mod.list)) {
      model.df$call[i] <- paste0(mod.list[[i]]$call$formula)[3]
      model.df$AIC[i] <- AIC(mod.list[[i]])
      model.df$BIC[i] <- BIC(mod.list[[i]])
      model.df$R2_marg[i] <- as.numeric(r2_nakagawa(mod.list[[i]])$R2_marginal)
      model.df$R2_cond[i] <- as.numeric(r2_nakagawa(mod.list[[i]])$R2_conditional)
      model.df$formula[i] <- as.character(mod.list[[i]][['call']])[2]
    }
  output <- list(model.df, mod.list)
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GENERALIZED LINEAR MIXED EFFECTS MODEL SELECTION (INTERACTIONS)

# NOTE: y.var must be a character vector
glmm.int.selection <- function(y.var, predictors, df, mod.family, zero.inflation, rem.str = '(1|species/plant_id)'){
  # FIRST: List models
  # Generate all possible combinations of two predictors for interaction
  interaction_combinations <- combn(predictors, 2, simplify = FALSE)
  
  # Create an empty list to store model results
  mod.list <- list()
  
  # Loop through the interaction combinations
  for (combo in interaction_combinations) {
    # Create the formula for the model (e.g., outcome ~ predictor1 * predictor2)
    formula <- as.formula(paste(y.var, "~", paste(combo, collapse = "*"), "+", rem.str))
    
    # Fit the mixed-effects model
    model <- glmmTMB(formula, data = df, family = mod.family, ziformula = as.formula(paste0('~', zero.inflation, sep = '')))
    
    # Store the model result in the list
    mod.list[[paste(combo, collapse = "_")]] <- model
  }
  
  # CREATE DATAFRAME
  model.df <- data.frame(model.id = c(1:length(mod.list)),
                         call = c(rep(NA, length(mod.list))),
                         multicollinearity = c(rep(NA, length(mod.list))),
                         AIC = c(rep(NA, length(mod.list))),
                         BIC = c(rep(NA, length(mod.list))),
                         R2_marg = c(rep(NA, length(mod.list))),
                         R2_cond = c(rep(NA, length(mod.list))),
                         formula = c(rep(NA, length(mod.list))))
  
  
  # THIRD: CALCULATE AIC, BIC, MALLOWS CP
  for(i in 1:length(mod.list)) {
    model.df$call[i] <- paste0(mod.list[[i]]$call$formula)[3]
    model.df$AIC[i] <- AIC(mod.list[[i]])
    model.df$BIC[i] <- BIC(mod.list[[i]])
    model.df$R2_marg[i] <- as.numeric(r2_nakagawa(mod.list[[i]])$R2_marginal)
    model.df$R2_cond[i] <- as.numeric(r2_nakagawa(mod.list[[i]])$R2_conditional)
    model.df$formula[i] <- as.character(mod.list[[i]][['call']])[2]
  }
  output <- list(model.df, mod.list)
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GENERALIZED LINEAR MIXED EFFECTS MODEL SELECTION (INTERACTIONS AND FIXED EFFECTS)

# NOTE: y.var must be a character vector
glmm.int.fixed.selection <- function(y.var, interaction, predictors, mod.family, zero.inflation, df, rem.str = '(1|species/plant_id)'){
  # FIRST: List models
  # Create an empty list to store model results
  mod.list <- list()
  call.vec <- c()
  
  # Remove terms from interaction term from predictors list so they are not duplicated in the function call; interaction terms are automatically included as single predictors as well as first-order interactions in lmer, so there is no need to have instances with lmer(y ~ a*b + a + b)
  predictors <- predictors %>% 
    select(-interaction)
  
  # Loop through the interaction combinations
  for(i in 1:ncol(predictors)){
    call <- colnames(predictors) %>%
      combinations(n = ncol(predictors), r = i, repeats.allowed = F) %>% 
      apply(1, paste0, collapse = ' + ') # all possible combinations of models from 2 - # of predictors we're interested in (note: cannot include only 1 predictor, as this breaks multicollinearity for loop below; also, we're probably not interested in a model with only one explanatory variable)
    for(j in (1+length(call.vec)):(length(call)+length(call.vec))){ # adding linear mixed effects model to mod.list
      mod.list[[j]] <- glmmTMB(as.formula(paste(y.var, '~',  paste(interaction, collapse = "*"), '+', call[j-length(call.vec)], '+', rem.str, sep = '')), data = df, family = mod.family, ziformula = as.formula(paste0('~', zero.inflation, sep = '')))
    }
    call.vec <- append(call.vec, call) # to index where to put model into model list
  }
  
  # CREATE DATAFRAME
  model.df <- data.frame(model.id = c(1:length(mod.list)),
                         call = c(rep(NA, length(mod.list))),
                         multicollinearity = c(rep(NA, length(mod.list))),
                         AIC = c(rep(NA, length(mod.list))),
                         BIC = c(rep(NA, length(mod.list))),
                         R2_marg = c(rep(NA, length(mod.list))),
                         R2_cond = c(rep(NA, length(mod.list))),
                         formula = c(rep(NA, length(mod.list))))
  
  
  # THIRD: CALCULATE AIC, BIC, MALLOWS CP
  for(i in 1:length(mod.list)) {
    model.df$call[i] <- paste0(mod.list[[i]]$call$formula)[3]
    model.df$AIC[i] <- AIC(mod.list[[i]])
    model.df$BIC[i] <- BIC(mod.list[[i]])
    model.df$R2_marg[i] <- as.numeric(r2_nakagawa(mod.list[[i]])$R2_marginal)
    model.df$R2_cond[i] <- as.numeric(r2_nakagawa(mod.list[[i]])$R2_conditional)
    model.df$formula[i] <- as.character(mod.list[[i]][['call']])[2]
  }
  output <- list(model.df, mod.list)
  return(output)
}

