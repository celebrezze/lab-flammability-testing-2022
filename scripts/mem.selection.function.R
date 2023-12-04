
# mallows' cp function is required:
mallows.cp <- function(model, k, n) {
  return(AIC(model) + 2*k / (n - k - 1))
}

# NOTE: y.var and rem.str must be a character vector, predictors must be a dataframe
mem.selection <- function(y.var, predictors, df, rem.str = '(1|plant_id)'){
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
      mod.list[[j]] <- lmer(as.formula(paste(y.var, '~', call[j-length(call.vec)], '+', rem.str, sep = '')), data = df)
    }
    call.vec <- append(call.vec, call) # to index where to put model into model list
  }
  
  # CREATE DATAFRAME
  model.df <- data.frame(model.id = c(1:length(mod.list)),
                         call = c(rep(NA, length(mod.list))),
                         multicollinearity = c(rep(NA, length(mod.list))),
                         AIC = c(rep(NA, length(mod.list))),
                         BIC = c(rep(NA, length(mod.list))),
                         CP = c(rep(NA, length(mod.list))))
  
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
      model.df$CP[i] <- mallows.cp(mod.list[[i]], k = length(mod.list[[i]]@beta - 1), n = nrow(df))
    }
  }
  output <- list(model.df, mod.list)
  return(output)
}

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