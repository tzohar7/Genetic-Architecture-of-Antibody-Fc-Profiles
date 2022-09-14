require('lme4')
require('optimx')

# Data input and processing
input <- read.csv(paste0("Main_Dataframe.csv", sep=""), header = TRUE)

col_len <- ncol(input)
maindata <- input[,10:col_len]

# Rank-based inverse normal transformation (INT)
for (i in 1:228) {
  maindata[,i] <- qnorm((rank(maindata[,i],na.last="keep")-0.5)/sum(!is.na(maindata[,i])))
}

maindata$MZ <- (input[['Zygosity']] == 'MZ')*1
maindata$DZ <- (input[['Zygosity']] == 'DZ')*1

# NOTE: Demographic data must be requested and this assumes age is current dataframe. Age must be removed from code if running with no demographic data.
maindata['Age'] <- input['Age']
maindata$Age <- qnorm((rank(maindata$Age,na.last="keep")-0.5)/sum(!is.na(maindata$Age)))

maindata['T1'] <- sqrt(0.5)*maindata['DZ']
maindata['T2'] <- maindata['MZ'] + sqrt(0.5)*maindata['DZ']

maindata['Pair'] <- as.integer(input$Pair)

#maindata['Twin'] <- 1:dim(maindata)[1]
maindata['Twin'] <- as.integer(input$SampleID)

#maindata$Twin_ID <- as.character(maindata$Twin_ID)
#maindata['Twin'] <- paste(maindata$Pair_Number, maindata$Twin_ID, sep = "_")

#maindata <- subset(maindata, select = -c(MZ,DZ,Twin_ID))
#names(maindata) <- sub("^Pair_Number$", "Pair", names(maindata))

# Storage dataframe and Setup
df <- data.frame(matrix(ncol = 8, nrow = 228))
x <- c("Phenotype", "H2", "Model", "AIC", 'A', 'C', 'E', 'dAIC')
colnames(df) <- x

model_str <- c('ACE','AE','CE')

# Linear Mixed Model Fitting for each feature
for (j in 1:228) {

  mydata <- maindata[,c(j,231, 234:237)]
  phenoname <- names(mydata)[1]
  names(mydata)[1] <- "y"

  ################################################ ACE MODEL #######################################################

  lmod_ACE <- lFormula(y ~ 1 + Age + (0 + T1|Twin) + (0 + T2|Pair) + (1|Pair), data=mydata,
                   REML=FALSE, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
                                                 check.nobs.vs.nRE="ignore"))

  devfun_ACE <- do.call(mkLmerDevfun, lmod_ACE)

  optfunc <- function(theta) {
    tvec <- c(theta[1],theta[1],theta[2])
    devfun_ACE(tvec)
  }

  opt_ACE <- nloptwrap(par=c(0.5,0.5),optfunc,lower=c(-Inf,-Inf),upper=c(Inf,Inf))
  m_ACE <- mkMerMod(environment(devfun_ACE), opt_ACE, lmod_ACE$reTrms, fr = lmod_ACE$fr)
  ACE_res <- summary(m_ACE)
  ACE_AIC <- ACE_res$AICtab[1]

  sig2A_ACE <- ACE_res$varcor[1]$Twin[1]
  sig2C_ACE <- ACE_res$varcor[3]$Pair.1[1]
  sig2E_ACE <- ACE_res$sigma^2
  h2_ACE <- sig2A_ACE/(sig2A_ACE+sig2C_ACE+sig2E_ACE)

  ############################################## AE MODEL ##########################################################

  lmod_AE <- lFormula(y ~ 1 + Age + (0 + T1|Twin) + (0 + T2|Pair), data=mydata,
                   REML=FALSE, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
                                                   check.nobs.vs.nRE="ignore"))

  devfun_AE <- do.call(mkLmerDevfun, lmod_AE)

  optfunc <- function(theta) {
    tvec <- c(theta,theta)
    devfun_AE(tvec)
  }

  opt_AE <- nloptwrap(par=0.5,optfunc,lower=-Inf,upper=Inf)
  m_AE <- mkMerMod(environment(devfun_AE), opt_AE, lmod_AE$reTrms, fr = lmod_AE$fr)
  AE_res <- summary(m_AE)
  AE_AIC <- AE_res$AICtab[1]

  sig2A_AE <- AE_res$varcor[1]$Twin[1]
  sig2C_AE <- 0
  sig2E_AE <- AE_res$sigma^2
  h2_AE <- sig2A_AE/(sig2A_AE+sig2C_AE+sig2E_AE)

  ############################################## CE MODEL ##########################################################

  lmod_CE <- lFormula(y ~  1 + Age + (1|Pair), data=mydata,
                   REML=FALSE, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
                                                 check.nobs.vs.nRE="ignore"))

  devfun_CE <- do.call(mkLmerDevfun, lmod_CE)

  optfunc <- function(theta) {
    tvec <- theta
    devfun_CE(tvec)
  }

  opt_CE <- nloptwrap(par=0.5,optfunc,lower=-Inf,upper=Inf)
  m_CE <- mkMerMod(environment(devfun_CE), opt_CE, lmod_CE$reTrms, fr = lmod_CE$fr)
  CE_res <- summary(m_CE)
  CE_AIC <- CE_res$AICtab[1]

  sig2A_CE <- 0
  sig2C_CE <- CE_res$varcor[1]$Pair[1]
  sig2E_CE <- CE_res$sigma^2
  h2_CE <- 0

  ######################## CHOOSING BEST MODEL AND GET HERITABILITY ########################

  AICs <- c(ACE_AIC, AE_AIC, CE_AIC)
  chosen_model <- which.min(AICs)
  h2res <- c(h2_ACE, h2_AE, h2_CE)

  # Opitional: for ACE only results
  chosen_model <- 1

  finalAIC <- AICs[chosen_model]
  finalmodelstr <- model_str[chosen_model]
  finalh2 <- h2res[chosen_model]
  dAIC <- CE_AIC - finalAIC

  ace_var <- rbind(c(sig2A_ACE, sig2C_ACE, sig2E_ACE), c(sig2A_AE, sig2C_AE, sig2E_AE), c(sig2A_CE, sig2C_CE, sig2E_CE))

  df[j,] <- c(phenoname,finalh2,finalmodelstr,finalAIC, ace_var[chosen_model,], dAIC)
  print(j/228)
}

write.csv(df,"Heritabilty_Results.csv", row.names = FALSE)
