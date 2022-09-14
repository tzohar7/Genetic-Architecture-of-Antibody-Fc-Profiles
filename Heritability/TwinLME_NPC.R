require('lme4')
require('optimx')

# Data input and processing
input <- read.csv(paste0("Main_Dataframe.csv", sep=""), header = TRUE)

nperm <- 1000

col_len <- ncol(input)
maindata <- input[,10:col_len]

# Rank-based inverse normal transformation (INT)
for (i in 1:228) {
  maindata[,i] <- qnorm((rank(maindata[,i],na.last="keep")-0.5)/sum(!is.na(maindata[,i])))
}

names(maindata) <- gsub("Strep", "Pneumo", names(maindata))
maindata <- maindata[,!sapply("EBOV", grepl, names(maindata))]

maindata$MZ <- (input[['Zygosity']] == 'MZ')*1
maindata$DZ <- (input[['Zygosity']] == 'DZ')*1

# NOTE: Demographic data must be requested and this assumes age is current dataframe. Age must be removed from code if running with no demographic data.
maindata['Age'] <- input['Age']
maindata$Age <- scale(maindata$Age)

maindata['T1'] <- sqrt(0.5)*maindata['DZ']
maindata['T2'] <- maindata['MZ'] + sqrt(0.5)*maindata['DZ']

maindata['Pair'] <- as.integer(input$Pair)
maindata['Twin'] <- as.integer(input$SampleID)

feature_list <- as.list(unique(sapply(strsplit(names(maindata[,1:227]),"_"), `[`, 1)))
antigen_list <- as.list(unique(sapply(strsplit(names(maindata[,1:227]),"_"), `[`, 2)))

antigen_list[[14]] <- c("Pre.RSV", "RSV")
antigen_list[[15]] <- c("Post.RSV", "RSV")
antigen_list[[19]] <- c("H1N1.CA", "Flu")
antigen_list[[20]] <- c("H3N2.TX", "Flu")
antigen_list <- antigen_list[-c(21,22)]

subfeatures <- append(feature_list,antigen_list)

# Adding higher classes
subfeatures <- append(subfeatures,c("IgG1", "FcR2A", "ADCD", "IgG1", "IgG1",
                                    "IgG1", "IgG1"))

# Titers, FcRs, Titers
subfeatures[[34]] <- c("IgG1", "IgG2", "IgG3", "IgG4")
subfeatures[[35]] <- c("FcR2A", "FcR2b", "FcR3AV", "FcR3b")
subfeatures[[36]] <- c("ADCD", "ADNP", "ADCP")

#subfeatures[[32]] <- c("IgG1", "IgG2", "IgG3", "IgG4")
#subfeatures[[33]] <- c("FcR2A", "FcR2b", "FcR3AV", "FcR3b")
#subfeatures[[34]] <- c("ADCD", "ADNP", "ADCP")

# Bacterial vs Viral
subfeatures[[37]] <- c("Ptox", "Ttox", "PPD", "Hib", "Pneumo", "Dip")
subfeatures[[38]] <- c("HBV", "Measles", "Mumps", "Rubella", "Noro", "Polio", "Post.RSV",
                       "Pre.RSV", "RSV", "VZV", "EBV", "H1N1.CA", "H3N2.TX", "Flu")
                       
# Viral Vaccines Vs Not Vaccinated
subfeatures[[39]] <- c("HBV", "Measles", "Mumps", "Rubella", "Polio", "H1N1.CA", "H3N2.TX", "Flu")
subfeatures[[40]] <- c("Noro", "Post.RSV", "Pre.RSV", "RSV", "EBV", "VZV")

output_file <- paste0("TwinLME_NPC_Results.csv", sep="")

res <- try({read.csv(output_file, header = TRUE)})

if (class(res) == "try-error") {
  df_perm <- data.frame(matrix(ncol = 2 + nperm, nrow = length(subfeatures)))
  x <- c("Phenotype", "Global", paste(1:nperm))
  colnames(df_perm) <- x
  df_perm$Phenotype <- as.character(df_perm$Phenotype)
  starting_val <- 1
} else {
  df_perm <- res
  df_perm$Phenotype <- as.character(df_perm$Phenotype)
  x <- c("Phenotype", "Global", paste(1:nperm))
  colnames(df_perm) <- x
  starting_val <- min(which(apply(!is.na(res), 1, all) == FALSE))
}

# Make Phenotype names
for (i in 1:length(subfeatures)) {
  subfeaturename <- subfeatures[[i]][1]
  df_perm[i,1] <- subfeaturename
}

print("Beginning at:")
print(starting_val)

####################################################################################################
# NPC LME Algorithm
###################################################################################################

for (i in starting_val:length(subfeatures)) {

  idx <- which(rowSums(sapply(subfeatures[[i]], grepl, names(maindata))) == 1)
  mydata <- maindata[,c(idx ,230, 233:236)]
  #mydata <- maindata[,c(idx ,210, 213:216)]

  ####################################################################################################
  # NORMAL TESTING
  ####################################################################################################

  # for features in subfeatures
  H2norm <- list()
  for (k in 1:(length(mydata)-5)) {

    mydatan <- mydata[,c(k ,(length(mydata)-4):length(mydata))]

    names(mydatan)[1] <- "y"

    ########### ACE MODEL ###########
    lmod_ACE <- lFormula(y ~ 1 + Age + (0 + T1|Twin) + (0 + T2|Pair) + (1|Pair), data=mydatan,
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

    H2norm <- rbind(H2norm, h2_ACE)
  }

  global_stat <- mean(unlist(H2norm))

  ####################################################################################################
  # PERMUATION TESTING
  ####################################################################################################

  # Permuting
  idx <- which(rowSums(sapply(subfeatures[[i]], grepl, names(maindata))) == 1)
  mydata <- maindata[,c(idx, 230, 229, 228, 235, 236)
  #mydata <- maindata[,c(idx, 210, 209, 208, 215, 216)]

  # Initate permuations for subfeature
  perm_global_stats <- list()
  j <- 0
  while (j < nperm) {

    res <- try({
                 # Permute

                 pdata <- mydata

                 pdata[,c('MZ', 'DZ')] <- pdata[rep(sample((seq(1,dim(maindata)[1],2))), each=2), c('MZ', 'DZ')]
                 pdata['T1'] <- sqrt(0.5)*pdata['DZ']
                 pdata['T2'] <- pdata['MZ'] + sqrt(0.5)*pdata['DZ']

                 pdata[,1:(length(pdata)-5)] <- pdata[, 1:(length(pdata)-5)]

                 # Permuation for features in subfeature
                 H2sperm <- list()
                      for (k in 1:(length(pdata)-7)) {

                        pidata <- pdata[,c(k ,(length(pdata)-6):length(pdata))]
                        names(pidata)[1] <- "y"

                        ########### ACE MODEL ###########
                        lmod_ACE <- lFormula(y ~ 1 + Age + (0 + T1|Twin) + (0 + T2|Pair) + (1|Pair), data=pidata,
                                             REML=FALSE, control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                                              check.nobs.vs.rankZ = "ignore",
                                                                             check.nobs.vs.nRE="ignore"))
                        devfun_ACE <- do.call(mkLmerDevfun, lmod_ACE)

                        optfunc <- function(theta) {
                          tvec <- c(theta[1],theta[1],theta[2])
                          devfun_ACE(tvec)
                        }

                        opt_ACEp <- nloptwrap(par=c(0.5,0.5),optfunc,lower=c(-Inf,-Inf),upper=c(Inf,Inf))
                        m_ACEp <- mkMerMod(environment(devfun_ACE), opt_ACEp, lmod_ACE$reTrms, fr = lmod_ACE$fr)
                        ACE_resp <- summary(m_ACEp)
                        ACE_AICp <- ACE_resp$AICtab[1]

                        sig2A_ACE <- ACE_resp$varcor[1]$Twin[1]
                        sig2C_ACE <- ACE_resp$varcor[3]$Pair.1[1]
                        sig2E_ACE <- ACE_resp$sigma^2

                        h2_ACEp <- sig2A_ACE/(sig2A_ACE+sig2C_ACE+sig2E_ACE)

                        H2sperm <- rbind(H2sperm, h2_ACEp)
                      }
               }, silent = FALSE)

    if (class(res) != "try-error") {

      perm_global_stats <- rbind(perm_global_stats, mean(unlist(H2sperm), na.rm=T))

      j <- j + 1

      print((j/nperm)*100)

    }

  }

  global_pvalue <- sum(perm_global_stats >= global_stat)/nperm

  df_perm[i,2] <- global_stat
  df_perm[i,3:dim(df_perm)[2]] <- perm_global_stats

  print(paste('Total: ', i/length(subfeatures)))

  write.csv(df_perm,"TwinLME_NPC_Results.csv", row.names = FALSE)

}
