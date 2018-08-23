##########################################################################
# Code for analyzing a simulated dataset to illustrate the weighted and 
# replicated longitudinal analysis methods for SMART trials proposed by 
# Lu, Nahum-Shani, Kasari, Lynch, Oslin, Pelham, Fabiano & Almirall (2016).
# Code by John J. Dziak and Jamie Yap, based on sample code by Xi Lu.
##########################################################################
# Choose options:
WorkingStructureOption <- "independence";  # Try choosing "independence", "AR-1", or "unstructured"
# to get a covariance matrix that treats observations as independent, autoregressive of order
# one, or having a freely estimated structure within each replication of a subject.; 
WeightsEstimationOption <- "known";        # Try choosing "known" or "estimated";
##########################################################################
# Load the data and prepare to start the analysis;
library(geepack);  # Load a package by S. Hojsgaard, U. Halekoh, & J. Yan for doing  
#  generalized estimating equations (Liang & Zeger, 1986);
DataWideFormat <- read.table("TutorialDatasetWideFormat.txt",header=TRUE,na.strings=".");
nTimes <- 3; 
##########################################################################
# Create the known weights based on the design.  
DataWideFormat$DesignWeight <- NA;
DataWideFormat$DesignWeight[which(DataWideFormat$R==0)] <- 4;
DataWideFormat$DesignWeight[which(DataWideFormat$R==1)] <- 2;
# Recall that DataWideFormat$R tells whether a subject is considered a responder or not. 
###########################################################################
# One option is to use the known weights as the weights for the analysis.  The other
# option is to estimate weights which might provide better performance.
if (WeightsEstimationOption=="known") {
  DataWideFormat$FinalWeight <- DataWideFormat$DesignWeight;
}
if (WeightsEstimationOption=="estimated") {
  DataWideFormat$A1DummyCoded <- 1*(DataWideFormat$A1==+1);
  logisticModel1 <- glm(formula=A1DummyCoded ~ CenteredAge + Male + Education + HDDays + Y1, 
                        # It is not yet known whether Y1 is needed here.  
                        # The choice of other covariates is up to the user.;
                        family=binomial,
                        data=DataWideFormat);
  # This model predicts stage 1 treatment assignment from baseline covariates,
  # potentially including baseline outcome Y1.
  DataWideFormat$p1 <- fitted(logisticModel1);
  DataWideFormat$EstimatedWeight1 <- DataWideFormat$A1DummyCoded/DataWideFormat$p1 +
    (1-DataWideFormat$A1DummyCoded)/(1-DataWideFormat$p1);
  # The estimated weight at stage 1 is the inverse estimated probability
  # of receiving the treatment which the participant did in fact receive 
  # at stage 1.  This is analogous to an inverse propensity weight.
  DataWideFormat$A2DummyCoded <- 1*(DataWideFormat$A2==+1);
  DataWideFormat$A2DummyCoded[which(DataWideFormat$R==1)] <- NA;
  # responders were not re-randomized, so A2 is neither -1 nor +1;
  logisticModel2 <- glm(formula=A2DummyCoded ~ CenteredAge + Male + Education + HDDays + Y1 + Y2, 
                        # It is not yet known whether Y1 is needed here.  
                        # The choice of other covariates is up to the user.;
                        family=binomial,
                        na.action=na.exclude,
                        data=DataWideFormat); 
  # Now we predict stage 2 treatment assignment.
  DataWideFormat$p2 <- fitted(logisticModel2);
  DataWideFormat$EstimatedWeight2 <- 1;     
  who.responded <- which(DataWideFormat$R==0);
  DataWideFormat$EstimatedWeight2[who.responded] <- 
    DataWideFormat$A2DummyCoded[who.responded]/DataWideFormat$p2[who.responded] +
    (1-DataWideFormat$A2DummyCoded[who.responded])/(1-DataWideFormat$p2[who.responded]);
  # Responders have a stage 2 weight of 1;  nonresponders have a stage 2 weight which is their
  # estimated probability of getting the treatment they did in fact get.
  DataWideFormat$FinalWeight <- DataWideFormat$EstimatedWeight1 * DataWideFormat$EstimatedWeight2;
  # The overall estimated weight is the product of the estimated weights from the two randomizations.; 
}
###########################################################################
# Reorganize the data from "wide" form to "long" form.
Time1 <- data.frame(cbind(DataWideFormat[,c("SubjectID","A1","R","A2","DesignWeight","FinalWeight")], 
                          Time=1, Y=DataWideFormat$Y1));
Time2 <- data.frame(cbind(DataWideFormat[,c("SubjectID","A1","R","A2","DesignWeight","FinalWeight")], 
                          Time=2, Y=DataWideFormat$Y2));
Time3 <- data.frame(cbind(DataWideFormat[,c("SubjectID","A1","R","A2","DesignWeight","FinalWeight")], 
                          Time=3, Y=DataWideFormat$Y3));
DataLongFormat <- rbind(Time1,Time2,Time3);
DataLongFormat <- DataLongFormat[order(DataLongFormat$SubjectID, DataLongFormat$Time),];
###########################################################################
# Recode time by creating the variables that tell how much time has been spent in  
# Stage One (before Time 2) or in Stage Two (after Time 2).
DataLongFormat$S1 <- NA;
DataLongFormat$S1[which(DataLongFormat$Time==1)] <- 0;
DataLongFormat$S1[which(DataLongFormat$Time==2)] <- 1;
DataLongFormat$S1[which(DataLongFormat$Time==3)] <- 1;
DataLongFormat$S2 <- NA;
DataLongFormat$S2[which(DataLongFormat$Time==1)] <- 0;
DataLongFormat$S2[which(DataLongFormat$Time==2)] <- 0;
DataLongFormat$S2[which(DataLongFormat$Time==3)] <- 1;
stopifnot(all(DataLongFormat$S1+DataLongFormat$S2==DataLongFormat$Time-1)); 
# If S1 and S2 are coded correctly, S1 + S2 should equal Time - 1, that is, 
# the total time since the beginning of treatment.
##############################################################
# Create "replications." 
# Replicate people who got A1=+1 and responded with R=1.  They were not
# re-randomized, so we have to replicate their data to inform both of
# the dynamic treatment regimens which they could have received had they
# been re-randomized.  It is somewhat as if we are creating two clones of
# each of them and counting each in a different treatment, because in
# reality it is indistinguishable which one they received.
##############################################################
DataLongFormat$wave <- DataLongFormat$Time;
RowsToReplicate <- DataLongFormat[which(DataLongFormat$R==1),];
RowsNotToReplicate <- DataLongFormat[which(DataLongFormat$R==0),];
PlusOnePseudodata <- RowsToReplicate;
PlusOnePseudodata$A2 <- 1;
MinusOnePseudodata <- RowsToReplicate;
MinusOnePseudodata$A2 <- -1;
MinusOnePseudodata$wave <- MinusOnePseudodata$Time + nTimes;  
# We keep the same subject ID to show that we don't really have all those
# new participants.  So we have to distinguish the new observations somehow,
# and so we treat them as new waves of data on the same person.  Although
# it seems very ad-hoc, this method has been shown to be valid.
##############################################################
# Create the final analysis dataset including replicates.
##############################################################
DataForAnalysis <- rbind(PlusOnePseudodata, MinusOnePseudodata, RowsNotToReplicate);
DataForAnalysis <- DataForAnalysis[order(DataForAnalysis$SubjectID,DataForAnalysis$wave),] 
# This sorts by the variables SubjectID and wave;
##############################################################
#  Do analysis with GEE and the weighted and replicated data, 
# using (at least initially) the working independence structure.
############################################################## 
GEEIndependent <- geeglm(formula = Y ~   S1 +
                           S2 +
                           S1:A1 + 
                           S2:A1 + 
                           S2:A2 +
                           S2:A1:A2,  
                         id=SubjectID,
                         weights = FinalWeight,    
                         data=DataForAnalysis,
                         corstr = "independence");
##############################################################
#  Do analysis with the non-independence structure, if desired.
############################################################## 
if (WorkingStructureOption=="independence") {
  FinalGEEOutput <- GEEIndependent;
  BlockWorkCorr <- diag(rep(1,nTimes)); # identity matrix;
}
if (WorkingStructureOption=="AR-1") {
  if(sum(is.na(DataLongFormat$Y))>0) {
    stop("This version of the analysis code assumes there is no missing data.")
  }
  # Estimate the AR-1 correlation parameter using the method of moments, by
  # finding the average cross-product of residuals of adjacent observations
  # and then dividing it by the average squared residuals of observations.
  # That is, cor(Y[t],Y[t+1]) = cov(Y[t],Y[t+1]) / var(Y).  This assumes
  # homoskedasticity (equal variance across time) for Y.  However, more 
  # elaborate code could work around this assumption.
  residualsByWave <- list(GEEIndependent$residuals[which(DataForAnalysis$wave==1)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==2)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==3)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==4)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==5)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==6)]);
  weightsByWave <- list(GEEIndependent$weights[which(DataForAnalysis$wave==1)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==2)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==3)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==4)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==5)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==6)]);
  idsByWave <- list(GEEIndependent$id[which(DataForAnalysis$wave==1)],
                    GEEIndependent$id[which(DataForAnalysis$wave==2)],
                    GEEIndependent$id[which(DataForAnalysis$wave==3)],
                    GEEIndependent$id[which(DataForAnalysis$wave==4)],
                    GEEIndependent$id[which(DataForAnalysis$wave==5)],
                    GEEIndependent$id[which(DataForAnalysis$wave==6)]);
  stopifnot(all(idsByWave[[1]]==idsByWave[[2]]) &
              all(idsByWave[[1]]==idsByWave[[3]]) &
              all(idsByWave[[4]]==idsByWave[[5]]) &
              all(idsByWave[[4]]==idsByWave[[6]]));  # double check that residuals are not being mismatched;
  stopifnot(all(weightsByWave[[1]]==weightsByWave[[2]]) &
              all(weightsByWave[[1]]==weightsByWave[[3]]) &
              all(weightsByWave[[4]]==weightsByWave[[5]]) &
              all(weightsByWave[[4]]==weightsByWave[[6]]));  # double check that weights are equal across wave within replicant as we assume;
  averageSquaredResidual <- weighted.mean(x=c(residualsByWave[[1]]^2,
                                              residualsByWave[[2]]^2,
                                              residualsByWave[[3]]^2,
                                              residualsByWave[[4]]^2,
                                              residualsByWave[[5]]^2,
                                              residualsByWave[[6]]^2),
                                          w=c(weightsByWave[[1]],
                                              weightsByWave[[2]],
                                              weightsByWave[[3]],
                                              weightsByWave[[4]],
                                              weightsByWave[[5]],
                                              weightsByWave[[6]]));
  averageCrossProductResidual <- weighted.mean(x=c(residualsByWave[[1]]*residualsByWave[[2]],
                                                   residualsByWave[[2]]*residualsByWave[[3]],
                                                   residualsByWave[[4]]*residualsByWave[[5]],
                                                   residualsByWave[[5]]*residualsByWave[[6]]),
                                               w=c(weightsByWave[[1]],
                                                   weightsByWave[[2]],
                                                   weightsByWave[[4]],
                                                   weightsByWave[[5]]));
  correlationEstimate <- averageCrossProductResidual / averageSquaredResidual;
  BlockWorkCorr <- matrix(0,nTimes,nTimes);
  for (thisRow in 1:nrow(BlockWorkCorr)) {
    for (thisColumn in 1:ncol(BlockWorkCorr)) {
      BlockWorkCorr[thisRow,thisColumn] <- correlationEstimate^abs(thisRow-thisColumn);
    }
  } 
  WorkCorr <- rbind(cbind(BlockWorkCorr,0*BlockWorkCorr),cbind(0*BlockWorkCorr,BlockWorkCorr)); 
  WorkCorrAsZCor <- fixed2Zcor(cor.fixed=WorkCorr, 
                               id=DataForAnalysis$SubjectID,  
                               waves=DataForAnalysis$wave); 
  FinalGEEOutput <- geeglm(formula = Y ~   S1 +
                             S2 +
                             S1:A1 + 
                             S2:A1 + 
                             S2:A2 +
                             S2:A1:A2,  
                           id=SubjectID,
                           weights = FinalWeight,    
                           data=DataForAnalysis,
                           corstr = "fixed",
                           zcor=WorkCorrAsZCor);   
}
if (WorkingStructureOption=="unstructured") {
  if(sum(is.na(DataLongFormat$Y))>0) {
    stop("This version of the analysis code assumes there is no missing data.")
  }
  # We have three correlation parameters to estimate using the method of 
  # moments here.  Specifically, they are cor(Y[1],Y[2]), cor(Y[1],Y[3]),
  # and cor(Y[2],Y[3])
  # As in the AR-1 case, we assume var(Y[1])=var(Y[2])=var(Y[3]).
  # 
  residualsByWave <- list(GEEIndependent$residuals[which(DataForAnalysis$wave==1)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==2)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==3)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==4)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==5)],
                          GEEIndependent$residuals[which(DataForAnalysis$wave==6)]);
  weightsByWave <- list(GEEIndependent$weights[which(DataForAnalysis$wave==1)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==2)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==3)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==4)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==5)],
                        GEEIndependent$weights[which(DataForAnalysis$wave==6)]);
  idsByWave <- list(GEEIndependent$id[which(DataForAnalysis$wave==1)],
                    GEEIndependent$id[which(DataForAnalysis$wave==2)],
                    GEEIndependent$id[which(DataForAnalysis$wave==3)],
                    GEEIndependent$id[which(DataForAnalysis$wave==4)],
                    GEEIndependent$id[which(DataForAnalysis$wave==5)],
                    GEEIndependent$id[which(DataForAnalysis$wave==6)]);
  stopifnot(all(idsByWave[[1]]==idsByWave[[2]]) &
              all(idsByWave[[1]]==idsByWave[[3]]) &
              all(idsByWave[[4]]==idsByWave[[5]]) &
              all(idsByWave[[4]]==idsByWave[[6]]));  # double check that residuals are not being mismatched;
  stopifnot(all(weightsByWave[[1]]==weightsByWave[[2]]) &
              all(weightsByWave[[1]]==weightsByWave[[3]]) &
              all(weightsByWave[[4]]==weightsByWave[[5]]) &
              all(weightsByWave[[4]]==weightsByWave[[6]]));  # double check that weights are equal across wave within replicant as we assume;
  averageSquaredResidual <- weighted.mean(x=c(residualsByWave[[1]]^2,
                                              residualsByWave[[2]]^2,
                                              residualsByWave[[3]]^2,
                                              residualsByWave[[4]]^2,
                                              residualsByWave[[5]]^2,
                                              residualsByWave[[6]]^2),
                                          w=c(weightsByWave[[1]],
                                              weightsByWave[[2]],
                                              weightsByWave[[3]],
                                              weightsByWave[[4]],
                                              weightsByWave[[5]],
                                              weightsByWave[[6]]));
  averageCrossProductResidual12 <- weighted.mean(x=c(residualsByWave[[1]]*residualsByWave[[2]],
                                                     residualsByWave[[4]]*residualsByWave[[5]]),
                                                 w=c(weightsByWave[[1]],
                                                     weightsByWave[[4]]));
  averageCrossProductResidual13 <- weighted.mean(x=c(residualsByWave[[1]]*residualsByWave[[3]],
                                                     residualsByWave[[4]]*residualsByWave[[6]]),
                                                 w=c(weightsByWave[[1]],
                                                     weightsByWave[[4]]));
  averageCrossProductResidual23 <- weighted.mean(x=c(residualsByWave[[2]]*residualsByWave[[3]],
                                                     residualsByWave[[5]]*residualsByWave[[6]]),
                                                 w=c(weightsByWave[[2]],
                                                     weightsByWave[[5]]));
  correlationEstimate12 <- averageCrossProductResidual12 / averageSquaredResidual;
  correlationEstimate13 <- averageCrossProductResidual13 / averageSquaredResidual;
  correlationEstimate23 <- averageCrossProductResidual23 / averageSquaredResidual;
  BlockWorkCorr <- matrix(c(1, correlationEstimate12, correlationEstimate13,
                            correlationEstimate12, 1, correlationEstimate23,
                            correlationEstimate13, correlationEstimate23, 1),
                          nrow=3,ncol=3,
                          byrow = TRUE);
  WorkCorr <- rbind(cbind(BlockWorkCorr,0*BlockWorkCorr),cbind(0*BlockWorkCorr,BlockWorkCorr)); 
  WorkCorrAsZCor <- fixed2Zcor(cor.fixed=WorkCorr, 
                               id=DataForAnalysis$SubjectID,  
                               waves=DataForAnalysis$wave); 
  FinalGEEOutput <- geeglm(formula = Y ~   S1 +
                             S2 +
                             S1:A1 + 
                             S2:A1 + 
                             S2:A2 +
                             S2:A1:A2,  
                           id=SubjectID,
                           weights = FinalWeight,    
                           data=DataForAnalysis,
                           corstr = "fixed",
                           zcor=WorkCorrAsZCor);   
  GEECoefficients <- coef(FinalGEEOutput);
  GEECovarianceMatrix <- FinalGEEOutput$geese$vbeta;
}
GEECoefficients <- coef(FinalGEEOutput);
GEECovarianceMatrix <- FinalGEEOutput$geese$vbeta;
##############################################################
# Convert the regression coefficients estimates into the desired linear contrast estimates. 
############################################################## 
ContrastCoefficients <- matrix(c(
  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,   # Time 1 Mean,   Any Regimen
  1.00,  1.00,  0.00,  1.00,  0.00,  0.00,  0.00,   # Time 2 Mean,   ++ or +-
  1.00,  1.00,  0.00, -1.00,  0.00,  0.00,  0.00,   # Time 2 Mean,  -+ or --
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,   # Time 3 Mean,   ++
  1.00,  1.00,  1.00,  1.00,  1.00, -1.00, -1.00,   # Time 3 Mean,   +-
  1.00,  1.00,  1.00, -1.00, -1.00,  1.00, -1.00,   # Time 3 Mean,  -+
  1.00,  1.00,  1.00, -1.00, -1.00, -1.00,  1.00,   # Time 3 Mean,  --
  0.00,  0.00,  0.00,  2.00,  0.00,  0.00,  0.00,   # Time 2 Mean,   +  versus - 
  0.00,  0.00,  0.00,  0.00,  0.00,  2.00,  2.00,   # Time 3 Mean,   ++ versus +-
  0.00,  0.00,  0.00,  2.00,  2.00,  0.00,  2.00,   # Time 3 Mean,   ++ versus -+
  0.00,  0.00,  0.00,  2.00,  2.00,  2.00,  0.00,   # Time 3 Mean,   ++ versus --
  0.00,  0.00,  0.00,  2.00,  2.00, -2.00,  0.00,   # Time 3 Mean,   +- versus -+
  0.00,  0.00,  0.00,  2.00,  2.00,  0.00, -2.00,   # Time 3 Mean,   +- versus --
  0.00,  0.00,  0.00,  0.00,  0.00,  2.00, -2.00,   # Time 3 Mean,  -+ versus --
  0.00,  1.00,  0.00,  1.00,  0.00,  0.00,  0.00,   # Stage 1 Slope,   ++ or +-
  0.00,  1.00,  0.00, -1.00,  0.00,  0.00,  0.00,   # Stage 1 Slope,  -+ or --
  0.00,  0.00,  1.00,  0.00,  1.00,  1.00,  1.00,   # Stage 2 Slope,   ++
  0.00,  0.00,  1.00,  0.00,  1.00, -1.00, -1.00,   # Stage 2 Slope,   +-
  0.00,  0.00,  1.00,  0.00, -1.00,  1.00, -1.00,   # Stage 2 Slope,  -+
  0.00,  0.00,  1.00,  0.00, -1.00, -1.00,  1.00,   # Stage 2 Slope,  --
  0.00,  0.00,  0.00,  0.00,  0.00,  2.00,  2.00,   # Stage 2 Slope,   ++ vs. +-
  0.00,  0.00,  0.00,  0.00,  2.00,  0.00,  2.00,   # Stage 2 Slope,   ++ vs. -+
  0.00,  0.00,  0.00,  0.00,  2.00,  2.00,  0.00,   # Stage 2 Slope,   ++ vs. --
  0.00,  0.00,  0.00,  0.00,  2.00, -2.00,  0.00,   # Stage 2 Slope,   +- vs. -+
  0.00,  0.00,  0.00,  0.00,  2.00,  0.00, -2.00,   # Stage 2 Slope,   +- vs. --
  0.00,  0.00,  0.00,  0.00,  0.00,  2.00, -2.00,   # Stage 2 Slope,  -+ vs. --
  2.00,  1.50,  0.50,  1.50,  0.50,  0.50,  0.50,   # Area under Curve,   ++
  2.00,  1.50,  0.50,  1.50,  0.50, -0.50, -0.50,   # Area under Curve,   +-
  2.00,  1.50,  0.50, -1.50, -0.50,  0.50, -0.50,   # Area under Curve,  -+
  2.00,  1.50,  0.50, -1.50, -0.50, -0.50,  0.50,   # Area under Curve,  --
  0.00,  0.00,  0.00,  0.00,  0.00, -1.00,  1.00,   # Area,   ++ vs. +-
  0.00,  0.00,  0.00,  3.00,  1.00,  0.00,  1.00,   # Area,   ++ vs. -+
  0.00,  0.00,  0.00,  3.00,  1.00,  1.00,  0.00,   # Area,   ++ vs. --
  0.00,  0.00,  0.00,  3.00,  1.00, -1.00,  0.00,   # Area,   +- vs. -+
  0.00,  0.00,  0.00,  3.00,  1.00,  0.00, -1.00,   # Area,   +- vs. --
  0.00,  0.00,  0.00,  3.00,  1.00,  0.00, -1.00,   # Area,   +- vs. --
  1.00,  0.75,  0.25,  0.75,  0.25,  0.25,  0.25,   # Average value,   ++
  1.00,  0.75,  0.25,  0.75,  0.25, -0.25, -0.25,   # Average value,   +-
  1.00,  0.75,  0.25, -0.75, -0.25,  0.25, -0.25,   # Average value,  -+
  1.00,  0.75,  0.25, -0.75, -0.25, -0.25,  0.25,   # Average value,  --
  0.00,  0.00,  0.00,  0.00,  0.00, -0.50,  0.50,   # Average,   ++ vs. +-
  0.00,  0.00,  0.00,  1.50,  0.50,  0.00,  0.50,   # Average,   ++ vs. -+
  0.00,  0.00,  0.00,  1.50,  0.50,  0.50,  0.00,   # Average,   ++ vs. --
  0.00,  0.00,  0.00,  1.50,  0.50, -0.50,  0.00,   # Average,   +- vs. -+
  0.00,  0.00,  0.00,  1.50,  0.50,  0.00, -0.50,   # Average,   +- vs. --
  0.00,  0.00,  0.00,  1.50,  0.50,  0.00, -0.50,   # Average,   +- vs. --
  0.00,  0.00,  0.00,  0.00,  2.00,  0.00,  2.00,   # Delayed Effect,   ++ vs. -+
  0.00,  0.00,  0.00,  0.00,  2.00,  2.00,  0.00,   # Delayed Effect,   ++ vs. --
  0.00,  0.00,  0.00,  0.00,  2.00, -2.00,  0.00,   # Delayed Effect,   +- vs. -+
  0.00,  0.00,  0.00,  0.00,  2.00,  0.00, -2.00,   # Delayed Effect,   +- vs. --
  0.00,  0.00,  0.00,  0.00,  2.00,  0.00,  0.00),   # Ave. Delayed Eff.,  + vs -  
  byrow=TRUE,ncol=7);                               
rownames(ContrastCoefficients) <- c(
  "Time 1 Mean,  Any Regimen",
  "Time 2 Mean,  ++ or +-",
  "Time 2 Mean,  -+ or --",
  "Time 3 Mean,  ++",
  "Time 3 Mean,  +-",
  "Time 3 Mean,  -+",
  "Time 3 Mean,  --",
  "Time 2 Mean,  +  versus - ",
  "Time 3 Mean,  ++ versus +-",
  "Time 3 Mean,  ++ versus -+",
  "Time 3 Mean,  ++ versus --",
  "Time 3 Mean,  +- versus -+",
  "Time 3 Mean,  +- versus --",
  "Time 3 Mean,  -+ versus --",
  "Stage 1 Slope,  ++ or +-",
  "Stage 1 Slope,  -+ or --",
  "Stage 2 Slope,  ++",
  "Stage 2 Slope,  +-",
  "Stage 2 Slope,  -+",
  "Stage 2 Slope,  --",
  "Stage 2 Slope,  ++ vs. +-",
  "Stage 2 Slope,  ++ vs. -+",
  "Stage 2 Slope,  ++ vs. --",
  "Stage 2 Slope,  +- vs. -+",
  "Stage 2 Slope,  +- vs. --",
  "Stage 2 Slope,  -+ vs. --",
  "Area under Curve,  ++",
  "Area under Curve,  +-",
  "Area under Curve,  -+",
  "Area under Curve,  --",
  "Area,  ++ vs. +-",
  "Area,  ++ vs. -+",
  "Area,  ++ vs. --",
  "Area,  +- vs. -+",
  "Area,  +- vs. --",
  "Area,  +- vs. --",
  "Average value,  ++",
  "Average value,  +-",
  "Average value,  -+",
  "Average value,  --",
  "Average,  ++ vs. +-",
  "Average,  ++ vs. -+",
  "Average,  ++ vs. --",
  "Average,  +- vs. -+",
  "Average,  +- vs. --",
  "Average,  +- vs. --",
  "Delayed Effect,  ++ vs. -+",
  "Delayed Effect,  ++ vs. --",
  "Delayed Effect,  +- vs. -+",
  "Delayed Effect,  +- vs. --",
  "Ave. Delayed Eff.,  + vs -");
colnames(ContrastCoefficients) <- names(coef(FinalGEEOutput));
ContrastEstimates <- as.vector(ContrastCoefficients%*%GEECoefficients);  # Note the matrix multiplication;
if (WeightsEstimationOption=="estimated") {
  # Calculate the score functions for the logistic regressions, using the original data:
  nsub <- nrow(DataWideFormat);
  Stage1LogisticResiduals <- DataWideFormat$A1DummyCoded - DataWideFormat$p1;
  scoreLogistic1 <- model.matrix(logisticModel1) * as.vector(Stage1LogisticResiduals);
  Stage2LogisticResiduals <- DataWideFormat$A2DummyCoded - DataWideFormat$p2;
  scoreLogistic2 <- matrix(0,nsub,ncol(model.matrix(logisticModel2)));
  scoreLogistic2[which(!is.na(Stage2LogisticResiduals)),] <- model.matrix(logisticModel2) * as.vector(Stage2LogisticResiduals[which(!is.na(Stage2LogisticResiduals))]); 
  # Finish calculating the standard errors, using the replicated data. 
  predictors <- model.matrix(FinalGEEOutput);
  resid <- FinalGEEOutput$residuals;
  nobs <- nrow(predictors);
  fitted <- FinalGEEOutput$fitted.values;
  fitted[which(fitted<.0000000000001)] <- .0000000000001; 
  fitted[which(fitted>.9999999999999)] <- .9999999999999; 
  stopifnot(length(resid)==nobs); # try to make sure nothing has gone wrong;
  meat <- matrix(0,ncol(predictors),ncol(predictors));
  bread <- matrix(0,ncol(predictors),ncol(predictors));
  # "I" and "J" are the notation used  in Xi Lu's code for the two matrices
  # we call "meat" and "bread," respectively. 
  # They don't represent the usual I and J matrices in matrix algebra which 
  #	are the identity (diagonal) and all-ones matrix respectively.  
  # Instead, they are the two matrices in Theorem 1.2 of the supplemental
  # material of Lu (2016).  The matrix "meat" will be an empirical covariance
  # estimator of the score function of the GEE equation, minus a correction 
  # for uncertainty in estimating the weights.  The matrix "bread" will be
  # the naive information matrix, that is, its inverse would be the GEE 
  # model-based covariance estimate of the GEE coefficients.  The sandwich
  # matrix t(bread) * meat * bread will be a better covariance estimate of 
  # the GEE coefficients.
  U <- matrix(0,nsub,ncol(predictors));
  # The rows of U will be set to equal the score function for each subject (derivative of that subject's
  # log-pseudo-likelihood contribution with respect to the GEE coefficients) in the FinalGEEOutput model,
  # with each subject in the original data set having one row. */   
  indicesNextSubject <- 1:nTimes;  
  for (i in 1:nsub) {
    indicesThisSubject <- indicesNextSubject;
    # The code as written assumes that no data are missing, and that the data are sorted by
    # person and then time, with nTimes rows for nonresponders and 2*nTimes rows for responders.
    weightThisSubject <- DataForAnalysis$DesignWeight[indicesThisSubject[1]];
    # The weights are the same for all observations within the 
    # subject, so we just read the first one.
    residualsThisSubject <- resid[indicesThisSubject];
    predictorsThisSubject <- predictors[indicesThisSubject,] ;  
    fittedValuesThisSubject <- fitted[indicesThisSubject];
    GlmWeightThisSubject <- fittedValuesThisSubject*(1-fittedValuesThisSubject); 
    errorThisSubject <- weightThisSubject * t(predictorsThisSubject) %*% 
      diag(as.vector(sqrt(GlmWeightThisSubject)))%*% 
      solve(BlockWorkCorr) %*% 
      diag(as.vector(1/sqrt(GlmWeightThisSubject)))%*%
      residualsThisSubject;
    # errorThisSubject is a k by 1 vector representing this subject's score function.
    # where k=ncol(predictors). ; 
    if (indicesThisSubject[nTimes] == nobs) {
      # This is the case where the current observations are for the last individual,
      # so no more observations are forthcoming. 
      U[i,] <- t(errorThisSubject);
      thisBread <- weightThisSubject * t(predictorsThisSubject)%*%
        diag(as.vector(sqrt(GlmWeightThisSubject)))%*% 
        solve(BlockWorkCorr) %*% 
        diag(as.vector(1/sqrt(GlmWeightThisSubject)))%*%
        predictorsThisSubject; 
      bread <- bread + thisBread;
    } else {
      if (DataForAnalysis$wave[[indicesThisSubject[nTimes]+1]] == 1) {
        # This is the case where the next observations are from an actual new individual. */
        U[i,] <- t(errorThisSubject);
        thisBread <- weightThisSubject * t(predictorsThisSubject)%*%
          diag(as.vector(sqrt(GlmWeightThisSubject)))%*% 
          solve(BlockWorkCorr) %*% 
          diag(as.vector(1/sqrt(GlmWeightThisSubject)))%*%
          predictorsThisSubject; 
        bread <- bread + thisBread;
        indicesNextSubject <- (indicesThisSubject[nTimes]+1):(indicesThisSubject[nTimes]+nTimes); 
      }
      if (DataForAnalysis$wave[[indicesThisSubject[nTimes]+1]] == nTimes+1) {
        # This is the case where the next observations are a replicate of this individual.;
        indicesReplicate <- (indicesThisSubject[nTimes]+1):(indicesThisSubject[nTimes]+nTimes);
        residReplicate <- resid[indicesReplicate];
        predictorsReplicate <- predictors[indicesReplicate,];
        errorReplicate <- weightThisSubject * t(predictorsReplicate)%*%
          diag(as.vector(sqrt(GlmWeightThisSubject)))%*% 
          solve(BlockWorkCorr) %*% 
          diag(as.vector(1/sqrt(GlmWeightThisSubject)))%*%
          residReplicate;
        U[i,] <- t(errorThisSubject + errorReplicate);
        thisBread <- weightThisSubject * t(predictorsThisSubject)%*%
          diag(as.vector(sqrt(GlmWeightThisSubject)))%*% 
          solve(BlockWorkCorr) %*% 
          diag(as.vector(1/sqrt(GlmWeightThisSubject)))%*%
          predictorsThisSubject + 
          weightThisSubject * t(predictorsReplicate)%*%
          diag(as.vector(sqrt(GlmWeightThisSubject)))%*% 
          solve(BlockWorkCorr) %*% 
          diag(as.vector(1/sqrt(GlmWeightThisSubject)))%*%
          predictorsReplicate;
        bread <- bread + thisBread;
        indicesNextSubject <- (indicesReplicate[nTimes]+1):(indicesReplicate[nTimes]+nTimes);
      }
      if ((DataForAnalysis$wave[indicesThisSubject[nTimes]+1] != 1) & 
          (DataForAnalysis$wave[indicesThisSubject[nTimes]+1] != nTimes+1) ) {
        print("Unexpected error in code or dataset;");
        print(indicesThisSubject); 
      }
    }
  }
  bread <- bread/nsub;  # re-express as an average, not a sum, of contributions from participants;
  AdjustedSampleSizeForEmpiricalCovariance <- nsub-ncol(predictors); 
                        # some analysts would just use 1/nsub here, 
  # and they are asymptotically equivalent, but this is slightly more conservative in adjusting
  # somewhat for overfitting;
  # If weights were known, then meat would be the empirical covariance matrix 
  # (1/AdjustedSampleSizeForEmpiricalCovariance)*t(U)%*%U.  However, using estimated
  # weights can reduce the covariates, especially if the covariates used to estimate the weights 
  # are good predictors.  To reflect this reduction in error variance, Lu et al (2016, Theorem 1.2) 
  # shows that the covariance estimator should be corrected.  U will be replaced by Uproj, the
  # projection of U onto the space defined by the "hat" or projection matrix of the regression  
  # of the weights on the predictors of the weights.
  scores <- cbind(scoreLogistic1, scoreLogistic2);
  hat <- scores%*%solve(t(scores)%*%scores)%*%t(scores)
  Uproj <- U - hat%*%U;
  meat <- (1/AdjustedSampleSizeForEmpiricalCovariance)*t(Uproj)%*%Uproj;
  invBread <- solve(bread);   
  GEECovarianceMatrix <- (1/nsub)* invBread %*% meat %*% invBread; 
  # This matrix is analogous to the sandwich covariance matrix in standard GEE.
  # See Theorem I.2 in the supplemental material for Lu et al (2016).
}
ContrastCovarianceMatrix <- ContrastCoefficients%*%GEECovarianceMatrix%*%t(ContrastCoefficients);
# This is Cramer's delta method, a simple application of Taylor linearization.
ContrastStdErrors <- sqrt(diag(ContrastCovarianceMatrix)); 
ContrastZRatios <- ContrastEstimates / ContrastStdErrors;
ContrastPValues <- 2*(1-pnorm(abs(ContrastZRatios)));
ContrastUpper95 <- ContrastEstimates + 1.96*ContrastStdErrors;
ContrastLower95 <- ContrastEstimates - 1.96*ContrastStdErrors;
print( cbind(ContrastEstimates, ContrastStdErrors,ContrastZRatios,ContrastPValues,ContrastLower95,ContrastUpper95));