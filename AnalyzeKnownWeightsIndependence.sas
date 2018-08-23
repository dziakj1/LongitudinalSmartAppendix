PROC IMPORT OUT= TutorialDatasetWideFormat
            DATAFILE= "C:\Users\jjd264\Documents\LongitudinalSmart\TutorialAppendix-UnknownWeights-Revised\TutorialDatasetWideFormat.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

DATA TutorialDatasetWideFormat;
	SET TutorialDatasetWideFormat;
    /* Calculate weights: */
    WeightRandomization1 = 2;
    IF (R = 0) THEN WeightRandomization2 = 2; ELSE WeightRandomization2 = 1;
    Weight = WeightRandomization1 * WeightRandomization2;
RUN;

DATA TutorialDatasetLongFormat;
	SET TutorialDatasetWideFormat;
	Time = 1;
	Y = Y1;
	OUTPUT;
	Time = 2;
	Y = Y2;
	OUTPUT;
	Time = 3;
	Y = Y3;
	OUTPUT;
	KEEP SubjectID Weight Time A1 R A2 Y;
RUN;


DATA TutorialDatasetLongFormat;
    SET TutorialDatasetLongFormat;
    /* Calculate recoded time: */
    IF Time > 2 THEN S1 = 1; ELSE S1 = Time - 1; 
            /*S1 is time since the beginning of the study, truncated above at 1;  
           that is, the total amount of time spent so far in Stage 1. */
    IF Time > 2 THEN S2 = Time - 2; ELSE S2 = 0;
            /*S1 is time since the second randomization, truncated below at 0;  
           that is, the total amount of time spent so far in Stage 2. */ 
    /* Note:  S1 + S2 = Time - 1, that is, time since the beginning of the study. */ 
RUN;

/* Do replication: */
DATA PlusOneReplicates;
    SET TutorialDatasetLongFormat;
    WHERE R = 1;
    A2 = +1;
    Wave = Time;  /* We make a new copy of Time here, named Wave, because we are 
                     going to count replicates as different waves. */
RUN;
DATA MinusOneReplicates;
    SET TutorialDatasetLongFormat;
    WHERE R = 1;
    A2 = -1;
    Wave = Time + 3;  /* See above. */
RUN;
DATA RowsNotToReplicate;
    SET TutorialDatasetLongFormat;
    WHERE R = 0;
    Wave = Time;
RUN;
DATA DataForAnalysis;
    SET PlusOneReplicates MinusOneReplicates RowsNotToReplicate;
RUN;
PROC SORT DATA=DataForAnalysis;
    BY SubjectID Wave;
RUN;
 
/* Perform analysis using Generalized Estimating Equations */
PROC GENMOD DATA=DataForAnalysis;
    CLASS SubjectID;
    MODEL Y = S1 S2 S1*A1 S2*A1 S2*A2 S2*A1*A2;
    REPEATED SUBJECT=SubjectID;
    WEIGHT Weight;
    /* Time-Specific Means: */
    /* In the model statement above, we assumed that 
            E(Y) = B0 + B1*S1 + B2*S2 + B3*S1*A1 + B4*S2*A1 + B5*S2*A2 + B6*S2*A1*A2.
       At time 1, S1=0 and S2=0, so E(Y) = B0.
       At time 2, S1=1 and S2=0, so E(Y) = B0 + B1 + B3*A1.
       At time 3, S1=1 and S2=1, so E(Y) = B0 + B1 + B2 + B3*A1 + B4*A1 + B5*A2 + B6*A1*A2. 
       Substituting the different possible values of A1 and A2 into these equations, we find that the 
       coefficients representing the time-specific means are as given below. */ 
    ESTIMATE "Time 1 Mean, Any Regimen"   INTERCEPT 1  S1 0    S2 0    S1*A1  0    S2*A1  0    S2*A2  0    S2*A1*A2  0;
    ESTIMATE "Time 2 Mean, ++ or +-"      INTERCEPT 1  S1 1    S2 0    S1*A1 +1    S2*A1  0    S2*A2  0    S2*A1*A2  0;
    ESTIMATE "Time 2 Mean, -+ or --"      INTERCEPT 1  S1 1    S2 0    S1*A1 -1    S2*A1  0    S2*A2  0    S2*A1*A2  0;
    ESTIMATE "Time 3 Mean, ++"            INTERCEPT 1  S1 1    S2 1    S1*A1 +1    S2*A1 +1    S2*A2 +1    S2*A1*A2 +1; 
    ESTIMATE "Time 3 Mean, +-"            INTERCEPT 1  S1 1    S2 1    S1*A1 +1    S2*A1 +1    S2*A2 -1    S2*A1*A2 -1; 
    ESTIMATE "Time 3 Mean, -+"            INTERCEPT 1  S1 1    S2 1    S1*A1 -1    S2*A1 -1    S2*A2 +1    S2*A1*A2 -1; 
    ESTIMATE "Time 3 Mean, --"            INTERCEPT 1  S1 1    S2 1    S1*A1 -1    S2*A1 -1    S2*A2 -1    S2*A1*A2 +1;  
    /* Contrasts in Time-Specific Means: */
    /* These are found by subtracting corresponding coefficients representing each pair of time-specific means. */
    ESTIMATE "Time 2 Mean, +  versus - "  INTERCEPT 0  S1 0    S2 0    S1*A1 +2    S2*A1  0    S2*A2  0    S2*A1*A2  0;  
    ESTIMATE "Time 3 Mean, ++ versus +-"  INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1  0    S2*A2 +2    S2*A1*A2 +2;  
    ESTIMATE "Time 3 Mean, ++ versus -+"  INTERCEPT 0  S1 0    S2 0    S1*A1 +2    S2*A1 +2    S2*A2  0    S2*A1*A2 +2; 
    ESTIMATE "Time 3 Mean, ++ versus --"  INTERCEPT 0  S1 0    S2 0    S1*A1 +2    S2*A1 +2    S2*A2 +2    S2*A1*A2  0; 
    ESTIMATE "Time 3 Mean, +- versus -+"  INTERCEPT 0  S1 0    S2 0    S1*A1 +2    S2*A1 +2    S2*A2 -2    S2*A1*A2  0; 
    ESTIMATE "Time 3 Mean, +- versus --"  INTERCEPT 0  S1 0    S2 0    S1*A1 +2    S2*A1 +2    S2*A2  0    S2*A1*A2 -2; 
    ESTIMATE "Time 3 Mean, -+ versus --"  INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1  0    S2*A2 +2    S2*A1*A2 -2; 
    /* Stage-Specific Slopes: */
    /* These are also found by subtracting corresponding entries for the appropriate pair of time-specific means,
       noting that because the stages are of length 1, the stage T slope is the time T+1 mean minus the time T 
       mean. */
    ESTIMATE "Stage 1 Slope, ++ or +-"    INTERCEPT 0  S1 1    S2 0    S1*A1 +1    S2*A1  0    S2*A2  0    S2*A1*A2  0;
    ESTIMATE "Stage 1 Slope, -+ or --"    INTERCEPT 0  S1 1    S2 0    S1*A1 -1    S2*A1  0    S2*A2  0    S2*A1*A2  0;
    ESTIMATE "Stage 2 Slope, ++"          INTERCEPT 0  S1 0    S2 1    S1*A1  0    S2*A1 +1    S2*A2 +1    S2*A1*A2 +1; 
    ESTIMATE "Stage 2 Slope, +-"          INTERCEPT 0  S1 0    S2 1    S1*A1  0    S2*A1 +1    S2*A2 -1    S2*A1*A2 -1; 
    ESTIMATE "Stage 2 Slope, -+"          INTERCEPT 0  S1 0    S2 1    S1*A1  0    S2*A1 -1    S2*A2 +1    S2*A1*A2 -1; 
    ESTIMATE "Stage 2 Slope, --"          INTERCEPT 0  S1 0    S2 1    S1*A1  0    S2*A1 -1    S2*A2 -1    S2*A1*A2 +1; 
    /* Contrasts in stage-specific slopes: */
    /* Contrasts in stage 1 slopes are omitted here because they the same as contrasts in time 2 means, provided 
       that time 1 means are assumed to be the same.  Contrasts in stage 2 slopes must be interpreted with great
       caution, because the means at the end of stage 1 are *not* necessarily the same, and so a high stage 2 slope could 
       represent either a high time-3 mean or a low time-2 mean, or both. */
    ESTIMATE "Stage 2 Slope, ++ vs. +-"   INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1  0    S2*A2 +2    S2*A1*A2 +2;
    ESTIMATE "Stage 2 Slope, ++ vs. -+"   INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1 +2    S2*A2  0    S2*A1*A2 +2;
    ESTIMATE "Stage 2 Slope, ++ vs. --"   INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1 +2    S2*A2 +2    S2*A1*A2  0;
    ESTIMATE "Stage 2 Slope, +- vs. -+"   INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1 +2    S2*A2 -2    S2*A1*A2  0;
    ESTIMATE "Stage 2 Slope, +- vs. --"   INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1 +2    S2*A2  0    S2*A1*A2 -2;
    ESTIMATE "Stage 2 Slope, -+ vs. --"   INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1  0    S2*A2 +2    S2*A1*A2 -2;
    /* Areas under the curve: */
    /* The area under the curve between time 1 and time 2 is approximated as the area of a trapezoid with base
       length 1, left-side height E(Y1), and right-side height E(Y2).  The area under the curve between time 2
       and time 3 is approximated as the area of a trapezoid with base 1, left-side height E(Y2), and right-
       side height E(Y3).  Thus, the total area is (E(Y1)+E(Y2))/2 + (E(Y2)+E(Y3))/2 = 0.5*E(Y1) + E(Y2) + 0.5*E(Y3).
       The coefficients representing E(Y1), E(Y2), and E(Y3) for each intervention are given above, and the
       resulting linear combinations of them follow by linear algebra or by substituting in and simplifying 
       by hand. */   
    ESTIMATE "Area under Curve, ++"       INTERCEPT 2  S1 1.5  S2 0.5  S1*A1 +1.5  S2*A1 +0.5  S2*A2 +0.5  S2*A1*A2 +0.5;
    ESTIMATE "Area under Curve, +-"       INTERCEPT 2  S1 1.5  S2 0.5  S1*A1 +1.5  S2*A1 +0.5  S2*A2 -0.5  S2*A1*A2 -0.5;
    ESTIMATE "Area under Curve, -+"       INTERCEPT 2  S1 1.5  S2 0.5  S1*A1 -1.5  S2*A1 -0.5  S2*A2 +0.5  S2*A1*A2 -0.5;
    ESTIMATE "Area under Curve, --"       INTERCEPT 2  S1 1.5  S2 0.5  S1*A1 -1.5  S2*A1 -0.5  S2*A2 -0.5  S2*A1*A2 +0.5;    
    /* Contrasts in areas under the curve: */
    /* These are the pairwise differences in the corresponding coefficients of the pairs of rows above.*/
    ESTIMATE "Area, ++ vs. +-"            INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1  0    S2*A2 -1    S2*A1*A2 +1;
    ESTIMATE "Area, ++ vs. -+"            INTERCEPT 0  S1 0    S2 0    S1*A1  3    S2*A1  1    S2*A2  0    S2*A1*A2 +1;
    ESTIMATE "Area, ++ vs. --"            INTERCEPT 0  S1 0    S2 0    S1*A1  3    S2*A1  1    S2*A2 +1    S2*A1*A2  0;
    ESTIMATE "Area, +- vs. -+"            INTERCEPT 0  S1 0    S2 0    S1*A1  3    S2*A1  1    S2*A2 -1    S2*A1*A2  0;
    ESTIMATE "Area, +- vs. --"            INTERCEPT 0  S1 0    S2 0    S1*A1  3    S2*A1  1    S2*A2  0    S2*A1*A2 -1;
    ESTIMATE "Area, +- vs. --"            INTERCEPT 0  S1 0    S2 0    S1*A1  3    S2*A1  1    S2*A2  0    S2*A1*A2 -1;
	/* Average values under the curve (=area under curve/time): */
	/* Our example is two time units (months) long, that is, from time one to time three.  Therefore, the
	   average expected response value over this time period is the area under the curve divided by two.  */
    ESTIMATE "Average value, ++"          INTERCEPT 1  S1 0.75 S2 0.25 S1*A1 +0.75 S2*A1 +0.25 S2*A2 +0.25 S2*A1*A2 +0.25;
    ESTIMATE "Average value, +-"          INTERCEPT 1  S1 0.75 S2 0.25 S1*A1 +0.75 S2*A1 +0.25 S2*A2 -0.25 S2*A1*A2 -0.25;
    ESTIMATE "Average value, -+"          INTERCEPT 1  S1 0.75 S2 0.25 S1*A1 -0.75 S2*A1 -0.25 S2*A2 +0.25 S2*A1*A2 -0.25;
    ESTIMATE "Average value, --"          INTERCEPT 1  S1 0.75 S2 0.25 S1*A1 -0.75 S2*A1 -0.25 S2*A2 -0.25 S2*A1*A2 +0.25;    
	/* Contrasts in average values under the curve: */
	/* Note that the p-values for these contrasts will be the same as the contrasts in areas under the curve,
	   because they are only rescaled versions of the same quantity.  However, the estimates will of course be rescaled. */
    ESTIMATE "Average value, ++ vs. +-"   INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1  0    S2*A2 -0.5  S2*A1*A2 +0.5;
    ESTIMATE "Average value, ++ vs. -+"   INTERCEPT 0  S1 0    S2 0    S1*A1  1.5  S2*A1  0.5  S2*A2  0    S2*A1*A2 +0.5;
    ESTIMATE "Average value, ++ vs. --"   INTERCEPT 0  S1 0    S2 0    S1*A1  1.5  S2*A1  0.5  S2*A2 +0.5  S2*A1*A2  0;
    ESTIMATE "Average value, +- vs. -+"   INTERCEPT 0  S1 0    S2 0    S1*A1  1.5  S2*A1  0.5  S2*A2 -0.5  S2*A1*A2  0;
    ESTIMATE "Average value, +- vs. --"   INTERCEPT 0  S1 0    S2 0    S1*A1  1.5  S2*A1  0.5  S2*A2  0    S2*A1*A2 -0.5;
    ESTIMATE "Average value, +- vs. --"   INTERCEPT 0  S1 0    S2 0    S1*A1  1.5  S2*A1  0.5  S2*A2  0    S2*A1*A2 -0.5;
	/* Delayed effects: */
    /* The delayed effect of the stage-1 treatment in a regimen, relative to the stage-1 treatment in another
       regimen, is the difference between the contrast of the time-3 means and the contrast of the time-2 means 
       for these regimens.  This is not calculated here if the stage-1 treatment is the same between regimens, 
       because it then simplifies to being the same as the contrast of the time-3 means.  The delayed effect
       is easiest to understand when the A2 assignments for the pair of regimen are the same but the A1 
       assignments differ (that is, for comparing ++ to -+ or comparing +- and --) but are also defined in 
       the other cases.*/
    ESTIMATE "Delayed Effect, ++ vs. -+"  INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1 +2    S2*A2  0    S2*A1*A2 +2;
    ESTIMATE "Delayed Effect, ++ vs. --"  INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1 +2    S2*A2 +2    S2*A1*A2  0;    
    ESTIMATE "Delayed Effect, +- vs. -+"  INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1 +2    S2*A2 -2    S2*A1*A2  0;  
    ESTIMATE "Delayed Effect, +- vs. --"  INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1 +2    S2*A2  0    S2*A1*A2 -2;
    ESTIMATE "Ave. Delayed Eff., + vs -"  INTERCEPT 0  S1 0    S2 0    S1*A1  0    S2*A1 +2    S2*A2  0    S2*A1*A2  0;
RUN;
 
