C     Probability distribution parameter estimation program.          
C     Based on work of Kite.
C     Modified  to use IMSL routines (attached at bottom) for goodness of
C     fit tests (April 1991).
C
      PROGRAM KITE   

C Programming
C     The main program of FA1 is KITE.  The purpose of KITE is  to  gather
C     the  information  and  send  it ot the appointed option.  Each option is
C     placed in a seperate subroutine.  The "option"  subroutines  call  other 
C     "supporting"  subroutines  or  functions.   The  following  list  is the
C     programs, subroutines and functions in order of their appearance in FA1. 
C     Common  variables  and  distinguised  variables  are  discused after the 
C     program listing.
C 
C 
C     KITE,main  program:  Contains all the interactive componments.  When the 
C          data  is  read  in,  it  will  check  non-negativity  and  non-zero
C          constraints.   If  data  is  acceptable, KITE will sort the data in 
C          descending order.  Each sorted datum  is  then  assigned  a  sample
C          probability  of  precedence  using  equation  7.   The  sum, sum of
C          squares and sum of cubes of the data and of the  natural  logarithm 
C          of the data is calculated.  Arrays which will contain the estimated
C          parameters and the G-O-F test statistics  are  all  initialized  to 
C          zero.   The  selected  return  periods are place in the appropriate 
C          arrays (TT = 2, 5, 10, 20, 50, 100).  All information  is  sent  to 
C          the designated option subroutine. 
C 
C     GENR,option subroutine: This subroutine is called by  the  main  program
C          KITE.  Subroutine GENR first calls all the other option subroutines 
C          in the order STAT, TNMR, LN2R, LN3R, T1ER,  PT3R,  AND  LP3R.   The
C          estimated  parameters  from  all  the  option  subroutines are then 
C          tabulated as well as all the G-O-F test statistics.  The G-O-F test 
C          statistics  are  ranked  according  to  the best fit and matched to
C          their corresponding distribution.
C 
C     STAT,option  subroutine:  This  subroutine is called by the main program
C          KITE and by the option  subroutine  GENR.   Subroutine  STAT  first
C          prints  out  the sorted data along with their sample probability of 
C          precedence.  It then prints  out  the  unbiased  and  biased  (when
C          appropriate)  forms  of  common sample statistics; sum, square sum,
C          cubed sum, mean, median, variance, standard deviation,  coefficient
C          of    variation,   skewness   coefficient,   coeffient   of   exess
C          (kutosis - 3.0).  These statistics are calculated for the data  and 
C          the natural logarithm of the data. 
C 
C     TNMR,option subroutine: This subroutine is called by  the  main  program
C          KITE and the option subroutine GENR.  Subroutine TNMR estimates and 
C          tests the parameters for the normal distribution.   The  method  of
C          moments  and  the  maximum likelihood proceedure for estimating the
C          parameters for the normal distribution use the same equations.
C 
C     LN2R,option  subroutine:  This  subroutine is called by the main program
C          KITE and the option subroutine GENR.  Subroutine LN2R estimates and 
C          tests  the parameters for the two-parameter lognormal distribution. 
C          The method of moments and maximum likelihood procedure are used  to 
C          estimate the parameters.
C 
C     LN3R,option subroutine: This subroutine is called by  the  main  program
C          KITE and the option subroutine GENR.  Subroutine LN3R estimates and 
C          tests   the   parameters   for   the   three-parameter    lognormal 
C          distribution.    The  method  of  moments  and  maximum  likelihood 
C          procedure are used to estimate the parameters.
C 
C     T1ER,option  subroutine:  This  subroutine is called by the main program
C          KITE and the option subroutine GENR.  Subroutine T1ER estimates and 
C          tests the parameters for the type I extremal (Gumbel) distribution. 
C 
C          The method of moments and maximum likelihood procedure are used  to 
C          estimate the parameters.
C 
C     PT3R,option subroutine: This subroutine is called by  the  main  program
C          KITE and the option subroutine GENR.  Subroutine PT3R estimates and 
C          test the parameters for the Pearson  type  III  distribution.   The 
C          method  of  moments  and  maximum  likelihood procedure are used to 
C          estimate the parameters.
C 
C 
C     LP3R,option  subroutine:  This  subroutine is called by the main program
C          KITE and the option subroutine GENR.  Subroutine LP3R estimates and 
C          tests  the parameters for the log-Pearson type III distribution.  A 
C          direct and indirect method of moments is used, as well as a maximum
C          likelihood procedure, to estimate the parameters.
C 
C      GOF,supportive subroutine: This subroutine  is  called  by  the  option 
C          subroutines TNMR, LN2R, TN3R, T1ER, PT3R, AND LP3R.  Subroutine GOF
C          calculates the chi-square and  K-S  G-O-F  statistics.   The  data, 
C          class  intervals  and theoretical cumulative data is passed to GOF.
C          The chi-square statistic,  the  degrees  of  freedom,  the  1-ALPHA
C          confidence  level  and  the  K-S  statistic  are passed back to the
C          calling routine. 
C 
C     ORDER,supportive  subroutine:  This  subroutine  is called by the option 
C          subroutine GENR.  Subroutine ORDER ranks the G-O-F test statistics. 
C          The  statistics  are  passed  to ORDER as well as a flag indicating
C          descending  or  ascending  order.   The  ranked  data   and   thier
C          distribution character codes are passed back. 
C 
C     SNDV,supportive  function:  This  function  is  used   in   the   option 
C          subroutines TNMR, LN2R, LN3R, PT3R AND LP3R.  Function SNDV is used
C          to approximate the standard normal deviate given a  probability  of
C          precedence (Abramowitz,1965). 
C 
C     SORT,supportive subroutine:  This  subroutine  is  called  by  the  main 
C          program  KITE.   Subroutine SORT ranks the data in descending order 
C          using a "shell" type method.  The data and the number  of  data  is 
C          passed  to  the  subroutine.  The ranked data is passed back to the
C          calling routine via the origional data array.
C 
C     TINV,supportive   function:  This  function  is  called  by  the  option 
C          subroutine TNMR, LN2R, LN3R, PT3R, AND LP3R.  Function TINV is used 
C          to  approximate  the  probability  of  precedence  given a standard 
C          normal deviate (Abramowitz, 1965). 
C 
C 
C     VARIABLES 
C         Six common statements are passed between the main  program  and  the
C     option subroutines.
C 
C           . COMMON /ONE/ N,X(1500),P(1500),TT(6)  
C           . COMMON /TWO/ TITLE
C           . COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,s4,sl4 
C           . COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)
C           . COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)
C           . COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6) 
C 
C     The arrays designated by (3,6) are used to indicate simular variable for 
C     the different distributions.  For variable (I,J) these designations are
C 
C                      I   J     DISTRITUTION     ESTIMATING  
C                                                 PROCEDURE 
C 
C                      1   1       NORMAL          MM & ML
C 
C                      1   2     2-P LOGNORMAL        MM
C                      2   2     2-P LOGNORMAL        ML
C 
C                      1   3     3-P LOGNORMAL        MM
C                      2   3     3-P LOGNORMAL        ML
C 
C                      1   4     T 1 EXTREMAL         MM
C                      2   4     T 1 EXTREMAL         ML
C 
C 
C                      1   5      PEARSON T 3         MM
C                      2   5      PEARSON T 3         ML
C 
C                      1   6    LOG-PEARSON T 3   MM(DIRECT)  
C                      2   6    LOG-PEARSON T 3   MM(INDIRECT)
C                      3   6    LOG-PEARSON T 3       ML
C 
C     There  are  also  a  number  of variables that are designated to contain 
C     certain values only.  All these variable are given below.
C
C     A = estimated parameter of a distribution. 
C     ALPHA = probability of preceding the chi-square G-O-F statistic. 
C     APAR(3,6) = array containing the estimated A parameters. 
C     B = estimated parameter of a distribution. 
C     BPAR(3,6) = array containing the estimated B parameters. 
C     C = estimated parameter of a distribution. 
C     CHI(3,6) = array containing the chi-square G-O-G statistics. 
C     CL(302) = array containing the lower value of designated class limits. 
C     CON(3,6) = array containing the 1-ALPHA confidence levels.
C     CPAR(3,6) = array containing the C parameters.
C     CPR = the present cumulative probability level. 
C     DF = degrees of freedom for the chi-square G-O-F statistic = N-PAR-1. 
C     DN = K-S G-O-F statistic.
C     FN1 = character variable containing the data input file name.
C     FN2 = character variable containing the FA1 output. 
C     I = counter variable for do loops. 
C     N = number of data. 
C     NCL = number of class intervals.
C     PAR = real variable containing the number of parameters being estimated. 
C     PP(1500) = sample probability of precedence corresponding to data array. 
C     PR = probability of a class interval.
C     RKS(3,6) = array containing the K-S G-O-F statistics. 
C     S = sum of the partial standard error statistics. 
C     S1 = sum of the data. 
C     S2 = sum of the square of the data.
C     S3 = sum of the cube of the data.
C     SL1 = sum of the natural logarithm of the data.
C     SL2 = sum of the squares of the natural logarithm of the data. 
C     SL3 = sum of the cube of the natural logarithm of the data. 
C     SE = standard error statistic.
C     SPACE = character variable containing a single blank (" ").
C     ST(6) = array containing standard error for the respective return period.
C     STE(3,6) = array containing the standard error statistics. 
C     T = standard normal deviate (may by used otherwise)
C     TITLE = character variable containing the title of the computer run. 
C     TT(6) = array containing the six selected return periods. 
C     X(1500) = array containing the data.
C     X2 = chi-square statistic.
C     XG = skewness coefficient for the data.
C     XM = sample mean of the data.
C     XN = real variable containing the number of data.
C     XSD = unbiased sample standard deviation of the data. 
C     XT(6) = event corresponding to the appropriate return period. 
C     XV = unbiased sample variance of the data. 
C     XX(1500) = theoretical sorted data from the distribution. 
C     YM = sample mean of the natural logarithm of the data.
C     YSD = unbiased sample standard deviation of the natural log of the data. 
C     YV = unbiased sample variance of the natural logarithm of the data. 
C 
C
C-------------------------------------------------------------------------
C
C The following is a comparison between the way that Kites original program
C computed parameters, and the way that Baxter modified them (to make them
C consistent) in this modified version of Kite's program.
C
C     KITE'S
C     METHOD OF MOMENTS PARAMETER ESTIMATION
C 
C     TWO PARAMETER LOGNORMAL 
C 
C     THE ENTIRE PROCEDURE IS PERFORMED ON THE RAW DATA.
C     THE PARAMETERS, THE MEAN, VARIANCE AND SKEW, ARE CALCULATED IN THE
C     USUAL MANNER. THE COEFFICIENT OF SKEW IS BIASED AND INCORPORATES
C     AN UNBIASED VARIANCE. THE FREQUENCY FACTOR IS CALCULATED WITH 
C     EQUATION 6-29. THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS
C     CALCULATED WITH EQUATION 6-31. THE STANDARD ERROR OF THIS ESTIMATE
C     IS CALCULATED WITH EQUATION 6-40. 
C 
C     THREE PARAMETER LOGNORMAL 
C 
C     THE COEFFICIENT OF SKEW IS BIASED AND INCORPORATES A BIASED 
C     VARIANCE. THE VARIANCE IS THEN CORRECTED FOR BIAS AND THIS CORRECT
C     FORM IS USED THROUGHOUT THE REST OF THE ROUTINE. THE PARAMETERS ARE
C     CALCULATED WITH EQUATIONS 7-12, 7-16 AND 7-17. THE FREQUENCY FACTOR
C     IS CALCULATED WITH EQUATION 7-32. THE EVENT MAGNITUDE AT A
C     GIVEN RETURN PERIOD IS CALCULATED WITH EQUATION 7-30. THE 
C     STANDARD ERROR OF THIS ESTIMATE IS CALCULATED WITH EQUATION 7-37. 
C 
C     PEARSON TYPE THREE
C     THE COEFFICIENT OF SKEW IS UNBIASED (VIA EQUATION 9-22) AND INCORPATES
C     A BIASED VARIANCE. THE VARIANCE IS THEN CORRECTED FOR BIAS AND THIS
C     FORM IS USED THROUGHOUT THE REST OF THE ROUTINE. THE PARAMETERS 
C     ARE CALCULATED WITH EQUATIONS 9-23, 9-24 AND 9-25. THE FREQUENCY
C     FACTOR IS CALCULATED WITH EQUATION 9-52. THE EQUATION IN THE BOOK 
C     IN ERROR, SINCE THE 
C     LAST TERM IN THIS EQUATION SHOULD READ -(1/3)*(SKEW/6)**5.
C     THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED WITH 
C     EQUATION 3-31.  THE STANDARD ERROR OF THIS ESTIMATE IS CALCULATED 
C     WITH EQUATION 9-56. 
C 
C     LOGPEARSON TYPE THREE 
C 
C     DIRECT APPLICATION OF LOGPEARSON TYPE THREE 
C     THE SKEW IS BIASED AND INCORPORATES A BIASED VARIANCE.
C     A BIASED VARIANCE IS USED IN THIS ROUTINE UNTIL AFTER THE COEFFICIENT
C     OF SKEW IS CALCULATED FOR THE INDIRECT METHOD, AFTER WHICH THE UNBIASED
C     FORM IS USED. THE MOMENTS ARE CALCULATED WITH EQUATIONS 10-14,
C     10-15 AND 10-16. THE PARAMETERS ARE CALCULATED WITH EQUATIONS 10-1
C     10-13 AND 10-8. THE FREQUENCY FACTOR IS CALCULATED WITH 
C     EQUATION 9-52. THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS
C     CALCULATED WITH EQUATION 10-21, WHERE THE RIGHT HAND SIDE OF THE
C     EQUATION IS RAISED TO THE EXPONENTIAL FUNCTION, E. THE STANDARD ERROR
C     OF THIS ESTIMATE IS CALCULATED WITH EQUATION 10-22 VIA 9-56.
C 
C     APPLICATION OF PEARSON TYPE THREE TO THE LOGARITHM OF THE DATA
C     THE NATURAL LOGARITHM OF EACH DATAPOINT IS TAKEN. THE COEFFICIENT 
C     OF SKEW IS UNBIASED (VIA EQUATION 9-22) AND INCORPORATES A BIASED 
C     VARIANCE. THE PARAMETERS ARE CALCULATED WITH EQUATIONS 9-23, 9-24 
C     AND 9-25. THE FREQUENCY FACTOR IS CALCULATED WITH EQUATION 9-52.
C     THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED
C     WITH EQUATION 10-21, WHERE THE RIGHT HAND SIDE IS RAISED TO THE 
C     EXPONENTIAL FUNCTION, E. THE STANDARD ERROR OF THIS ESTIMATE
C     IS CALCULATED WITH EQUATION 10-22 VIA 9-56. 
C 
C     EXTREMAL TYPE ONE 
C 
C     THE COEFFICIENT OF SKEW IS BIASED AND INCORPORATES AN UNBIASED
C     VARIANCE. THE UNBIASED VARIANCE IS USED THROUGHOUT THE ENTIRE 
C     ROUTINE. THE PARAMETERS ARE CALCULATED WITH EQUATIONS 8-21 AND 8-2
C     THE FREQUENCY FACTOR IS CALCULATED WITH EQUATION 8-54, WHICH
C     INCORPORATES THE Y CALCULATED IN EQUATION 8-46. THE EVENT MAGNITUDE
C     AT A GIVEN RETURN PERIOD IS CALCULATED WITH EQUATION 8-53. THE
C     STANDARD ERROR OF THIS ESTIMATE IS CALCULATED WITH EQUATION 8-61. 
C 
C     MAXIMUM LIKELIHOOD PARAMETER ESTIMATES
C 
C     TWO PARAMETER LOGNORMAL 
C 
C     THE NATURAL LOGARITHM IS TAKEN OF EACH POINT. THE SKEW IS BIASED ARE
C     INCORPORATES AN UNBIASED VARIANCE. THE FREQUENCY FACTOR IS CALCULATED
C     WITH EQUATION 6-29. THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS
C     CALCULATED WITH EQUATION 6-43. THE STANDARD ERROR OF THIS ESTIMATE
C     IS CALCULATED WITH EQUATION 6-58. THIS APPROACH IS A MODIFIED 
C     VERSION OF THE ORGINAL PROGRAM. THE ORGINAL USED THE METHOD 
C     OF MOMENTS PARAMETERS TO CALCULATE THE EVENT MAGNITUDE AT A 
C     GIVEN RETURN PERIOD AND THE STANDARD ERROR OF THIS ESTIMATE.
C     THE EQUATIONS USED TO CALCULATE THE EVENT MAGNITUDE AT A GIVEN
C     RETURN PERIOD AND THE STANDARD ERROR OF THIS ESTIMATE HAVE
C     ALSO BEEN CHANGED. THE ORGINAL VERSION USED EQUATIONS 
C     6-31 AND 6-62 RESPECTIVELY. 
C 
C     THREE PARAMETER LOGNORMAL 
C 
C     THE SKEW IS BIASED AND INCORPORATES AN UNBIASED VARIANCE. 
C     THE PARAMETERS ARE CALCULATED WITH EQUATIONS 7-19, 7-20, AND
C     7-21. THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED
C     WITH EQUATION 7-29. THE STANDARD ERROR OF THIS ESTIMATE IS CALCULATED
C     WITH EQUATION 7-57. 
C 
C     PEARSON TYPE THREE
C 
C     THE MOMENTS ARE CALCULATED WITH EQUATIONS 10-14, 10-15 AND 10-16. 
C     THE PARAMETERS ARE CALCULATED WITH EQUATIONS 9-32, 9-33 AND 
C     9-28 . THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS
C     CALCULATED WITH EQUATION 9-49. THE STANDARD ERROR OF THIS ESTIMATE
C     IS CALCULATED WITH EQUATION 9-57.
C 
C     LOGPEARSON TYPE THREE 
C 
C     THE MOMENTS ARE CALCULATED WITH EQUATIONS 10-14, 10-15 AND 10-16. 
C     THE PARAMETERS ARE CALCULATED WITH EQUATIONS 9-32, 9-33 AND 
C     9-28 . THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS
C     CALCULATED WITH EQUATION 9-49, WHERE THE RIGHT HAND SIDE IS RAISED
C     TO THE EXPONENTIAL FUNCTION, E. THE STANDARD ERROR OF THIS
C     ESTIMATE IS CALCULATED WITH EQUATION 10-22 VIA EQUATION 9-57. 
C 
C     TYPE ONE EXTREAMAL
C 
C     THE PARAMETERS ARE CALCULATED WITH EQUATIONS 8-32 AND 8-40. THE 
C     EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED WITH EQUATION
C     8-64 (THE SAME AS BAXTER'S METHOD OF MOMENTS ESTIMATE). THE 
C     STANDARD ERROR OF THIS ESTIMATE IS CALCULATED WITH EQUATION 8-77. 
C 
C     BAXTER'S
C     IN ALL CASES THE SKEW IS UNBIASED WITH EQUATION 3-38 IN 
C     STATISTICAL METHODS IN HYDROLOGY BY C.T. HAAN.
C 
C     METHOD OF MOMENTS 
C 
C     NORMAL
C     THE SKEW IS UNBIASED AND INCORPORATES AN UNBIASED VARIANCE. 
C     THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED WITH 
C     EQUATION 5-11. THE STANDARD ERROR OF THIS ESTIMATE IS CALCULATED
C     WITH EQUATION 5-19. 
C 
C     TWO PARAMETER LOGNORMAL 
C 
C     THE COEFFICIENT OF SKEW IS UNBIASED AND INCORPORATES AN 
C     UNBIASED VARIANCE. THE PARAMETERS ARE CALCULATED WITH EQUATIONS 
C     6-30, AND 6-31 IN STATISTICAL METHODS IN HYDROLOGY BY C.T. HAAN 
C     THE FREQUENCY FACTOR IS CALCULATED WITH 6-29. THE EVENT MAGNITUDE 
C     AT A GIVEN RETURN PERIOD IS CALCULATED WITH EQUATION 6-30. THE
C     STANDARD ERROR OF THIS ESTIMATE IS CALCULATED WITH EQUATION 6-40. 
C 
C     THREE PARAMETER LOGNORMAL 
C 
C     THE COEFFICIENT OF SKEW IS UNBIASED AND INCORPORATES AN 
C     UNBIASED VARIANCE. THE PARAMETERS ARE CALCULATED WITH EQUATIONS 
C     7-12, 7-16 AND 7-17. THE FREQUENCY FACTOR IS CALCULATED WITH
C     EQUATION 7-32. THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS
C     CALCULATED WITH EQUATION 7-29. THE STANDARD ERROR OF THIS 
C     ESTIMATE IS CALCULATED WITH EQUATION 7-37.
C 
C     PEARSON TYPE THREE
C 
C     THE COEFFICIENT OF SKEW IS UNBIASED AND INCORPORATES AN UNBIASED
C     VARIANCE. THE PARAMETERS ARE CALCULATED WITH EQUATIONS 9-23, 9-24 
C     AND 9-25. THE FREQUENCY FACTOR IS CALCULATED WITH 9-52. THE EVENT 
C     MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED WITH EITHER EQUATION
C     9-45 OR 3-31. THE STANDARD ERROR OF THIS ESTIMATE IS CALCULATED WITH
C     EQUATION 9-56.
C 
C     LOGPEARSON TYPE THREE 
C 
C     DIRECT APPLICATION OF LOGPEARSON TYPE THREE 
C     THE COEFFICIENT OF SKEW IS UNBIASED AND INCORPORATES AN 
C     UNBIASED VARIANCE. THE PARAMETERS ARE CALCULATED WITH EQUATIONS 
C     10-8, 10-12, AND 10-13. THE FREQUENCY FACTOR IS CALCULATED WITH 
C     EQUATION 9-52. THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS
C     CALCULATED WITH EITHER EQUATION 9-45 OR 10-21 , WHERE THE RIGHT HAND
C     SIDE OF THE EQUATION IS RAISED TO THE EXPONENTIAL FUNCTION, E.G.
C     THE STANDARD ERROR OF THIS ESTIMATE IS CALCULATED WITH EQUATION 
C     10-22 VIA 9-56. 
C 
C     APPLICATION OF PEARSON TYPE THREE TO THE LOGARITHM OF DATA. 
C     THE NATURAL LOGARITHM IS TAKED OF EACH DATA POINT. THE COEFFICIENT
C     OF SKEW IS UNBIASED AND INCORPORATES AN UNBIASED VARIANCE. THE
C     PARAMETERS ARE CALCULATED WITH EQUATIONS 9-23, 9-24 AND 9-25. 
C     THE FREQUENCY FACTOR IS CALCULATED WITH EQUATION 9-52. THE
C     EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED WITH EITHER
C     EQUATION 9-45 OR 10-21, WHERE THE RIGHT HAND SIDE IS RAISED TO THE
C     EXPONENTIAL FUNCTION, E. THE STANDARD ERROR OF THIS ESTIMATE
C     IS CALCULATED WITH EQUATION 10-22 VIA 9-56. 
C 
C     TYPE ONE EXTREMAL 
C 
C     THE PARAMETERS ARE CALCULATED WITH EQUATIONS 8-21 AND 8-22. 
C     THE FREQUENCY FACTOR IS CALCULATED WITH EQUATION 8-57. THE
C     EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED WITH 
C     EQUATION 8-55. THE STANDARD ERROR OF THIS ESTIMATE IS CALCULATED
C     WITH EQUATION 8-61. 
C 
C     MAXIMUM LIKELIHOOD PARAMETER ESTIMATES
C 
C     NORMAL
C 
C     SAME AS THE METHOD OF MOMENTS ESTIMATE. 
C 
C     TWO PARAMETER LOGNORMAL 
C 
C     THE COEFFICIENT OF SKEW IS UNBIASED AND INCORPORATES AN UNBIASED
C     VARIANCE.THE PARAMETERS,THE MEAN AND VARIANCE, ARE CALCULATED 
C     IN THE USUAL MANNER. THE FREQUENCY FACTOR IS CALCULATED WITH
C     EQUATION 6-29. THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS
C     CALCULATED WITH EQUATION 6-43. THE STANDARD ERROR OF THIS ESTIMATE
C     IS CALCULATED WITH EQUATION 6-62. BAXTER'S RESULTS FOR THE STANDARD
C     ERROR DIFFER FROM KITE'S BECAUSE BAXTER USES THE VARIANCE 
C     OF THE RAW DATA. THIS VARIANCE DIFFERS FROM THE VARIANCE
C     OBTAINED BY BACK CALCULATING FROM THE PARAMETERS (THE MEAN AND
C     VARIANCE OF THE LOGGED DATA). THIS BACK CALCULATION COULD BE
C     PERFORMED WITH EQUATION 6-14. 
C 
C     THREE PARAMETER LOGNORMAL 
C 
C     THE COEFFICIENT OF SKEW IS UNBIASED AND INCORPORATES AN UNBIASED
C     VARIANCE. THE PARAMETERS ARE CALCULATED WITH EQUATIONS 7-19, 7-20 
C     AND 7-21. THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS 
C     CALCULATED WITH EQUATION 7-29. THE STANDARD ERROR OF THIS ESTIMATE
C     IS CALCULATED WITH EQUATION 7-57. 
C 
C     PEARSON TYPE THREE
C 
C     THE PARAMETERS ARE CALCULATED WITH EQUATIONS 9-32 AND 9-33. 
C     THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED WITH 
C     EQUATION 9-45 OR 9-49. THE STANDARD ERROR OF THIS ESTIMATE IS 
C     CALCULATED WITH EQUATION 9-57.
C 
C     LOGPEARSON TYPE THREE.
C     THE PARAMETERS ARE CALCULATED WITH EQUATIONS 9-32, 9-33 AND 
C     9-28 . THE EVENT MAGNITUDE AT A GIVEN RETURN PERIOD IS
C     CALCULATED WITH EITHER EQUATION 9-45 OR 9-49, WHERE THE RIGHT HAND
C     SIDE IS RAISED TO THE EXPONENTIAL FUNCTION, E. THE STANDARD 
C     ERROR OF THIS ESTIMATE IS CALCULATED WITH EQUATION 10-22 VIA 9-57.
C 
C     EXTREMAL TYPE ONE 
C 
C     THE PARAMETERS ARE CALCULATED WITH 8-32 AND 8-40. THE EVENT 
C     MAGNITUDE AT A GIVEN RETURN PERIOD IS CALCULATED WITH EQUATION
C     8-64. THE STANDARD ERROR OF THIS ESTIMATE IS CALCULATED WITH EQUATION
C     8-78. 
C 
      COMMON /ONE/ N,X(1500),PP(1500),TT(6)        
      COMMON /TWO/ TITLE                  
      COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,S4,SL4     
      COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)   
      COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)                  
      COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6),PX(1500)                
      CHARACTER TITLE*72                  
      CHARACTER ANS*1,LINE*71,FN1*60,FN2*60        
C           
C *** MODIFIED FORTRAN 5 VERSION OF THE KITE PROGRAMS                   
C *** MAIN PROGRAM - SELECTION OF INCOMING DATA    
C           
C *** DATA ENTRY EITHER BY FILE OF TERMINAL INPUT  
1     CONTINUE       
C     PRINT*,'IS THE DATA CONTAINED IN A LOCAL FILE OR DO YOU WISH'     
C     PRINT*,'TO GENERATE A DATA FILE (L=LOCAL, G=GENERATE)'            
C     READ 100,ANS   
      PRINT 2100,'Name of the data file? '
2100  FORMAT(A,$)
      READ 100,FN1   
      OPEN(12,FILE=FN1,STATUS='OLD')      
      REWIND 12      
      PRINT 2100,'Number of observations? '
      READ*,N        
      READ(12,*)(X(I),I=1,N)              
      close(12)
C     IF(ANS.EQ.'L')THEN                  
C         READ(2,*)(X(I),I=1,N)           
C     ELSE  
C         PRINT*,'ENTER YOUR DATA.  SEPERATE EACH VALUE BY A COMMA,'    
C         PRINT*,'BLANK, OR CARRIAGE RETURN'       
C         READ*,(X(I),I=1,N)              
C         WRITE(2,*)(X(II),II=1,N)        
C         PRINT*,'THANKS, YOUR DATA IS NOW IN LOCAL FILE ',FN1          
C     END IF
C           
C *** SORT DATA IN DESCENDING ORDER - X(1) LARGEST VALUE                
      CALL SORT (X,N)
C           
C *** SET UP PROBABILITY FOR PLOTTING POSITION     
C *** CHECK FOR ZERO OR NEGATIVE VALUES IN DATA    
C *** FORM SUMMATIONS OF DATA, SQUARED-DATA AND CUBED-DATA              
C *** FORM SUMMATIONS OF LN(DATA), SQUARED-LN(DATA) AND CUBED-LN(DATA)  
2     CONTINUE       
      XN=N*1.        
      S1=0.0
      S2=0.0
      S3=0.0
      S4=0.0
      SL1=0.0        
      SL2=0.0        
      SL3=0.0        
      SL4=0.0 
      DO 20 I=1,N    
          XI=I*1.    
          PP(I)=1.-XI/(XN+1.)             
          IF(X(I).EQ.0)THEN               
              print 101                   
              STOP   
          ELSE IF(X(I).LT.0)THEN          
              print 102
              STOP   
          END IF     
          S1=S1+X(I) 
          S2=S2+X(I)*X(I)                 
          S3=S3+X(I)*X(I)*X(I)            
          S4=S4+(X(I)**4) 
          Y=ALOG(X(I))                    
          SL1=SL1+Y  
          SL2=SL2+Y*Y
          SL3=SL3+Y*Y*Y                   
          SL4=SL4+(Y**4)
20    CONTINUE       
C           
C *** INITIALIZE /FOR/ AND /FIV/ VARIABLES
      DO 30 I=1,3    
          DO 31 II=1,6                    
              APAR(I,II)=0.0              
              BPAR(I,II)=0.0              
              CPAR(I,II)=0.0              
              CHI(I,II)=0.0               
              CON(I,II)=0.0               
              STE(I,II)=0.0               
              RKS(I,II)=0.0               
31        CONTINUE   
30    CONTINUE       
C           
C *** ESTABLISH SELECTIVE RETURN PERIODS  
      TT(1)=2.       
      TT(2)=5.       
      TT(3)=10.      
      TT(4)=20.      
      TT(5)=50.      
      TT(6)=100.     
C           
C *** INPUT TITLE FOR PROGRAM RUN         
      print*,' '
      PRINT 100,'Title of this run (71 characters maximum):'            
      READ 100,TITLE 
      print*,' '
C           
C *** ESTABLISH OUTPUT, EITHER FILE OR TERMINAL    
c     PRINT*,'DO YOU WISH THE OUTPUT TO GO TO THE TERMINAL OR A'        
c     PRINT*,'LOCAL FILE (T=TERMINAL, L=LOCAL FILE)'                    
c     READ 100,ANS   
c     IF(ANS.EQ.'L')THEN                  
c         PRINT*,'OUTPUT FILE NAME'       
c         READ 100,FN2                    
c     ELSE  
c         FN2='OUTPUT'                    
c     ENDIF 
      print 2100,'Output file name? '
      read 100,fn2
      OPEN(2,FILE=FN2)                    
C           
C *** ESTABLISHED DESIRED RUN OPTION      
3     PRINT 103      
      PRINT 2100,'Option desired for this run? '    
      READ*,NOP      
      IF(NOP.EQ.1)THEN                    
          CALL GENR  
      ELSE IF(NOP.EQ.2)THEN               
          CALL STAT  
      ELSE IF(NOP.EQ.3)THEN               
          CALL TNMR  
      ELSE IF(NOP.EQ.4)THEN               
          CALL LN2R  
      ELSE IF(NOP.EQ.5)THEN               
          CALL LN3R  
      ELSE IF(NOP.EQ.6)THEN               
          CALL T1ER  
      ELSE IF(NOP.EQ.7)THEN               
          CALL PT3R  
      ELSE IF(NOP.EQ.8)THEN               
          CALL LP3R  
      ELSE  
          PRINT*     
          PRINT*,'*** ERROR IN OPTION ASSIGNMENT, ',                    
     +           'MUST BE BETWEEN 1 AND 8'
          PRINT*     
          GOTO 3     
      END IF
      CLOSE(2)       
C           
C *** DECISION TO MADE ANOTHER PROGRAM RUN
      print*,' '
      PRINT 2100,'Do you wish another run (yes or no)? '                
      READ 100,ANS   
      IF(ANS.EQ.'N'.or.ANS.EQ.'n')STOP    
      PRINT 2100 ,
     *      'Are you going to use a different data set (yes or no)? '   
      READ 100,ANS   
      IF(ANS.EQ.'N'.or.ANS.EQ.'n')GO TO 2 
      GO TO 1        
C           
C *** FORMAT STATEMENTS                   
C           
100   FORMAT(A)      
101   FORMAT(///'*** ERROR IN DATA ***'/  
     +'THERE IS A ZERO VALUE IN YOUR DATA SERIES.  THIS PROGRAM IS'/    
     +'UNABLE TO PROCESS ZERO VALUES. YOUR OPTIONS ARE:'/               
     +5X,'1) ADD 1.0 TO ALL DATA.'/       
     +5X,'2) ADD SMALL POSITIVE VALUE (SUCH AS 0.1, 0.01, 0.001, ETC)'/ 
     +5X,'   TO ALL DATA.'/               
     +5X,'3) SUBSTITUTE 1.0 IN PLACE OF ALL ZERO VALUES.'/              
     +5X,'4) SUBSTITUTE SMALL POSITIVE VALUES FOR ALL ZERO VALUES.'/    
     +5X,'5) DELETE ALL ZERO OBSERVATIONS.'//)     
102   FORMAT(///'*** ERROR IN DATA ***'/  
     +'YOUR DATA CONTAINS A NEGATIVE VALUE.  THIS PROGRAM IS UNABLE '/  
     +'TO PROCESS NEGATIVE VALUES'//)     
103   FORMAT(/,'Program options are',//,  
     +       10X,'(1)  generalized frequency analysis'/                 
     +       10X,'(2)  data analysis'/    
     +       10X,'(3)  normal distribution'/       
     +       10X,'(4)  two-parameter lognormal distribution'/           
     +       10X,'(5)  three-parameter lognormal distribution'/         
     +       10X,'(6)  type 1 extremal distribution'/                   
     +       10X,'(7)  pearson type three distribution'/                
     +       10X,'(8)  log-pearson type three distribution'/)           
C           
      END   
C1
C*GENR******************************************************************
C 
      SUBROUTINE GENR
      COMMON /ONE/ N,X(1500),PP(1500),TT(6)        
      COMMON /TWO/ TITLE                  
      COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,S4,SL4     
      COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)   
      COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)                  
      COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6),PX(1500)                
      CHARACTER TITLE*72,C1(12)*7,C2(12)*7,PROC*12 
      REAL SS1(12),SS2(12)                
C           
C *** CALLS SUBROUTINE STAT,TNMR,LN2R,LN3R,T1ER,PT3R,LP3R               
C *** SUMMARIZES THE RESULTS OF ALL OPTIONS        
C           
      CALL STAT      
      CALL TNMR      
      CALL LN2R      
      CALL LN3R      
      CALL T1ER      
      CALL PT3R      
      CALL LP3R      
C           
C *** WRITE TITLE AND HEADING             
      WRITE(2,200)TITLE                   
C           
C *** WRITE PARAMETER HEADINGS            
      WRITE(2,210)   
      Z=0.0 
C           
C *** NORMAL PARAMETERS                   
      WRITE(2,220)APAR(1,1),BPAR(1,1)     
C           
C *** 2-PAR LOGNORMAL PARAMETERS          
      WRITE(2,230)   
      WRITE(2,291)APAR(1,2),BPAR(1,2)     
      WRITE(2,292)APAR(2,2),BPAR(2,2)     
C           
C *** 3-PAR LOGNORMAL PARAMETERS          
      WRITE(2,240)   
      IF(STE(1,3).LT.Z)THEN               
          WRITE(2,296)'MOMENT      '      
      ELSE  
          WRITE(2,291)APAR(1,3),BPAR(1,3),CPAR(1,3)
      END IF
      IF(STE(2,3).LT.Z)THEN               
          WRITE(2,296)'LIKELIHOOD  '      
      ELSE  
          WRITE(2,292)APAR(2,3),BPAR(2,3),CPAR(2,3)
      END IF
C           
C *** TYPE I EXTREMAL PARAMETERS          
      WRITE(2,250)   
      WRITE(2,291)APAR(1,4),BPAR(1,4)     
      IF(STE(2,4).LT.Z)THEN               
          WRITE(2,295)'LIKELIHOOD  '      
      ELSE  
          WRITE(2,292)APAR(2,4),BPAR(2,4) 
      END IF
C           
C *** PEARSON TYPE III PARAMETERS         
      WRITE(2,260)   
      WRITE(2,291)APAR(1,5),BPAR(1,5),CPAR(1,5)    
      IF(STE(2,5).LT.Z)THEN               
          WRITE(2,296)'LIKELIHOOD  '      
      ELSE  
          WRITE(2,292)APAR(2,5),BPAR(2,5),CPAR(2,5)
      END IF
C           
C *** LOG-PEARSON TYPE III PARAMETERS     
      WRITE(2,270)   
      IF(STE(1,6).LT.Z)THEN               
          WRITE(2,296)'MOMENT (D)  '      
      ELSE  
          WRITE(2,293)APAR(1,6),BPAR(1,6),CPAR(1,6)
      END IF
      WRITE(2,294)APAR(2,6),BPAR(2,6),CPAR(2,6)    
      IF(STE(3,6).LT.Z)THEN               
          WRITE(2,296)'LIKELIHOOD  '      
      ELSE  
          WRITE(2,292)APAR(3,6),BPAR(3,6),CPAR(3,6)
      END IF
C           
      WRITE(2,999)   
C           
C *** WRITE ORDERED STATISTIC'S HEADING   
      WRITE(2,300)TITLE                   
      Z=0.0 
C           
C *** WRITE(STANDARD ERROR AND K-S STATISTICS      
      WRITE(2,999)   
      CALL ORDER (STE,-1,C1,SS1)          
      CALL ORDER (RKS,-1,C2,SS2)          
      WRITE(2,310)   
      DO 10 I=1,12   
          IF(SS1(I).LT.Z.AND.SS2(I).LT.Z)THEN      
              WRITE(2,391)C1(I),C2(I)     
          ELSE IF(SS1(I).LT.Z)THEN        
              WRITE(2,392)C1(I),C2(I),SS2(I)       
          ELSE IF(SS2(I).LT.Z)THEN        
              WRITE(2,393)C1(I),SS1(I),C2(I)       
          ELSE       
              WRITE(2,394)C1(I),SS1(I),C2(I),SS2(I)
          END IF     
10    CONTINUE       
C           
C *** WRITE CHI-SQUARE AND (1-ALPHA) STATISTICS    
      WRITE(2,999)   
      CALL ORDER (CHI,-1,C1,SS1)          
      CALL ORDER (CON,+1,C2,SS2)          
      WRITE(2,320)   
      DO 20 I=1,12   
          IF(SS1(I).LT.Z.AND.SS2(I).LT.Z)THEN      
              WRITE(2,391)C1(I),C2(I)     
          ELSE IF(SS1(I).LT.Z)THEN        
              WRITE(2,392)C1(I),C2(I),SS2(I)       
          ELSE IF(SS2(I).LT.Z)THEN        
              WRITE(2,393)C1(I),SS1(I),C2(I)       
          ELSE       
              WRITE(2,394)C1(I),SS1(I),C2(I),SS2(I)
          END IF     
20    CONTINUE       
      WRITE(2,998)   
C           
      RETURN
C           
C *** FORMAT STATEMENTS                   
200   FORMAT(1H1/29X,'- DATA SUMMARY -'// 
     +1X,'TITLE: ',A72//////////)         
210   FORMAT(/44X,'PARAMETERS'/           
     +5X,'DISTRIBUTION',15X,'A',16X,'B',16X,'C')   
220   FORMAT(/1X,'NORMAL'/                
     +11X,'MOMENT',4X,3(5X,E12.5))        
230   FORMAT(/1X,'2-PAR LOGNORMAL')       
240   FORMAT(/1X,'3-PAR LOGNORMAL')       
250   FORMAT(/1X,'TYPE I EXTREMAL')       
260   FORMAT(/1X,'PEARSON TYPE III')      
270   FORMAT(/1X,'LOG-PEARSON TYPE III')  
291   FORMAT(11X,'MOMENT',4X,3(5X,E12.5)) 
292   FORMAT(11X,'LIKELIHOOD',3(5X,E12.5))
293   FORMAT(11X,'MOMENT (D)',5X,E12.5,2(5X,E12.5))
294   FORMAT(11X,'MOMENT (ID)',4X,E12.5,2(5X,E12.5))                    
295   FORMAT(11X,A12,2(5X,12('-')))       
296   FORMAT(11X,A12,3(5X,12('-')))       
300   FORMAT(1H1/25X,'ORDERED G-O-F STATISTICS'/1X,'TITLE; ',A72///)    
310   FORMAT(11X,'STANDARD ERROR',21X,'KOLMOGOROV-SMIRNOV'/             
     +5X,'DISTRIBUTION',5X,'STATISTIC',   
     +11X,'DISTIBUTION',5X,'STATISTIC')   
320   FORMAT(/43X,'CHI-SQUARED SIGNIFICANCE'/      
     +13X,'CHI-SQUARE',25X,'(1-ALPHA)'/   
     +5X,'DISTRIBUTION',5X,'STATISTIC',   
     +11X,'DISTRIBUTION',5X,'STATISTIC'/) 
391   FORMAT(8X,A7,5X,12('-'),13X,A7,5X,12('-'))   
392   FORMAT(8X,A7,5X,12('-'),13X,A7,5X,E12.5)     
393   FORMAT(8X,A7,5X,E12.5,13X,A7,5X,12('-'))     
394   FORMAT(8X,A7,5X,E12.5,13X,A7,5X,E12.5)       
998   FORMAT(//1X,25X,'LEGEND'///         
     +1X,'NRM = NORMAL DISTRIBUTION',20X,'MM = METHOD OF MOMENTS'/      
     +1X,'LN2 = 2 PAR LOGNORMAL DISTRIBUTION',11X, 
     +'M1 = METHOD OF MOMENTS (DIRECT), LPR ONLY'/ 
     +1X,'LN3 = 3 PAR LOGNORMAL DISTRIBUTION',11X, 
     +'M2 = METHOD OF MOMENTS (INDIRECT), LPR ONLY'/                    
     +1X,'T1E = TYPE I EXTREMAL DISTRIBUTION',11X, 
     +'LK = MAXIMUM LIKELIHOOD PROCEDURE'/
     +1X,'PRS = PEARSON TYPE 3 DISTRIBUTION'/      
     +1X,'LPR = LOG-PEARSON TYPE 3 DISTRIBUTION'///)                    
999   FORMAT(///)    
C           
      END   
C1
C*STAT******************************************************************
C 
      SUBROUTINE STAT
      COMMON /ONE/ N,X(1500),P(1500),TT(6)
      COMMON /TWO/ TITLE                  
      COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,S4,SL4     
      COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)   
      COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)                  
      COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6),PX(1500)                
      CHARACTER TITLE*72,SPACE*1          
      REAL Y(1500)   
C           
C *** OPTION 1       
C *** SORTS THE DATA, X(1) BEING THE LARGEST       
C *** PERFORMS BIASED AND UNBIASED STATISTICS ON THE DATA AND LN(DATA)  
C ***         1) SUM 
C ***         2) SQUARE SUM               
C ***         3) CUBED SUM                
C ***         4) MEAN
C ***         5) MEDIAN                   
C ***         6) VARIANCE                 
C ***         7) STANDARD DEVIATION       
C ***         8) COEFFICIENT OF VARIATION 
C ***         9) SKEWNESS COEFFICIENT     
C ***        10) COEFFICIENT OF EXCESS (KURTOSIS-3)
C           
C *** WRITE TITLE AND HEADINGS            
      WRITE(2,200)TITLE                   
      WRITE(2,210)   
C           
C *** WRITE OUT RANKED DATA               
      L=N   
      SPACE=' '      
1     IF(L-3.GT.0)THEN                    
          WRITE(2,220)(SPACE,P(I),X(I),I=L,L-3,-1) 
      ELSE  
          WRITE(2,220)(SPACE,P(I),X(I),I=L,1,-1)   
          GOTO 2     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 1                    
C           
2     CONTINUE       
      WRITE(2,230)N  
C           
      XN=N*1.        
C           
C *** MEAN  
      XM=S1/XN       
      YM=SL1/XN      
C           
C *** VARIANCE       
      XVB=(S2-S1*S1/XN)/XN                
      XVUB=XN*XVB/(XN-1.)                 
      YVB=(SL2-SL1*SL1/XN)/XN             
      YVUB=XN*YVB/(XN-1.)                 
C           
C *** STANDARD DEVIATION                  
      XSDB=SQRT(XVB) 
      XSDUB=SQRT(XVUB)                    
      YSDB=SQRT(YVB) 
      YSDUB=SQRT(YVUB)                    
C           
C *** COEFFICIENT OF VARIATION            
      XCVB=XSDB/XM   
      XCVUB=XSDUB/XM 
      YCVB=YSDB/YM   
      YCVUB=YSDUB/YM 
C           
C *** SKEW  
      XSK1=(S3)-(3.*S1*S2/XN)+(2.*S1*S1*S1/(XN*XN))
      XSKB=XSK1/XN   
      XSKUB=XSK1*XN/(XN-1.)/(XN-2.)       
      YSK1=(SL3)-(3.*SL1*SL2/XN)+(2.*SL1*SL1*SL1/(XN*XN))               
      YSKB=YSK1/XN   
      YSKUB=YSK1*XN/(XN-1.)/(XN-2.)       
C           
C *** SKEWNESS COEFFICIENT                
      XCSB=XSKB/(XSDB**3)                 
      XCSUB=XSKUB/(XSDUB**3)              
      YCSB=YSKB/(YSDB**3)                 
      YCSUB=YSKUB/(YSDUB**3)              
C           
C *** COEFFICIENT OF EXCESS (PEAKEDNESS)  
      S=(S4-(4.*(S1/N)*S3)+(6.*S2*(S1/N)**2)-(3.*N*(S1/N)**4))/XN       
      SL=(SL4-(4.*(SL1/N)*SL3)+(6.*SL2*(SL1/N)**2)-(3.*N*(SL1/N)**4))/XN
      XCEB=S/(XSDB**4)-3.                 
      XCEUB=(S*XN**3)/((XN-1.)*(XN-2.)*(XN-3.)*(XSDUB**4))-3.           
      YCEB=SL/(YSDB**4)-3.                
      YCEUB=(SL*XN**3)/((XN-1.)*(XN-2.)*(XN-3.)*(YSDUB**4))-3.          
C           
C *** MEDIAN VALUE   
      IF((N/2).EQ.(XN/2.))THEN            
          J=N/2      
          XMD=(X(J)+X(J+1))/2.            
          YMD=(ALOG(X(J))+ALOG(X(J+1)))/2.
      ELSE  
          J=(N/2)+1  
          XMD=X(J)   
          YMD=ALOG(X(J))                  
      END IF
C           
      WRITE(2,240)S1,S2,S3,XM,XMD,XVB,XVUB,XSDB,XSDUB,XCVB,XCVUB,       
     +             XSKB,XSKUB,XCSB,XCSUB,XCEB,XCEUB
C           
      WRITE(2,250)   
C           
      WRITE(2,240)SL1,SL2,SL3,YM,YMD,YVB,YVUB,YSDB,YSDUB,YCVB,YCVUB,    
     +             YSKB,YSKUB,YCSB,YCSUB,YCEB,YCEUB
C           
      WRITE(2,299)   
C           
      RETURN
C           
C *** FORMAT STATEMENTS                   
C           
200   FORMAT(1H1/30X,'- GENERAL STATISTICS -'//    
     +1X,'TITLE: ',A72)                   
210   FORMAT(//1X,'RANKED DATA WITH THE PROBABILITY OF PRECEDENCE ',    
     +'PARENTHESIS.'/
     +6X,'NOTE: PROBABILITY OF PRECEDENCE = RANK/(N+1).'/)              
220   FORMAT(4(A1,'(',F5.4,')',E11.5,2X)) 
230   FORMAT(/4X,'NUMBER OF DATA =',I5//  
     +29X,'DATA STATISTICS'//             
     +13X,'STATISTIC',19X,'BIASED',11X,'UNBIASED'/ 
     +54X,'(*LEAST BIASED)')              
240   FORMAT(/       
     +7X,'SUM',31X,'------',9X,E12.5/     
     +7X,'SQUARED SUM',23X,'------',9X,E12.5/      
     +7X,'CUBED SUM',25X,'------',9X,E12.5/        
     +7X,'MEAN',30X,'------',9X,E12.5/    
     +7X,'MEDIAN',28X,'------',9X,E12.5/  
     +7X,'VARIANCE',22X,E12.5,7X,E12.5/   
     +7X,'STANDARD DEVIATION',12X,E12.5,7X,E12.5,' *'/                  
     +7X,'COEFFICIENT OF VARIATION',6X,E12.5,7X,E12.5,' *'/             
     +7X,'SKEW',26X,E12.5,7X,E12.5/       
     +7X,'SKEWNESS COEFFICIENT',10X,E12.5,7X,E12.5/
     +7X,'COEFFICIENT OF EXCESS',9X,E12.5,7X,E12.5,'*'/                 
     +10X,'(KUTOSIS-3 ; PEAKEDNESS)')     
250   FORMAT(///27X,'LN(DATA) STATISTICS'//        
     +13X,'STATISTIC',19X,'BIASED',11X,'UNBIASED'/ 
     +54X,'(*LEAST BIASED)')              
299   FORMAT(///)    
C           
      END   
C1
C*TNMR******************************************************************
C 
      SUBROUTINE TNMR
      COMMON /ONE/ N,X(1500),P(1500),TT(6)
      COMMON /TWO/ TITLE                  
      COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,S4,SL4     
      COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)   
      COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)                  
      COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6),PX(1500)                
      CHARACTER TITLE*72,SPACE*1          
C           
C *** OPTION 2       
C *** ESTIMATES PARAMETERS AND TESTS G-O-F FOR NORMAL DISTRIBUTION      
C           
C *** WRITE TITLE AND HEADINGS            
      WRITE(2,200)TITLE                   
      WRITE(2,210)   
C           
C *** CALCULAT MEAN AND STANDARD DEVIATION OF DATA 
      XN=N*1.        
      PAR=2.
      XM=S1/XN       
      XV=(S2-S1*S1/XN)/(XN-1.)            
      XSD=SQRT(XV)   
C           
C *** MOMENT AND LIKELIHOOD ESTIMATION    
      A=XM  
      B=XSD 
      WRITE(2,220)A,B
      APAR(1,1)=A    
      BPAR(1,1)=B    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR DETERMINATION          
      S=0.0 
      DO 20 I=1,N    
          XI=I*1.    
          T=SNDV(P(I))                    
          XX(I)=A+T*B
          S=S+(X(I)-XX(I))**2             
20    CONTINUE       
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
1     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 2     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 1                    
C           
C *** CHI-SQUARE AND K-S GOODNESS OF FIT TESTS     
2     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 29 I=1,N    
          T=(X(I)-A)/B                    
          PX(I)=TINV(T)                   
29    CONTINUE       
      CL(1)=-9999999999.                  
      DO 30 I=2,NCL  
          CPR=CPR+PR 
          T=SNDV(CPR)
          CL(I)=A+T*B
30    CONTINUE       
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(1,1)=X2    
      CON(1,1)=1.-ALPHA                   
      STE(1,1)=SE    
      RKS(1,1)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
      DO 40 I=1,6    
          T=SNDV(1.-1./TT(I))             
          XT(I)=A+B*T
          ST(I)=SQRT( (B*B/XN)*(1.+T*T/2.) )       
40    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
C           
      RETURN
C           
C *** FORMAT STATEMENTS                   
C           
200   FORMAT(1H1/32X,'- NORMAL -'//       
     +1X,'TITLE: ',A72)                   
210   FORMAT(/20X,36('=')/20X,'I',34X,'I'/
     +20X,'I',31X,'2',2X,'I'/             
     +20X,'I',13X,'1',11X,'-(X-A)',3X,'I'/
     +20X,'I',2X,'P(X)= ------------ EXP -------',2X,'I'/               
     +20X,'I',8X,'B*SQRT(2*PI)',10X,'2',3X,'I'/    
     +20X,'I',27X,'2*B',4X,'I'/           
     +20X,'I',34X,'I'/20X,36('='))        
220   FORMAT(//30X,'A =',E12.5/30X,'B =',E12.5/)   
230   FORMAT(/1X,'THEORETICAL SORTED EVENTS'/      
     +6X,'NOTE: PROBABILITY IN PARENTHESIS = RANK/(N+1)'/)              
240   FORMAT(4(A1,'(',F5.4,')',E11.5,2X)) 
250   FORMAT(/6X,'STANDARD ERROR STATISTIC =',E12.5/                    
     +6X,'KOLMOGOROV-SMIRNOV STATISTIC =',F6.4/    
     +6X,'CHI-SQUARE G-O-F STATISTICS: '/ 
     +11X,'CHI-SQUARE=',F7.3/             
     +11X,'DEGREES OF FREEDOM=',F3.0/     
     +11X,'SIGNIFICANT AT 1-ALPHA=',F6.4) 
260   FORMAT(//1X,'SELECTED RETURN-PERIOD EVENTS'//
     +1X,'T,PERIOD',3X,'2',11X,'5',10X,'10',10X,'20',10X,'50',9X,'100'//
     +1X,'X',3X,6E12.5/2X,'T'/            
     +1X,'S',3X,6E12.5/2X,'T'//1X,'WHERE:'/        
     +7X,'X  = RETURN PERIOD FLOW FROM ESTIMATED DISTRIBUTION'/9X,'T'/  
     +7X,'S  = STANDARD ERROR ESTIMATE OF FLOW AT THAT RETURN PERIOD'/  
     +8X,'T')        
299   FORMAT(///)    
C           
      END   
C1
C*LN2R******************************************************************
C 
      SUBROUTINE LN2R
      COMMON /ONE/ N,X(1500),P(1500),TT(6)
      COMMON /TWO/ TITLE                  
      COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,S4,SL4     
      COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)   
      COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)                  
      COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6),PX(1500)                
      CHARACTER SPACE*1,TITLE*72          
C           
C *** OPTION 3       
C *** ESTIMATES PARAMETERS AND TESTS G-O-F FOR     
C *** TWO-PARAMETER LOGNORMAL DISTRIBUTION, BOTH   
C *** METHOD OF MOMENTS AND MAXIMUM LIKELIHOOD PROCEDURE.               
C           
C *** CALCULATE MEAN AND STANDARD DEVIATION OF DATA AND LN(DATA)        
      XN=N*1.        
      PAR=2.
      XM=S1/XN       
      YM=SL1/XN      
      XV=(S2-S1*S1/XN)/(XN-1.)            
      YV=(SL2-SL1*SL1/XN)/(XN-1.)         
      XSD=SQRT(XV)   
      YSD=SQRT(YV)   
C           
C *** WRITE TITLE AND HEADINGS FOR METHOD OF MOMENTS                    
      WRITE(2,200)TITLE                   
      WRITE(2,210)   
C           
C *** METHOD OF MOMENTS PARAMETER ESTIMATES        
      B2=ALOG((XSD*XSD)/(XM*XM)+1.)       
      B=SQRT(B2)     
      A=ALOG(XM)-B2/2.                    
      WRITE(2,220)A,B
      APAR(1,2)=A    
      BPAR(1,2)=B    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR MOMENT ESTIMATION  
      S=0   
      DO 20 I=1,N    
          T=SNDV(P(I))                    
          XX(I)=EXP(A+T*B)                
          S=S+(X(I)-XX(I))**2             
20    CONTINUE       
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
1     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 2     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 1                    
C           
C *** CHI-SQUARE AND K-S G-O-F TESTS FOR MOMENT ESTIMATION              
2     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 29 I=1,N    
          T=(ALOG(X(I))-A)/B              
          PX(I)=TINV(T)                   
29    CONTINUE       
      CL(1)=-9999999999.                  
      DO 30 I=2,NCL  
          CPR=CPR+PR 
          T=SNDV(CPR)
          CL(I)=EXP(A+T*B)                
30    CONTINUE       
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(1,2)=X2    
      CON(1,2)=1.-ALPHA                   
      STE(1,2)=SE    
      RKS(1,2)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR METHOD OF MOMENT ESTIMATION     
      Z=SQRT(EXP(YSD*YSD)-1.)             
      RK1=ALOG(1.+Z*Z)                    
      DO 40 I=1,6    
          T=SNDV(1.-1./TT(I))             
          XT(I)=EXP(A+T*B)                
          RK=(EXP(T*SQRT(RK1)-RK1/2.)-1.)/Z        
          DELTA=SQRT(1.+(Z**3+3.*Z)*RK+   
     +    (Z**8+6.*Z**6+15.*Z**4+16.*Z**2+2.)*RK*RK/4.)                 
          ST(I)=DELTA*XSD/(SQRT(XN))      
40    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
C           
C *** WRITE TITLE AND HEADINGS FOR MAXIMUM LIKELIHOOD PROCEDURE         
      WRITE(2,201)TITLE                   
      WRITE(2,210)   
C           
C *** MAXIMUM LIKELIHOOD PROCEDURE        
      A=YM  
      B=YSD 
      WRITE(2,220)A,B
      APAR(2,2)=A    
      BPAR(2,2)=B    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR LIKELIHOOD         
      S=0   
      DO 50 I=1,N    
          T=SNDV(P(I))                    
          XX(I)=EXP(A+T*B)                
          S=S+(X(I)-XX(I))**2             
50    CONTINUE       
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
3     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 4     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 3                    
C           
C *** CHI-SQUARE AND K-S G-O-F FOR LIKELIHOOD PROCEDURE                 
4     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 59 I=1,N    
          T=(ALOG(X(I))-A)/B              
          PX(I)=TINV(T)                   
59    CONTINUE       
      CL(1)=-9999999999.                  
      DO 60 I=2,NCL  
          CPR=CPR+PR 
          T=SNDV(CPR)
          CL(I)=EXP(A+T*B)                
60    CONTINUE       
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(2,2)=X2    
      CON(2,2)=1.-ALPHA                   
      STE(2,2)=SE    
      RKS(2,2)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR LIKELIHOOD PROCEDURE            
      Z=SQRT(EXP(YSD*YSD)-1.)             
      RK1=ALOG(1.+Z*Z)                    
      DO 70 I=1,6    
          T=SNDV(1.-1./TT(I))             
          XT(I)=EXP(A+T*B)                
          RK=(EXP(T*SQRT(RK1)-RK1/2.)-1.)/Z        
          DELTA=SQRT(RK1*(1.+RK*Z)**2*(1.+T*T/2.)/Z/Z)                  
          ST(I)=DELTA*XSD/(SQRT(XN))      
70    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
C           
      RETURN
C *** FORMAT STATEMENTS                   
C           
200   FORMAT(1H1/24X,'- TWO-PARAMETER LOGNORMAL -'//                    
     +29X,'METHOD OF MOMENTS'//           
     +1X,'TITLE: ',A72)                   
201   FORMAT(1H1/24X,'- TWO-PARAMETER LOGNORMAL -'//                    
     +24X,'MAXIMUM LIKELIHOOD PROCEDURE'//
     +1X,'TITLE: ',A72)                   
210   FORMAT(//      
     +16X,42('=')/16X,'I',40X,'I'/        
     +16X,'I',37X,'2',2X,'I'/             
     +16X,'I',14X,'1',12X,'-(LN(X)-A)',3X,'I'/     
     +16X,'I',2X,'P(X)= -------------- EXP -----------',2X,'I'/         
     +16X,'I',8X,'X*B*SQRT(2*PI)',11X,'2',6X,'I'/  
     +16X,'I',30X,'2*B',7X,'I'/           
     +16X,'I',40X,'I'/16X,42('='))        
220   FORMAT(//30X,'A =',E12.5/30X,'B =',E12.5/)   
230   FORMAT(/1X,'THEORETICAL SORTED EVENTS'/      
     +6X,'NOTE: PROBABILITY IN PARENTHESIS = RANK/(N+1)'/)              
240   FORMAT(4(A1,'(',F5.4,')',E11.5,2X)) 
250   FORMAT(/6X,'STANDARD ERROR STATISTIC =',E12.5/                    
     +6X,'KOLMOGOROV-SMIRNOV STATISTIC =',F6.4/    
     +6X,'CHI-SQUARE G-O-F STATISTICS: '/ 
     +11X,'CHI-SQUARE=',F7.3/             
     +11X,'DEGREES OF FREEDOM=',F3.0/     
     +11X,'SIGNIFICANT AT 1-ALPHA=',F6.4) 
260   FORMAT(//1X,'SELECTED RETURN-PERIOD EVENTS'//
     +1X,'T,PERIOD',3X,'2',11X,'5',10X,'10',10X,'20',10X,'50',9X,'100'//
     +1X,'X',3X,6E12.5/2X,'T'/            
     +1X,'S',3X,6E12.5/2X,'T'//1X,'WHERE:'/        
     +7X,'X  = RETURN PERIOD FLOW FROM ESTIMATED DISTRIBUTION'/9X,'T'/  
     +7X,'S  = STANDARD ERROR ESTIMATE OF FLOW AT THAT RETURN PERIOD'/  
     +8X,'T')        
299   FORMAT(///)    
C           
      END   
C1
C*LN3R******************************************************************
C 
      SUBROUTINE LN3R
      COMMON /ONE/ N,X(1500),P(1500),TT(6)
      COMMON /TWO/ TITLE                  
      COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,S4,SL4     
      COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)   
      COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)                  
      COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6),PX(1500)                
      CHARACTER SPACE*1,TITLE*72          
      REAL M2,M3,M4,M5,M6,MY,K            
C           
C *** OPTION 4       
C *** ESTIMATES PARAMETERS AND TESTS G-O-F FOR     
C *** THREE-PARAMETER LOGNORMAL DISTRIBUTION, BOTH 
C *** METHOD OF MOMENTS AND MAXIMUM LIKELIHOOD PROCEDURE.               
C           
C *** CALCULATE MEAN, STAND DEV AND COEFF OF SKEW FOR DATA AND LN(DATA) 
      XN=N*1.        
      PAR=3.
      XM=S1/XN       
      XV=(S2-S1*S1/XN)/(XN-1.)            
      XSD=SQRT(XV)   
      XSK1=(S3)-(3.*S1*S2/XN)+(2.*S1*S1*S1/(XN*XN))
      XSK=XSK1*XN/((XN-1.)*(XN-2.))       
      XG=XSK/(XSD**3)
C           
C *** WRITE TITLE AND HEADINGS FOR METHOD OF MOMENTS                    
      WRITE(2,200)TITLE                   
      WRITE(2,210)   
C           
C *** METHOD OF MOMENTS PARAMETER ESTIMATES        
      Z1=XSD/XM      
      W=(SQRT(XG*XG+4)-XG)/2.             
      Z2=(1.-W**(2./3.))/(W**(1./3.))     
      IF(XG.LE.0)THEN
          WRITE(2,221)XG                  
          STE(1,3)=-1
          RKS(1,3)=-1
          CHI(1,3)=-1
          CON(1,3)=-1
          GOTO 7     
      END IF
      A=ALOG(XSD/Z2)-0.5*ALOG(Z2*Z2+1.)   
      B=SQRT(ALOG(Z2*Z2+1.))              
      C=XM-XSD/Z2    
      WRITE(2,220)A,B,C                   
      APAR(1,3)=A    
      BPAR(1,3)=B    
      CPAR(1,3)=C    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR MOMENT ESTIMATION  
      S=0   
      DO 20 I=1,N    
          T=SNDV(P(I))                    
          XX(I)=C+EXP(A+T*B)              
          S=S+(X(I)-XX(I))**2             
20    CONTINUE       
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
1     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 2     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 1                    
C           
C *** CHI-SQUARE AND K-S G-O-F TESTS FOR MOMENT ESTIMATION              
2     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 29 I=1,N    
          T=(ALOG(X(I)-C)-A)/B            
          PX(I)=TINV(T)                   
29    CONTINUE       
      CL(1)=-9999999999.                  
      DO 30 I=2,NCL  
          CPR=CPR+PR 
          T=SNDV(CPR)
          CL(I)=EXP(A+T*B)                
30    CONTINUE       
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(1,3)=X2    
      CON(1,3)=1.-ALPHA                   
      STE(1,3)=SE    
      RKS(1,3)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR METHOD OF MOMENT ESTIMATION     
      SY=(ALOG(Z2**2+1.0))**0.5           
      SY2=SY**2      
      MY=ALOG(XSD/Z2)-0.5*ALOG(Z2**2+1.0) 
      E=EXP(SY2)     
      EA=EXP(2.0*SY2)
      EB=EXP(2.5*SY2)
      EC=EXP(3.0*SY2)
      ED=EXP(4.0*SY2)
      EF=EXP(6.0*SY2)
      EG=EXP(10.0*SY2)                    
      EH=EXP(15.0*SY2)                    
      EI=EXP(4.0*MY) 
      EJ=EXP(5.0*MY) 
      EK=EXP(6.0*MY) 
      EL=(EXP(SY2)-1.0)**2                
      M2=XSD*XSD     
      M3=XG*(XSD**3) 
      M4=EA*EI*EL*(ED+2.0*EC+3.0*EA-3.0)  
      M5=EB*EJ*(EG-5.0*EF+10.0*EC-10.0*E+4.0)      
      M6=EC*EK*(EH-6.0*EG+15.0*EF-20.0*EC+15.0*E-5.0)                   
      VM1=M2/XN      
      VM2=(M4-M2**2)/XN                   
      VM3=(M6-M3**2-6.0*M4*M2+9.0*M2**3)/XN        
      CM1M2=M3/XN    
      CM1M3=(M4-3.0*M2**2)/XN             
      CM2M3=(M5-4.0*M3*M2)/XN             
      DO 40 I=1,6    
          T=SNDV(1.-1./TT(I))             
          DXDM1=1.0  
          DWDG=-0.5+XG/(2.0*(XG**2+4.0)**0.5)      
          DZ2DW=(-1./3.)*(W**(-4./3.)+W**(-2./3.)) 
          D1=ALOG(Z2**2+1.0)              
          D2=EXP((SQRT(D1))*T-D1/2)       
          D3=(2.0*Z2)/(1.0+Z2**2)         
          D4=T/(2.0*Z2*SQRT(D1))          
          D5=1.0/Z2  
          D6=1.0/(2.0*Z2**3)              
          DKDZ2=D3*(D2*(D4-D5-D6)+D6+D5/2.)        
          K=(D2-1.0)/Z2                   
          DKDG=DKDZ2*DZ2DW*DWDG           
          DXDM2=(1.0/(2.0*SQRT(M2)))*(K-3.0*XG*DKDG)                    
          DXDM3=DKDG/M2                   
          ST(I)=SQRT((DXDM1**2)*VM1+(DXDM2**2)*VM2+(DXDM3**2)*VM3+2.*   
     &    DXDM1*DXDM2*CM1M2+2.*DXDM1*DXDM3*CM1M3+2.*DXDM2*DXDM3*CM2M3)  
          XT(I)=C+EXP(A+T*B)              
40     CONTINUE      
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
7     WRITE(2,299)   
C           
C *** WRITE TITLE AND HEADINGS FOR MAXIMUM LIKELIHOOD PROCEDURE         
      WRITE(2,201)TITLE                   
      WRITE(2,210)   
C           
C *** MAXIMUM LIKELIHOOD PROCEDURE        
      XMIN=10000000 
      DO 21 I=1,N 
        IF(X(I) .LT. XMIN)XMIN=X(I) 
21    CONTINUE
      AML=XMIN*.8 
      ICOUNT=0       
6     ICOUNT=ICOUNT+1
      AA=0.0
      BB=0.0
      CC=0.0
      D=0.0 
      E=0.0 
      F=0.0 
      DO 50 I=1,N    
          AA=AA+ALOG(X(I)-AML)            
          BB=BB+(ALOG(X(I)-AML))**2       
          CC=CC+1.0/((X(I)-AML))          
          D=D+1.0/((X(I)-AML)**2)         
          E=E+(1.0/((X(I)-AML)))*ALOG(X(I)-AML)    
          F=F+(1.0/((X(I)-AML)**2))*ALOG(X(I)-AML) 
50    CONTINUE       
      G=(BB/XN)-(AA/XN)**2-(AA/XN)        
      H=(-2.*E/XN)+(2.*AA/XN)*(CC/XN)+(CC/XN)      
      FCN=CC*G+E     
      FPN=CC*H+D*G+F-D                    
      AS=AML-(FCN/FPN)                    
      DELTA=ABS(0.00001*AS)               
      IF (ABS(AS-AML).LT.DELTA) GOTO 8    
      IF (ICOUNT.GT.25) GOTO 9            
      AML=AS
      GOTO 6
8     CONTINUE       
      A=AA/XN        
      B=SQRT((BB-AA*AA/XN)/(XN-1))        
      C=AS  
      WRITE(2,220)A,B,C                   
      APAR(2,3)=A    
      BPAR(2,3)=B    
      CPAR(2,3)=C    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR LIKELIHOOD         
      S=0   
      DO 60 I=1,N    
          T=SNDV(P(I))                    
          XX(I)=C+EXP(A+T*B)              
          S=S+(X(I)-XX(I))**2             
60    CONTINUE       
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
3     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 4     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 3                    
C           
C *** CHI-SQUARE AND K-S G-O-F FOR LIKELIHOOD PROCEDURE                 
4     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 69 I=1,N    
          T=(ALOG(X(I)-C)-A)/B            
          PX(I)=TINV(T)                   
69    CONTINUE       
      CL(1)=-9999999999.                  
      DO 70 I=2,NCL  
          CPR=CPR+PR 
          T=SNDV(CPR)
          CL(I)=C+EXP(A+T*B)              
70    CONTINUE       
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(2,3)=X2    
      CON(2,3)=1.-ALPHA                   
      STE(2,3)=SE    
      RKS(2,3)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR LIKELIHOOD PROCEDURE            
      VAR=B*B        
      AA=EXP(VAR-2.0*A)                   
      BB=EXP(2.0*VAR-2.0*A)               
      CC=EXP(VAR/2.0-A)                   
      D1=(VAR+1.0)/(2.0*VAR)              
      D2=1.0/(2.0*VAR)                    
      D=D1*BB-D2*AA-AA                    
      E=1.0/(XN*D)   
      VA=E*0.5       
      VMU=(VAR*E)*(D1*BB-AA)              
      VVAR=VAR*E*((VAR+1.0)*BB-AA)        
      CAMU=CC*E/2.0  
      CAVAR=VAR*E*CC 
      CMUVAR=VAR*E*AA
      CAMU=-CAMU     
      CMUVAR=-CMUVAR 
      DO 80 I=1,6    
          T=SNDV(1.-1./TT(I))             
          Z=EXP(A+T*B)                    
          VX=VA+(VMU*Z**2)+(T*Z*CAVAR/B  )+(2.0*Z*CAMU)+                
     &    (T*Z**2*CMUVAR/B  )+(T**2*Z**2*VVAR/(4.*VAR))                 
          XT(I)=C+Z  
          ST(I)=SQRT(VX)                  
80     CONTINUE      
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
C           
      RETURN
C           
9     WRITE(2,222)   
      STE(2,3)=-1    
      RKS(2,3)=-1    
      CHI(2,3)=-1    
      CON(2,3)=-1    
      WRITE(2,299)   
      RETURN
C           
C *** FORMAT STATEMENTS                   
C           
200   FORMAT(1H1/23X,'- THREE-PARAMETER LOGNORMAL -'//                  
     +29X,'METHOD OF MOMENTS'//           
     +1X,'TITLE: ',A72)                   
201   FORMAT(1H1/23X,'- THREE-PARAMETER LOGNORMAL -'//                  
     +24X,'MAXIMUM LIKELIHOOD PROCEDURE'//
     +1X,'TITLE: ',A72)                   
210   FORMAT(//      
     +13X,48('=')/13X,'I',46X,'I'/        
     +13X,'I',43X,'2',2X,'I'/             
     +13X,'I',16X,'1',14X,'-(LN(X-C)-A)',3X,'I'/   
     +13X,'I',2X,'P(X)= ------------------ EXP -------------',2X,'I'/   
     +13X,'I',8X,'(X-C)*B*SQRT(2*PI)',12X,'2',7X,'I'/                   
     +13X,'I',35X,'2*B',8X,'I'/           
     +13X,'I',46X,'I'/13X,48('='))        
220   FORMAT(//30X,'A =',E12.5/30X,'B =',E12.5/30X,'C=',E12.5/)         
221   FORMAT(///11X,'*** UNABLE TO PERFORM THE MOMENT ESTIMATION ***'/  
     +1X,10X,'*** DUE TO THE NEGATIVE SKEW (=',E12.5,')')               
222   FORMAT(///1X,'MAXIMUM LIKELIHOOD PROCEDURE DID NOT CONVERGE')     
230   FORMAT(/1X,'THEORETICAL SORTED EVENTS'/      
     +6X,'NOTE: PROBABILITY IN PARENTHESIS = RANK/(N+1)'/)              
240   FORMAT(4(A1,'(',F5.4,')',E11.5,2X)) 
250   FORMAT(/6X,'STANDARD ERROR STATISTIC =',E12.5/                    
     +6X,'KOLMOGOROV-SMIRNOV STATISTIC =',F6.4/    
     +6X,'CHI-SQUARE G-O-F STATISTICS: '/ 
     +11X,'CHI-SQUARE=',F7.3/             
     +11X,'DEGREES OF FREEDOM=',F3.0/     
     +11X,'SIGNIFICANT AT 1-ALPHA=',F6.4) 
260   FORMAT(//1X,'SELECTED RETURN-PERIOD EVENTS'//
     +1X,'T,PERIOD',3X,'2',11X,'5',10X,'10',10X,'20',10X,'50',9X,'100'//
     +1X,'X',3X,6E12.5/2X,'T'/            
     +1X,'S',3X,6E12.5/2X,'T'//1X,'WHERE:'/        
     +7X,'X  = RETURN PERIOD FLOW FROM ESTIMATED DISTRIBUTION'/9X,'T'/  
     +7X,'S  = STANDARD ERROR ESTIMATE OF FLOW AT THAT RETURN PERIOD'/  
     +8X,'T')        
299   FORMAT(///)    
C           
      END   
C1
C*T1ER******************************************************************
C 
      SUBROUTINE T1ER
      COMMON /ONE/ N,X(1500),P(1500),TT(6)
      COMMON /TWO/ TITLE                  
      COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,S4,SL4     
      COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)   
      COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)                  
      COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6),PX(1500)                
      CHARACTER SPACE*1,TITLE*72          
C           
C *** OPTION 5       
C *** ESTIMATES PARAMETERS AND TESTS G-O-F FOR     
C *** EXTREMAL TYPE ONE (GUMBEL) DISTRIBUTION, BOTH
C *** METHOD OF MOMENTS AND MAXIMUM LIKELIHOOD PROCEDURE.               
C           
C *** CALCULATE MEAN AND STANDARD DEVIATION OF DATA
      XN=N*1.        
      PAR=2.
      XM=S1/XN       
      XV=(S2-S1*S1/XN)/(XN-1.)            
      XSD=SQRT(XV)   
C           
C *** WRITE TITLE AND HEADINGS FOR METHOD OF MOMENTS                    
      WRITE(2,200)TITLE                   
      WRITE(2,210)   
C           
C *** METHOD OF MOMENTS PARAMETER ESTIMATES        
      A=1.282550/XSD 
      B=XM-0.4500532*XSD                  
      WRITE(2,220)A,B
      APAR(1,4)=A    
      BPAR(1,4)=B    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR MOMENT ESTIMATION  
      S=0.0 
      DO 20 I=1,N    
          T=1./(1.-P(I))                  
          XX(I)=B-(1./A)*ALOG(-1.*ALOG(1.-1./T))   
          S=S+(X(I)-XX(I))**2             
20    CONTINUE       
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
1     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 2     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 1                    
C           
C *** CHI-SQUARE AND K-S G-O-F TESTS FOR MOMENT ESTIMATION              
2     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 29 I=1,N    
29    PX(I)=EXP(-1*EXP(-1*A*(X(I)-B)))    
      CL(1)=-9999999999.                  
      DO 30 I=2,NCL  
          CPR=CPR+PR 
          T=1./(1.-CPR)                   
          CL(I)=B-(1./A)*ALOG(-1.*ALOG(1.-1./T))   
30    CONTINUE       
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(1,4)=X2    
      CON(1,4)=1.-ALPHA                   
      STE(1,4)=SE    
      RKS(1,4)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR METHOD OF MOMENT ESTIMATION     
      DO 40 I=1,6    
          T=1./(1./TT(I))                 
          XT(I)=B-(1./A)*ALOG(-1.*ALOG(1.-1./T))   
          RK=-1.*(0.45+0.7797*ALOG(-1.*ALOG(1.-1./T)))                  
          ST(I)=XSD*SQRT((1.+1.139547093*RK+1.100000027*RK*RK)/XN)      
40    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
C           
C *** WRITE TITLE AND HEADINGS FOR MAXIMUM LIKELIHOOD PROCEDURE         
      WRITE(2,201)TITLE                   
      WRITE(2,210)   
C           
C *** MAXIMUM LIKELIHOOD PROCEDURE        
      AML=A 
5     ICOUNT=ICOUNT+1
      AA=1./(AML*AML)
      BB=XM-1./AML   
      C=0.0 
      D=0.0 
      E=0.0 
      DO 45 I=1,N    
          TEMP=EXP(-AML*X(I))             
          C=C+TEMP   
          D=D+TEMP*X(I)                   
          E=E+TEMP*X(I)*X(I)              
45    CONTINUE       
      FCN=D-BB*C     
      FPN=BB*D-E-AA*C
      AS=AML-(FCN/FPN)                    
      DELTA=ABS(0.0000001*AS)             
      IF(ABS(AS-AML).LT.DELTA)GOTO 6      
      IF(ICOUNT.GT.25)GOTO 9              
      AML=AS
      GOTO 5
6     CONTINUE       
      A=AS  
      B=(1./AS)*ALOG(XN/C)                
      WRITE(2,220)A,B
      APAR(2,4)=A    
      BPAR(2,4)=B    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR LIKELIHOOD         
8     S=0.0 
      DO 50 I=1,N    
          T=1./(1/1.-P(I))                
          XX(I)=B-(1./A)*ALOG(-1.*ALOG(1.-1./T))   
          S=S+(X(I)-XX(I))**2             
50    CONTINUE       
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
3     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 4     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 3                    
C           
C *** CHI-SQUARE AND K-S G-O-F FOR LIKELIHOOD PROCEDURE                 
4     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 59 I=1,N    
59    PX(I)=EXP(-1.*EXP(-1.*A*(X(I)-B)))  
      CL(1)=-9999999999.                  
      DO 60 I=2,NCL  
          CPR=CPR+PR 
          T=1./(1.-CPR)                   
          CL(I)=B-(1./A)*ALOG(-1.*ALOG(1.-1./T))   
60    CONTINUE       
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(2,4)=X2    
      CON(2,4)=1.-ALPHA                   
      STE(2,4)=SE    
      RKS(2,4)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR LIKELIHOOD PROCEDURE            
      DO 70 I=1,6    
          T=1./(1./TT(I))                 
          XT(I)=B-(1./A)*ALOG(-1.*ALOG(1.-1./T))   
          GE=0.5772157                    
          PSD6=1.644934                   
          YT=-1.*ALOG(-1.*ALOG((T-1.)/T)) 
          DELTA=SQRT(1.+(1.-GE+YT)**2/PSD6)        
          ST(I)=DELTA/A/SQRT(XN)          
70    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
      RETURN
C           
9     WRITE(2,221)   
      STE(2,4)=-1    
      RKS(2,4)=-1    
      CHI(2,4)=-1    
      CON(2,4)=-1    
      WRITE(2,299)   
      RETURN
C           
C *** FORMAT STATEMENTS                   
C           
200   FORMAT(1H1/29X,'- TYPE I EXTREMAL -'//       
     +29X,'METHOD OF MOMENTS'//           
     +1X,'TITLE: ',A72)                   
201   FORMAT(1H1/29X,'- TYPE I EXTREMAL -'//       
     +24X,'MAXIMUM LIKELIHOOD PROCEDURE'//
     +1X,'TITLE: ',A72)                   
210   FORMAT(//      
     +14X,45('=')/14X,'I',43X,'I'/        
     +14X,'I',2X,'P(X)= A*EXP( -A*(X-B) - EXP(-A*(X-B)) )',2X,'I'/      
     +14X,'I',43X,'I'/14X,45('='))        
220   FORMAT(//30X,'A =',E12.5/30X,'B =',E12.5/)   
221   FORMAT(///1X,'MAXIMUM LIKELIHOOD PROCEDURE DID NOT CONVERGE')     
230   FORMAT(/1X,'THEORETICAL SORTED EVENTS'/      
     +6X,'NOTE: PROBABILITY IN PARENTHESIS = RANK/(N+1)'/)              
240   FORMAT(4(A1,'(',F5.4,')',E11.5,2X)) 
250   FORMAT(/6X,'STANDARD ERROR STATISTIC =',E12.5/                    
     +6X,'KOLMOGOROV-SMIRNOV STATISTIC =',F6.4/    
     +6X,'CHI-SQUARE G-O-F STATISTICS: '/ 
     +11X,'CHI-SQUARE=',F7.3/             
     +11X,'DEGREES OF FREEDOM=',F3.0/     
     +11X,'SIGNIFICANT AT 1-ALPHA=',F6.4) 
260   FORMAT(//1X,'SELECTED RETURN-PERIOD EVENTS'//
     +1X,'T,PERIOD',3X,'2',11X,'5',10X,'10',10X,'20',10X,'50',9X,'100'//
     +1X,'X',3X,6E12.5/2X,'T'/            
     +1X,'S',3X,6E12.5/2X,'T'//1X,'WHERE:'/        
     +7X,'X  = RETURN PERIOD FLOW FROM ESTIMATED DISTRIBUTION'/9X,'T'/  
     +7X,'S  = STANDARD ERROR ESTIMATE OF FLOW AT THAT RETURN PERIOD'/  
     +8X,'T')        
299   FORMAT(///)    
C           
      END   
C1
C*PT3R******************************************************************
C 
      SUBROUTINE PT3R
      COMMON /ONE/ N,X(1500),P(1500),TT(6)
      COMMON /TWO/ TITLE                  
      COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,S4,SL4     
      COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)   
      COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)                  
      COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6),PX(1500)                
      CHARACTER SPACE*1,TITLE*72          
C           
C *** OPTION 6       
C *** ESTIMATES PARAMETERS AND TESTS G-O-F FOR     
C *** PEARSON TYPE III, BOTH              
C *** METHOD OF MOMENTS AND MAXIMUM LIKELIHOOD PROCEDURE.               
C           
C *** CALCULATE MEAN, STAN DEV AND COEFF OF SKEW OF DATA                
      XN=N*1.        
      PAR=3.
      XM=S1/XN       
      XV=(S2-S1*S1/XN)/(XN-1.)            
      XSD=SQRT(XV)   
      XSK1=(S3)-(3.*S1*S2/XN)+(2.*S1*S1*S1/(XN*XN))
      XSK=XSK1*XN/((XN-1.)*(XN-2.))       
      XG=XSK/(XSD**3)
C     XG=XG1*SQRT(XN*(XN-1.))*(1.+8.5/XN)/(XN-2.)  
C           
C *** WRITE TITLE AND HEADINGS FOR METHOD OF MOMENTS                    
      WRITE(2,200)TITLE                   
      WRITE(2,210)   
C           
C *** METHOD OF MOMENTS PARAMETER ESTIMATES        
      B=(2./XG)**2   
      A=XSD/SQRT(B)  
      C=XM-XSD*SQRT(B)                    
      WRITE(2,220)A,B,C                   
      APAR(1,5)=A    
      BPAR(1,5)=B    
      CPAR(1,5)=C    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR MOMENT ESTIMATION  
      S=0   
      V=2.*B
      IF(V.GE.(0.5).AND.V.LE.(200000.))THEN        
          DO 20 I=1,N
              CALL MDCHI (P(I),V,X2,IER)  
              XX(I)=X2*A/2.+C             
              S=S+(X(I)-XX(I))**2         
20        CONTINUE   
      ELSE  
          DO 21 I=1,N
              T=SNDV(P(I))                
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              XX(I)=C+A*E**3              
              S=S+(X(I)-XX(I))**2         
21        CONTINUE   
      END IF
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
1     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 2     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 1                    
C           
C *** CHI-SQUARE AND K-S G-O-F TESTS FOR MOMENT ESTIMATION              
2     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 29 I=1,N    
      IF(X(I).LT.C)THEN                   
          E=((C-X(I))/A)**(1./3.)         
          E=-E       
      ELSE IF(X(I).EQ.C)THEN              
          E=0.0      
      ELSE  
          E=((X(I)-C)/A)**(1./3.)         
      END IF
          T=(E-B**(1./3.)+1./(9.*B**(2./3.)))*3.*B**(1./6.)             
          PX(I)=TINV(T)                   
29    CONTINUE       
      CL(1)=-9999999999.                  
      V=2.*B
      IF(V.GE.(0.5).AND.V.LE.(200000.))THEN        
          DO 30 I=2,NCL                   
             CPR=CPR+PR                   
             CALL MDCHI(CPR,V,X2,IER)     
             CL(I)=X2*A/2.+C              
30       CONTINUE    
      ELSE  
          DO 31 I=2,NCL                   
              CPR=CPR+PR                  
              T=SNDV(CPR)                 
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              CL(I)=C+A*E**3              
31        CONTINUE   
      END IF
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(1,5)=X2    
      CON(1,5)=1.-ALPHA                   
      STE(1,5)=SE    
      RKS(1,5)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR METHOD OF MOMENT ESTIMATION     
      V=2.*B
      DO 40 I=1,6    
          PR=1.-1./TT(I)                  
          T=SNDV(PR) 
          T1=T       
          T2=(T*T-1.)/6.                  
          T3=2.*(T*T*T-6.*T)/(6*6*6)      
          T4=(T*T-1.)/(6*6*6)             
          T5=T/(6.**4)                    
          T6=2./(6.**6)                   
          RK=T1+T2*XG+T3*XG*XG-T4*XG**3+T5*XG**4-T6*XG**5               
          SLOPE=T2+T3*2.*XG-T4*3.*XG*XG+T5*4.*XG**3-T6*5.*XG**4         
          T7=(1.+0.75*XG*XG)*(0.5*RK*RK)  
          T8=RK*XG   
          T9=6.*(1.+.25*XG*XG)*SLOPE      
          T10=SLOPE*(1.+1.25*XG*XG)+(XG*RK/2.)     
          DELTA=T7+T8+T9*T10              
          IF(V.GE.(0.5).AND.V.LE.(200000.))THEN    
              CALL MDCHI(PR,V,X2,IER)     
              XT(I)=X2*A/2.+C             
          ELSE       
              XT(I)=XM+RK*XSD             
          END IF     
          ST(I)=XSD*SQRT(DELTA/XN)        
40    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
C           
C *** WRITE TITLE AND HEADINGS FOR MAXIMUM LIKELIHOOD PROCEDURE         
      WRITE(2,201)TITLE                   
      WRITE(2,210)   
C           
C *** MAXIMUM LIKELIHOOD PROCEDURE        
      IT=1  
      ICOUNT=0       
      GML=X(N)*0.99  
5     ICOUNT=ICOUNT+1
      AA=0.0
      BB=0.0
      CC=0.0
      R=0.0 
      DO 45 I=1,N    
      AA=AA+1.0/(X(I)-GML)                
      BB=BB+(X(I)-GML)                    
      IF(X(I).LE.GML)GOTO 8               
      CC=CC+ALOG(X(I)-GML)                
      R=R+1.0/((X(I)-GML)**2)             
45    CONTINUE       
      BETA=AA/(AA-(XN**2)/BB)             
      ALPHA=BB/(XN*BETA)                  
      D=BETA+2.      
      PSI=ALOG(D)-(1./(2.*D))-(1./(12.*D**2))+(1./(120.*D**4))-         
     &    (1./(252.*D**6))-(1./(BETA+1.))-(1./BETA)
      FCN=-XN*PSI+CC-XN*ALOG(ALPHA)       
      TRI=(1./D)+(1./(2.*D**2))+(1./(6.*D**3))-(1./(30.*D**5))+         
     &    (1./(42.*D**7))-(1./(30.*D**9))+(1./((BETA+1.)**2))+          
     &    (1./(BETA**2))                  
      VV=AA-(XN**2)/BB                    
      U=AA  
      W=(BB/XN)-(XN/AA)                   
      DU=R  
      DV=R-(XN**3)/(BB**2)                
      DW=-1.+(XN*R)/(AA**2)               
      FPN=-XN*TRI*((VV*DU-U*DV)/(VV**2))-AA-XN*DW/W
      AS=GML-(FCN/FPN)                    
      DELTA=ABS(0.00000001*AS)            
      IF (ABS(AS-GML).LT.DELTA) GOTO 6    
      IF (ICOUNT.GT.25) GOTO 8            
      GML=AS
      GOTO 5
8     CONTINUE       
      IF(IT.EQ.3)GOTO 9                   
      IF(IT.EQ.2)THEN
          IT=3       
          ICOUNT=0   
          GML=X(N)*0.5                    
          GOTO 5     
      END IF
      IF(IT.EQ.1)THEN
          IT=2       
          ICOUNT=0   
          GML=C      
          GOTO 5     
      END IF
6     CONTINUE       
      A=ALPHA        
      B=BETA
      C=AS  
      WRITE(2,220)A,B,C                   
      APAR(2,5)=A    
      BPAR(2,5)=B    
      CPAR(2,5)=C    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR LIKELIHOOD         
      S=0.0 
      V=2.*B
      IF(V.GE.(0.5).AND.V.LE.(200000.))THEN        
          DO 50 I=1,N
              CALL MDCHI (P(I),V,X2,IER)  
              XX(I)=X2*A/2.+C             
              S=S+(X(I)-XX(I))**2         
50        CONTINUE   
      ELSE  
          DO 51 I=1,N
              T=SNDV(P(I))                
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              XX(I)=C+A*E**3              
              S=S+(X(I)-XX(I))**2         
51        CONTINUE   
      END IF
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
3     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 4     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 3                    
C           
C *** CHI-SQUARE AND K-S G-O-F FOR LIKELIHOOD PROCEDURE                 
4     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 59 I=1,N    
      IF(X(I).LT.C)THEN                   
          E=((C-X(I))/A)**(1./3.)         
          E=-E       
      ELSE IF(X(I).EQ.C)THEN              
          E=0.0      
      ELSE  
          E=((X(I)-C)/A)**(1./3.)         
      END IF
          T=(E-B**(1./3.)+1./(9.*B**(2./3.)))*3.*B**(1./6.)             
          PX(I)=TINV(T)                   
59    CONTINUE       
      CL(1)=-9999999999.                  
      V=2.*B
      IF(V.GE.(0.5).AND.V.LE.(200000.))THEN        
          DO 60 I=2,NCL                   
             CPR=CPR+PR                   
             CALL MDCHI(CPR,V,X2,IER)     
             CL(I)=X2*A/2.+C              
60       CONTINUE    
      ELSE  
          DO 61 I=2,NCL                   
              CPR=CPR+PR                  
              T=SNDV(CPR)                 
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              CL(I)=C+A*E**3              
61        CONTINUE   
      END IF
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(2,5)=X2    
      CON(2,5)=1.-ALPHA                   
      STE(2,5)=SE    
      RKS(2,5)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR LIKELIHOOD PROCEDURE            
      ALPHA=A        
      BETA=B
      GAMMA=C        
      D=BETA+2.0     
      TRI=(1./D)+(1./(2.*D**2))+(1./(6.*D**3))-(1./(30.*D**5))+         
     &    (1./(42.*D**7))-(1./(30.*D**9))+(1./((BETA+1.)**2))+          
     &    (1./(BETA**2))                  
      H=(BETA-2.)*ALPHA**4                
      PP=2.*TRI-(2.*BETA-3.)/((BETA-1.)**2)        
      DET=PP/H       
      VARA=(1./(XN*(ALPHA**2)*DET))*((TRI/(BETA-2.))-1./((BETA-1.)**2)) 
      VARB=2./(XN*DET*(BETA-2.)*ALPHA**4) 
      VARG=(BETA*TRI-1.)/(XN*DET*ALPHA**2)
      COVAB=(-1./(XN*DET*ALPHA**3))*((1./(BETA-2.))-(1./(BETA-1.)))     
      COVAG=(1./(XN*DET*ALPHA**2))*((1./(BETA-1.))-TRI)                 
      COVBG=(-1./(XN*DET*ALPHA**3))*((BETA/(BETA-1.))-1.)               
      V=2.*B
      DO 70 I=1,6    
          PR=1.-1./TT(I)                  
          T=SNDV(PR) 
          E=BETA**(1./3.)-1./(9.*BETA**(2./3.))+T/(3.*BETA**(1./6.))    
          F=1./(3.*BETA**(2./3.))+2./(27.*BETA**(5./3.))-               
     &      T/(18.*BETA**(7./6.))         
          DXDA=E**3  
          DXDB=3.*ALPHA*E**2*F            
          DXDG=1.0   
          IF(V.GE.(0.5).AND.V.LE.(200000.))THEN    
              CALL MDCHI(PR,V,X2,IER)     
              XT(I)=X2*A/2.+C             
          ELSE       
              XT(I)=C+A*E**3              
          END IF     
          ST(I)=SQRT(VARA*DXDA**2+VARB*DXDB**2+VARG*DXDG**2+2.*DXDA*    
     &    DXDB*COVAB+2.*DXDA*DXDG*COVAG+2.*DXDB*DXDG*COVBG)             
70    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
C           
      RETURN
C           
9     WRITE(2,221)   
      STE(2,5)=-1    
      RKS(2,5)=-1    
      CHI(2,5)=-1    
      CON(2,5)=-1    
      WRITE(2,299)   
      RETURN
C           
C *** FORMAT STATEMENTS                   
C           
200   FORMAT(1H1/28X,'- PEARSON TYPE III -'//      
     +29X,'METHOD OF MOMENTS'//           
     +1X,'TITLE: ',A72)                   
201   FORMAT(1H1/28X,'- PEARSON TYPE III -'//      
     +24X,'MAXIMUM LIKELIHOOD PROCEDURE'//
     +1X,'TITLE: ',A72)                   
210   FORMAT(//      
     +16X,43('=')/16X,'I',41X,'I'/        
     +16X,'I',22X,'(B-1)',14X,'I'/        
     +16X,'I',11X,'1',6X,'(X-C)',9X,'-(X-C)',3X,'I'/                    
     +16X,'I',2X,'P(X)=',10('-'),1X,9('-'),' EXP ',7('-'),2X,'I'/       
     +16X,'I',7X,'A*GAMMA(B)',3X,'(B-1)',10X,'A',5X,'I'/                
     +16X,'I',19X,'A',21X,'I'/            
     +16X,'I',41X,'I'/16X,43('='))        
220   FORMAT(//30X,'A =',E12.5/30X,'B =',E12.5/30X,'C =',E12.5/)        
221   FORMAT(///1X,'MAXIMUM LIKELIHOOD DOES NOT CONVERGE WITHIN THE',   
     +' ALLOWED STEPS (25)')              
230   FORMAT(/1X,'THEORETICAL SORTED EVENTS'/      
     +6X,'NOTE: PROBABILITY IN PARENTHESIS = RANK/(N+1)'/)              
240   FORMAT(4(A1,'(',F5.4,')',E11.5,2X)) 
250   FORMAT(/6X,'STANDARD ERROR STATISTIC =',E12.5/                    
     +6X,'KOLMOGOROV-SMIRNOV STATISTIC =',F6.4/    
     +6X,'CHI-SQUARE G-O-F STATISTICS: '/ 
     +11X,'CHI-SQUARE=',F7.3/             
     +11X,'DEGREES OF FREEDOM=',F3.0/     
     +11X,'SIGNIFICANT AT 1-ALPHA=',F6.4) 
260   FORMAT(//1X,'SELECTED RETURN-PERIOD EVENTS'//
     +1X,'T,PERIOD',3X,'2',11X,'5',10X,'10',10X,'20',10X,'50',9X,'100'//
     +1X,'X',3X,6E12.5/2X,'T'/            
     +1X,'S',3X,6E12.5/2X,'T'//1X,'WHERE:'/        
     +7X,'X  = RETURN PERIOD FLOW FROM ESTIMATED DISTRIBUTION'/9X,'T'/  
     +7X,'S  = STANDARD ERROR ESTIMATE OF FLOW AT THAT RETURN PERIOD'/  
     +8X,'T')        
299   FORMAT(///)    
C           
      END   
C1
C*LP3R******************************************************************
C 
      SUBROUTINE LP3R
      COMMON /ONE/ N,X(1500),P(1500),TT(6)
      COMMON /TWO/ TITLE                  
      COMMON /THR/ S1,S2,S3,SL1,SL2,SL3,S4,SL4     
      COMMON /FOR/ APAR(3,6),BPAR(3,6),CPAR(3,6)   
      COMMON /FIV/ CHI(3,6),CON(3,6),STE(3,6),RKS(3,6)                  
      COMMON /WKA/ XX(1500),CL(302),XT(6),ST(6),PX(1500)                
      CHARACTER SPACE*1,TITLE*72          
      REAL Y(1500),L1,L2,L3,NSX,M1,M2,M3  
C           
C *** OPTION 7       
C *** ESTIMATES PARAMETERS AND TESTS G-O-F FOR     
C *** LOG-PEARSON TYPE III, BOTH          
C *** METHOD OF MOMENTS AND MAXIMUM LIKELIHOOD PROCEDURE.               
C           
C *** CALCULATE MEAN, STAN DEV AND COEFF OF SKEW OF DATA AND LN(DATA)   
      XN=N*1.        
      PAR=3.
      XM=S1/XN       
      YM=SL1/XN      
      XV=(S2-S1*S1/XN)/(XN-1.)            
      YV=(SL2-SL1*SL1/XN)/(XN-1.)         
      XSD=SQRT(XV)   
      YSD=SQRT(YV)   
      XSK1=(S3)-(3.*S1*S2/XN)+(2.*S1*S1*S1/(XN*XN))
      YSK1=(SL3)-(3.*SL1*SL2/XN)+(2.*SL1*SL1*SL1/(XN*XN))               
      XSK=XSK1*XN/((XN-1.)*(XN-2.))       
      YSK=YSK1*XN/((XN-1.)*(XN-2.))       
      XG=XSK/(XSD**3)
      YG=YSK/(YSD**3)
C     XG=XG1*SQRT(XN*(XN-1.))*(1.+8.5/XN)/(XN-2.)  
C     YG=YG1*SQRT(XN*(XN-1.))*(1.+8.5/XN)/(XN-2.)  
C           
C *** WRITE TITLE AND HEADINGS FOR METHOD OF MOMENTS (DIRECT)           
      WRITE(2,200)TITLE                   
      WRITE(2,210)   
C           
C *** METHOD OF MOMENT (DIRECT) PARAMETER ESTIMATES
      L1=S1/XN       
      L2=S2/XN       
      L3=S3/XN       
      BB=(ALOG(L3)-3.*ALOG(L1))/(ALOG(L2)-2.*ALOG(L1))                  
      CC=1./(BB-3.)  
      IF(BB.GT.3.5.AND.BB.LE.6.0)THEN     
          AA=-.23019+1.65262*CC+.20911*CC*CC-.04557*CC*CC*CC            
      ELSE IF(BB.GT.3.0.AND.BB.LE.3.5)THEN
          AA=-.47157+1.99955*CC           
      ELSE  
          WRITE(2,221)L1,L2,L3,BB         
          STE(1,6)=-1
          RKS(1,6)=-1
          CHI(1,6)=-1
          CON(1,6)=-1
          GOTO 19    
      END IF
      A=1./(AA+3.)   
      B=(ALOG(L2)-2.*ALOG(L1))/(2.*ALOG(1.-A)-ALOG(1-2.*A))             
      C=ALOG(L1)+B*ALOG(1.-A)             
      WRITE(2,220)A,B,C                   
      APAR(1,6)=A    
      BPAR(1,6)=B    
      CPAR(1,6)=C    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR MOMENT (D) ESTIMATI
      S=0   
      V=2.*B
      IF(V.GE.(0.5).AND.V.LE.(200000.))THEN        
          DO 10 I=1,N
              CALL MDCHI (P(I),V,X2,IER)  
              XX(I)=EXP(X2*A/2.+C)        
              S=S+(X(I)-XX(I))**2         
10        CONTINUE   
      ELSE  
          DO 12 I=1,N
              T=SNDV(P(I))                
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              XX(I)=EXP(C+A*E**3)         
              S=S+(X(I)-XX(I))**2         
12        CONTINUE   
      END IF
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
1     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 2     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 1                    
C           
C *** CHI-SQUARE AND K-S G-O-F TESTS FOR MOMENT (DIRECT) ESTIMATION     
2     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 14 I=1,N    
      IF(ALOG(X(I)).LT.C)THEN             
          E=((C-ALOG(X(I)))/A)**(1./3.)   
          E=-E       
      ELSE IF(ALOG(X(I)).EQ.C)THEN        
          E=0.0      
      ELSE  
          E=((ALOG(X(I))-C)/A)**(1./3.)   
      END IF
          T=(E-B**(1./3.)+1./(9.*B**(2./3.)))*3.*B**(1./6.)             
          PX(I)=TINV(T)                   
14    CONTINUE       
      CL(1)=-9999999999.                  
      V=2.*B
      IF(V.GE.(0.5).AND.V.LE.(200000.))THEN        
          DO 15 I=2,NCL                   
              CPR=CPR+PR                  
              CALL MDCHI(CPR,V,X2,IER)    
              CL(I)=EXP(X2*A/2.+C)        
15        CONTINUE   
      ELSE  
          DO 16 I=2,NCL                   
              CPR=CPR+PR                  
              T=SNDV(CPR)                 
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              CL(I)=EXP(C+A*E**3)         
16        CONTINUE   
      END IF
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(1,6)=X2    
      CON(1,6)=1.-ALPHA                   
      STE(1,6)=SE    
      RKS(1,6)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR METHOD OF MOMENT (DIRECT) ESTIMATION     
      V=2.*B
      DO 20 I=1,6    
          PR=1.-1./TT(I)                  
          T=SNDV(1.-1./TT(I))             
          T1=T       
          T2=(T*T-1.)/6.                  
          T3=2.*(T*T*T-6.*T)/(6*6*6)      
          T4=(T*T-1.)/(6*6*6)             
          T5=T/(6.**4)                    
          T6=2./(6.**6)                   
          RK=T1+T2*XG+T3*XG*XG-T4*XG**3+T5*XG**4-T6*XG**5               
          SLOPE=T2+T3*2.*XG-T4*3.*XG*XG+T5*4.*XG**3-T6*5.*XG**4         
          T7=1.0     
          T8=RK*XG   
          T9=(1.+.75*XG*XG)*(RK*RK/2.)    
          T10=3.*SLOPE*RK*(XG+.25*XG*XG*XG)        
          T11=3.*(SLOPE*SLOPE)*(2.+3.*XG*XG+(5./8.)*XG**4)              
          DELTA=T7+T8+T9+T10+T11          
          M1=C+A*B   
          M2=B*A*A   
          IF(V.GE.(0.5).AND.V.LE.(200000.))THEN    
              CALL MDCHI(PR,V,X2,IER)     
              XT(I)=EXP(X2*A/2.+C)        
          ELSE       
              XT(I)=EXP(M1+RK*SQRT(M2))   
          END IF     
          SX=SQRT(M2*DELTA/XN)            
          PSX=XT(I)*(EXP(SX)-1.)          
          NSX=-XT(I)*(EXP(-SX)-1.)        
          ST(I)=(PSX+NSX)/2.              
20    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
19    WRITE(2,299)   
C           
C *** WRITE TITLE AND HEADINGS FOR METHOD OF MOMENTS (INDIRECT)         
      WRITE(2,201)TITLE                   
      WRITE(2,210)   
C           
C *** METHOD OF MOMENTS PARAMETER ESTIMATES (INDIRECT)                  
      B=(2./YG)**2   
      A=YSD/SQRT(B)  
      C=YM-YSD*SQRT(B)                    
      WRITE(2,220)A,B,C                   
      APAR(2,6)=A    
      BPAR(2,6)=B    
      CPAR(2,6)=C    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR MOMENT (ID) ESTIMAT
      S=0   
      V=2.*B
      IF(V.GE.(0.5).AND.V.GE.(200000.))THEN        
          DO 25 I=1,N
              CALL MDCHI (P(I),V,X2,IER)  
              XX(I)=EXP(X2*A/2.+C)        
              S=S+(X(I)-XX(I))**2         
25        CONTINUE   
      ELSE  
          DO 26 I=1,N
              T=SNDV(P(I))                
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              XX(I)=EXP(C+A*E**3)         
              S=S+(X(I)-XX(I))**2         
26        CONTINUE   
      END IF
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
3     IF(L-3.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 4     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 3                    
C           
C *** CHI-SQUARE AND K-S G-O-F TESTS FOR MOMENT (INDIRECT) ESTIMATION   
4     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 29 I=1,N    
      IF(ALOG(X(I)).LT.C)THEN             
          E=((C-ALOG(X(I)))/A)**(1./3.)   
          E=-E       
      ELSE IF(ALOG(X(I)).EQ.C)THEN        
          E=0.0      
      ELSE  
          E=((ALOG(X(I))-C)/A)**(1./3.)   
      END IF
          T=(E-B**(1./3.)+1./(9.*B**(2./3.)))*3.*B**(1./6.)             
          PX(I)=TINV(T)                   
29    CONTINUE       
      CL(1)=-9999999999.                  
      V=2.*B
      IF(V.GE.(0.5).AND.V.LE.(200000.))THEN        
          DO 30 I=2,NCL                   
              CPR=CPR+PR                  
              CALL MDCHI(CPR,V,X2,IER)    
              CL(I)=EXP(X2*A/2.+C)        
30        CONTINUE   
      ELSE  
          DO 31 I=2,NCL                   
              CPR=CPR+PR                  
              T=SNDV(CPR)                 
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              CL(I)=EXP(C+A*E**3)         
31        CONTINUE   
      END IF
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(2,6)=X2    
      CON(2,6)=1.-ALPHA                   
      STE(2,6)=SE    
      RKS(2,6)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR METHOD OF MOMENT (INDIRECT) ESTIMATION   
      V=2.*B
      DO 35 I=1,6    
          PR=1.-1./TT(I)                  
          T=SNDV(PR) 
          T1=T       
          T2=(T*T-1.)/6.                  
          T3=2.*(T*T*T-6.*T)/(6*6*6)      
          T4=(T*T-1.)/(6*6*6)             
          T5=T/(6.**4)                    
          T6=2./(6.**6)                   
          RK=T1+T2*YG+T3*YG*YG-T4*YG**3+T5*YG**4-T6*YG**5               
          SLOPE=T2+T3*2.*YG-T4*3.*YG*YG+T5*4.*YG**3-T6*5.*YG**4         
          T7=1.0     
          T8=RK*YG   
          T9=(1.+.75*YG*YG)*(RK*RK/2.)    
          T10=3.*SLOPE*RK*(YG+.25*YG*YG*YG)        
          T11=3.*(SLOPE*SLOPE)*(2.+3.*YG*YG+(5./8.)*YG**4)              
          DELTA=T7+T8+T9+T10+T11          
          M2=YSD*YSD 
          IF(V.GE.(0.5).AND.V.LE.(200000.))THEN    
              CALL MDCHI(PR,V,X2,IER)     
              XT(I)=EXP(X2*A/2.+C)        
          ELSE       
              XT(I)=EXP(YM+RK*YSD)        
          END IF     
          SX=SQRT(M2*DELTA/XN)            
          PSX=XT(I)*(EXP(SX)-1.)          
          NSX=-XT(I)*(EXP(-SX)-1.)        
          ST(I)=(PSX+NSX)/2.              
35    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
C           
C *** WRITE TITLE AND HEADINGS FOR MAXIMUM LIKELIHOOD PROCEDURE         
      WRITE(2,202)TITLE                   
      WRITE(2,210)   
C           
C *** MAXIMUM LIKELIHOOD PROCEDURE        
      IT=1  
      ICOUNT=0       
      DO 40 I=1,N    
40    Y(I)=ALOG(X(I))
      GML=Y(N)*0.99  
5     ICOUNT=ICOUNT+1
      AA=0.0
      BB=0.0
      CC=0.0
      R=0.0 
      DO 45 I=1,N    
          IF(Y(I).LE.GML)GOTO 6           
          AA=AA+1.0/(Y(I)-GML)            
          BB=BB+(Y(I)-GML)                
          CC=CC+ALOG(Y(I)-GML)            
          R=R+1.0/((Y(I)-GML)**2)         
45    CONTINUE       
      BETA=1./(1.-(XN*XN)/(BB*AA))        
      ALPHA=(BB/XN)-(XN/AA)               
      D=BETA+2.      
      PSI=ALOG(D)-(1./(2.*D))-(1./(12.*D**2))+(1./(120.*D**4))-         
     &    (1./(252.*D**6))-(1./(BETA+1.))-(1./BETA)
      FCN=-XN*PSI+CC-XN*ALOG(ALPHA)       
      TRI=(1./D)+(1./(2.*D**2))+(1./(6.*D**3))-(1./(30.*D**5))+         
     &    (1./(42.*D**7))-(1./(30.*D**9))+(1./((BETA+1.)**2))+          
     &    (1./(BETA**2))                  
      VV=AA-(XN**2)/BB                    
      U=AA  
      W=(BB/XN)-(XN/AA)                   
      DU=R  
      DV=R-(XN**3)/(BB**2)                
      DW=-1.+(XN*R)/(AA**2)               
      FPN=-XN*TRI*((VV*DU-U*DV)/(VV**2))-AA-XN*DW/W
      AS=GML-(FCN/FPN)                    
      DELTA=ABS(0.00000001*AS)            
      IF (ABS(AS-GML).LT.DELTA) GOTO 7    
      IF (ICOUNT.GT.25) GOTO 6            
      GML=AS
      GOTO 5
6     CONTINUE       
      IF(IT.EQ.3)GOTO 11                  
      IF(IT.EQ.2)THEN
          IT=3       
          ICOUNT=0   
          GML=X(N)*0.5                    
          GOTO 5     
      END IF
      IF(IT.EQ.1)THEN
          IT=2       
          ICOUNT=0   
          GML=C      
          GOTO 5     
      END IF
7     CONTINUE       
      A=ALPHA        
      B=BETA
      C=AS  
      WRITE(2,220)A,B,C                   
      APAR(3,6)=A    
      BPAR(3,6)=B    
      CPAR(3,6)=C    
C           
C *** THEORETICAL SORTED DATA AND STANDARD ERROR FOR LIKELIHOOD         
      S=0.0 
      V=2.*B
      IF(V.GE.(0.5).AND.V.LE.(200000.))THEN        
          DO 50 I=1,N
              CALL MDCHI (P(I),V,X2,IER)  
              XX(I)=EXP(X2*A/2.+C)        
              S=S+(X(I)-XX(I))**2         
50        CONTINUE   
      ELSE  
          DO 51 I=1,N
              T=SNDV(P(I))                
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              XX(I)=EXP(C+A*E**3)         
              S=S+(X(I)-XX(I))**2         
51        CONTINUE   
      END IF
      SE=SQRT(S/(XN-PAR))                 
      WRITE(2,230)   
      L=N   
      SPACE=' '      
8     IF(L-8.GT.0)THEN                    
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,L-3,-1)
      ELSE  
          WRITE(2,240)(SPACE,P(I),XX(I),I=L,1,-1)  
          GOTO 9     
      END IF
      L=L-4 
      IF(L.GT.0)GOTO 8                    
C           
C *** CHI-SQUARE AND K-S G-O-F FOR LIKELIHOOD PROCEDURE                 
9     NCL=INT(1.+3.3*ALOG10(XN))          
      PR=1./NCL      
      CPR=0.0        
      DO 54 I=1,N    
      IF(ALOG(X(I)).LT.C)THEN             
          E=((C-ALOG(X(I)))/A)**(1./3.)   
          E=-E       
      ELSE IF(ALOG(X(I)).EQ.C)THEN        
          E=0.0      
      ELSE  
          E=((ALOG(X(I))-C)/A)**(1./3.)   
      END IF
          T=(E-B**(1./3.)+1./(9.*B**(2./3.)))*3.*B**(1./6.)             
          PX(I)=TINV(T)                   
54    CONTINUE       
      CL(1)=-9999999999.                  
      V=2.*B
      IF(V.GE.(0.5).AND.V.LE.(200000.))THEN        
          DO 55 I=2,NCL                   
              CPR=CPR+PR                  
              CALL MDCHI(CPR,V,X2,IER)    
              CL(I)=EXP(X2*A/2.+C)        
55        CONTINUE   
      ELSE  
          DO 56 I=2,NCL                   
              CPR=CPR+PR                  
              T=SNDV(CPR)                 
              E=B**(1./3.)-1./(9.*B**(2./3.))+T/(3.*B**(1./6.))         
              CL(I)=EXP(C+A*E**3)         
56        CONTINUE   
      END IF
      CL(NCL+1)=+9999999999.              
      CALL GOF(X,N,CL,PX,NCL,2,X2,DF,ALPHA,DN)     
      WRITE(2,250)SE,DN,X2,DF,1.-ALPHA    
      CHI(3,6)=X2    
      CON(3,6)=1.-ALPHA                   
      STE(3,6)=SE    
      RKS(3,6)=DN    
C           
C *** SELECTED RETURN PERIOD FLOWS AND THEIR STANDARD ERROR             
C *** FOR LIKELIHOOD PROCEDURE            
      ALPHA=A        
      BETA=B
      GAMMA=C        
      D=BETA+2.0     
      TRI=(1./D)+(1./(2.*D**2))+(1./(6.*D**3))-(1./(30.*D**5))+         
     &    (1./(42.*D**7))-(1./(30.*D**9))+(1./((BETA+1.)**2))+          
     &    (1./(BETA**2))                  
      H=(BETA-2.)*ALPHA**4                
      PP=2.*TRI-(2.*BETA-3.)/((BETA-1.)**2)        
      DET=PP/H       
      VARA=(1./(XN*(ALPHA**2)*DET))*((TRI/(BETA-2.))-1./((BETA-1.)**2)) 
      VARB=2./(XN*DET*(BETA-2.)*ALPHA**4) 
      VARG=(BETA*TRI-1.)/(XN*DET*ALPHA**2)
      COVAB=(-1./(XN*DET*ALPHA**3))*((1./(BETA-2.))-(1./(BETA-1.)))     
      COVAG=(1./(XN*DET*ALPHA**2))*((1./(BETA-1.))-TRI)                 
      COVBG=(-1./(XN*DET*ALPHA**3))*((BETA/(BETA-1.))-1.)               
      V=2.*B
      DO 60 I=1,6    
          PR=1.-1./TT(I)                  
          T=SNDV(PR) 
          E=BETA**(1./3.)-1./(9.*BETA**(2./3.))+T/(3.*BETA**(1./6.))    
          F=1./(3.*BETA**(2./3.))+2./(27.*BETA**(5./3.))-               
     &      T/(18.*BETA**(7./6.))         
          DXDA=E**3  
          DXDB=3.*ALPHA*E**2*F            
          DXDG=1.0   
          IF(V.GE.(0.5).AND.V.LE.(200000.))THEN    
              CALL MDCHI (PR,V,X2,IER)    
              XT(I)=EXP(X2*A/2.+C)        
          ELSE       
              XT(I)=EXP(C+A*E**3)         
          END IF     
          SX=SQRT(VARA*DXDA**2+VARB*DXDB**2+VARG*DXDG**2+2.*DXDA*       
     &    DXDB*COVAB+2.*DXDA*DXDG*COVAG+2.*DXDB*DXDG*COVBG)             
          PSX=XT(I)*(EXP(SX)-1.)          
          NSX=-XT(I)*(EXP(-SX)-1.)        
          ST(I)=(PSX+NSX)/2.              
60    CONTINUE       
      WRITE(2,260)(XT(I),I=1,6),(ST(II),II=1,6)    
C           
      WRITE(2,299)   
      RETURN
C           
11    WRITE(2,222)   
      STE(3,6)=-1    
      RKS(3,6)=-1    
      CHI(3,6)=-1    
      CON(3,6)=-1    
      WRITE(2,299)   
      RETURN
C           
C *** FORMAT STATEMENTS                   
C           
200   FORMAT(1H1/26X,'- LOG-PEARSON TYPE III -'//  
     +25X,'METHOD OF MOMENTS (DIRECT)'//  
     +1X,'TITLE: ',A72)                   
201   FORMAT(1H1/26X,'- LOG-PEARSON TYPE III -'//  
     +24X,'METHOD OF MOMENTS (INDIRECT)'//
     +1X,'TITLE: ',A72)                   
202   FORMAT(1H1/26X,'- LOG-PEARSON TYPE III -'//  
     +24X,'MAXIMUM LIKELIHOOD PROCEDURE'//
     +1X,'TITLE: ',A72)                   
210   FORMAT(//      
     +11X,52('=')/11X,'I',50X,'I'/        
     +11X,'I',28X,'(B-1)',17X,'I'/        
     +11X,'I',12X,'1',7X,'(LN(X)-C)',9X,'-(LN(X)-C)',2X,'I'/            
     +11X,'I',2X,'P(X)=',12('-'),1X,13('-'),' EXP ',10('-'),2X,'I'/     
     +11X,'I',7X,'A*X*GAMMA(B)',7X,'(B-1)',12X,'A',6X,'I'/              
     +11X,'I',25X,'A',24X,'I'/            
     +11X,'I',50X,'I'/11X,52('='))        
220   FORMAT(//30X,'A =',E12.5/30X,'B =',E12.5/30X,'C =',E12.5/)        
221   FORMAT(///1X,'METHOD OF MOMENTS (DIRECT) CAN NOT BE APPLIED'/     
     +20X,'L1=',E12.5/20X,'L2=',E12.5/20X,'L3=',E12.5/                  
     +6X,'THEREFORE ==> B=',E12.5//       
     +1X,'FOR THIS METHOD TO BE USED, B MUST BE WITHIN THE RANGE'/      
     +/20X,'3.0 < B <= 6.0')              
222   FORMAT(///1X,'MAXIMUM LIKELIHOOD DOES NOT CONVERGE WITHIN THE',   
     +' ALLOWED STEPS (25)')              
230   FORMAT(/1X,'THEORETICAL SORTED EVENTS'/      
     +6X,'NOTE: PROBABILITY IN PARENTHESIS = RANK/(N+1)'/)              
240   FORMAT(4(A1,'(',F5.4,')',E11.5,2X)) 
250   FORMAT(/6X,'STANDARD ERROR STATISTIC =',E12.5/                    
     +6X,'KOLMOGOROV-SMIRNOV STATISTIC =',F6.4/    
     +6X,'CHI-SQUARE G-O-F STATISTICS: '/ 
     +11X,'CHI-SQUARE=',F7.3/             
     +11X,'DEGREES OF FREEDOM=',F3.0/     
     +11X,'SIGNIFICANT AT 1-ALPHA=',F6.4) 
260   FORMAT(//1X,'SELECTED RETURN-PERIOD EVENTS'//
     +1X,'T,PERIOD',3X,'2',11X,'5',10X,'10',10X,'20',10X,'50',9X,'100'//
     +1X,'X',3X,6E12.5/2X,'T'/            
     +1X,'S',3X,6E12.5/2X,'T'//1X,'WHERE:'/        
     +7X,'X  = RETURN PERIOD FLOW FROM ESTIMATED DISTRIBUTION'/9X,'T'/  
     +7X,'S  = STANDARD ERROR ESTIMATE OF FLOW AT THAT RETURN PERIOD'/  
     +8X,'T')        
299   FORMAT(///)    
C           
      END   
C1
C*GOF*******************************************************************
C 
      SUBROUTINE GOF(X,N,CL,PX,NCL,NP,X2,DF,ALPHA,DN)                   
      REAL X(*),CL(302),PX(*)             
C
C     Modification to use IMSL routine MDCH which has been modified to
C     have double precision arguments.
      
      doubleprecision P
C           
C *** RUNS CHI-SQUARE AND KOLMOGOROV-SMIRNOV GOODNESS-OF-FIT TESTS      
C *** FOR ESTABLISED CLASS HISTOGRAM.     
C           
C *** CHI-SQUARE TEST
      XN=N  
      XNCL=NCL*1.    
      PR=1./XNCL     
      EF=XN*PR       
      X2=0.0
C           
      DO 10 I=1,NCL  
          XL=CL(I)   
          XH=CL(I+1) 
          OF=0.0     
          DO 20 II=1,N                    
              IF(X(II).GE.XL.AND.X(II).LT.XH)OF=OF+1.                   
20        CONTINUE   
          X2=X2+(OF-EF)**2/EF             
10    CONTINUE       
C           
      DF=XNCL-NP-1.  
      CALL MDCH (DBLE(X2),DBLE(DF),P,IER)             
      ALPHA=P        
C           
C *** KOLMOGOROV-SMIRNOV TEST             
      DN=PX(N)       
      DO 30 K=N,1,-1 
          SN=(N-K+1)/XN                   
          D1=ABS(PX(K)-SN)                
          IF(K.NE.1)D2=ABS(PX(K-1)-SN)    
          IF(D1.GT.DN)DN=D1               
          IF(D2.GT.DN)DN=D2               
30    CONTINUE       
C           
      RETURN
      END   
C1
C*ORDER*****************************************************************
C 
      SUBROUTINE ORDER (STAT,IO,CH,S)     
      INTEGER IO     
      REAL STAT(3,6),S(12)                
      CHARACTER CH(12)*7,CT*7             
C           
C *** ORDERS THE G-O-F STATISTICS         
C *** IO=+1, DESCENDING ORDER             
C *** IO=-1, ASCENDING ORDER              
C           
      CH(1)='NRM(MM)'
      S(1)=STAT(1,1) 
C           
      CH(2)='LN2(MM)'
      S(2)=STAT(1,2) 
      CH(3)='LN2(LK)'
      S(3)=STAT(2,2) 
C           
      CH(4)='LN3(MM)'
      S(4)=STAT(1,3) 
      CH(5)='LN3(LK)'
      S(5)=STAT(2,3) 
C           
      CH(6)='T1E(MM)'
      S(6)=STAT(1,4) 
      CH(7)='T1E(LK)'
      S(7)=STAT(2,4) 
C           
      CH(8)='PRS(MM)'
      S(8)=STAT(1,5) 
      CH(9)='PRS(LK)'
      S(9)=STAT(2,5) 
C           
      CH(10)='LPR(M1)'                    
      S(10)=STAT(1,6)
      CH(11)='LPR(M2)'                    
      S(11)=STAT(2,6)
      CH(12)='LPR(LK)'                    
      S(12)=STAT(3,6)
C           
      Z=0.0 
      N=12  
      IL=1  
      LU=N  
C           
1     IU=N-1
      IT=+1 
      DO 10 I=IL,IU,+1                    
          IF(S(I).LT.S(I+1))THEN          
              IT=-1  
              CT=CH(I+1)                  
              ST=S(I+1)                   
              CH(I+1)=CH(I)               
              S(I+1)=S(I)                 
              CH(I)=CT                    
              S(I)=ST
          END IF     
10    CONTINUE       
C           
      IF(IT.GT.0)GOTO 3                   
C           
      IL=IL+1        
      IT=+1 
      DO 20 I=IU,IL,-1                    
          IF(S(I).GT.S(I-1))THEN          
              IT=-1  
              CT=CH(I-1)                  
              ST=S(I-1)                   
              CH(I-1)=CH(I)               
              S(I-1)=S(I)                 
              CH(I)=CT                    
              S(I)=ST
          END IF     
20    CONTINUE       
C           
      IF(IT.GT.0)GOTO 3                   
      GOTO 1
C           
3     CONTINUE       
      IF(IO.GT.0)RETURN                   
      IL=1  
      IH=N  
      DO 30 I=1,N/2  
          CT=CH(IL)  
          ST=S(IL)   
          CH(IL)=CH(IH)                   
          S(IL)=S(IH)
          CH(IH)=CT  
          S(IH)=ST   
          IL=IL+1    
          IH=IH-1    
30    CONTINUE       
C           
      DO 40 I=1,N    
          IF(S(1).GE.Z)GOTO 4             
          CT=CH(1)   
          ST=S(1)    
          DO 41 II=1,N-1                  
              CH(II)=CH(II+1)             
              S(II)=S(II+1)               
41        CONTINUE   
          CH(N)=CT   
          S(N)=ST    
40    CONTINUE       
C           
4     CONTINUE       
      RETURN
      END   
C1
C*SNDV******************************************************************
C 
      FUNCTION SNDV (P)                   
      REAL P
C           
C *** GIVES THE STANDARD NORMAL DEVIATE SUCH THAT THE PROBABILITY       
C *** OF PRECEDING THE STANDARD NORMAL DEVIATE IS P = GIVEN PROBABILITY 
C           
      PR=P  
      IF (P.GT.0.5) PR=1.0-P              
C           
      C0=2.515517    
      C1=0.802853    
      C2=0.010328    
      D1=1.432788    
      D2=0.189269    
      D3=0.001308    
C           
      W=SQRT(ALOG(1.0/PR/PR))             
      SNDV=W-(C0+C1*W+C2*W*W)/(1.0+D1*W+D2*W*W+D3*W*W*W)                
      IF (P.GT.0.5) SNDV=(-1.0)*SNDV      
C           
      SNDV=(-1.0)*SNDV                    
C           
      RETURN
      END   
C1
C*SORT******************************************************************
C 
      SUBROUTINE SORT (X,N)               
      INTEGER N      
      REAL X(*)      
C           
C *** PUTS THE ARRAY X(I) INTO DESCENDING ORDER WITH X(1) LARGEST.      
C           
      L=N   
10    CONTINUE       
      L=L/2 
      IF(L.GE.1)THEN 
         M=N-L       
         DO 30 I=1,M 
             J=I     
20           J1=J+L  
             IF(X(J).LT.X(J1))THEN        
                 TX=X(J)                  
                 X(J)=X(J1)               
                 X(J1)=TX                 
                 J=J-L                    
                 IF(J.GT.0)GOTO 20        
             END IF  
30       CONTINUE    
         GOTO 10     
      END IF
C           
      RETURN
      END   
C1
C*TINV******************************************************************
C 
      FUNCTION TINV(T)                    
      REAL TINV      
C           
C *** FINDS THE PROBABILITY OF PRECEEDING A STANDARD NORMAL DEVIATE     
C           
      D1=.0498673470 
      D2=.0211410061 
      D3=.0032776263 
      D4=.0000380036 
      D5=.0000488906 
      D6=.0000053830 
C           
      X=ABS(T)       
      PT=1.-0.5*(1. + D1*X + D2*X**2 + D3*X**3 + D4*X**4 + D5*X**5 +    
     + D6*X**6)**(-16)                    
C           
      TINV=PT        
      IF(T.LT.0)TINV=1.-PT                
C           
      RETURN
      END   
C
C   This set of IMSL subroutines has been adapted for use with the FA1
C   package.  The routines are used to find the chi squared and inverse
C   chi-squared cdf values.  All internal math has been changed to double
C   precision, but the calls into MDCH and MDCHI still pass in single
C   precision values.
C
C
C   IMSL ROUTINE NAME   - MERRC=ERFC      
C           
C-----------------------------------------------------------------------
C           
C   COMPUTER            - CDC/SINGLE      
C           
C   LATEST REVISION     - JUNE 1, 1981    
C           
C   PURPOSE             - EVALUATE THE COMPLEMENTED ERROR FUNCTION      
C           
C   USAGE               - RESULT = ERFC(Y)
C           
C   ARGUMENTS    Y      - INPUT ARGUMENT OF THE COMPLEMENTED ERROR      
C      FUNCTION.     
C                ERFC   - OUTPUT VALUE OF THE COMPLEMENTED ERROR        
C      FUNCTION.     
C           
C   PRECISION/HARDWARE  - SINGLE/ALL      
C    NOTE - ERFC MAY NOT BE SUPPLIED BY IMSL IF IT 
C      RESIDES IN THE MATHEMATICAL SUBPROGRAM      
C      LIBRARY SUPPLIED BY THE MANUFACTURER.       
C           
C   REQD. IMSL ROUTINES - NONE REQUIRED   
C           
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C      CONVENTIONS IS AVAILABLE IN THE MANUAL      
C      INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C           
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C           
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C      APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C      EXPRESSED OR IMPLIED, IS APPLICABLE.        
C           
C-----------------------------------------------------------------------
C           
      doubleprecision FUNCTION ERFC(Y)  
C         SPECIFICATIONS FOR ARGUMENTS         
      doubleprecision    Y   
C         SPECIFICATIONS FOR LOCAL VARIABLES   
      INTEGER            ISW,I                 
      DIMENSION          P(5),Q(3),P1(8),Q1(7),P2(5),Q2(4)              
      doubleprecision    P,Q,P1,Q1,P2,Q2,XMIN,XLARGE,SSQPI,X,           
     * RES,XSQ,XNUM,XDEN,XI,xbig  
C         COEFFICIENTS FOR 0.0 .LE. Y .LT.     
C         .477        
      DATA               P(1)/-.44422647396874/,      
     1 P(2)/10.731707253648/,
     2 P(3)/15.915606197771/,
     3 P(4)/374.81624081284/,
     4 P(5)/2.5612422994823d-02/    
      DATA               Q(1)/17.903143558843/,
     1 Q(2)/124.82892031581/,
     2 Q(3)/332.17224470532/ 
C         COEFFICIENTS FOR .477 .LE. Y         
C         .LE. 4.0    
      DATA               P1(1)/7.2117582508831/,      
     1 P1(2)/43.162227222057/,      
     2 P1(3)/152.98928504694/,      
     3 P1(4)/339.32081673434/,      
     4 P1(5)/451.91895371187/,      
     5 P1(6)/300.45926102016/,      
     6 P1(7)/-1.3686485738272d-07/, 
     7 P1(8)/.56419551747897/
      DATA               Q1(1)/77.000152935229/,      
     1 Q1(2)/277.58544474399/,      
     2 Q1(3)/638.98026446563/,      
     3 Q1(4)/931.35409485061/,      
     4 Q1(5)/790.95092532790/,      
     5 Q1(6)/300.45926095698/,      
     6 Q1(7)/12.782727319629/
C         COEFFICIENTS FOR 4.0 .LT. Y          
      DATA               P2(1)/-.22695659353969/,     
     1 P2(2)/-4.9473091062325d-02/, 
     2 P2(3)/-2.9961070770354d-03/, 
     3 P2(4)/-2.2319245973418d-02/, 
     4 P2(5)/-2.7866130860965d-01/  
      DATA               Q2(1)/1.0516751070679/,      
     1 Q2(2)/.19130892610783/,      
     2 Q2(3)/1.0620923052847d-02/,  
     3 Q2(4)/1.9873320181714/
C         CONSTANTS   
      DATA               XMIN/1.0d-8/,XLARGE/5.6875d0/
C         ERFC(XBIG) .APPROX. SETAP            
      DATA               XBIG/25.90625/        
      DATA               SSQPI/.56418958354776/
C         FIRST EXECUTABLE STATEMENT           
      X = Y           
      ISW = 1         
      IF (X.GE.0.0d0) GO TO 5
      ISW = -1        
      X = -X          
    5 IF (X.LT..477d0) GO TO 10                
      IF (X.LE.4.0d0) GO TO 30                 
      IF (ISW .GT. 0) GO TO 40                 
      IF (X.LT.XLARGE) GO TO 45                
      RES = 2.0d0     
      GO TO 65        
C         ABS(Y) .LT. .477, EVALUATE           
C         APPROXIMATION FOR ERFC               
   10 IF (X.LT.XMIN) GO TO 20
      XSQ = X*X       
      XNUM = P(5)     
      DO 15 I = 1,4   
         XNUM = XNUM*XSQ+P(I)
   15 CONTINUE        
      XDEN = ((Q(1)+XSQ)*XSQ+Q(2))*XSQ+Q(3)    
      RES = X*XNUM/XDEN      
      GO TO 25        
   20 RES = X*P(4)/Q(3)      
   25 IF (ISW.EQ.-1) RES = -RES                
      RES = 1.0d0-RES 
      GO TO 65        
C         .477 .LE. ABS(Y) .LE. 4.0            
C         EVALUATE APPROXIMATION FOR ERFC      
   30 XSQ = X*X       
      XNUM = P1(7)*X+P1(8)   
      XDEN = X+Q1(7)  
      DO 35 I=1,6     
         XNUM = XNUM*X+P1(I) 
         XDEN = XDEN*X+Q1(I) 
   35 CONTINUE        
      RES = XNUM/XDEN 
      GO TO 55        
C         4.0 .LT. ABS(Y), EVALUATE            
C         MINIMAX APPROXIMATION FOR ERFC       
   40 IF (X.GT.XBIG) GO TO 60
   45 XSQ = X*X       
      XI = 1.0d0/XSQ  
      XNUM = P2(4)*XI+P2(5)  
      XDEN = XI+Q2(4) 
      DO 50 I = 1,3   
         XNUM = XNUM*XI+P2(I)
         XDEN = XDEN*XI+Q2(I)
   50 CONTINUE        
      RES = (SSQPI+XI*XNUM/XDEN)/X             
   55 RES = RES*EXP(-XSQ)    
      IF (ISW.EQ.-1) RES = 2.0d0-RES           
      GO TO 65        
   60 RES = 0.0d0     
   65 ERFC = RES      
      RETURN          
      END             
C   IMSL ROUTINE NAME   - MERFI       
C   
C-----------------------------------------------------------------------
C   
C   COMPUTER            - CDC/SINGLE  
C   
C   LATEST REVISION     - JANUARY 1, 1978 
C   
C   PURPOSE             - INVERSE ERROR FUNCTION   
C   
C   USAGE               - CALL MERFI (P,Y,IER)     
C   
C   ARGUMENTS    P      - INPUT VALUE IN THE EXCLUSIVE RANGE (-1.0,1.0) 
C                Y      - OUTPUT VALUE OF THE INVERSE ERROR FUNCTION    
C                IER    - ERROR PARAMETER (OUTPUT) 
C    TERMINAL ERROR  
C      IER = 129 INDICATES P LIES OUTSIDE THE LEGAL
C        RANGE. PLUS OR MINUS MACHINE INFINITY IS  
C        GIVEN AS THE RESULT (SIGN IS THE SIGN OF  
C        THE FUNCTION VALUE OF THE NEAREST LEGAL   
C        ARGUMENT).  
C   
C   PRECISION/HARDWARE  - SINGLE/ALL  
C   
C   REQD. IMSL ROUTINES - UERTST,UGETIO   
C   
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C      CONVENTIONS IS AVAILABLE IN THE MANUAL      
C      INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C   
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C   
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C      APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C      EXPRESSED OR IMPLIED, IS APPLICABLE.        
C   
C-----------------------------------------------------------------------
C   
      SUBROUTINE MERFI (P,Y,IER)      
C          SPECIFICATIONS FOR ARGUMENTS         
      implicit doubleprecision(a-h,o-z)
      doubleprecision    P,Y          
      INTEGER            IER          
C          SPECIFICATIONS FOR LOCAL VARIABLES   
      doubleprecision    A(65),H1,H2,H3,H4,XINF,X,SIGMA,Z,W,X3,X4,X5,   
     *                   X6,B         
      INTEGER            N,IPP,L,LB2  
      DATA               A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8),A(9),  
     *                   A(10),A(11),A(12),A(13),A(14),A(15),A(16),     
     *                   A(17),A(18),A(19),A(20),A(21),A(22),A(23)      
     *                   /.99288537662,.12046751614,                    
     *                   .16078199342D-01,.26867044372D-02,             
     *                   .49963473024D-03,.98898218599D-04,             
     *                   .20391812764D-04,.4327271618D-05,              
     *                   .938081413D-06,.206734721D-06,                 
     *                   .46159699D-07,.10416680D-07,                   
     *                   .2371501D-08,.543928D-09, 
     *                   .125549D-09,.29138D-10,
     *                   .6795D-11,.1591D-11,   
     *                   .374D-12,.88D-13,      
     *                   .21D-13,.5D-14,        
     *                   .1D-14/      
      DATA               A(24),A(25),A(26),A(27),A(28),A(29),A(30),     
     *                   A(31),A(32),A(33),A(34),A(35),A(36),A(37),     
     *                   A(38),A(39),A(40)      
     *                   /.91215880342D00,-.16266281868D-01,            
     *                   .43355647295D-03,.21443857007D-03,             
     *                   .2625751076D-05,-.302109105D-05,               
     *                   -.12406062D-07,.62406609D-07,                  
     *                   -.540125D-09,-.142328D-08,
     *                   .34384D-10,.33585D-10, 
     *                   -.1458D-11,-.81D-12,   
     *                   .53D-13,.2D-13,        
     *                   -.2D-14/     
      DATA               A(41),A(42),A(43),A(44),A(45),A(46),A(47),     
     *                   A(48),A(49),A(50),A(51),A(52),A(53),A(54),     
     *                   A(55),A(56),A(57),A(58),A(59),A(60),A(61),     
     *                   A(62),A(63),A(64),A(65)
     *                   /.95667970902,-.023107004309,                  
     *                   -.43742360975D-02,-.57650342265D-03,           
     *                   -.10961022307D-04,.25108547025D-04,            
     *                   .10562336068D-04,.275441233D-05,               
     *                   .432484498D-06,-.20530337D-07,                 
     *                   -.43891537D-07,-.1768401D-07,                  
     *                   -.3991289D-08,-.186932D-09,                    
     *                   .272923D-09,.132817D-09,  
     *                   .31834D-10,.1670D-11,  
     *                   -.2036D-11,-.965D-12,  
     *                   -.22D-12,-.1D-13,      
     *                   .13D-13,.6D-14,        
     *                   .1D-14/      
      DATA               H1,H2,H3,H4/-1.5488130424,
     *                   2.5654901231,-.55945763133,                    
     *                   2.2879157163/
      DATA               XINF/1.0d30/                    
CFIRST EXECUTABLE STATEMENT           
      X = P   
      IER = 0 
      SIGMA = dSIGN(1.d0,X)              
      IF (.NOT.(X.GT.-1.d0.AND.X.LT.1.d0)) GO TO 35 
      Z = dABS(X) 
      IF(Z.GT..8d0) GO TO 20            
      W = Z*Z/.32-1.                  
      N = 22  
      IPP = 1 
      L = 1   
    5 LB2 = 1 
      X3 = 1. 
      X4 = W  
      X6 = A(IPP)
   10 X6 = X6 + A(IPP+LB2) * X4       
      X5 = X4 * W * 2.-X3             
      X3 = X4 
      X4 = X5 
      LB2 = LB2 + 1                   
      IF (LB2 .LE. N) GO TO 10        
      GO TO (15,30),L                 
   15 Y = Z * X6 * SIGMA              
      GO TO 9005 
   20 B = SQRT(-dLOG(1.-Z*Z))         
      IF (Z .GT. .9975d0) GO TO 25      
      W = H1*B+H2
      IPP = 24
      L = 2   
      N = 16  
      GO TO 5 
   25 W = H3 * B + H4                 
      IPP = 41
      N = 24  
      L = 2   
      GO TO 5 
   30 Y = B * X6 * SIGMA              
      GO TO 9005 
   35 Y = SIGMA*XINF                  
      IER = 129  
 9000 CONTINUE
C     CALL UERTST(IER,6HMERFI )       
 9005 RETURN  
      END     
C   IMSL ROUTINE NAME   - MGAMA=GAMMA       
C              
C-----------------------------------------------------------------------
C              
C   COMPUTER            - CDC/SINGLE        
C              
C   LATEST REVISION     - JUNE 1, 1982      
C              
C   PURPOSE             - EVALUATE THE GAMMA FUNCTION                   
C              
C   USAGE               - RESULT = GAMMA(X) 
C              
C   ARGUMENTS    X      - INPUT ARGUMENT.   
C      GAMMA IS SET TO MACHINE INFINITY, WITH THE  
C      PROPER SIGN, WHENEVER  
C      X IS ZERO,    
C      X IS A NEGATIVE INTEGER,                   
C      ABS(X) .LE. XMIN, OR  
C      ABS(X) .GE. XMAX.     
C      XMIN IS OF THE ORDER OF 10**(-39) AND      
C      XMAX IS AT LEAST 34.8. THE EXACT VALUES    
C      OF XMIN AND XMAX MAY ALLOW LARGER RANGES   
C      FOR X ON SOME COMPUTERS.                   
C      SEE THE PROGRAMMING NOTES IN THE MANUAL    
C      FOR THE EXACT VALUES. 
C                GAMMA  - OUTPUT SINGLE PRECISION VALUE OF THE GAMMA    
C      FUNCTION.       
C              
C   PRECISION/HARDWARE  - SINGLE/ALL        
C    NOTE - GAMMA MAY NOT BE SUPPLIED BY IMSL IF   
C      IT RESIDES IN THE MATHEMATICAL SUBPROGRAM   
C      LIBRARY SUPPLIED BY THE MANUFACTURER.       
C              
C   REQD. IMSL ROUTINES - UERTST,UGETIO     
C              
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C      CONVENTIONS IS AVAILABLE IN THE MANUAL      
C      INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C              
C   REMARKS      AN ERROR MESSAGE PRINTED BY UERTST FROM GAMMA SHOULD   
C                BE INTERPRETED AS FOLLOWS  
C                IER    - ERROR INDICATOR   
C    TERMINAL ERROR    
C      IER = 129 INDICATES THAT THE ABSOLUTE VALUE 
C   OF THE INPUT ARGUMENT X WAS SPECIFIED 
C   GREATER THAN OR EQUAL TO XMAX. GAMMA  
C   IS SET TO MACHINE INFINITY.           
C      IER = 130 INDICATES THAT THE INPUT ARGUMENT 
C   X WAS SPECIFIED AS ZERO OR A NEGATIVE 
C   INTEGER OR THAT THE ABSOLUTE VALUE OF 
C   INPUT ARGUMENT WAS LESS THAN OR EQUAL 
C   XMIN. GAMMA IS SET TO MACHINE INFINITY
C   WITH THE PROPER SIGN FOR THE GAMMA    
C   FUNCTION. IF X IS ZERO OR AN EVEN     
C   NEGATIVE INTEGER, GAMMA HAS A NEGATIVE
C   SIGN. OTHERWISE IT HAS A POSITIVE SIGN
C              
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C              
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C      APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C      EXPRESSED OR IMPLIED, IS APPLICABLE.        
C              
C-----------------------------------------------------------------------
C              
      doubleprecision FUNCTION GAMMA (X)               
      implicit doubleprecision (a-h,o-z)
C     SPECIFICATIONS FOR ARGUMENTS          
      doubleprecision    X                  
C     SPECIFICATIONS FOR LOCAL VARIABLES    
      doubleprecision    P,Q,P4,BIG1,PI,XMIN,XINF,XSIGN,Y,T,R,A,TOP,    
     *                   DEN,B              
      INTEGER            I,IEND,IEND1,IEND2,IER,J        
      LOGICAL            MFLAG              
      DIMENSION          P(7),Q(6),P4(5)    
C     COEFFICIENTS FOR MINIMAX              
C     APPROXIMATION TO GAMMA(X),            
C     2.0 .LE. X .LE. 3.0                   
      DATA               P(1)/3.4109112397125d+01/,      
     1                   P(2)/ -4.8341273405598d+01/,    
     2                   P(3)/4.3005887829594d+02/,      
     3                   P(4)/-5.5688734338586d+01/,     
     4                   P(5)/2.0585220673644d+03/,      
     5                   P(6)/7.7192407739800d-01/,      
     6                   P(7)/-3.1721064346240d+00/      
      DATA               Q(1)/2.4455138217658d+02/,      
     1                   Q(2)/-1.0174768492818d+03/,     
     2                   Q(3)/1.1615998272754d+03/,      
     3                   Q(4)/2.0512896777972d+03/,      
     4                   Q(5)/6.8080353498091d-01/,      
     5                   Q(6)/-2.5386729086746d+01/      
C     COEFFICIENTS FOR MINIMAX              
C     APPROXIMATION TO LN(GAMMA(X)),        
C     12.0 .LE. X      
      DATA               P4(1)/9.1893853320467d-01/,     
     1                   P4(2)/8.3333333333267d-02/,     
     2                   P4(3)/-2.7777776505151d-03/,    
     3                   P4(4)/7.9358449768d-04/,        
     4                   P4(5)/-5.82399983d-04/          
      DATA               IEND/7/,IEND1/6/,IEND2/5/       
      DATA               XINF/.123456789e30 /            
      DATA               PI/3.1415926535898/
C     GAMMA(XMIN) .APPROX. XINF             
C     GAMMA(BIG1) .APPROX. XINF             
      DATA               XMIN/0.0/          
      DATA               BIG1/177.803/      
C     FIRST EXECUTABLE STATEMENT            
      IER = 0  
      MFLAG = .FALSE.  
      T = X    
      IF (ABS(T).GT.XMIN) GO TO 5           
      IER = 130
      GAMMA = XINF     
      IF (T.LE.0.0d0) GAMMA = -XINF           
      GO TO 9000       
    5 IF (ABS(T).LT.BIG1) GO TO 10          
      IER = 129
      GAMMA = XINF     
      GO TO 9000       
   10 IF (T.GT.0.0d0) GO TO 25                
C     ARGUMENT IS NEGATIVE                  
      MFLAG = .TRUE.   
      T = -T   
      R = dINT(T)      
      XSIGN = 1.0      
      IF (DMOD(R,2.0d0) .EQ. 0.0) XSIGN = -1. 
      R = T-R  
      IF (R.NE.0.0d0) GO TO 20                
      IER = 130
      GAMMA = XINF     
      IF (XSIGN.EQ.-1.0d0) GAMMA = -XINF      
      GO TO 9000       
C     ARGUMENT IS NOT A NEGATIVE INTEGER    
   20 R = PI/SIN(R*PI)*XSIGN                
      T = T+1.0
C     EVALUATE APPROXIMATION FOR GAMMA(T)   
C       T .GT. XMIN    
   25 IF (T.GT.12.0d0) GO TO 60               
      I = T    
      A = 1.0  
      IF (I.GT.2) GO TO 40                  
      I = I+1  
      GO TO (30,35,50),I                    
C     0.0 .LT. T .LT. 1.0                   
   30 A = A/(T*(T+1.0))
      T = T+2.0
      GO TO 50 
C     1.0 .LE. T .LT. 2.0                   
   35 A = A/T  
      T = T+1.0
      GO TO 50 
C     3.0 .LE. T .LE. 12.0                  
   40 DO 45 J=3,I      
         T = T-1.0     
         A = A*T       
   45 CONTINUE 
C     2.0 .LE. T .LE. 3.0                   
   50 TOP = P(IEND1)*T+P(IEND)              
      DEN = T+Q(IEND1) 
      DO 55 J=1,IEND2  
         TOP = TOP*T+P(J)                   
         DEN = DEN*T+Q(J)                   
   55 CONTINUE 
      Y = (TOP/DEN)*A  
      IF (MFLAG) Y = R/Y                    
      GAMMA = Y
      GO TO 9005       
C     T .GT. 12.0      
   60 TOP = DLOG(T)    
      TOP = (T-1.5)*(TOP-1.0)+TOP-1.5       
      T = 1.0/T
      B = T*T  
      Y = (((P4(5)*B+P4(4))*B+P4(3))*B+P4(2))*T+P4(1)+TOP               
      Y = DEXP(Y)       
      IF (MFLAG) Y = R/Y                    
      GAMMA = Y
      GO TO 9005       
 9000 CONTINUE 
C     CALL UERTST(-IER,6H MGAMA)            
C     CALL UERTST(IER,6HGAMMA )             
 9005 RETURN   
      END      
C   IMSL ROUTINE NAME   - MDCH          
C       
C-----------------------------------------------------------------------
C       
C   COMPUTER            - CDC/SINGLE    
C       
C   LATEST REVISION     - JUNE 1, 1980  
C       
C   PURPOSE             - CHI-SQUARED PROBABILITY DISTRIBUTION FUNCTION 
C       
C   USAGE               - CALL MDCH (CS,DF,P,IER)  
C       
C   ARGUMENTS    CS     - INPUT VALUE FOR WHICH THE PROBABILITY IS      
C      COMPUTED. CS MUST BE GREATER THAN OR EQUAL  
C      TO ZERO.    
C                DF     - INPUT NUMBER OF DEGREES OF FREEDOM OF THE     
C      CHI-SQUARED DISTRIBUTION. DF MUST BE GREATER
C      THAN OR EQUAL TO .5 AND LESS THAN OR EQUAL  
C      TO 200,000. 
C                P      - OUTPUT PROBABILITY THAT A RANDOM VARIABLE     
C      WHICH FOLLOWS THE CHI-SQUARED DISTRIBUTION  
C      WITH DF DEGREES OF FREEDOM IS LESS THAN OR  
C      EQUAL TO CS.
C                IER    - ERROR PARAMETER. (OUTPUT)
C    TERMINAL ERROR
C      IER = 129 INDICATES THAT CS OR DF WAS       
C        SPECIFIED INCORRECTLY.                    
C    WARNING ERROR 
C      IER = 34 INDICATES THAT THE NORMAL PDF      
C        WOULD HAVE PRODUCED AN UNDERFLOW.         
C       
C   PRECISION/HARDWARE  - SINGLE/ALL    
C       
C   REQD. IMSL ROUTINES - H32/MDNOR,MERRC=ERFC,MGAMAD=DGAMMA,UERTST,    
C      UGETIO      
C  - H36,H48,H60/MDNOR,MERRC=ERFC,MGAMA=GAMMA,     
C      UERTST,UGETIO 
C       
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C      CONVENTIONS IS AVAILABLE IN THE MANUAL      
C      INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C       
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C       
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C      APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C      EXPRESSED OR IMPLIED, IS APPLICABLE.        
C       
C-----------------------------------------------------------------------
C       
      SUBROUTINE MDCH (CS,DF,P,IER)     
C  SPECIFICATIONS FOR ARGUMENTS         
C*    Note that this routine now uses only doubleprecision arguments
      integer ier
      doubleprecision    CS,DF,P        
C  SPECIFICATIONS FOR LOCAL VARIABLES   
      doubleprecision    A,Z,DGAM,EPS,W,W1,B,Z1,ONE,HALF,RINFM,X,C      
      doubleprecision    GAMMA,THRD,PT2,FUNC       
      doubleprecision    thrten       
      DATA               EPS/1.d-14/,HALF/5.d-1/   
      DATA               THRTEN/13.d0/,ONE/1.d0/   
      DATA               THRD/.33333333333333d0/   
      DATA               PT2/.22222222222222d0/    
      DATA               RINFM/1.d20/       
      FUNC(W,A,Z)=W*DEXP(A*DLOG(Z)-Z)   
C  FIRST EXECUTABLE STATEMENT           
C  TEST FOR INVALID INPUT VALUES        
      IF (DF .GE. .5 .AND. DF .LE. 2.d5 .AND. CS .GE. 0.0) GO TO 5      
      IER=129      
      P=RINFM      
      GO TO 9000   
    5 IER=0        
C  SET P=0. IF CS IS LESS THAN OR       
C  EQUAL TO 10.**(-12)                  
      IF (CS .GT. 1.d-12) GO TO 15      
   10 P=0.0        
      GO TO 9005   
   15 IF(DF.LE.100.d0) GO TO 20           
C  USE NORMAL DISTRIBUTION APPROXIMATION
C  FOR LARGE DEGREES OF FREEDOM         
      IF(CS.LT.2.0d0) GO TO 10            
      X=((CS/DF)**THRD-(ONE-PT2/DF))/SQRT(PT2/DF)  
      IF (X .GT. 5.0d0) GO TO 50          
      IF (X .LT. -11.313708d0) GO TO 55   
      CALL MDNOR (X,P)                  
      GO TO 9005   
C  INITIALIZATION FOR CALCULATION USING 
C  INCOMPLETE GAMMA FUNCTION            
   20 IF (CS .GT. 200.d0) GO TO 50        
      A=HALF*DF    
      Z=HALF*CS    
      C = A        
      DGAM = GAMMA(C)                   
      W=DMAX1(HALF*A,THRTEN)            
      IF (Z .GE. W) GO TO 35            
      IF (DF .GT. 25. .AND. CS .LT. 2.) GO TO 10   
C  CALCULATE USING EQUATION NO. 6.5.29  
      W=ONE/(DGAM*A)                    
      W1=W
         DO 25 I=1,50                   
         B=I       
         W1=W1*Z/(A+B)                  
         IF (W1 .LE. EPS*W) GO TO 30    
         W=W+W1    
   25    CONTINUE  
   30 P=FUNC(W,A,Z)
      GO TO 9005   
C  CALCULATE USING EQUATION NO. 6.5.32  
   35 Z1=ONE/Z     
      B=A-ONE      
      W1=B*Z1      
      W=ONE+W1     
         DO 40 I=2,50                   
         B=B-ONE   
         W1=W1*B*Z1
         IF (W1 .LE. EPS*W) GO TO 45    
         W=W+W1    
   40    CONTINUE  
   45 W=Z1*FUNC(W,A,Z)                  
      P=ONE-W/DGAM 
      GO TO 9005   
   50 P=1.0        
      GO TO 9005   
C  WARNING ERROR - UNDERFLOW WOULD HAVE 
C  OCCURRED        
   55 P=0.0        
      IER=34       
 9000 CONTINUE     
C     CALL UERTST (IER,6HMDCH  )        
c9005 RETURN       
9005  continue
      zqxP=P
      return
      END 
C   IMSL ROUTINE NAME   - MDCHI           
C           
C-----------------------------------------------------------------------
C           
C   COMPUTER            - CDC/SINGLE      
C           
C   LATEST REVISION     - JANUARY 1, 1978 
C           
C   PURPOSE             - INVERSE CHI-SQUARED PROBABILITY DISTRIBUTION  
C      FUNCTION      
C           
C   USAGE               - CALL MDCHI (P,DF,X,IER)  
C           
C   ARGUMENTS    P      - INPUT PROBABILITY IN THE EXCLUSIVE RANGE      
C      (0,1)         
C                DF     - INPUT NUMBER OF DEGREES OF FREEDOM. DF MUST BE
C      IN THE EXCLUSIVE RANGE (.5,200000.).        
C                X      - OUTPUT CHI-SQUARED VALUE, SUCH THAT A RANDOM  
C      VARIABLE, DISTRIBUTED AS CHI-SQUARED WITH   
C      DF DEGREES OF FREEDOM, WILL BE LESS THAN OR 
C      EQUAL TO X WITH PROBABILITY P.              
C                IER    - ERROR PARAMETER. (OUTPUT)
C    TERMINAL ERROR  
C      IER = 129 INDICATES THAT THE BOUNDS WHICH   
C        ENCLOSED P COULD NOT BE FOUND WITHIN 20   
C        (NCT) ITERATIONS     
C      IER = 130 INDICATES THAT AN ERROR OCCURRED  
C        IN IMSL ROUTINE MDNRIS (P IS NOT A VALID  
C        VALUE)      
C      IER = 131 INDICATES THAT AN ERROR OCCURRED  
C        IN IMSL ROUTINE MDCH 
C      IER = 132 INDICATES THAT THE VALUE X COULD  
C        NOT BE FOUND WITHIN 50 (ITMAX) ITERATIONS,
C        SO THAT THE ABSOLUTE VALUE OF P1-P WAS    
C        LESS THAN OR EQUAL TO EPS. (P1 IS THE     
C        CALCULATED PROBABILITY AT X, EPS = .0001) 
C           
C   PRECISION/HARDWARE  - SINGLE/ALL      
C           
C   REQD. IMSL ROUTINES - H32/MDCH,MDNOR,MDNRIS,MERFI,MERRC=ERFC,       
C      MGAMAD=DGAMMA,UERTST,UGETIO                 
C  - H36,H48,H60/MDCH,MDNOR,MDNRIS,MERFI,          
C      MERRC=ERFC,MGAMA=GAMMA,UERTST,UGETIO        
C           
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C      CONVENTIONS IS AVAILABLE IN THE MANUAL      
C      INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C           
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C           
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C      APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C      EXPRESSED OR IMPLIED, IS APPLICABLE.        
C        
C-----------------------------------------------------------------------
C        
      SUBROUTINE MDCHI   (zqxP,zqxDF,zqxX,IER)              
      implicit doubleprecision (a-h,o-z)
      real zqxP,zqxDF,zqxX
C        
      DATA               EPS/.0001/,ITMAX/50/,NCT/20/,NSIG/5/           
      P=zqxP
      DF=zqxDF
      
C             FIRST EXECUTABLE STATEMENT           
      IC = 0                  
      IER = 0                 
C             ESTIMATE STARTING X                  
      CALL MDNRIS (P,XP,IER)  
      IF (IER .NE. 0) GO TO 80
      DFF = .22222222222222   
      DFF = DFF/DF            
      D = dSQRT(DFF)           
      X  = DF*(1. -DFF + XP*D)**3                  
C             IS THE CASE ASYMPTOTICALLY NORMAL    
      IF (DF .GE. 40.) GO TO 9005                  
C             FIND BOUNDS (IN X) WHICH ENCLOSE P   
      NCNT = 0                
      IST = 0                 
      ISW = 0                 
      DX = X * .125           
    5 IF (X) 10,15,20         
   10 X = 0.0                 
      DX = -DX                
      GO TO 20                
   15 DX = .1                 
   20 CALL MDCH (X,DF,P1,IER) 
      DX = DX + DX            
      IF (IER .NE. 0) GO TO 85
      CSS = X                 
      NCNT = NCNT + 1         
      IF (NCNT .GT. NCT) GO TO 90                  
      IF (P1 - P) 25,9005,30  
   25 X = X + DX              
      ISW = 1                 
      IF (IST .EQ. 0) GO TO 5 
      GO TO 35                
   30 X = X - DX              
      IST = 1                 
      IF (ISW .EQ. 0) GO TO 5 
      XR = CSS                
      XL = X                  
      GO TO 40                
   35 XL = CSS                
      XR = X                  
C             PREPARE FOR ITERATION TO FIND X      
   40 EPSP = 10.**(-NSIG)     
      IF (XL .LT. 0.) XL = 0.0
      CALL MDCH (XL,DF,P1,IER)
      IF (IER .NE. 0) GO TO 85
      FXL = P1 - P            
      CALL MDCH (XR,DF,P1,IER)
      IF (IER .NE. 0) GO TO 85
      FXR = P1-P              
      IF (FXL*FXR .NE. 0.0d0) GO TO 45               
      X = XR                  
      IF (FXL .EQ. 0.0) X = XL
      GO TO 9005              
   45 IF (DF .LE. 2. .OR. P .GT. .98) GO TO 50     
C             REGULA FALSI METHOD                  
      X = XL + FXL*(XR-XL)/(FXL-FXR)               
      GO TO 55                
C             BISECTION METHOD
   50 X = (XL+XR) * .5        
   55 CALL MDCH (X,DF,P1,IER) 
      IF (IER .NE. 0) GO TO 85
      FCS = P1-P              
      IF (dABS(FCS) .GT. EPS) GO TO 60              
      GO TO 9005              
   60 IF (FCS * FXL .GT. 0.0d0) GO TO 65             
      XR = X                  
      FXR = FCS               
      GO TO 70                
   65 XL = X                  
      FXL = FCS               
   70 IF (XR-XL .GT. EPSP*dABS(XR)) GO TO 75        
      GO TO 9005              
   75 IC = IC+1               
      IF (IC .LE. ITMAX) GO TO 45                  
      IER = 132               
      GO TO 9000              
C             ERROR RETURNED FROM MDNRIS           
   80 IER = 130               
      GO TO 9000              
C             ERROR RETURNED FROM MDCH             
   85 IER = 131               
      GO TO 9000              
   90 IER = 129               
 9000 CONTINUE                
C     CALL UERTST (IER,6HMDCHI )                   
c9005 RETURN                  
9005  continue
      zqxX=X
      return
      END
C   IMSL ROUTINE NAME   - MDNOR           
C        
C-----------------------------------------------------------------------
C        
C   COMPUTER            - CDC/SINGLE               
C        
C   LATEST REVISION     - JANUARY 1, 1978          
C        
C   PURPOSE             - NORMAL OR GAUSSIAN PROBABILITY DISTRIBUTION   
C      FUNCTION               
C        
C   USAGE               - CALL MDNOR (Y,P)         
C        
C   ARGUMENTS    Y      - INPUT VALUE AT WHICH FUNCTION IS TO BE        
C      EVALUATED.             
C                P      - OUTPUT PROBABILITY THAT A RANDOM VARIABLE     
C      HAVING A NORMAL (0,1) DISTRIBUTION WILL BE  
C      LESS THAN OR EQUAL TO Y.                    
C        
C   PRECISION/HARDWARE  - SINGLE/ALL               
C        
C   REQD. IMSL ROUTINES - MERRC=ERFC               
C        
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C      CONVENTIONS IS AVAILABLE IN THE MANUAL      
C      INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C        
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C        
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C      APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C      EXPRESSED OR IMPLIED, IS APPLICABLE.        
C        
C-----------------------------------------------------------------------
C        
      SUBROUTINE MDNOR  (Y,P) 
C             SPECIFICATIONS FOR ARGUMENTS         
      doubleprecision    P,Y  
C             SPECIFICATIONS FOR LOCAL VARIABLES   
      doubleprecision erfc
      doubleprecision    SQR1D2                    
      DATA               SQR1D2/.70710678118655/   
C             FIRST EXECUTABLE STATEMENT           
      P = .5 * ERFC(-Y*SQR1D2)
      RETURN                  
      END
C   IMSL ROUTINE NAME   - MDNRIS                    
C           
C-----------------------------------------------------------------------
C           
C   COMPUTER            - CDC/SINGLE                
C           
C   LATEST REVISION     - JUNE 1, 1981              
C           
C   PURPOSE             - INVERSE STANDARD NORMAL (GAUSSIAN)            
C       PROBABILITY DISTRIBUTION FUNCTION           
C           
C   USAGE               - CALL MDNRIS (P,Y,IER)     
C           
C   ARGUMENTS    P      - INPUT VALUE IN THE EXCLUSIVE RANGE (0.0,1.0)  
C                Y      - OUTPUT VALUE OF THE INVERSE NORMAL (0,1)      
C       PROBABILITY DISTRIBUTION FUNCTION           
C                IER    - ERROR PARAMETER (OUTPUT)  
C     TERMINAL ERROR            
C       IER = 129 INDICATES P LIES OUTSIDE THE LEGAL
C         RANGE. PLUS OR MINUS MACHINE INFINITY IS  
C         GIVEN AS THE RESULT (SIGN IS THE SIGN OF  
C         THE FUNCTION VALUE OF THE NEAREST LEGAL   
C         ARGUMENT).            
C           
C   PRECISION/HARDWARE  - SINGLE/ALL                
C           
C   REQD. IMSL ROUTINES - MERFI,UERTST,UGETIO       
C           
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
C       CONVENTIONS IS AVAILABLE IN THE MANUAL      
C       INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
C           
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
C           
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
C       APPLIED TO THIS CODE. NO OTHER WARRANTY,    
C       EXPRESSED OR IMPLIED, IS APPLICABLE.        
C           
C-----------------------------------------------------------------------
C           
      SUBROUTINE MDNRIS (P,Y,IER)                   
      implicit doubleprecision(a-h,o-z)
C              SPECIFICATIONS FOR ARGUMENTS         
      doubleprecision    P,Y    
      INTEGER            IER    
C              SPECIFICATIONS FOR LOCAL VARIABLES   
      doubleprecision    D(25),A,B,EPS,W,H3,H4,X3,X4,X5,X6              
      doubleprecision    SIGMA,SQRT2,X,XINF         
      DATA               XINF/1.0d30/
      DATA               SQRT2/1.4142135623731/     
      DATA               EPS/1.0D-12/ 
      DATA               D(1)/.95667970902049/      
      DATA               D(2)/-.02310700430907/     
      DATA               D(3)/-.00437423609751/     
      DATA               D(4)/-.00057650342265/     
      DATA               D(5)/-.00001096102231/     
      DATA               D(6)/.00002510854703/      
      DATA               D(7)/.00001056233607/      
      DATA               D(8)/.00000275441233/      
      DATA               D(9)/.00000043248450/      
      DATA               D(10)/-.00000002053034/    
      DATA               D(11)/-.00000004389154/    
      DATA               D(12)/-.00000001768401/    
      DATA               D(13)/-.00000000399129/    
      DATA               D(14)/-.00000000018693/    
      DATA               D(15)/.00000000027292/     
      DATA               D(16)/.00000000013282/     
      DATA               D(17)/.00000000003183/     
      DATA               D(18)/.00000000000167/     
      DATA               D(19)/-.00000000000204/    
      DATA               D(20)/-.00000000000097/    
      DATA               D(21)/-.00000000000022/    
      DATA               D(22)/-.00000000000001/    
      DATA               D(23)/.00000000000001/     
      DATA               D(24)/.00000000000001/     
      DATA               D(25)/.00000000000000/     
      DATA               H3/-.55945763132983/       
      DATA               H4/2.2879157162634/        
C              FIRST EXECUTABLE STATEMENT           
      IER = 0                   
      IF (P .GT. 0.0 .AND. P .LT. 1.0) GO TO 5      
      IER = 129                 
      SIGMA = DSIGN(1.0D0,P)       
      Y = SIGMA * XINF          
      GO TO 9000                
    5 IF(P.LE.EPS) GO TO 10     
      X = 1.0 -(P + P)          
      CALL MERFI (X,Y,IER)      
      Y = -SQRT2 * Y            
      GO TO 9005                
C              P TOO SMALL, COMPUTE Y DIRECTLY      
   10 A = P+P                   
      B = DSQRT(-DLOG(A+(A-A*A)))
      W = H3 * B + H4           
      X3 = 1.0                  
      X4 = W
      X6 = D(1)                 
      DO 15 I=2,25              
         X6 = X6 + D(I) * X4    
         X5 = X4 * W * 2. - X3  
         X3 = X4                
         X4 = X5                
   15 CONTINUE                  
      Y = B*X6                  
      Y = -SQRT2 * Y            
      GO TO 9005                
 9000 CONTINUE                  
c     CALL UERTST(IER,6HMDNRIS) 
 9005 RETURN
      END