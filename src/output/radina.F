c
c Copyright (c) 1998-2012 The OPIUM Group
c
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
c
      SUBROUTINE RADINA (R,F,M,NP,H,P,INV,IEN)                          
C RADIAL INTEGRAL ON LOGARITHMIC GRID- CORRECTED TRAPEZOIDAL RULE       
C THIS METHOD IS PARTICULARLY CONVENIENT AND EFFECTIVE FOR DEF INT--    
C USES ENDPT CORRECTIONS IN TERMS OF DIFFERENCES UP TO FOURTH ORDER.    
C R = R(Y) WHERE Y IS UNIFORM VARIABLE; F = INPUT FUNCTION ON Y GRID;   
C M = POWER OF R USED IF MOMENTS ARE DESIRED (OTHERWISE=0); NP=NUMBER   
C OF GRID POINTS USED IN INTEGRAL; H = UNIFORM SPACING OF Y VARIABLE;   
C P = POWER OF R IN ASYMPTOTIC FORM OF F(R) NEAR R=0: F(R)=A*R**P.      
C ON RETURN, P IS REDEFINED AS VALUE OF DEFINITE INTEGRAL.              
C ASYMPTOTIC FORM ALLOWS ANALYTIC INTEGRAL FROM R=0 TO R=R(1).          
C                                             DOUG ALLAN                
C MODIFIED JAN 1990 BY ANDREW RAPPE.                                    
C NOW HIGHER ORDER INTEGRATIONS SCHEME IS USED THROUGHOUT RANGE OF      
C INTEGRATION.  SPECIAL FORMULA FOR ENDPOINTS IS REMOVED.               
C FOR FUNCTIONS WHICH ARE PIECEWISE CONTINUOUS, RADINA INTEGRATES       
C F OVER INV INTERVALS. IF INV>1, THE ARRAY IEN CONTAINS THE GRIDPT #'S 
C OF THE (INV-1) POINTS OF DISCONTINUITY.                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                         
#include "fortdim.h"
      DIMENSION R(npdm), F(npdm), IEN(10)                                      
C     WRITE (7,*) 'START RADINA WITH NP=',NP,'H=',H,'P=',P,'M=',M       
C STEP 1:  INTEGRATE FROM R=0 TO R=R(1) ANALYTICALLY.                   
      N = M+1                                                           
      S1 = (R(1)**N)*F(1)                                               
      P = S1/(H*(N+P))                                                  
C     WRITE (7,*) 'ANALYTIC PART',P,H,R(1),N,H                          
C STEP 2:  INTEGRATE ONE INTERVAL AT A TIME.                            
C          THE STANDARD METHOD, WHICH IS CORRECT FOR POLYNOMIALS UP TO  
C          R**(IST+1), NEEDS (K * IST)+1 POINTS IN AN INTERVAL.         
C          TO ADJUST FOR THIS,                                          
C          THE LAST IST+1 TO 2*IST-1                                    
C          POINTS ARE HANDLED SPECIALLY.                                
C THE NEXT LINES ARE SPECIFIC TO IST=4.                                 
      IST = 4                                                           
      A51 = 16.0E0/7.0E0                                                
      A52 = 6.0E0/7.0E0                                                 
      A53 = 28.0E0/45.0E0                                               
      A8A = 122663.0E0/403200.0E0/A53                                   
      A8B = 584245.0E0/403200.0E0/A53                                   
      A8C = 216087.0E0/403200.0E0/A53                                   
      A8D = 488205.0E0/403200.0E0/A53                                   
      A7A =  2952.0E0/10080.0E0/A53                                     
      A7B = 15552.0E0/10080.0E0/A53                                     
      A7C =  1944.0E0/10080.0E0/A53                                     
      A7D = 19584.0E0/10080.0E0/A53                                     
      A6A =  475.0E0/1440.0E0/A53                                       
      A6B = 1875.0E0/1440.0E0/A53                                       
      A6C = 1250.0E0/1440.0E0/A53                                       
      IF (INV.LE.0) INV = 1                                             
      IEN(INV) = FLOAT(NP)                                              
      LEND = 1                                                          
      DO 1 L = 1,INV                                                    
         IREND = INT(IEN(L))                                            
         NUMPM1 = IREND - LEND                                          
         IF (NUMPM1.LT.IST) THEN                                        
            WRITE (7,*) 'TERMINAL ERROR IN RADINA.DISCONTINUITIES       
     $ TOO CLOSE TOGETHER.  CANNOT INTEGRATE WITH HIGH ACCURACY.'       
            WRITE (7,*) 'INV,L,LEND,IREND= ',INV,L,LEND,IREND,IEN           
            STOP                                                        
            END IF                                                      
         NQ = NUMPM1/IST - 1                                            
         IREM = NUMPM1 - NQ * IST + 1                                   
C        WRITE (7,*) 'NQ,IREM=',NQ,IREM                                 
         SUM = R(LEND)**N * F(LEND) * (-0.5E0)                          
C        WRITE (7,*) 'SUM=',SUM                                         
         DO 2 J = 1,NQ                                                  
C           WRITE (7,*) 'WE ARE HERE',J                                 
            INDL = LEND + (J-1) * IST                                   
            SUM = SUM + R(INDL  )**N * F(INDL  )                        
     $                + R(INDL+1)**N * F(INDL+1) * A51                  
     $                + R(INDL+2)**N * F(INDL+2) * A52                  
     $                + R(INDL+3)**N * F(INDL+3) * A51                  
C           IF ((J/100)*100.EQ.J) WRITE (7,*) J,SUM                     
 2       CONTINUE                                                       
            INDL = LEND + NQ * IST                                      
C THE NEXT LINES ARE SPECIFIC TO IST=4.                                 
         IF (IREM.EQ.5) THEN                                            
c           write (7,*) 'irem=5'                                        
            SUM = SUM + R(INDL  )**N * F(INDL  )                        
     $                + R(INDL+1)**N * F(INDL+1) * A51                  
     $                + R(INDL+2)**N * F(INDL+2) * A52                  
     $                + R(INDL+3)**N * F(INDL+3) * A51                  
     $                + R(INDL+4)**N * F(INDL+4) * 0.5E0                
            END IF                                                      
         IF (IREM.EQ.6) THEN                                            
            SUM = SUM + R(INDL  )**N * F(INDL  ) * (0.5E0 + A6A)        
     $                + R(INDL+1)**N * F(INDL+1) * A6B                  
     $                + R(INDL+2)**N * F(INDL+2) * A6C                  
     $                + R(INDL+3)**N * F(INDL+3) * A6C                  
     $                + R(INDL+4)**N * F(INDL+4) * A6B                  
     $                + R(INDL+5)**N * F(INDL+5) * A6A                  
            END IF                                                      
         IF (IREM.EQ.7) THEN                                            
            SUM = SUM + R(INDL  )**N * F(INDL  ) * (0.5E0 + A7A)        
     $                + R(INDL+1)**N * F(INDL+1) * A7B                  
     $                + R(INDL+2)**N * F(INDL+2) * A7C                  
     $                + R(INDL+3)**N * F(INDL+3) * A7D                  
     $                + R(INDL+4)**N * F(INDL+4) * A7C                  
     $                + R(INDL+5)**N * F(INDL+5) * A7B                  
     $                + R(INDL+6)**N * F(INDL+6) * A7A                  
            END IF                                                      
         IF (IREM.EQ.8) THEN                                            
            SUM = SUM + R(INDL  )**N * F(INDL  ) * (0.5E0 + A8A)        
     $                + R(INDL+1)**N * F(INDL+1) * A8B                  
     $                + R(INDL+2)**N * F(INDL+2) * A8C                  
     $                + R(INDL+3)**N * F(INDL+3) * A8D                  
     $                + R(INDL+4)**N * F(INDL+4) * A8D                  
     $                + R(INDL+5)**N * F(INDL+5) * A8C                  
     $                + R(INDL+6)**N * F(INDL+6) * A8B                  
     $                + R(INDL+7)**N * F(INDL+7) * A8A                  
            END IF                                                      
         IF (IREM.GT.8.OR.IREM.LT.5) THEN                               
            WRITE (7,*) 'TERMINAL ERROR IN RADINA.  IREM WRONG.'        
            STOP                                                        
            END IF                                                      
         P = P + SUM * A53                                              
         LEND = IREND                                                   
 1    CONTINUE                                                          
C      ENDPT  = (23.75*(     S1       +F(NP  )*(R(NP  )**N))          
C     1        + 95.10*(F(2)*(R(2)**N)+F(NP-1)*(R(NP-1)**N))          
C     2        + 55.20*(F(3)*(R(3)**N)+F(NP-2)*(R(NP-2)**N))          
C     3        + 79.30*(F(4)*(R(4)**N)+F(NP-3)*(R(NP-3)**N))          
C     4        + 70.65*(F(5)*(R(5)**N)+F(NP-4)*(R(NP-4)**N)))/ 72.  
C      I2 = NP - 5                                                      
C      DO 100 I=6,I2                                                    
C      P = P + (R(I)**N)*F(I)                                           
C  100 CONTINUE                                                         
C      P = (P + ENDPT ) * H                                             
      P = P * H                                                         
      RETURN                                                            
      END                                                               
