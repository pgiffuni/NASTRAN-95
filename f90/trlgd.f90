SUBROUTINE trlgd(fct,fco,ap,as,ad,ah,  &
        ppo,pso,pdo,pdt,pht,iflag1,scr1,iflag      )
     
!     THE PURPOSE OF THIS SUBROUTINE IS TO COMPUTE LOAD FACTORS
!         BOTH AT APPLIED TIMES(T) AND OUTPUT TIMES(O).
 
!     INPUTS (6)
 
!       FCT --MATRIX OF TIME FUNCTIONS--ALL TIMES
!       FCO --MATRIX OF TIME FUNCTIONS--OUTPUT TIMES
!       IFLAG =-1   IMPLIES FCT = FCO AND ONLY FCT EXISTS
!         NOTE THAT ALSO IMPLIES PDO = PDT
!       AP,AS,AD,AH  ARE TRANSFORMATION MATRICIES TO P,S,D,AND H SET RES
 
!     OUTPUTS (5)
!       PPO,PSO,PDO  LOADS AT OUTPUT TIMES(ANY MAY NOT EXIST)
!       PDT,PHT      LOADS AT ALL TIMES (ANY MAY NOT EXIST)
 
!     SCR1          SCRATCH FILE FOR MPYAD
 
!     IFLAG1 =-1  IMPLIES THAT AP = AD
 
!     FCT MAY BE FCO, IN CASE OF EQUALITY THE T FILES WILL EXIST
 
 
 INTEGER, INTENT(IN OUT)                  :: fct
 INTEGER, INTENT(IN OUT)                  :: fco
 INTEGER, INTENT(IN OUT)                  :: ap
 INTEGER, INTENT(IN)                      :: as
 INTEGER, INTENT(IN OUT)                  :: ad
 INTEGER, INTENT(IN OUT)                  :: ah
 INTEGER, INTENT(IN)                      :: ppo
 INTEGER, INTENT(IN)                      :: pso
 INTEGER, INTENT(IN)                      :: pdo
 INTEGER, INTENT(IN)                      :: pdt
 INTEGER, INTENT(IN)                      :: pht
 INTEGER, INTENT(IN OUT)                  :: iflag1
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN OUT)                  :: iflag
 INTEGER :: mcb(7) , trnsp,SIGN,prec
 
 COMMON /system/iskip(54),iprec
 
 SIGN = +1
 trnsp = 0
 prec = iprec
 
!     FORM PPO
 
 mcb(1) = ppo
 CALL rdtrl(mcb)
 IF(mcb(1) <= 0) GO TO 10
 CALL ssg2b(ap,fco,0,ppo,trnsp,prec,SIGN,scr1)
 mcb(1) = ppo
 CALL rdtrl(mcb)
 mcb(2) = mcb(2) -1
 CALL wrttrl(mcb)
 10 CONTINUE
 
!     FORM   PSO
 
 mcb(1) = pso
 CALL rdtrl(mcb)
 IF (mcb(1) <= 0) GO TO 20
 mcb(1) = as
 CALL rdtrl(mcb)
 IF (mcb(2) <= 0)  GO TO 20
 CALL ssg2b(as,fco,0,pso,trnsp,prec,SIGN,scr1)
 mcb(1) = pso
 CALL rdtrl(mcb)
 mcb(2) = mcb(2) -1
 CALL wrttrl(mcb)
 20 CONTINUE
 
!     BUILD PDO
 
 IF(iflag1 == -1) GO TO 30
 mcb(1) = pdo
 CALL rdtrl(mcb)
 IF (mcb(1) <= 0) GO TO 30
 CALL ssg2b(ad,fco,0,pdo,trnsp,prec,SIGN,scr1)
 mcb(1) = pdo
 CALL rdtrl(mcb)
 mcb(2) = mcb(2) -1
 CALL wrttrl(mcb)
 30 CONTINUE
 
!     BUILD PDT
 
 mcb(1) = pdt
 CALL rdtrl(mcb)
 IF (mcb(1) <= 0) GO TO 40
 CALL ssg2b(ad,fct,0,pdt,trnsp,prec,SIGN,scr1)
 40 CONTINUE
 
!     BUILD PHT
 
 mcb(1) = pht
 CALL rdtrl(mcb)
 IF(mcb(1) <= 0) GO TO 50
 CALL ssg2b(ah,fct,0,pht,trnsp,prec,SIGN,scr1)
 50 CONTINUE
 RETURN
END SUBROUTINE trlgd
