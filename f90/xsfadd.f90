SUBROUTINE xsfadd
     
!XSFABD
 
!     REVISED  8/89 BY G.C./UNISYS
!          1.  THE ORDER OF COMM AND XFIAT IN /XSFA1/ ARE REVERSED IN
!              THIS ROUTINE AND IN THE FOLLOWING 7 SUBROUTINES -
!              XCLEAN, XDPH, XPOLCK, XPUNP, XPURGE, XSFA AND XSOSGN.
!              ANY INCREASE IN SIZE OF XFIAT CAN THEREFORE BE MADE
!              EASILY THROUGH OUT THESE GROUP OF ROUTINES BY JUST
!              CHANGING THE XFIAT DIMENSION HERE.
!          2.  IN THIS GROUP OF ROUTINES, THE ARRAY XFIAT IN /XSFA1/ IS
!              RENAMED TO XFIAT, NOT TO BE CONFUSED WITH THE XFIAT ARRAY
!              IN /XFIAT/
!          3.  ENTN1 MUST EQUAL ICFIAT, THE 24TH WORD OF /SYSTEM/
!              HOWEVER, XSFA AND XPURGE ROUTINES INITIALIZE ENTN1 AGAIN
!              TO ICFIAT, JUST TO BE SURE.
!          4.  THE DIMENSION OF XFIAT SHOULD BE 800 WHEN ENTN1 = 8, OR
!              1100 WHEN ENTN1 IS 11
 
 INTEGER :: almsk,apndmk,comm,cursno,entn1,entn2,entn3,  &
     entn4,flag,fnx,rmsk,rxmsk,s,scornt,sos,tapmsk, thcrmk,xfiat,zap
!WKBR COMMON /XSFA1 / MF(401),SOS(1501),COMM(20),XFIAT(1100)
 COMMON /xsfa1 / mf(401),sos(1501),comm(20),xfiat(1320)
 EQUIVALENCE            (comm (1),almsk ),(comm (2),apndmk),  &
     (comm (3),cursno),(comm (4),entn1 ),(comm (5),entn2 ),  &
     (comm (6),entn3 ),(comm (7),entn4 ),(comm (8),flag  ),  &
     (comm (9),fnx   ),(comm(10),lmsk  ),(comm(11),lxmsk ),  &
     (comm(12),macsft),(comm(13),rmsk  ),(comm(14),rxmsk ),  &
     (comm(15),s     ),(comm(16),scornt),(comm(17),tapmsk),  &
     (comm(18),thcrmk),(comm(19),zap   )
 entn1 = 11
 entn2 = 3
 entn3 = 4
 entn4 = 3
 flag  = 0
 DO  i = 1, 1320
   xfiat(i) = 0
 END DO
 tapmsk = 32768
!            TAPMSK = O 000000100000  = Z 00008000
 apndmk = 1073741824
!            APNDMK = O 010000000000  = Z 40000000
 rmsk   = 32767
!            RMSK   = O 000000077777  = Z 00007FFF
 rxmsk  = 65535
!            RXMSK  = O 000000177777  = Z 0000FFFF
 lmsk   = 1073676288
!            LMSK   = O 007777600000  = Z 3FFF0000
 lxmsk  = 2147418112
!            LXMSK  = O 017777600000  = Z 7FFF0000
 scornt = 1073708992
!            SCORNT = O 007777677700  = Z 3FFF7FC0
 zap    = 32767
!            ZAP    = O 000000077777  = Z 00007FFF
END SUBROUTINE xsfadd
