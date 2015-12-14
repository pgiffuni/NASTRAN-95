SUBROUTINE qparam
     
!     PARAM PERFORMS THE FOLLOWING OPERATIONS ON PARAMETERS--
!      1. OUT = IN1 .AND. IN2
!      2. OUT = IN1 .OR . IN2
!      3. OUT = IN1   +   IN2
!      4. OUT = IN1   -   IN2
!      5. OUT = IN1   *   IN2
!      6. OUT = IN1   /   IN2
!      7. OUT = .NOT. IN1
!      8. OUT = IN1  .IMP. IN2
!      9. STORE VALUE OF OUT IN VPS.
!     10. OUT = VALUE OF PRECISION CELL FROM /SYSTEM/
!     11. OUT = CURRENT TIME
!     12. OUT = TIME TO GO
!     13. OUT = SYSTEM(IN1) = IN2
!     14. OUT = SYSTEM(25) WITH BITS IN1 THRU IN2 TURNED ON OR OFF.
!     15. OUT = SYSTEM CELL IN1.
!     16. SAVE AND RESTORES SENSE SWITCHES
!     17. SETS SENSE SWITCHES
!     18. SAVE AND RESTORES SYSTEM CELLS
!     19. OUT = -1 IF IN1 .EQ. IN2, OUT = +1 OTHERWISE.
!     20. OUT = -1 IF IN1 .GT. IN2, OUT = +1 OTHERWISE.
!     21. OUT = -1 IF IN1 .LT. IN2, OUT = +1 OTHERWISE.
!     22. OUT = -1 IF IN1 .LE. IN2, OUT = +1 OTHERWISE.
!     23. OUT = -1 IF IN1 .GE. IN2, OUT = +1 OTHERWISE.
!     24. OUT = -1 IF IN1 .NE. IN2, OUT = +1 OTHERWISE.
!     25. UNDEFINED.
!     26. UNDEFINED.
!     27. UNDEFINED.
!     28. UNDEFINED.
!     29. UNDEFINED.
!     30. UNDEFINED.
 
 EXTERNAL        lshift,orf,andf
 INTEGER :: switch,off,orf,xorf,op,opcode,out,outtap,andf,vps, oscar
 DIMENSION       opcode(30),switch(2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / op(2),out,in1,in2
 COMMON /system/ ksystm(80)
 COMMON /oscent/ oscar(16)
 COMMON /xvps  / vps(1)
 EQUIVALENCE     (ksystm( 2),outtap),(ksystm(23),lsystm),  &
     (ksystm(55),iprec ),(ksystm(79),switch(1))
 DATA    opcode/ 4HAND ,4HOR  ,4HADD ,4HSUB ,4HMPY  &
     , 4HDIV ,4HNOT ,4HIMPL,4HNOP ,4HPREC , 4HKLOC,4HTMTO,4HSYST,4HDIAG,4HSYSR  &
     , 4HSSSR,4HSSST,4HSTSR,4HEQ  ,4HGT , 4HLT  ,4HLE  ,4HGE  ,4HNE  ,4H****  &
     , 4H****,4H****,4H****,4H****,4H**** /
 DATA    off   / 4HOFF /
 
!     BRANCH ON OPERATION CODE.
 
 DO  i = 1,30
   IF (op(1) == opcode(i)) GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90,100,  &
       110,120,130,140,150,160,170,180,190,200,  &
       210,220,230,240,250,260,270,280,290,300), i
 END DO
 GO TO 990
 
!     .AND.
 
 10 out = -1
 IF (in1 >= 0 .OR. in2 >= 0) out = +1
 GO TO 900
 
!     .OR.
 
 20 out = +1
 IF (in1 < 0 .OR . in2 < 0) out = -1
 GO TO 900
 
!     ADD
 
 30 out = in1 + in2
 GO TO 900
 
!     SUB
 
 40 out = in1 - in2
 GO TO 900
 
!     MPY
 
 50 out = in1*in2
 GO TO 900
 
!     DIV
 
 60 out = in1/in2
 GO TO 900
 
!     NOT
 
 70 out = -in1
 GO TO 900
 
!     IMPLY
 
 80 out = +1
 IF (in1 >= 0 .OR. in2 < 0) out = -1
 GO TO 900
 
!     NOP
 
 90 GO TO 900
 
!     PROVIDE PRECISION FROM /SYSTEM/.
 
 100 out = iprec
 GO TO 900
 
!     PROVIDE CURRENT TIME
 
 110 CALL klock (out)
 GO TO 900
 
!     PROVIDE TIME-TO-GO
 
 120 CALL tmtogo (out)
 GO TO 900
 
!     MODIFY SYSTEM CELL.
 
 130 out = in2
 ksystm(in1) = in2
 IF (in1 <= 0 .OR. in1 > lsystm) WRITE (outtap,135) uwm,in1
 135 FORMAT (a25,' 2317, PARAM HAS STORED OUTSIDE DEFINED RANGE OF ',  &
     'COMMON BLOCK /SYSTEM/.', /32X,'INDEX VALUE =',i20)
 GO TO 900
 
!     TURN DIAG SWITCH ON OR OFF.
 
 140 IF (in2 < in1) in2 = in1
 DO  i = in1,in2
   IF (i > 31) GO TO 142
   out = lshift(1,i-1)
   switch(1) = orf(switch(1),out)
   IF (op(2) == off) switch(1) = switch(1) - out
   CYCLE
   142 out = i - 31
   out = lshift(1,out-1)
   switch(2) = orf(switch(2),out)
   IF (op(2) == off) switch(2) = switch(2) - out
   out = out + 31
 END DO
 out = switch(1)
 IF (i > 31) out = switch(2)
 GO TO 900
 
!     RETURN VALUE OF IN1-TH WORD OF /SYSTEM/.
 
 150 out = ksystm(in1)
 GO TO 900
 
!     SAVE OR RESTORE SSWITCH WORD
 
 160 IF (in1 <  0) GO TO 165
 IF (in1 > 31) GO TO 161
 out = switch(1)
 GO TO 900
 161 CONTINUE
 out = switch(2)
 GO TO 900
 165 IF (IABS(in1) > 31) GO TO 166
 switch(1) = out
 GO TO 900
 166 switch(2) = out
 GO TO 900
 
!     TURN SSWITCH ON OR OFF
 
 170 IF (out == 0) GO TO 900
 IF (out > 0) GO TO 175
 IF (IABS(out) > 31) GO TO 171
 mask = lshift(1,IABS(out)-1)
 switch(1) = xorf(mask,orf(mask,switch(1)))
 GO TO 900
 171 CONTINUE
 out  = out + 31
 mask = lshift(1,IABS(out)-1)
 switch(2) = xorf(mask,orf(mask,switch(2)))
 out  = out - 31
 GO TO 900
 175 CONTINUE
 IF (out > 31) GO TO 176
 switch(1) = orf(lshift(1,out-1),switch(1))
 GO TO 900
 176 CONTINUE
 out = out - 31
 switch(2) = orf(lshift(1,out-1),switch(2))
 out = out + 31
 GO TO 900
 
!     SAVE OR RESTORE A CELL OF SYSTEM
 
!     SAVE
 
 180 CONTINUE
 IF (in1 < 0) GO TO 185
 out = ksystm(in1)
 GO TO 900
 
!     RESTORE
 
 185 in1 = IABS(in1)
 ksystm(in1) = out
 GO TO 900
 
!     ARITHMETIC RELATIONAL OPERATORS.
 
 190 IF (in1-in2 == 0) THEN
   GO TO   192
 END IF
 191 out = +1
 GO TO 900
 192 out = -1
 GO TO 900
 200 IF (in1-in2 > 0) THEN
   GO TO   192
 ELSE
   GO TO   191
 END IF
 210 IF (in1-in2 < 0) THEN
   GO TO   192
 ELSE
   GO TO   191
 END IF
 220 IF (in1-in2 > 0) THEN
   GO TO   191
 ELSE
   GO TO   192
 END IF
 230 IF (in1-in2 < 0) THEN
   GO TO   191
 ELSE
   GO TO   192
 END IF
 240 IF (in1-in2 == 0) THEN
   GO TO   191
 ELSE
   GO TO   192
 END IF
 
!     UNDEFINED.
 
 250 GO TO 900
 
!     UNDEFINED.
 
 260 GO TO 900
 
!     UNDEFINED.
 
 270 GO TO 900
 
!     UNDEFINED.
 
 280 GO TO 900
 
!     UNDEFINED.
 
 290 GO TO 900
 
!     UNDEFINED.
 
 300 GO TO 900
 
!     SAVE OUT IN THE VPS.
 
 900 i = andf(oscar(16),65535)
 vps(i) = out
 RETURN
 
!     OPERATION CODE NOT DEFINED-- WRITE MESSAGE.
 
 990 WRITE  (outtap,998) ufm,op(1),op(2)
 998 FORMAT (a23,' 2024, OPERATION CODE ',2A4,' NOT DEFINED FOR ',  &
     'MODULE PARAM.')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE qparam
