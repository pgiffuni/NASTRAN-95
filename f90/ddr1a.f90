SUBROUTINE ddr1a (pd,k2dd,b2dd,mdd,vud,pad,frl,frqset,scr1,scr2,  &
        scr3,scr4,itype,scr5)
     
!     ROUTINE TO COMPUTE PAD FROM MODAL APPROXIMATION TO SYSTEM
 
 
 INTEGER, INTENT(IN)                      :: pd
 INTEGER, INTENT(IN)                      :: k2dd
 INTEGER, INTENT(IN)                      :: b2dd
 INTEGER, INTENT(IN OUT)                  :: mdd
 INTEGER, INTENT(IN OUT)                  :: vud
 INTEGER, INTENT(IN)                      :: pad
 INTEGER, INTENT(IN)                      :: frl
 INTEGER, INTENT(IN OUT)                  :: frqset
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER, INTENT(IN OUT)                  :: scr4
 INTEGER, INTENT(IN OUT)                  :: itype
 INTEGER, INTENT(IN)                      :: scr5
 INTEGER :: sysbuf,FILE,mcb(7),mcb1(7),mcb2(7),sr1,sr3,freq,  NAME(2)
 DIMENSION       iblk(60),b(2),imcb(21),ifile(3)
 COMMON /system/ sysbuf
 COMMON /condas/ consts(5)
 COMMON /zzzzzz/ core(1)
 COMMON /zntpkx/ a(4),ii,ieol,IEOR
 EQUIVALENCE     (mcb(1),imcb(1)),(mcb2(1),imcb(8)),  &
     (mcb1(1),imcb(15)),(consts(2),twopi)
 DATA    NAME  / 4HDDR1,4HA   /
 DATA    freq  / 4HFREQ       /
 
!     INITIALIZE + FIND OUT WHAT EXISTS
 
 sr1  = scr1
 sr3  = scr3
 ibuf = korsz(core) - sysbuf + 1
 nok2dd = 1
 mcb(1) = k2dd
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) nok2dd = -1
 nob2dd = 1
 mcb(1) = b2dd
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) nob2dd = -1
 mcb(1) = pd
 CALL rdtrl (mcb)
 
!     IS THIS FREQRES OR TRANSIENT
 
 IF (itype /= freq) GO TO 160
 
!     BRING IN FRL
 
 FILE = frl
 CALL OPEN  (*280,frl,core(ibuf),0)
 CALL fread (frl,0,-2,0)
 CALL READ  (*300,*20,frl,core(1),ibuf,0,nfreq)
 GO TO 310
 20 CALL CLOSE (frl,1)
 nload = mcb(2)/nfreq
 it = 3
 
!     BUILD  ACCELERATION AND VELOCITY IF NEEDED
 
 30 CALL gopen (vud,core(ibuf),0)
 
!     PUT  ACCELERATION VECTOR ON SCR1
 
 nz = ibuf - sysbuf
 CALL gopen (scr1,core(nz),1)
 CALL makmcb (mcb1,scr1,mcb(3),2,it)
 IF (nob2dd < 0) GO TO 40
 
!     PUT VELOCITY VECTOR ON SCR2
 
 nz = nz - sysbuf
 CALL gopen (scr2,core(nz),1)
 CALL makmcb (mcb2,scr2,mcb(3),2,it)
 40 IF (itype /= freq) GO TO 170
 
!     COMPUTE  VECTORS
 
 DO  i = 1,nfreq
   core(i) = core(i)*twopi
 END DO
 DO  j = 1,nload
   DO   i = 1,nfreq
     w  = core(i)
     w2 = -w*w
     CALL bldpk (3,3,scr1,iblk(1),1)
     IF (nob2dd < 0) GO TO 50
     CALL bldpk (3,3,scr2,iblk(21),1)
     50 CALL intpk (*80,vud,0,3,0)
     60 IF (ieol == 0) THEN
       GO TO    70
     ELSE
       GO TO    80
     END IF
     70 CALL  zntpki
     b(1) = w2*a(1)
     b(2) = w2*a(2)
     CALL bldpki (b(1),ii,scr1,iblk(1))
     IF (nob2dd < 0) GO TO 60
     b(1)  =-w*a(2)
     b(2)  = w*a(1)
     CALL bldpki (b(1),ii,scr2,iblk(21))
     GO TO 60
     
!     END OF COLUMN
     
     80 CALL bldpkn (scr1,iblk(1),mcb1(1))
     IF (nob2dd < 0) CYCLE
     CALL bldpkn (scr2,iblk(21),mcb2(1))
   END DO
 END DO
 110 CALL CLOSE  (scr1,1)
 CALL CLOSE  (vud,1)
 CALL wrttrl (mcb1(1))
 IF (nob2dd < 0) GO TO 120
 CALL CLOSE  (scr2,1)
 CALL wrttrl (mcb2(1))
 
!     MULTIPLY OUT
 
 120 IF (nob2dd < 0 .AND. nok2dd < 0) sr3 = pad
 CALL ssg2b (mdd,scr1,pd,sr3,0,1,0,scr4)
 IF (nok2dd < 0) GO TO 130
 
!     MULTIPLY  IN K2DD
 
 IF (nob2dd < 0) sr1 = pad
 CALL ssg2b (k2dd,scr5,sr3,sr1,0,1,0,scr4)
 GO TO 140
 
!     NO  K2DD
 
 130 sr1 = sr3
 
!     MULTIPLY IN B2DD
 
 140 IF (nob2dd < 0) GO TO 150
 CALL ssg2b (b2dd,scr2,sr1,pad,0,1,0,scr4)
 150 RETURN
 
!     TRANSIENT ANALYSIS
 
 160 nload = mcb(2)
 
!     PUT DISPLACEMENT ON SCR5,VELOCITY ON SCR2,ACCELERATION SCR1
 
 it = 1
 
!     PUT HEADERS ON FILES
 
 GO TO 30
 
!     PUT DISPLACEMENT ON SCR5
 
 170 FILE = scr5
 nz   = nz - sysbuf
 CALL gopen (scr5,core(nz),1)
 mcb(1) = scr5
 mcb(2) = 0
 mcb(4) = 2
 mcb(5) = 1
 ifile(1) = scr5
 ifile(2) = scr2
 ifile(3) = scr1
 mcb(6) = 0
 DO  kk = 1,nload
   CALL bldpk (1,1,scr5,iblk,1)
   IF (nob2dd < 0) GO TO 190
   CALL bldpk (1,1,scr2,iblk(21),1)
   190 CALL bldpk (1,1,scr1,iblk(41),1)
   DO  i = 1,3
     l =  i*7 - 6
     k = 20*i - 19
     FILE = ifile(i)
     
!     FWDREC OVER  UNNEEDED STUFF
     
     IF (i == 2 .AND. nob2dd < 0) GO TO 250
     CALL intpk (*240,vud,0,1,0)
     220 IF (ieol == 0) THEN
       GO TO   230
     ELSE
       GO TO   240
     END IF
     230 CALL zntpki
     CALL bldpki (a,ii,FILE,iblk(k))
     GO TO 220
     
!     END COLUMN
     
     240 CALL bldpkn (FILE,iblk(k),imcb(l))
     CYCLE
     250 CALL skprec (vud,1)
   END DO
 END DO
 
!     FINISH OFF
 
 CALL CLOSE  (scr5,1)
 CALL wrttrl (mcb)
 GO TO 110
 
!     ERROR MESAGES
 
 280 ip1 = -1
 290 CALL mesage (ip1,FILE,NAME)
 300 ip1 = -2
 GO TO 290
 310 ip1 = -8
 GO TO 290
END SUBROUTINE ddr1a
