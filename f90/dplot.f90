SUBROUTINE dplot
     
 IMPLICIT INTEGER (a-z)
 INTEGER :: tit(32),NAME(2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /plothd/ iused
 COMMON /system/ bufsiz ,nout
 COMMON /zzzzzz/ x(1)
 COMMON /BLANK / ngp    ,lsil   ,nsets  ,pltflg ,pltnum ,ngpset  ,  &
     nodef  ,skp1(3),pltpar ,gpsets ,elsets ,casecc  ,  &
     bgpdt  ,eqexin ,sil    ,pdef1  ,pdef2  ,s2      ,  &
     plotx  ,setd   ,ecpt   ,oes1   ,scr1   ,scr2    , scr3   ,scr4
 
!     NOTE THAT NSETS IS DMAP PARAMETER JUMPPLOT
!     IUSED IS USED IN PLOT AND HDPLOT
 
 DATA    inprew, rew / 0,1   /,  &
     tit   / 12*1H ,4HMESS,4HAGES,4H fro,4HM th,4HE pl,4HOT m,  &
     4HODUL,1HE   ,12*1H    /
 DATA    NAME  / 4HDPLO,4HT   /
 
!     FILE NAMES FOR UNDEFORMED PLOTS MAY BE
!     108  = USET (GPTLBL - SPC DEGREES OF FREEDOM)
!     109  = ECT  (ELELBL - PROPERTY IDS)
!     110  = ECPT
!          = EPT (UNDEFORMED PLOT ONLY, DMAP NUMBER 25 OR LESS)
!            EPT IS NEEDED FOR PSHELL CARDS IN ORDER TO PICK UP ANY
!            OFFSET FOR CTRIA3 AND CQUAD4 (IN COMECT)
 
 pltpar = 101
 gpsets = 102
 elsets = 103
 casecc = 104
 bgpdt  = 105
 eqexin = 106
 sil    = 107
 pdef1  = 108
 pdef2  = 109
 ecpt   = 110
 oes1   = 111
 oes1l  = 112
 onrgy1 = 113
 plotx  = 201
 scr1   = 301
 scr2   = 302
 scr3   = 303
 scr4   = 304
 nodef  = 0
 IF (ngp <= 0 .OR. lsil <= 0) GO TO 80
 CALL totape (2,x(1))
 
!     OUTPUT THE TITLE FOR MESSAGE FILE
!     THE LAST BUFFER IS BUFSIZ+1 FOR SUBROUTINE ELELBL
 
 buf = korsz(x) - 4*bufsiz
 IF (buf-4*bufsiz < 10) GO TO 85
 IF (nsets <= 0) GO TO 60
 CALL gopen (plotx,x(buf),rew)
 
!     COMMENTS FROM G.CHAN/UNISYS       11/90
!     NEXT 2 LINES ADD TIT HEADING TO THE 4TH LINE OF NASTRAN HEADERS
!     WHEN THE PLOTX FILE IS READ AND PRINTED BY PRTMSG MODULE.
!     THIS SHORTCUT TECHNIQUE IS NO WHERE DISCUSSED IN THE USER'S NOR
!     PROGRAMMER'S MAUNALS
 
 CALL WRITE (plotx,-4,1,0)
 CALL WRITE (plotx,tit,32,0)
 
!     READ THE SETID-S FROM -GPSETS- FILE.  SET NEGATIVE SETID-S THAT
!     HAVE NO ASSOCIATED GRIDS.  FIND FIRST DEFINED SET OR EXIT IF NONE
 
 buf = buf - bufsiz
 CALL gopen (gpsets,x(buf),inprew)
 CALL fread (gpsets,x,nsets,1)
 setd = 0
 x(nsets+1) = 1
 
 DO  i = 1,nsets
   CALL READ (*30,*60,gpsets,x(nsets+2),1,1,i1)
   IF (x(nsets+2) > 0) GO TO 40
   30 WRITE  (nout,31) uwm,x(nsets+1)
   31 FORMAT (a25,' 697, SET',i9,  &
       ' NOT DEFINED.  FIRST SET DEFINED WILL BE USED.')
   x(i) = -x(i)
   CYCLE
   40 IF (setd == 0) setd = i
 END DO
 CALL CLOSE (gpsets,rew)
 IF (setd /= 0) GO TO 70
 60 WRITE  (nout,61) ufm
 61 FORMAT (a23,' 698, NO SETS DEFINED FOR PLOTS')
 CALL mesage (-61,0,0)
 
!     PROCESS PLOT REQUESTS
 
 70 CALL gopen (pltpar,x(buf),inprew)
 i1 = 1
 i2 = i1  + nsets
 buf= buf - bufsiz
 CALL param (x(i1),x(i2),buf-nsets)
 CALL CLOSE (pltpar,rew)
 
!     SET JUMPPLOT NEGATIVE IF NO FUTHER REQUESTS
 
 IF (pltflg >= 0 .AND. nodef == 0) nsets = -1
 CALL clstab (plotx,rew)
 CALL CLOSE  (gpsets,rew)
 pltflg = -1
 80 RETURN
 
!     INSUFFICIENT CORE
 
 85 CALL mesage (-8,buf,NAME)
 nsets = -1
 pltflg= -1
 GO TO 80
END SUBROUTINE dplot
