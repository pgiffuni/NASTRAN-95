SUBROUTINE aspro (dmap,var,nvar,obits,sol)
     
    !     THIS CODE  PERFORMS THE  ROUTINE PROCESSING OF THE  DMAP ALTERS
    !     FOR ASDMAP.  KEY TABLES ARE-
 
    !      DMAP  -  RAW  18 WORD PER CARD BCD DATA ON INPUT, VARIABLE
    !               CHARACTERS ARE ADDED AND  FIELDS AND CARDS ARE DELETED
    !               DEPENDING ON USER INPUT IN  VAR(IABLE) AND OPTION FLAGS.
 
    !      VAR      CONTROL DATA AND USER INPUT DATA, 3 WORDS, BCD + DATA
 
    !      PTBS     POSITIONS-TO-BE-SET TABLE, CONTENTS-PER ENTRY
 
    !                    1   CARD NUMBER IN DMAP
    !                    2   FIRST CHARACTER OF MODIFIED FIELD
    !                    3   FIRST CHARACTER FOR ADDED VARIABLE
    !                    4   NUMBER OF VARIABLE CHARACTERS
    !                    5   KEY OF VARIABLE TO BE INSERTED
    !                    6   MATRIX OPTION FLAG , 1= K, 2=M, 4=P  ETC
    !                    7   OUTPUT CONTROL FLAG, AVOIDS SAME DATA BLOCK
    !                        OUTPUT FROM TWO MODULES
 
    !      OCT      OPTIONAL CARDS TABLE - EACH ENTRY =
    !                   DMAP CARD NO. , DELETE BITS ,  KEEP BITS
 
    !      OBITS -  BITS ARE ON FOR REQUIRED MATRICES  =  SUM OF NUMBERS
    !                   K=1 , M=2 , P=4 , PA=8 , B=16 , K4=32
 
 
    INTEGER, INTENT(OUT)                     :: dmap(18,1)
    INTEGER, INTENT(IN)                      :: var(3,200)
    INTEGER, INTENT(IN)                      :: nvar
    INTEGER, INTENT(IN OUT)                  :: obits
    INTEGER, INTENT(IN OUT)                  :: sol
    EXTERNAL        andf
    LOGICAL :: rmv,rmvall
    INTEGER :: andf,dbs(2,50), flag,ii(2),NAME(4),  &
        oct(3,50),ptbs(7,200), vword, alter,BLANK,ast,slas,oball,rfmask(40)
    CHARACTER (LEN=25) :: sfm
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm,uwm,uim,sfm
    COMMON /asdbd / irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs, iph,nph,idat(1)
    COMMON /zzzzzz/ sbd(2)
    COMMON /system/ idum1,iout,nogo
    EQUIVALENCE     (ndbs,sbd(1)),(dbs(1,1),sbd(2))
    DATA    alter / 4HALTE /,      BLANK  / 4H     /
    DATA    ast   / 4H*    /,      slas   / 4H/    /
    DATA    rfmask/ 65536,131072,262144,0,0,0,0,524288,1048576,31*0 /
    DATA    oball / 63     /
 
    rmvall = .true.
    nxdel  = 0
    nold   = 0
 
    !     DELETE CARDS USING OCT TABLE.
 
    IF (noct == 0) GO TO 45
    m = ioct - 1
    DO  j = 1,noct
        DO  i = 1,3
            m = m + 1
            oct(i,j) = idat(m)
        END DO
    END DO
    DO  i = 1,noct
        icd = oct(1,i)
        IF (oct(3,i) == 0) GO TO 20
        IF (andf(oct(3,i),obits)       == 0) GO TO 35
20      IF (andf(oct(2,i),rfmask(sol)) == 0) CYCLE
35      dmap(1,icd) = -1
    END DO
45  IF (nptbs == 0) GO TO 2000
    m = iptbs - 1
    DO  j = 1,nptbs
        DO  i = 1,7
            m = m + 1
            ptbs(i,j) = idat(m)
        END DO
    END DO
    DO  i = 1,nptbs
        icd = ptbs(1,i)
        IF (dmap(1,icd) == -1) CYCLE
        IF (icd == 0) CYCLE
        rmv = .false.
   
        !     CHECK IF  OPTION IS ON
   
        kopt = ptbs(6,i)
        IF (andf(kopt,obits) == 0) rmv = .true.
        IF (andf(kopt,oball) == 0) rmv = .false.
        IF (andf(kopt,rfmask(sol)) /= 0) rmv = .true.
        nchar = ptbs(4,i)
        nc    = 0
        flag  = 0
        vword = BLANK
        icol  = ptbs(3,i)
   
        !     FIND  VARIABLE  IF  REQUIRED
   
        IF (rmv) GO TO 300
        key = ptbs(5,i)
        i3  = nvar/3
        DO  j = 1,i3
            IF (var(1,j) == key) GO TO 70
            IF (key == j .AND. var(1,j) == alter) GO TO 450
        END DO
   
        !     VARIABLE HAS NOT BEEN SET
   
        vword = BLANK
        rmv   = .true.
        GO TO 300
   
        !     VARIABLE  IS FOUND , IT IS IN VAR(2,J) AND/OR VAR(3,J)
   
70      vword   = var(2,j)
        NAME(1) = var(2,j)
        NAME(2) = var(3,j)
   
        !     TEST FOR REAL OR INTEGER
   
        IF (vword ==  0) GO TO 300
        IF (vword == -1) GO TO 74
        IF (vword == -2) GO TO 3010
   
        !     WORD IS REAL (TEMPORARY ERROR)
   
        GO TO 75
   
        !     WORD IS AN INTEGER
   
74      NAME(1) = NAME (2)
        NAME(2) = 0
        flag    = 1
75      IF (ptbs(7,i) == 0) GO TO 500
        nc    = ptbs(3,i) - ptbs(2,i)
        ii(1) = BLANK
        ii(2) = BLANK
        IF (nc > 0) GO TO 80
        IF (nc < 0) GO TO 3010
        ii(1) = NAME(1)
        ii(2) = NAME(2)
        GO TO 100
   
        !     CONSTRUCT WHOLE DATA BLOCK NAME
   
80      CALL pull (dmap(1,icd),ii,ptbs(2,i),nc,0)
        CALL push (NAME,ii,nc+1,nchar,flag)
   
        !     CHECK OUTPUT DATA BLOCKS AGAINST PREVIOUS OUTPUTS
   
100     IF (ndbs == 0) GO TO 142
   
        DO  l = 1, ndbs
            IF (ii(1) == dbs(1,l) .AND. ii(2) == dbs(2,l)) GO TO 150
        END DO
142     IF (ptbs(7,i) > 0) GO TO 200
   
        !     VARIABLE IS OK , ADD NAME TO LIST
   
        ndbs = ndbs + 1
        dbs(1,ndbs) = ii(1)
        dbs(2,ndbs) = ii(2)
        GO TO 500
150     IF (ptbs(7,i) > 0) GO TO 500
   
        !     DATA BLOCK IS OUTPUT, REMOVE IF ALLREADY DEFINED.
   
200     rmv  =.true.
   
        !     REMOVE WHOLE  NAME HERE  , CHECK FOR PARAMETER
   
300     ii(1)   = 0
        NAME(1) = BLANK
        NAME(2) = BLANK
        NAME(3) = BLANK
        NAME(4) = BLANK
        flag    = 0
        CALL pull (dmap(1,icd),ii,ptbs(2,i),1,0)
        IF (ii(1) == slas .OR. ii(1) == ast) GO TO 500
        icol  = ptbs(2,i)
        nchar = nchar + ptbs(3,i) - ptbs(2,i)
        GO TO 500
   
        !     CHECK IF ALTER CARD, OUTPUT AS BCD AND TWO INTEGERS
   
450     dmap(1,icd) = alter
        dmap(2,icd) = var(2,j)
        dmap(3,icd) = var(3,j)
        rmvall = .false.
        nxdel  = 0
        IF (var(2,j) == 0) rmvall = .true.
        GO TO 910
   
        !     ADD VARIABLES TO BCD DMAP
   
500     CALL push (NAME,dmap(1,icd),icol,nchar,flag)
   
        IF (.NOT.rmv) rmvall = .false.
   
        !     IF ALL VARIABLES ARE REMOVED FROM ONE CARD, DELETE THE CARD
   
        nnew = ptbs(1,i+1)
        IF (icd == nnew) CYCLE
   
        !     NEXT  COMMAND GOES TO NEW CARD,  CHECK IF CONTINUATION
   
905     CALL pull (dmap(1,icd+1),ii,1,4,0)
   
        IF (ii(1) /= BLANK) GO TO 910
   
        !     CONTINUATION FOUND
   
        nxdel = nxdel + 1
        IF (nnew == icd+1) CYCLE
        icd = icd+1
        GO TO 905
   
        !     FINISHED WITH  LOGICAL CARD
   
910     IF (.NOT.rmvall) GO TO 940
        dmap(1,icd) = -1
        IF (nxdel <= 0) CYCLE
        DO  l = 1,nxdel
            j = icd-l
            dmap(1,j) = -1
        END DO
940     rmvall = .true.
   
        !     END OF LOOP ON VARIABLE CHARACTERS
   
        nxdel = 0
   
    END DO
 
    !     PROCESS CARDS TO BE DELETED FROM SEQUENCE
 
2000 ikeep = 0
    DO  icd = 1,nrdm
   
        IF (dmap(1,icd) == -1) CYCLE
   
        !     KEEP CARD
   
        ikeep = ikeep + 1
        DO  j = 1,18
            dmap(j,ikeep) = dmap(j,icd)
        END DO
    END DO
    nrdm = ikeep
    RETURN
3010 WRITE  (iout,3020) sfm,dmap(1,icd)
3020 FORMAT (a25,' 6010, ILLEGAL VARIABLE TO BE SET IN DMAP STATEMENT',  &
        3X,a4)
 
    nogo = 1

    RETURN
END SUBROUTINE aspro
