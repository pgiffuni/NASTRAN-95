SUBROUTINE ifp1s (list,istor,nlist)
     
    !     THIS ROUTINE FINDS ANY OVERLAPPING INTERVALS IN A SET LIST.
    !     IT WILL ALSO CHECK SINGLES
 
    INTEGER, INTENT(IN OUT)                  :: list(1)
    INTEGER, INTENT(OUT)                     :: istor(1)
    INTEGER, INTENT(IN OUT)                  :: nlist
    INTEGER :: otpe
 
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm,uwm
    COMMON /system/ sysbuf,otpe,inx(6),nlpp,inx1(2),line
 
    ipair = 0
    DO  i = 1,nlist
        IF (list(i) > 0) CYCLE
10      IF (ipair   /= 0) GO TO 40
20      l = 2*ipair + 1
        istor(l  ) = list(i-1)
        istor(l+1) = list(i  )
        ipair = ipair + 1
30      list(i  ) = 0
        list(i-1) = 0
        CYCLE
   
        !     PAIR FOUND - CHECK FOR OVERLAP
   
40      l  = 1
        in = IABS(list(i-1))
        iout = IABS(list(i  ))
        k  = 2*l - 1
50      IF (in >= istor(k)   .AND. in <= IABS(istor(k+1))  ) GO TO 60
        IF (iout >= istor(k) .AND. iout <= IABS(istor(k+1))) GO TO 60
        l  = l + 1
        IF (l <= ipair) GO TO 50
   
        !     STORE NEW PAIR IN LIST
   
        GO TO 20
   
        !     ERROR IN INTERVAL
   
60      list(i-1) = MIN0(in,istor(k))
        list(i  ) =-MAX0(iout,IABS(istor(k+1)))
        IF (list(i-1) == istor(k) .AND. list(i) == istor(k+1)) GO TO 30
        ix = IABS(istor(k+1))
        WRITE  (otpe,70) uwm,in,iout,istor(k),ix
70      FORMAT (a25,' 621, INTERVAL',i8,' THRU',i8,' OVERLAPS INTERVAL',  &
            i8,' THRU', i8,'. THE MAXIMUM INTERVAL WILL BE USED.')
        line = line + 3
        IF (line >= nlpp) CALL page
   
        !     REMOVE PAIR L FROM ISTOR
   
80      IF (l >= ipair) GO TO 90
        m = 2*l + 1
        k = 2*l - 1
        istor(k  ) = istor(m  )
        istor(k+1) = istor(m+1)
        l = l + 1
        GO TO 80
90      ipair = ipair - 1
        GO TO 10
    END DO
 
    !     ALL PAIRS PROCESSED - TRY SINGLES
 
    ising = 0
    m  = 2*ipair
    loop180:  DO  i = 1,nlist
        in = list(i)
        IF (ipair   == 0) GO TO 140
        IF (list(i) == 0) CYCLE loop180
   
        !     CHECK EACH PAIR
   
        l = 1
110     k = 2*l - 1
        IF (in >= istor(k) .AND. in <= IABS(istor(k+1))) GO TO 120
        l = l + 1
        IF (l-ipair > 0) THEN
            GO TO   140
        ELSE
            GO TO   110
        END IF
   
        !     ERROR -- PAIR CONTAINS SINGLE
   
120     in1 = IABS(istor(k+1))
        WRITE  (otpe,130) uwm,in,istor(k),in1
130     FORMAT (a25,' 619, SET MEMBER',i8,' BELONGS TO',i8,' THRU',i8)
        line = line + 3
        IF (line >= nlpp) CALL page
        CYCLE loop180
   
        !     CHECK FOR DUPLICATE SINGLES
   
140     IF (ising == 0) GO TO 170
        DO  k = 1,ising
            l = 2*ipair + k
            IF (in /= istor(l)) CYCLE
            WRITE  (otpe,150) uwm,in
150         FORMAT (a25,' 620, SET MEMBER',i8,' IS DUPLICATED IN SET LIST.')
            line = line + 3
            IF (line >= nlpp) CALL page
            CYCLE loop180
        END DO
170     m = m + 1
        ising = ising + 1
        istor(m) = in
    END DO loop180
 
    !     COPY GOOD STUFF INTO LIST
 
    DO  i = 1,m
        list(i) = istor(i)
    END DO
    nlist = m
 
    !     SORT LIST
 
    n1 = m - 1
    DO  i = 1,n1
        n2 = i + 1
        DO  k = n2,m
            IF (IABS(list(i))-IABS(list(k)) > 0) THEN
                GO TO   210
            ELSE
                GO TO   220
            END IF
     
            !     SWITCH
     
210         in = list(i)
            list(i) = list(k)
            list(k) = in
220     CONTINUE
        END DO
    END DO
 
    RETURN
END SUBROUTINE ifp1s
