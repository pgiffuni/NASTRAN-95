SUBROUTINE smcspl ( mcol, zi )
     
! SMCSPL RETRIEVES COLUMN "MCOL" FROM THE SPILL FILE.
! IF THIS COLUMN IS THE PIVOT COLUMN AND NO SPACE IS AVAILABLE, THEN
! IN-MEMORY DATA WILL BE WRITTEN TO THE SPILL FILE TO MAKE SPACE
! AVAILABLE FOR THE COLUMN DATA.  IF THE COLUMN IS NOT THE PIVOT
! COLUMN, THEN THE DATA IS READ INTO THE SPILL ARRAY IN OPEN CORE.
! WHEN A NEW PIVOT COLUMN IS DETERMINED, AN ANALYSIS IS DONE TO
! FREE UP MEMORY OF COLUMN DATA NO LONGER NEEDED.
 
 
 INTEGER, INTENT(IN)                      :: mcol
 INTEGER, INTENT(IN OUT)                  :: zi(10)
 INTEGER :: itemp(4)
 INCLUDE  'SMCOMX.COM'
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON / xmssg / ufm, uwm, uim, sfm
 
 mdir = (mcol-1)*4 + 1
 
! POSITION SPILL FILE TO CORRECT RECORD FOR THIS COLUMN AND READ DATA
 
 CALL filpos ( iscr1, zi( mdir+3 ) )
 CALL READ ( *7001, *7002, iscr1, zi( ispill ), 4, 0, 4 )
 mm2    = zi( ispill+1 )
 mterms = zi( ispill+3 )
 mwords = mm2 + mterms * ivwrds
!      PRINT *,' SMCSPL,MM2,MTERMS,MWORDS=',MM2,MTERMS,MWORDS
 CALL READ ( *7001, *7002, iscr1, zi( ispill+4 ),mwords,1,mwords )
 
! CHECK IF WE HAVE ALREADY SCANNED FOR UNNEEDED COLUMNS FOR THIS PIVOT
 
!      PRINT *,' SMCSPL,MEMLCK,KCOL=',MEMLCK,KCOL
 IF ( memlck == kcol ) GO TO 300
 memlck = kcol
 
! SCAN FOR COLUMNS NO LONGER NEEDED AND ADD THEM TO THE FREE CHAIN
 
 ifirst = 0
 DO  i = memcol1, kcol
   idir = (i-1)*4 + 1
   
! CHECK TO SEE IF THIS COLUMN NEEDED BY ANY SUBSEQUENT COLUMNS TO FOLLOW
   
   IF ( zi( idir + 2 ) >= kcol ) GO TO 199
   
! DATA NO LONGER NEEDED, IS DATA IN MEMORY IF SO FREE THE SPACE TO THE
! FREE CHAIN
   
   IF ( zi( idir     ) == 0    ) CYCLE
   
! DATA IS IN MEMORY, RETURN SPACE TO FREE CHAIN
! FIRST, CHECK IF A FREE CHAIN EXISTS
   
   IF ( memfre /= 0 ) GO TO 100
   
! FREE CHAIN DOES NOT EXISTS, MAKE THIS SPACE THE FREE CHAIN
   
   iidx   = zi( idir )
   memfre = iidx
   memlas = iidx
   zi( iidx )   = 0
   zi( iidx+1 ) = 0
   zi( idir   ) = 0
   CYCLE
   
! FREE CHAIN EXISTS, ADD THIS SPACE TO IT
   
   100   lidx   = memlas
   iidx   = zi( idir )
   memlas = iidx
   zi( lidx+1 ) = memlas
   zi( iidx   ) = lidx
   zi( iidx+1 ) = 0
   zi( idir   ) = 0
   CYCLE
   199   IF ( ifirst == 0 ) ifirst = i
 END DO
 memcol1 = ifirst
 
! CHECK IF THE FREE CHAIN IS EMPTY
 
 300   IF ( memfre == 0 ) GO TO 1000
 
! LOOP THROUGH FREE CHAIN TO FIND BLOCK LARGE ENOUGH FOR DATA
 
 iidx = memfre
 400   CONTINUE
 IF ( zi( iidx+2 ) >= (mwords+4) ) GO TO 500
 iidx = zi( iidx+1 )
 IF ( iidx /= 0 ) GO TO 400
 
! FREE CHAIN EXHAUSTED WITHOUT LARGE ENOUGH BLOCK, MUST CREATE SPACE
 
!      PRINT *,' SMCSPL GOING TO 1000 FROM 400'
 GO TO 1000
 
! SPACE FOUND, USE THIS FOR THE COLUMN DATA READ FROM THE SPILL FILE.
! RECONNECT FREE CHAIN WITHOUT THIS SPACE
 
 500   zi ( mdir ) = iidx
 iprev  = zi( iidx   )
 inext  = zi( iidx+1 )
!      PRINT *,' SMCSPL,AFTER 500,IPREV,INEXT=',IPREV,INEXT
 IF ( iprev /= 0 ) GO TO 510
 IF ( inext == 0 ) GO TO 505
 zi( inext ) = 0
 memfre      = inext
 GO TO 530
 505   memfre = 0
 GO TO 530
 510   IF ( inext == 0 ) GO TO 520
!      PRINT *,' SMCSPL,AFTER 510,INEXT,IPREV=',INEXT,IPREV
 zi( iprev+1 ) = inext
 zi( inext   ) = iprev
 GO TO 530
 520   zi( iprev+1 ) = 0
 memlas        = iprev
 
! MOVE DATA TO IN MEMORY LOCATION
 
 530   CONTINUE
 zi( mdir    ) = iidx
 zi( mdir+3  ) = 0
 zi( iidx    ) = mcol
 zi( iidx+1  ) = mm2
 zi( iidx+3  ) = mterms
 DO  j = 1, mwords
   zi( iidx+j+3 ) = zi (ispill+j+3 )
 END DO
 memcoln = mcol
!      PRINT *,' SMCSPL,A540,IIDX,ZI(1-5=',IIDX,(ZI(IIDX+KB),KB=0,4)
 GO TO 7777
 
! NO SPACE FOUND IN MEMORY FOR THIS DATA.
! CHECK IF COLUMN BEING REQUESTED IS THE PIVOT COLUMN
 
 1000 CONTINUE
!      PRINT *,' SMCSPL,MCOL,KCOL=',MCOL,KCOL
 IF ( mcol /= kcol ) GO TO 2000
 
! COLUMN REQUESTED IS THE PIVOT COLUMN, FIRST DETERMINE IF THERE
! ARE CONTIGUOUS BLOCKS IN THE FREE CHAIN THAT CAN BE MERGED TOGETHER
 
 IF ( memfre == 0 ) GO TO 1400
 index1 = memfre
 indext = memfre
 1100  CONTINUE
 index2 = zi( indext + 1 )
 IF ( index2 == 0 ) GO TO 1300
 
! COMPUTE THE LAST ADDRESS (PLUS 1) OF THIS FREE BLOCK AND COMPARE
! IT WITH THE BEGINNING OF BLOCK REFERENCED BY VARIABLE "INDEX1"
 
 iend   = index2 + zi( index2 + 2 )
 IF ( iend == index1 ) GO TO 1200
 
! BLOCK IS NOT CONTIGUOUS, GO AND TEST NEXT BLOCK IN CHAIN
 
 indext = index2
 GO TO 1100
 
! BLOCK IS CONTIGUOUS, MERGE THIS BLOCK AND THEN GO BACK TO
! TEST THE FREE CHAIN FOR SPACE FOR THE CURRENT PIVOT COLUMN.
!     EACH FREE CHAIN BLOCK HAS THE FOLLOWING FORMAT FOR THE FIRST 3
!     WORDS:
!             (1) = Pointer to previous block in chain
!             (2) = Pointer to next block in chain
!             (3) = Number of words in this block
!         (Note: Blocks are allocated from high memory to low:)
!              Memory Address N
!                     Block k
!                     Block k-1
!                        .
!                     Block 1
!              Memory Address N+M
 
 1200  CONTINUE
!      PRINT *,' SMCSPL,A1200,INDEX1,INDEX2=',INDEX1,INDEX2
 zi( index2+2 ) = zi( index1+2 ) + zi( index2+2 )
 
! RESET NEXT AND PREVIOUS POINTERS OF CHAIN BLOCKS
 
 indexp        = zi( index1 )
 zi( index2 )  = indexp
 IF ( indexp == 0 ) memfre = index2
 IF ( indexp /= 0 ) zi( indexp+1 ) = index2
!      PRINT *,' SMCSPL,A1200,MWORDS,ZI(INDEX1+2=',MWORDS,ZI(INDEX1+2)
 IF ( zi( index2+2 ) < (mwords+4) ) GO TO 1000
 iidx = index2
 GO TO 500
 
!  NO BLOCKS CONTIGUOUS WITH THIS BLOCK, GET NEXT BLOCK IN CHAIN
!  AND CHECK FOR CONTIGUOUS BLOCKS WITH IT.
 
 1300  CONTINUE
 index1 = zi( index1 + 1 )
 
!  FIRST CHECK THAT THERE IS ANOTHER BLOCK IN THE FREE CHAIN
 
 IF ( index1 == 0 ) GO TO 1400
 indext = memfre
 GO TO 1100
 1400  CONTINUE
 
! COLUMN REQUESTED IS THE PIVOT COLUMN, MUST FIND MEMORY TO READ
! THIS DATA INTO.  SEARCH FOR LAST COLUMN IN MEMORY WITH SUFFICIENT
! SPACE AND WRITE THAT COLUMN TO SPILL AND READ THE PIVOT COLUMN DATA
! INTO THE MEMORY THAT BECAME AVAILABLE.
 
 idir   = (memcoln-1) * 4 + 1
 kcolp1 = kcol + 1
 DO  i = memcoln, 1, -1
   IF ( i == kcol ) CYCLE
   idir = (i-1)*4 + 1
   
! CHECK TO SEE IF DATA ALREADY ON SPILL FILE
   
   IF ( zi( idir ) == 0 ) CYCLE
   
! DATA IS IN MEMORY, CHECK TO SEE IF ENOUGH SPACE
   
   imidx = zi( idir )
   IF ( zi( imidx+2 ) < (mwords+4) ) CYCLE
   
! SUFFICIENT SPACE, WRITE THIS COLUMN DATA TO THE SPILL FILE
! TO MAKE ROOM FOR THE PIVOTAL COLUMN DATA TO BE KEPT IN MEMORY.
! SKIP TO END OF FILE, BACKSPACE OVER EOF, CLOSE AND REOPEN FILE
! FOR WRITE WITH APPEND.
   
   CALL dssend( iscr1 )
   CALL skprec( iscr1, -1 )
   CALL CLOSE ( iscr1,  2 )
   CALL gopen ( iscr1, zi( ibuf2 ), 3 )
   im2    = zi( imidx+1 )
   iterms = zi( imidx+3 )
   length = im2 + iterms*ivwrds
   itemp( 1 ) = i
   itemp( 2 ) = im2
   itemp( 3 ) = 0
   itemp( 4 ) = iterms
   CALL WRITE ( iscr1, itemp , 4, 0 )
   CALL savpos( iscr1, kpos )
   CALL WRITE ( iscr1, zi( imidx+4 ), length, 1 )
   CALL CLOSE ( iscr1, 3 )
   CALL gopen ( iscr1, zi( ibuf2 ), 0 )
   
! SET DIRECTORY AND MOVE DATA INTO MEMORY LOCATION
   
!      PRINT *,' SMCSPL B1450,IMIDX,ISPILL=',IMIDX,ISPILL
   zi( idir    ) = 0
   zi( idir+3  ) = kpos
   zi( mdir    ) = imidx
   zi( mdir+3  ) = 0
   zi( imidx   ) = mcol
   zi( imidx+1 ) = mm2
   zi( imidx+3 ) = mterms
!      PRINT *,' SMCSPL,B1450,MCOL,MM2,MTERMS=',MCOL,MM2,MTERMS
   DO  j = 1, mwords
     zi( imidx+j+3 ) = zi (ispill+j+3 )
   END DO
   memcoln = mcol
!      PRINT *,' SMCSPL,A1450,ZI(1-5=',(ZI(IMIDX+KB),KB=0,4)
   GO TO 7777
 END DO
 
! NONE OF THE EXISTING IN-MEMORY ALLOCATIONS ARE LARGE ENOUGH.
! THEREFORE, MUST MERGE TWO TOGETHER TO TRY AND MAKE ENOUGH SPACE.
 
 loop1900:  DO  i = memcoln, 1, -1
   IF ( i == kcol ) CYCLE loop1900
   idir = ( i-1)*4 + 1
   IF ( zi( idir ) == 0 ) CYCLE loop1900
   imidx1  = zi( idir )
   ispace1 = zi( imidx1+2 )
!      PRINT *,' SMCSPL,B1800,IMIDX1,ISPACE1=',IMIDX1,ISPACE1
   iend1   = imidx1 + ispace1
   DO  j = memcoln, 1, -1
     IF ( j == kcol ) CYCLE
     IF ( j == i    ) CYCLE
     jdir = ( j-1 ) * 4 + 1
     IF ( zi( jdir ) == 0 ) CYCLE
     jmidx1  = zi( jdir )
     ispace2 = zi( jmidx1+2 )
!      PRINT *,' SMCSPL,I1800,JMIDX1,ISPACE2=',JMIDX1,ISPACE2
     iend2   = jmidx1 + ispace2
     IF ( IABS( imidx1-iend2 ) <= 4 ) GO TO 1700
     IF ( IABS( jmidx1-iend1 ) <= 4 ) GO TO 1700
     CYCLE
     
! COLUMNS J AND I HAVE CONTIGUOUS MEMORY, CHECK IF COMBINED SPACE IS
! LARGE ENOUGH FOR THIS COLUMN
     
     1700  itotal = ispace1 + ispace2
!      PRINT *,' SMCSPL,A1700,ISPACE1,ISPACE2,ITOTAL,MWORDS='
!     &,         ISPACE1,ISPACE2,ITOTAL,MWORDS
     IF ( itotal < (mwords+4) ) CYCLE loop1900
     
! SPACE IS LARGE ENOUGH, SO WRITE COLUMNS I AND J TO SPILL AND MERGE
! THE TWO AREAS TOGETHER.
! SKIP TO END OF FILE, BACKSPACE OVER EOF, CLOSE AND REOPEN FILE
! FOR WRITE WITH APPEND.
     
     CALL dssend ( iscr1 )
     CALL skprec ( iscr1, -1 )
     CALL CLOSE  ( iscr1,  2 )
     CALL gopen  ( iscr1, zi( ibuf2 ), 3 )
     
! WRITE COLUMN I TO SPILL FILE
     
     im2        = zi( imidx1+1 )
     iterms     = zi( imidx1+3 )
     ilen       = im2 + iterms*ivwrds
     itemp( 1 ) = i
     itemp( 2 ) = im2
     itemp( 3 ) = 0
     itemp( 4 ) = iterms
!      PRINT *,' SMCSPL WRITING COLUMN I=',I
     CALL WRITE ( iscr1, itemp, 4, 0 )
     CALL savpos( iscr1, kpos )
     CALL WRITE ( iscr1, zi( imidx1+4 ), ilen, 1 )
     
! RESET DIRECTORY FOR COLUMN I
     
     zi( idir   ) = 0
     zi( idir+3 ) = kpos
     
! WRITE COLUMN J TO THE SPILL FILE
     
     jm2    = zi( jmidx1+1 )
     jterms = zi( jmidx1+3 )
     jlen   = 4 + jm2 + jterms*ivwrds
     itemp( 1 ) = j
     itemp( 2 ) = jm2
     itemp( 3 ) = 0
     itemp( 4 ) = jterms
!      PRINT *,' SMCSPL,WRITING COLUMN J=',J
     CALL WRITE ( iscr1, itemp, 4, 0 )
     CALL savpos( iscr1, kpos )
     CALL WRITE ( iscr1, zi( jmidx1+4 ), jlen, 1 )
     
! RESET DIRECTORY FOR COLUMN J
     
     zi( jdir   ) = 0
     zi( jdir+3 ) = kpos
     CALL CLOSE ( iscr1, 3 )
     CALL gopen ( iscr1, zi( ibuf2 ), 0 )
     INDEX = jmidx1
     IF ( imidx1 < jmidx1 ) INDEX = imidx1
     
! MOVE DATA INTO MEMORY LOCATION
     
     PRINT*,' B1750,INDEX,ISPILL=',INDEX,ispill
     zi( INDEX   ) = mcol
     zi( INDEX+1 ) = mm2
     zi( INDEX+2 ) = itotal
     zi( INDEX+3 ) = mterms
     zi( mdir    ) = INDEX
     zi( mdir+3  ) = 0
     DO  k = 1, mwords
       zi( INDEX+k+3 ) = zi( ispill+k+3 )
     END DO
     memcoln = mcol
     GO TO 7777
   END DO
 END DO loop1900
 GO TO 7003
 
! NO SPACE FOUND AND COLUMN IS NOT THE PIVOTAL COLUMN, USE DATA
! FROM SPILL AREA
 
 2000  CONTINUE
 7777  CONTINUE
!      print *,' smcspl is returning, memfre=',memfre
!      ikb = memfre
!      do 9777 kk = 1, 100
!      if ( ikb .eq. 0 ) go to 9778
!      print *,' free block i,1-3=',kk,(zi(ikb+kb),kb=0,2)
!      ikb = zi( ikb+1 )
!9777  continue
!9778  continue
 RETURN
 7001  WRITE ( nout, 9001 ) ufm, kcol
9001  FORMAT(1X, a23,/,' UNEXPECTED END OF FILE FOR COLUMN ',i4  &
    ,' IN SUBROUTINE SMCSPL')
ierror = 3
GO TO 7070
7002  WRITE ( nout, 9002 ) ufm, kcol
9002  FORMAT(1X, a23,/,' UNEXPECTED END OF RECORD FOR COLUMN ',i4  &
    ,' IN SUBROUTINE SMCSPL')
ierror = 3
GO TO 7070
7003  WRITE ( nout, 9003 ) ufm, kcol
9003  FORMAT(1X,a23,/,' INSUFFICIENT CORE IN SUBROUTINE SMCSPL FOR'  &
    ,' SYMMETRIC DECOMPOSITION, COLUMN=',i6)
ierror = 1
GO TO 7070
7070  CALL smchlp
CALL mesage( -61, 0, 0 )
RETURN
END SUBROUTINE smcspl
