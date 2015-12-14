SUBROUTINE comugv
     
! FOR DDAM/EARTHQUAKE ANALYSES, COMBUGV COMBINES DISPLACEMENT
! COMPONENTS BY (1)ADDING THE COMPONENTS IN ABS VALUE AND (2)TAKING THE
! SQUARE ROOT OF THE SUMS OF THE SQUARES. AFTER THIS MODEULE, THE
! TWO OUTPUT DATA BLOCKS ARE N X NMODES, WHEREAS UGV IS N X (NMODES)(L)
! MODULE NRLSUM COMBINES STRESSES ACROSS MODES FOR EACH DIRECTION
! INDIVIDUALLY. THE OUTPUTS OF THIS MODULE HAVE THE DIRECTIONS
! COMBINED. BUT NRLSUM CAN WORK ON THEM (AFTER CASEGEN AND SDR2) BY
! SPECIFYING NDIR=1 IN THE DMAP STATEMENT FOR THOSE MODULES.
! THIS MODULE WILL ALSO COMBINE THE MAXIMUM RESPONSES ACROSS THE MODES
! BY USING SQRSS TO COME UP WITH ONE RESPONSE VECTOR. THEREFORE THIS
! MODULE COMBINES COMPONENTS TO GET MAXIMUM RESPONSES BY ADDING (UGVADD)
! AND BY SQRSS (UGVSQR). THEN IT TAKES EACH OF THESE AND TAKES SQRSS
! ACROSS THE MODES TO GET UGVADC AND UGVSQC, RESPECTIVELY.
! FINALLY, THE MODULE COMPUTES THE NRL SUMS FOR THE L DIRECTIONS
! TO USE CASEGEN,SDR2,ETC. ON UGVADC AND UGVSQC, IN CASEGEN,USE
! LMODES=NDIR=1 IN DMAP STATEMENT. FOR UGVNRL, JUST USE LMODES=1.
 
! COMBUGV UGV/UGVADD,UGVSQR,UGVADC,UGVSQC,UGVNRL/V,N,NMODES/V,N,NDIR $
 
 INTEGER :: buf1,buf2,buf3,ugv,ugvadd,ugvsqr,ugvadc,ugvsqc
 INTEGER :: ugvnrl
 INTEGER :: indb(2),oudb(2)
 DIMENSION nam(2),mcb(7),mcb1(7),mcb2(7)
 COMMON/unpakx/jout,iii,nnn,jncr
 COMMON/packx/iin,iout,ii,nn,incr
 COMMON/system/ibuf
 COMMON/BLANK/nmodes,ndir
 COMMON/zzzzzz/z(1)
 DATA ugv,ugvadd,ugvsqr,ugvadc,ugvsqc/101,201,202,203,204/
 DATA ugvnrl/205/
 DATA nam/4HCOMB,4HUGV /
 
! OPEN CORE AND BUFFERS
 
 lcore=korsz(z)
 buf1=lcore-ibuf+1
 buf2=buf1-ibuf
 buf3=buf2-ibuf
 lcore=buf3-1
 IF(lcore <= 0)GO TO 1008
 mcb(1)=ugv
 CALL rdtrl(mcb)
 ncol=mcb(2)
 nrow=mcb(3)
 IF(ncol /= nmodes*ndir)GO TO 1007
 IF(lcore < 4*nrow)GO TO 1008
 mcb1(1)=ugvadd
 mcb1(2)=0
 mcb1(3)=nrow
 mcb1(4)=2
 mcb1(5)=1
 mcb1(6)=0
 mcb1(7)=0
 mcb2(1)=ugvsqr
 mcb2(2)=0
 mcb2(3)=nrow
 mcb2(4)=2
 mcb2(5)=1
 mcb2(6)=0
 mcb2(7)=0
 
 jout=1
 iii=1
 nnn=nrow
 jncr=1
 iin=1
 iout=1
 ii=1
 nn=nrow
 incr=1
 
 CALL gopen(ugv,z(buf1),0)
 CALL gopen(ugvadd,z(buf2),1)
 CALL gopen(ugvsqr,z(buf3),1)
 
! UNPACK NDIR COLUMNS OF UGV WHICH CORRESPOND TO A SINGLE MODE
 
 nm1=nmodes-1
 nd1=ndir-1
 DO  i=1,nmodes
   
! POINTER TO PROPER MODE IN 1ST DIRECTION
   
   nskip=i-1
   IF(nskip == 0)GO TO 20
   DO  ll=1,nskip
     CALL fwdrec (*1002,ugv)
   END DO
   
! UNPACK VECTOR
   
   20 CALL unpack (*25,ugv,z(1))
   GO TO 40
   
   25 DO  j=1,nrow
     z(j)=0.
   END DO
   
! SKIP TO NEW DIRECTION, UNPACK, SKIP AND UNAPCK
   
   40 IF(nd1 == 0)GO TO 100
   DO  j=1,nd1
     IF(nm1 == 0)GO TO 50
     DO  jj=1,nm1
       CALL fwdrec (*1002,ugv)
     END DO
     
     50 jnrow = j*nrow
     CALL unpack (*55,ugv,z(jnrow+1))
     CYCLE
     55 DO  jj=1,nrow
       z(j*nrow+jj)=0.
     END DO
     
   END DO
   
! NOW PERFORM EACH OPERATION AND STORE INTO Z(3*NROW+1)
   
   DO  kk=1,nrow
     z(3*nrow+kk)=ABS(z(kk))+ABS(z(nrow+kk))+ABS(z(2*nrow+kk))
   END DO
   CALL pack(z(3*nrow+1),ugvadd,mcb1)
   
   DO  kk=1,nrow
     z(3*nrow+kk)=SQRT(z(kk)**2+z(nrow+kk)**2+z(2*nrow+kk)**2)
   END DO
   CALL pack(z(3*nrow+1),ugvsqr,mcb2)
   GO TO 110
   
! JUST ONE DIRECTION ON UGV- COPY TO DATA BLOCKS
   
   100 CALL pack(z(1),ugvadd,mcb1)
   CALL pack(z(1),ugvsqr,mcb2)
   
! DONE FOR THIS MODE - GET ANOTHER
   
   110 CALL REWIND(ugv)
   CALL fwdrec (*1002,ugv)
   
 END DO
 
 CALL CLOSE(ugvadd,1)
 CALL CLOSE(ugvsqr,1)
 CALL wrttrl(mcb1)
 CALL wrttrl(mcb2)
 
! NOW COMPUTE NRL SUMS FOR THE L DIRECTIONS
 
 mcb1(1)=ugvnrl
 mcb1(2)=0
 mcb1(3)=nrow
 mcb1(4)=2
 mcb1(5)=1
 mcb1(6)=0
 mcb1(7)=0
 CALL REWIND(ugv)
 CALL fwdrec (*1002,ugv)
 CALL gopen(ugvnrl,z(buf2),1)
 
 DO  nd=1,ndir
   
! SET UP VECTOR OF MAXIMUM DISPLACEMENT COMPONENTS AND VECTOR OF SUMS
   
   DO  i=1,nrow
     z(i)=0.
     z(2*nrow+i)=0.
   END DO
   
   DO  i=1,nmodes
     
     CALL unpack (*1220,ugv,z(nrow+1))
     
! COMPARE TO MAXIMUM COMPONENTS
     
     DO  j=1,nrow
       IF (ABS(z(nrow+j)) > z(j))z(j)=ABS(z(nrow+j))
       z(2*nrow+j)=z(2*nrow+j)+z(nrow+j)**2
     END DO
     
! GET ANOTHER DISPLACEMENT VECTOR CORRESPONDING TO ANOTHER MODE
     
   END DO
   
! SUBTRACT THE MAXIMA FROM THE SUMS
   
   DO  j=1,nrow
     z(2*nrow+j)=z(2*nrow+j)-z(j)**2
     
! TAKE SQUARE ROOT AND ADD IN THE MAXIMA
     
     z(2*nrow+j)=SQRT(z(2*nrow+j))+z(j)
   END DO
   
! PACK RESULTS ANG GET ANOTHER DIRECTION
   
   CALL pack(z(2*nrow+1),ugvnrl,mcb1)
 END DO
 
 CALL CLOSE(ugv,1)
 CALL CLOSE(ugvnrl,1)
 CALL wrttrl(mcb1)
 
! NOW LETS COMBINE RESPONSES OVER THE MODES USING SQRSS. DO FOR BOTH
! UGVADD AND UGVSQR. THE RESULT WILL BE ONE DISLPACEMENT VECTOR.
! (BOTH UGVADD AND UGVSQR ARE N X M ( M= NO. OF MODES)
 
 indb(1)=ugvadd
 indb(2)=ugvsqr
 oudb(1)=ugvadc
 oudb(2)=ugvsqc
 
 DO  i=1,2
   
   mcb(1)=indb(i)
   CALL rdtrl(mcb)
   ncol=mcb(2)
   nrow=mcb(3)
   mcb1(1)=oudb(i)
   mcb1(2)=0
   mcb1(3)=nrow
   mcb1(4)=2
   mcb1(5)=1
   mcb1(6)=0
   mcb1(7)=0
   IF(ncol /= nmodes)GO TO 1007
   
   CALL gopen(indb(i),z(buf1),0)
   CALL gopen(oudb(i),z(buf2),1)
   
   DO  j=1,nrow
     z(j)=0.
   END DO
   
! UNPACK THE COLUMNS OF INDB AND ACCUMULATE SUMS OF SQUARES
   
   DO  j=1,nmodes
     CALL unpack (*150,indb(i),z(nrow+1))
     
     DO  k=1,nrow
       z(k)=z(k)+z(nrow+k)**2
     END DO
     
   END DO
   
   DO  k=1,nrow
     z(k)=SQRT(z(k))
   END DO
   
   CALL pack(z(1),oudb(i),mcb1)
   
   CALL CLOSE(indb(i),1)
   CALL CLOSE(oudb(i),1)
   CALL wrttrl(mcb1)
   
 END DO
 
 RETURN
 
 1002 CALL mesage(-2,ugv,nam)
 1007 CALL mesage(-7,0,nam)
 1008 CALL mesage(-8,0,nam)
 RETURN
END SUBROUTINE comugv
