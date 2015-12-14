SUBROUTINE subcc
     
!     THIS ROUTINE WAS ORIGINALLY CALLED SUBD
 
 COMPLEX :: gusamp,sbkde1,sbkde2,  &
     f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,am5tt,am6,  &
     sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,am5t,  &
     ai,a,b,bsycon,alp,f1,am1,aln,blkapm,bkdel3,f1s,c1, c2p,c2n,  &
     c2,amtest,ft2,blam1,ft3,am2,sum1,sum2,f2,blam2,  &
     ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t,  &
     c1p,c1n,bkdel1,bkdel2,blkap1,arg,arg2,ft3tst,  &
     bc,bc2,bc3,bc4,bc5,ca1,ca2,ca3,ca4,  &
     clift,cmomt,pres1,pres2,pres3,pres4,qres4,fq7,  &
     fqa,fqb,ss,t1,t2,t3,t4,const,const2,const3,const4,  &
     const5,const6,c1a,c2a,cexp1,cexp2,cexp1a,cexp2a
 DIMENSION       pres1(21),pres2(21),pres3(21),pres4(21),qres4(21),  &
     sbkde1(201),sbkde2(201), sumsv1(201),sumsv2(201),svkl1(201),svkl2(201),  &
     xlsv1(21),xlsv2(21),xlsv3(21),xlsv4(21)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,ibbout
 COMMON /blk1  / scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 COMMON /blk2  / bsycon
 COMMON /blk3  / sbkde1,sbkde2,f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,  &
     am5tt,am6,sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,  &
     am5t,a,b,alp,f1,am1,aln,blkapm,bkdel3,f1s,c1,c2p,  &
     c2n,c2,amtest,ft2,blam1,ft3,am2,sum1,sum2,f2,  &
     blam2,ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t,c1p,c1n,  &
     bkdel1,bkdel2,blkap1,arg,arg2,ft3tst,bc,bc2,bc3,  &
     bc4,bc5,ca1,ca2,ca3,ca4,clift,cmomt,pres1,pres2,  &
     pres3,pres4,qres4,fqa,fqb,fq7
 COMMON /blk4  / i,r,y,a1,b1,c4,c5,gl,i6,i7,jl,nl,ri,rt,r5,sn,sp,  &
     xl,y1,amu,gam,idx,inx,nl2,rl1,rl2,rq1,rq2,xl1,  &
     alp1,alp2,gamn,gamp,iner,iout,redf,stag,step,  &
     amach,betnn,betnp,bkap1,xlsv1,xlsv2,xlsv3,xlsv4,  &
     alpamp,amoaxs,gusamp,disamp,pitaxs,pitcor
 
 am6 = 0.0
 f5  = 0.0
 am5 = 0.0
 s1  = sps + sns
 s2  = sigma - sps*del
 s3  = sps/(dstr**2)
 s4  = sns/dstr
 s5  = del*sns + sigma
 ss  = CEXP(-ai*sigma)
 DO  iout = 1,200
   IF (iout > i7) GO TO 240
   r5  = iout - 1
   rq1 = SQRT((r5*pi/sns)**2+scrk**2)
   rq2 =-rq1
   c4  = (rq1*s1+s2)/(2.0*pi)
   c5  = (rq2*s1+s2)/(2.0*pi)
   bc2 = bc/(2.0*svkl1(iout))*CEXP(-ai*(-s2)*(sps+3.0*sns)/  &
       (2.0*s1))/(2.0*pi*ai)
   bc3 = bc2*svkl1(iout)/svkl2(iout)
   bc4 = bc/(2.0*svkl1(iout))*CEXP(ai*(-s2)*(sns-sps)/ (2.0*s1))/(2.0*pi*ai)
   bc5 = bc4*svkl1(iout)/svkl2(iout)
   f5t = 0.0
   am5t= 0.0
   am5tt = 0.0
   DO  jl = 1,nl
     qres4(jl) = 0.0
   END DO
   DO  iner = 1,200
     r    = iner - 1
     gamp = 2.0*pi*r - s2
     gamn =-2.0*pi*r - s2
     c1p  = (gamp/dstr) - scrk
     c2p  = (gamp/dstr) + scrk
     alp  = gamp*s3 - s4*CSQRT(c1p)*CSQRT(c2p)
     bkdel3 = sbkde1(iner)
     IF (iner <= i6) GO TO 20
     CALL akapm (alp,bkdel3)
     sbkde1(iner) = bkdel3
     20 CONTINUE
     t1   = alp*sps-gamp
     t2   = alp*dstr**2-gamp*sps
     sum1 = sumsv1(iout)*CEXP(ai*t1)*bkdel3*t1/ (t2*svkl1(iout)*(alp-rq1))
     sum3 = sumsv2(iout)*CEXP(ai*t1)*bkdel3*t1/ (t2*svkl2(iout)*(alp-rq2))
     IF (iner == 1) GO TO 40
     c1n  = (gamn/dstr) - scrk
     c2n  = (gamn/dstr) + scrk
     aln  = gamn*s3 - s4*CSQRT(c1n)*CSQRT(c2n)
     bkdel3 = sbkde2(iner)
     IF (iner <= i6) GO TO 30
     CALL akapm (aln,bkdel3)
     sbkde2(iner) = bkdel3
     30 CONTINUE
     t1   = aln*sps - gamn
     t2   = aln*dstr**2 - gamn*sps
     sum2 = sumsv1(iout)*CEXP(ai*t1)*bkdel3*t1/ (t2*svkl1(iout)*(aln-rq1))
     sum4 = sumsv2(iout)*CEXP(ai*t1)*bkdel3*t1/ (t2*svkl2(iout)*(aln-rq2))
     40 CONTINUE
     IF (iner == 1) sum2 = 0.0
     IF (iner == 1) sum4 = 0.0
     c1p = CEXP(-ai*(alp-del)*sps)
     c2p = CEXP(-ai*(alp-del)*sns)
     c1n = CEXP(-ai*(aln-del)*sps)
     c2n = CEXP(-ai*(aln-del)*sns)
     f5t = f5t + (sum1+sum3)*ai*ss/(alp-del)*(c1p-c2p) +  &
         (sum2+sum4)*ss*ai/(aln-del)*(c1n-c2n)
     am5t= am5t + (sum1+sum3)*ss*(ai*sps*c1p/(alp-del) - ai*sns*c2p/  &
         (alp-del) + 1.0/((alp-del)**2)*(c1p-c2p) + ai*(2.0-sps)/  &
         (alp-del)*(c1p-c2p)) + (sum2+sum4)*ss*(ai*sps*c1n/(aln-del)  &
         - ai*sns*c2n/(aln-del) + 1.0/((aln-del)**2)*(c1n-c2n) +  &
         ai*(2.0-sps)/(aln-del)*(c1n-c2n))
     temp  = (sps-sns)/rl1
     const = (sum1+sum3)*ss
     const2= (sum2+sum4)*ss
     c1a   =-ai*(alp-del)
     c2a   =-ai*(aln-del)
     cexp1 = CEXP(c1a*sns)
     cexp2 = CEXP(c2a*sns)
     cexp1a= CEXP(c1a*temp)
     cexp2a= CEXP(c2a*temp)
     DO  jl = 1,nl
       qres4(jl) = qres4(jl) - (const*cexp1+const2*cexp2)
       cexp1 = cexp1*cexp1a
       cexp2 = cexp2*cexp2a
     END DO
     betnp = ( 2.0*r*pi-s5)/s1
     betnn = (-2.0*r*pi-s5)/s1
     c1p   = CEXP(-2.0*pi*r*ai*sns/s1)
     c2p   = CEXP(-2.0*pi*r*ai*sps/s1)
     c1n   = CEXP(2.0*pi*r*ai*sns/s1)
     c2n   = CEXP(2.0*pi*r*ai*sps/s1)
     t1    = CEXP(-ai*betnp*sps)
     t2    = CEXP(-ai*betnp*sns)
     t3    = CEXP(-ai*betnn*sps)
     t4    = CEXP(-ai*betnn*sns)
     ca1   = ai*ss/betnp*(t1-t2)
     ca2   = ai*ss/betnn*(t3-t4)
     ca3   = ss*(ai*sps/betnp*t1 - ai*sns*t2/betnp+(t1-t2)/  &
         betnp**2 + (2.0-sps)*ai/betnp*(t1-t2))
     ca4   = ss*(ai*sps*t3/betnn - ai*sns*t4/betnn+(t3-t4)/  &
         betnn**2 + (2.0-sps)*ai/betnn*(t3-t4))
     IF (iner > 1) GO TO 70
     f5t   = f5t - sumsv1(iout)*(bc2*c1p-bc4*c2p)/(r-c4)*ca1 -  &
         sumsv2(iout)*(bc3*c1p-bc5*c2p)/(r-c5)*ca1
     am5t  = am5t - sumsv1(iout)*(bc2*c1p-bc4*c2p)/(r-c4)*ca3 -  &
         sumsv2(iout)*(bc3*c1p-bc5*c2p)/(r-c5)*ca3
     temp  = (sps-sns)/rl1
     const = ss*sumsv1(iout)*(bc2*c1p-bc4*c2p)/(r-c4)
     const2= ss*sumsv2(iout)*(bc3*c1p-bc5*c2p)/(r-c5)
     c1a   =-ai*betnp
     cexp1 = CEXP(c1a*sns)
     cexp1a= CEXP(c1a*temp)
     DO  jl = 1,nl
       qres4(jl) = qres4(jl)+const*cexp1+const2*cexp1
       cexp1 = cexp1*cexp1a
     END DO
     GO TO 90
     70 CONTINUE
     f5t = f5t - sumsv1(iout)*((bc2*c1p-bc4*c2p)/(r-c4)*ca1 -  &
         (bc2*c1n-bc4*c2n)/(r+c4)*ca2) - sumsv2(iout)*  &
         ((bc3*c1p-bc5*c2p)/(r-c5)*ca1-(bc3*c1n-bc5*c2n)/(r+c5)*ca2)
     am5t= am5t - sumsv1(iout)*((bc2*c1p-bc4*c2p)/(r-c4)*ca3-(bc2*c1n-  &
         bc4*c2n)/(r+c4)*ca4)-sumsv2(iout)*((bc3*c1p-bc5*c2p)/  &
         (r-c5)*ca3-(bc3*c1n-bc5*c2n)/(r+c5)*ca4)
     temp   = (sps-sns)/rl1
     const  = (bc2*c1p-bc4*c2p)/(r-c4)
     const2 = (bc2*c1n-bc4*c2n)/(r+c4)
     const3 = (bc3*c1p-bc5*c2p)/(r-c5)
     const4 = (bc3*c1n-bc5*c2n)/(r+c5)
     const5 = ss*sumsv1(iout)
     const6 = ss*sumsv2(iout)
     c1a    =-ai*betnp
     c2a    =-ai*betnn
     cexp1  = CEXP(c1a*sns)
     cexp2  = CEXP(c2a*sns)
     cexp1a = CEXP(c1a*temp)
     cexp2a = CEXP(c2a*temp)
     DO  jl = 1,nl
       qres4(jl) = qres4(jl) + const5*(const*cexp1-const2*cexp2) +  &
           const6*(const3*cexp1-const4*cexp2)
       cexp1  = cexp1*cexp1a
       cexp2  = cexp2*cexp2a
     END DO
     90 CONTINUE
     IF (cabs((am5tt-am5t)/am5t) < 0.001) GO TO 110
     am5tt  = am5t
   END DO
   GO TO 200
   110 CONTINUE
   IF (iner  <=  i6) GO TO 120
   i6 = iner
   120 CONTINUE
   f5  = f5  + f5t
   am5 = am5 + am5t
   DO  jl = 1,nl
     pres4(jl) = pres4(jl) + qres4(jl)
   END DO
   alp1 = (2.0*pi*c4-del*sns-sigma)/s1
   alp2 = (2.0*pi*c5-del*sns-sigma)/s1
   t1   = 1.0 - CEXP(-2.0*pi*ai*c4)
   t2   = 1.0 - CEXP(-2.0*pi*ai*c5)
   c1p  = CEXP(-2.0*pi*ai*c4*sns/s1)/(t1)
   c2p  = CEXP( 2.0*pi*ai*c4*sns/s1)/(t1)
   c1n  = CEXP(-2.0*pi*ai*c5*sns/s1)/(t2)
   c2n  = CEXP( 2.0*pi*ai*c5*sns/s1)/(t2)
   t1   = CEXP(-ai*sps*alp1)
   t2   = CEXP(-ai*sns*alp1)
   t3   = CEXP(-ai*sps*alp2)
   t4   = CEXP(-ai*sns*alp2)
   ca1  = ai*ss/alp1*(t1-t2)
   ca2  = ai*ss/alp2*(t3-t4)
   ca3  = ss*(ai*sps*t1/alp1 - ai*sns*t2/alp1 + (t1-t2)/  &
       alp1**2 + (2.0-sps)*ai/alp1*(t1-t2))
   ca4  = ss*(ai*sps*t3/alp2 - ai*sns*t4/alp2 + (t3-t4)/  &
       alp2**2 + (2.0-sps)*ai/alp2*(t3-t4))
   f5   = f5 - 2.0*pi*ai*sumsv1(iout)*(bc2*c1p-bc4*c2p)*ca1 - 2.0*pi*  &
       ai*sumsv2(iout)*(bc3*c1n-bc5*c2n)*ca2
   am5  = am5 - 2.0*pi*ai*sumsv1(iout)*(bc2*c1p-bc4*c2p)*ca3 - 2.0*  &
       pi*ai*sumsv2(iout)*(bc3*c1n-bc5*c2n)*ca4
   temp = (sps-sns)/rl1
   const  = ss*2.0*pi*ai
   const2 = const*sumsv1(iout)*(bc2*c1p-bc4*c2p)
   const3 = const*sumsv2(iout)*(bc3*c1n-bc5*c2n)
   c1a    =-ai*alp1
   c2a    =-ai*alp2
   cexp1  = CEXP(c1a*sns)
   cexp2  = CEXP(c2a*sns)
   cexp1a = CEXP(c1a*temp)
   cexp2a = CEXP(c2a*temp)
   DO  jl = 1,nl
     pres4(jl) = pres4(jl)+const2*cexp1+const3*cexp2
     cexp1  = cexp1*cexp1a
     cexp2  = cexp2*cexp2a
   END DO
   IF (cabs((am5-am6)/am5) < 0.0009) GO TO 160
   am6 = am5
 END DO
 GO TO 220
 160 CONTINUE
 clift = f1 + f2 - f2p + f4 + f5
 cmomt = am1 + am2 - am2p + am4 + am5 - amoaxs*clift
 GO TO 270
 
 200 WRITE  (ibbout,210) ufm
 210 FORMAT (a23,' - AMG MODULE -SUBROUTINE SUBCC.  AM5T LOOP DID NOT',  &
     ' CONVERGE.')
 GO TO 260
 220 WRITE  (ibbout,230) ufm
 230 FORMAT (a23,' - AMG MODULE -SUBROUTINE SUBCC.  AM5 LOOP DID NOT',  &
     ' CONVERGE.')
 GO TO 260
 240 WRITE  (ibbout,250) ufm,i7
 250 FORMAT (a23,' - AMG MODULE -SUBROUTINE SUBCC.  OUTER LOOP OF AM5',  &
     ' EXCEEDED I7 (',i6,1H))
 260 CALL mesage (-61,0,0)
 270 CONTINUE
 RETURN
END SUBROUTINE subcc
