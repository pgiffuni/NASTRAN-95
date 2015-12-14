SUBROUTINE suba
     
!     UNSTEADY FLOW ANAYSIS OF A SUPERSONIC CASCADE
 
!     LIFT AND MOMENT COEFICIENT
 
 DIMENSION       pres1(21),pres2(21),pres3(21),pres4(21),qres4(21),  &
     sbkde1(201),sbkde2(201),sumsv1(201),sumsv2(201),  &
     svkl1(201),svkl2(201),xlsv1(21),xlsv2(21), xlsv3(21),xlsv4(21)
 COMPLEX :: sbkde1,sbkde2,f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,  &
     am5tt,am6,sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,  &
     am5t,ai,a,b,bsycon,alp,f1,am1,aln,blkapm,bkdel3,  &
     f1s,c1,c2p,c2n,c2,amtest,ft2,blam1,ft3,am2,sum1,  &
     sum2,f2,blam2,ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t,  &
     c1p,c1n,bkdel1,bkdel2,blkap1,arg,arg2,ft3tst,bc,  &
     bc2,bc3,bc4,bc5,ca1,ca2,ca3,ca4,clift,cmomt,  &
     pres1,pres2,pres3,pres4,qres4,fqa,fqb,t1,t2,t3,t4,  &
     gusamp,fq7,cexp3,cexp4,cexp5,const,c1a,c2a
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
 
 s1   = sps - sns
 s2   = sps*del - sigma
 s3   = sps/(dstr**2)
 s4   = sns/dstr
 s0   = 2.0 - sps + sns
 t1   = CEXP(-ai*sigma)
 t2   = CEXP(ai*sigma)
 a1   = 2.0*pi/s1
 b1   = s2/s1
 gam  = s2
 c1p  = gam/dstr - scrk
 c1n  = gam/dstr + scrk
 alp  = gam*s3 + s4*CSQRT(c1p)*CSQRT(c1n)
 bc   = -b1/alp*bsycon/SIN(pi*b1/a1)
 t3   = alp - del
 f1   = (alp-amu)/t3*ai*sns/(beta*(gam-alp*sps))
 arg2 = del
 CALL akapm (arg2,bkdel1)
 arg  = del - gl
 CALL akapm  (arg,bkdel2)
 CALL dlkapm (arg2,blkap1)
 inx  = 0
 CALL drkapm (alp,inx,blkapm)
 f1   = f1*bkdel1/blkapm*(-t3/(t3+gl)*a*ai*bkdel2/bkdel1 + b*blkap1+b/t3)
 f1s  = f1
 nl   = 10
 rl1  = nl - 1
 cexp3 = CEXP(-ai*t3/rl1*s1)
 pres1(1) = f1s
 nnl1 = nl - 1
 DO  jl = 1,nnl1
   pres1(jl+1) = pres1(jl)*cexp3
 END DO
 f1   = f1*ai/t3*(CEXP(-ai*t3*s1)-1.0)
 am1  = f1/(ai*t3)-f1s/(ai*t3)*s1*CEXP(-ai*t3*s1)
 amtest = 0.0
 fqb  = bkdel1/(beta*bc)*CEXP(ai*s2/2.0)* (-a*ai*bkdel2/bkdel1+b*blkap1)
 DO  i = 1,200
   r    = i
   gamp = 2.0*pi*r + s2
   gamn =-2.0*pi*r + s2
   c1p  = (gamp/dstr) - scrk
   c2p  = (gamp/dstr) + scrk
   alp  = gamp*s3 + s4*CSQRT(c1p)*CSQRT(c2p)
   t3   = alp - del
   idx  = i
   CALL drkapm (alp,idx,blkapm)
   c1   = (alp-amu)/t3*ai*sns/(beta*(gamp-alp*sps))*bkdel1/  &
       (blkapm)*(-t3/(t3+gl)*a*ai*bkdel2/bkdel1+b*blkap1+b/t3)
   c1n  = (gamn/dstr) - scrk
   c2n  = (gamn/dstr) + scrk
   aln  = gamn*s3 + s4*CSQRT(c1n)*CSQRT(c2n)
   t4   = aln - del
   idx  =-i
   CALL drkapm (aln,idx,blkapm)
   c2   = (aln-amu)/t4*ai*sns/(beta*(gamn-aln*sps))*bkdel1/  &
       (blkapm)*(-t4/(t4+gl)*a*ai*bkdel2/bkdel1+b*blkap1+b/t4)
   f1   = f1+c1*ai/t3*(CEXP(-ai*t3*s1)-1.0)+c2*ai/ t4*(CEXP(-ai*t4*s1)-1.0)
   am1  = am1+c1/(ai*t3)*(-s1*CEXP(-ai*t3*s1)+ai/  &
       t3*(CEXP(-ai*t3*s1)-1.0))+c2/(ai*t4)*  &
       (-s1*CEXP(-ai*t4*s1)+ai/t4*(CEXP(-ai*t4*s1)-1.0))
   c2a  = c2
   c1a  = c1
   aa   = s1/rl1
   cexp3 = CEXP(-ai*t3*aa)
   cexp4 = CEXP(-ai*t4*aa)
   temp  = 2.0*pi*r
   cexp5 = CEXP(ai*(sigma-sns*del)/s1*aa)
   const = 4.0*fqb/temp
   pres1(1) = pres1(1) + c1 + c2
   DO  jl = 1,nnl1
     const = const*cexp5
     c1a   = c1a*cexp3
     c2a   = c2a*cexp4
     pres1(jl+1) = pres1(jl+1) + c1a + c2a
     pres1(jl+1) = pres1(jl+1) + const*SIN(temp*jl/rl1)
   END DO
   IF (cabs((am1-amtest)/am1) < 0.0005) GO TO 45
   amtest = am1
 END DO
 GO TO 9992
 9992 WRITE  (ibbout,3005) ufm
 3005 FORMAT (a23,' FROM AMG MODULE. AM1 LOOP IN SUBROUTINE SUBA DID ',  &
     'NOT CONVERGE.')
 CALL mesage (-61,0,0)
 45 CONTINUE
 aa    = s1/rl1
 cexp3 = CEXP(ai*(sigma-sns*del)/rl1)
 const = fqb
 temp  = 2.0*aa/(sps-sns)
 pres1(1) = pres1(1) - fqb
 DO  jl = 1,nnl1
   const = const*cexp3
   pres1(jl+1) = pres1(jl+1) - const*(1.0-jl*temp)
 END DO
 y    = 0.0
 y1   = sns
 arg  = del - gl
 CALL alamda (arg,y,blam1)
 CALL alamda (arg,y1,blam2)
 CALL akappa (arg,bkap1)
 ft2  = a*ai*(del-gl-amu)*blam1/bkap1
 ft2t = a*ai*(del-gl-amu)*blam2/bkap1
 arg  = del
 CALL alamda (arg,y,blam1)
 CALL alamda (arg,y1,blam2)
 CALL akappa (arg,bkap1)
 gam  = SQRT(del**2-scrk**2)
 s5   = SIN(sns*gam)
 s6   = COS(sns*gam)
 c1   =-1.0/(beta*gam*s5)
 c1t  = c1*(ai*sps*t2*s6-sns*del/gam*t2*s5)-blam2/bkap1*del/gam*(s5  &
     +gam*sns*s6)/(gam*s5)
 c1   = c1*(arg/gam*sns*s5+ai*sps*t2)-blam1/bkap1*del/(gam*s5)*(s5/  &
     gam+sns*s6)
 ft3  =-b*(blam1/bkap1+(del-amu)*c1)
 ft3t =-b*(blam2/bkap1+(del-amu)*c1t)
 IF (gl == 0.0) GO TO 50
 f2   = ft2*(CEXP(2.0*ai*gl)-CEXP(ai*gl*s1))/(ai*gl)+  &
     ft3*s0+b*ai*(del-amu)*blam1/bkap1*(4.0-s1**2)/2.0
 am2  = ft2*(2.0*CEXP(2.0*ai*gl)/(ai*gl)-s1/(ai*gl)*CEXP(gl*ai*s1)+  &
     (CEXP(2.0*ai*gl)-CEXP(ai*s1*gl))/gl**2)+ft3*(4.0-s1**2)/2.0  &
     +b*ai*(del-amu)*blam1/bkap1*(8.0-s1**3)/3.0
 f2p  = ft2t*t1*CEXP(ai*gl*sns)/(ai*gl)*(CEXP(2.0*ai*gl)-  &
     CEXP(ai*gl*s1))+ft3t*t1*s0+b*ai*(del-amu)*t1*blam2/ bkap1*(s0**2/2.0+sps*s0)
 am2p = ft2t*t1*(CEXP(ai*gl*sps)/(ai*gl)*s0*CEXP(ai*gl*s0)+  &
     CEXP(ai*gl*sps)/(gl**2)*(CEXP(ai*gl*s0)-1.0))+  &
     ft3t*t1*s0**2/2.0+b*ai*(del-amu)*t1*blam2/bkap1*(s0**3/3.0+ sps*s0**2/2.0)
 GO TO 55
 50 CONTINUE
 f2   = ft2*s0+ft3*s0+b*ai*(del-amu)*blam1/bkap1*(4.-s1**2)/2.
 am2  = ft2*(4.0-s1**2)/2.0+ft3*(4.0-s1**2)/2.0+b*ai*(del-amu)*  &
     blam1/bkap1*(8.0-s1**3)/3.0
 f2p  = ft2t*t1*s0+ft3t*t1*s0+b*ai*(del-amu)*t1*blam2/bkap1*(s0**2  &
     /2.0+sps*s0)
 am2p = ft2t*t1*s0**2/2.0+ft3t*t1*s0**2/2.0+b*ai*(del-amu)*t1*blam2  &
     /bkap1*(s0**3/3.0+sps*s0**2/2.0)
 55 CONTINUE
 nl2  = 20
 rl2  = nl2 - 1
 aa   = sps - sns
 const = b*ai*(del-amu)*blam1/bkap1
 temp = s0/rl2
 c1a  = ai*gl
 cexp3 = CEXP(c1a*aa)
 cexp4 = CEXP(c1a*temp)
 DO  jl = 1,nl2
   xl = aa + temp*(jl-1)
   pres2(jl) = ft2*cexp3 + ft3+const*xl
   cexp3 = cexp3*cexp4
 END DO
 CALL subbb
 RETURN
END SUBROUTINE suba
