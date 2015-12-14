SUBROUTINE hdsket (x,y,z,np,nc)
     
!     THIS SUBROUTINE SETS UP PEN MOTION INDICATORS.
 
 
 REAL, INTENT(IN OUT)                     :: x(1)
 REAL, INTENT(IN OUT)                     :: y(1)
 REAL, INTENT(IN OUT)                     :: z(1)
 INTEGER, INTENT(IN OUT)                  :: np
 INTEGER, INTENT(IN OUT)                  :: nc
 INTEGER :: xcc,x1skt,y1skt,z1skt,x21,y21,z21,xe,ye,xu,yu, xi,yi,zi,di,w,iz(1)
 
 COMMON /hdptrs/ xdum,xcc,xasolv,yasolv,zasolv,x1skt,y1skt,z1skt,  &
     zcoef1,zcoef,icount,irct,x21,y21,z21,iia,xe,ye,  &
     xu,yu,xi,yi,zi,di,ibeg,iend,ict,icct,w
 COMMON /zzzzzz/ rz(1)
 COMMON /hdsc  / scx,yaw,rol,pit,lz,vp,jjj,icore
 EQUIVALENCE     (iz(1),rz(1))
 
 l   = np
 li  = np
 IF (l <= 2) GO TO 50
 lx  = 1
 npx = np
 1 npx = npx-1
 i   = lx
 DO  m = i,npx
   rx = 0
   a  = x(m+1) - x(m)
   b  = y(m+1) - y(m)
   c  = z(m+1) - z(m)
   IF (a /= 0.) CYCLE
   IF (b /= 0.) CYCLE
   IF (c /= 0.) CYCLE
   ix = m
   ix1skt = npx
   DO  mx = ix,ix1skt
     x(mx) = x(mx+1)
     y(mx) = y(mx+1)
     z(mx) = z(mx+1)
   END DO
   rx = 1
   lx = m
   IF (lx == npx) EXIT
   GO TO 1
 END DO
 10 CONTINUE
 IF (rx == 1.) npx = npx - 1
 np = npx + 1
 li = np
 IF (np <= 2) GO TO 50
 ix = 0
 m1 = 0
 m  = 1
 is = np - 1
 20 CONTINUE
 m  = m  + ix
 m1 = m1 + ix + 1
 IF (m-1 == li) GO TO 70
 
!     SEARCH FOR MATCHING COORDINATES.
 
 DO  j = m,is
   t = x(j+1) - x(m)
   u = z(j+1) - z(m)
   v = y(j+1) - y(m)
   IF (t /= 0.) CYCLE
   IF (v /= 0.) CYCLE
   IF (u /= 0.) CYCLE
   np = np + 1
   
!     MATCH FOUND.....STORE COORDINATES AND SET SWITCH TO LIFT PEN
!     AND/OR END SET.
   
   ix = j + 2 - m
   ix1skt = j - is + 1
   DO  ik = 1,ix
     rz(x1skt+m1-2+ik) = x(m-1+ik)
     rz(y1skt+m1-2+ik) = y(m-1+ik)
     rz(z1skt+m1-2+ik) = z(m-1+ik)
   END DO
   rz(z1skt-1+m1+ix) = -ISIGN(1,ix1skt)*9999.
   GO TO 20
 END DO
 50 CONTINUE
 DO  j = 1,li
   rz(x1skt-1+j) = x(j)
   rz(y1skt-1+j) = y(j)
   rz(z1skt-1+j) = z(j)
 END DO
 np = np + 1
 rz(z1skt-1+np) = -9999.
 70 CONTINUE
 CALL hdlin (rz(x1skt),rz(y1skt),rz(z1skt),np,nc,  &
     rz(xcc),iz(icount),iz(irct),rz(x21),rz(y21),rz(z21),  &
     iz(iia),rz(xe),rz(ye),rz(xu),rz(yu),rz(xi),rz(yi),rz(zi),  &
     rz(di),iz(ibeg),iz(iend),iz(ict),iz(icct),  &
     iz(w),iz(w),rz(w),rz(w),iz(w),iz(w),iz(w),rz(w),rz(w),rz(w),  &
     rz(w),rz(w),rz(w),rz(w),iz(w),iz(w),rz(w),rz(w),rz(w),rz(w), iz(w),iz(w))
 np = l
 
!     RESET VALUE FOR MAXIMUM NUMBER OF EDGES IF ARGUMENT IS COMPLETED.
 
 IF (vp > 0.) lz = lz/5
 RETURN
END SUBROUTINE hdsket
