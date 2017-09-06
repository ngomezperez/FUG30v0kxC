!================================================
! This subroutine calculates the right-hand sides
! of the hydrodynamic variables
!================================================
SUBROUTINE hydro_deriv()
  use types
  use comms
  use params
  use trans_params
  use dimens
  use control
  use constants
  use fft
  use dbydz
  use arrays
  use timing

  implicit none

!-------------------
! Local declarations
!-------------------
  integer       :: i,j,k
  real(kind=dp) :: z1, z2, z3, z4, rr, tt 
  real(kind=dp) :: u0, v0, w0
  real(kind=dp) :: temp1, temp2, temp3
  real(kind=dp) :: divu
  real(kind=dp) :: w1, w2, w3, w4

!-------------------------------
! Set wk1 to dr/dx, wk2 to dr/dy 
!-------------------------------
  call d1bydx(r, wk1 )
  call d1bydy(r, wk2 )

!-----------------
! Measure the time
!-----------------
  call findtime(time_tmp)

!--------------------------------------
! rdot = - u dr/dx - v dr/dy
! ex and ey are used as workspace
!--------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        rdot(i,j,k,3)  = - u(i,j,k)*wk1(i,j,k) - v(i,j,k)*wk2(i,j,k)
        ex(i,j,k) = -t(i,j,k)*wk1(i,j,k)
        ey(i,j,k) = -t(i,j,k)*wk2(i,j,k)
      end do
    end do
  end do

!------------------
! Increment t1_time
!------------------
  call findtime(time_tmp2)
  t1_time = t1_time + (time_tmp2-time_tmp)

!-------------------------------------------------------------------
! Set udot to tx, wk1 to txx, etc. Initially using udot as workspace
!-------------------------------------------------------------------
  call d2bydx( t, udot(1,1,1,3), wk1)
  call d2bydy( t, vdot(1,1,1,3), wk2)
  call d1upbydz( t, wdot(1,1,1,3), 5 )
  call d1bydz( t, wk4,  5 )
  call d2bydz( t, wk3,  0 )

!-------------
! Measure time
!-------------
  call findtime(time_tmp)
  call ini_trans_params()
!--------------------------------------------------
! Next find grad T for the U.grad T in the heat eqn
! and delsquared T for the head  eqn. 
! Also add viscous heating here
!--------------------------------------------------
! NB tdot is in fact rho tdot
!--------------------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        u0 = u(i,j,k)
        v0 = v(i,j,k)
        w0 = w(i,j,k)
        tdot(i,j,k,3)  = -r(i,j,k)*( u0*udot(i,j,k,3) +                  &
                                     v0*vdot(i,j,k,3) +                  &
                                     w0*wdot(i,j,k,3)      )             &
!=======================
!NGP Aug 31  
! adding the multiplying kappaZ(k) for varying thermal cond
!======================
                         + gck*kappaZ(k)*				 &
			     ( wk1(i,j,k) + wk2(i,j,k) + wk3(i,j,k) )    &
!=======================
!NGP Aug 31
! adding the term for varying thermal cond dkappa/dz*dTemp/dz
!======================
			 +gck*dkappadZ(k)*wdot(i,j,k,3)			 &
                         +fzckg*jsq(i,j,k)
        udot(i,j,k,3) = ex(i,j,k) - r(i,j,k)*udot(i,j,k,3)
        vdot(i,j,k,3) = ey(i,j,k) - r(i,j,k)*vdot(i,j,k,3)
        wdot(i,j,k,3) = -r(i,j,k)*wk4(i,j,k)
        ex(i,j,k) = 0.0d0
        ey(i,j,k) = 0.0d0
      end do
    end do
  end do

!------------------
! Increment t2_time
!------------------
  call findtime(time_tmp2)
  t2_time = t2_time + (time_tmp2 - time_tmp)

!------------------------------------------------------------
! First calculate derivatives
! wk1 = du/dx, wk2 = dv/dy, wk3 = dw/dz (upwind), wk4 = dw/dz
!------------------------------------------------------------
  call d1bydx(u, wk1)
  call d1bydy(v, wk2)
  call d1upbydz(w, wk3, 1)
  call d1bydz(w, wk4, 1)

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!------------------------------------------------------------
! Now find Ux, Vy and Wz for div U.
! Add the symmetric viscous heating terms to tdot.
! Need fo find Wz on boundaries here also
!------------------------------------------------------------
! cubic interpolation for dwdz on boundaries, az + bz^2 +cz^3
!------------------------------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        z1 = wk1(i,j,k)
        z2 = wk2(i,j,k)
        z3 = wk3(i,j,k)
        z4 = wk4(i,j,k)
        tt = t(i,j,k)
        rr = r(i,j,k)
        divu = z1 + z2 + z4
        rdot(i,j,k,3) = rdot(i,j,k,3) - rr*divu
        tdot(i,j,k,3) = tdot(i,j,k,3) - gm1*tt*divu*rr                 &
                      + 2*sckg*( z1*z1 + z2*z2 + z4*z4 - third*divu*divu)
        temp1 = fo2*bx(i,j,k)*bx(i,j,k)
        temp2 = fo2*by(i,j,k)*by(i,j,k)
        temp3 = fo2*bz(i,j,k)*bz(i,j,k)
        divu  = scko3*divu
        u0    = u(i,j,k)
        v0    = v(i,j,k)
        w0    = w(i,j,k)
        udot(i,j,k,3) = udot(i,j,k,3) - rr*u0*z1
        vdot(i,j,k,3) = vdot(i,j,k,3) - rr*v0*z2
        wdot(i,j,k,3) = wdot(i,j,k,3) - rr*w0*z3
        wk1(i,j,k) = divu + temp1 - temp2 - temp3
        wk2(i,j,k) = divu + temp2 - temp1 - temp3
        wk3(i,j,k) = scko3*z1 + scko3*z2 + temp3 - temp1 - temp2
        wk4(i,j,k) = 4.0d0*scko3*z4
      end do
    end do
  end do

!------------------
! Increment t3_time
!------------------
  call findtime(time_tmp2)
  t3_time = t3_time + (time_tmp2 - time_tmp)

!------------------------------------------
! Calculate grad terms in momentum equation
!------------------------------------------
! NB poldot etc used as workspace
!------------------------------------------
    call d1bydx( wk1, ex )
    call d1bydy( wk2, ey )
    call d1bydz( wk3, ez,  2 )
    call d1bydz( wk4, wk3, 5 )

!----------------------------------------------------------
! Find cross derivative Uy, Vx, for viscous heat and set up
! asymmetric gradient terms in momentum equations.
!----------------------------------------------------------
    call d1bydx( v, wk2 )
    call d1bydy( u, wk1 )

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!----------------
! Add derivatives
!----------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        udot(i,j,k,3) = udot(i,j,k,3) + ex(i,j,k)
        vdot(i,j,k,3) = vdot(i,j,k,3) + ey(i,j,k)
        wdot(i,j,k,3) = wdot(i,j,k,3) + ez(i,j,k) + wk3(i,j,k)
        ex(i,j,k) = 0.0d0
        ey(i,j,k) = 0.0d0
        ez(i,j,k) = 0.0d0
        rr = r(i,j,k)
        w1 = wk1(i,j,k)
        w2 = wk2(i,j,k)
        udot(i,j,k,3) = udot(i,j,k,3) - rr*v(i,j,k)*w1
        vdot(i,j,k,3) = vdot(i,j,k,3) - rr*u(i,j,k)*w2
        wk3(i,j,k) = f*by(i,j,k)*bx(i,j,k)
        tdot(i,j,k,3) = tdot(i,j,k,3) + sckg*( w1 + w2 )*( w1 + w2 )
      end do
    end do
  end do

!------------------
! Increment t4_time
!------------------
  call findtime(time_tmp2)
  t4_time = t4_time + (time_tmp2 - time_tmp)

!----------------------
! Calculate derivatives
!----------------------
  call d1bydx(wk3, wk2 )
  call d1bydy(wk3, wk1 )
  call d1bydx(w, wk3 )
  call d1upbydz(u, wk4, 2 )

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!---------------------------------------------------------------
! Do the 'cross' gradient terms set up above. Then find Uz for
! U eqn. Wx for W eqn. Funny order is because of viscous heating 
! term Uz*wx.
!---------------------------------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        udot(i,j,k,3) = udot( i,j, k,3) + wk1(i,j,k)
        vdot(i,j,k,3) = vdot( i,j, k,3) + wk2(i,j,k)
        rr = r(i,j,k)
        w3 = wk3(i,j,k)
        w4 = wk4(i,j,k)
        udot(i,j,k,3) = udot( i,j, k,3) - rr*w(i,j,k)*w4 
        wdot(i,j,k,3) = wdot( i,j, k,3) - rr*u(i,j,k)*w3
        temp1 = f*bz(i,j,k)*bx(i,j,k)
        wk1(i,j,k) = temp1 
        wk2(i,j,k) = temp1
        tdot(i,j,k,3) = tdot(i,j,k,3) + sckg*(w3+w4)*(w3+w4)
      end do
    end do
  end do

!------------------
! Increment t5_time
!------------------
  call findtime(time_tmp2)
  t5_time = t5_time + (time_tmp2 - time_tmp)

!----------------------
! Calculate derivatives
!----------------------
   call d1bydz(wk1, wk3, 5 )
   call d1bydx(wk2, wk4 )
   call d1bydy(w, wk1 )
   call d1upbydz(v, wk2, 2 )

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!-----------------------------------------------
! Do Vz, Wy and  remaining viscous heating terms
! Add gravity term to wdot.
! Add coriolis terms to udot, vdot and wdot
!-----------------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        rr = r(i,j,k)
        udot(i,j,k,3) = udot(i,j,k,3) + wk3(i,j,k)
        wdot(i,j,k,3) = wdot(i,j,k,3) + wk4(i,j,k) + thmplus1*rr
        w1 = wk1(i,j,k)
        w2 = wk2(i,j,k)
        vdot(i,j,k,3) = vdot(i,j,k,3) - rr*w(i,j,k)*w2
        wdot(i,j,k,3) = wdot(i,j,k,3) - rr*v(i,j,k)*w1 
        temp1 = f*bz(i,j,k)*by(i,j,k)
        wk3(i,j,k) = temp1
        wk4(i,j,k) = temp1 
        tdot(i,j,k,3) = tdot(i,j,k,3) + sckg*(w1 + w2)*(w1 + w2)
        udot(i,j,k,3) = udot(i,j,k,3)              & 
                      - rr*cor*( w(i,j,k)*cpsi + v(i,j,k)*spsi )    
        vdot(i,j,k,3) = vdot(i,j,k,3)              & 
                      - rr*cor*( - u(i,j,k)*spsi )    
        wdot(i,j,k,3) = wdot(i,j,k,3)              & 
                      - rr*cor*( - u(i,j,k)*cpsi )    
      end do
    end do
  end do

!------------------
! Increment t6_time
!------------------
  call findtime(time_tmp2)
  t6_time = t6_time + (time_tmp2 - time_tmp)

!----------------------
! Calculate derivatives
!----------------------
   call d1bydy( wk3, wk1 )
   call d1bydz( wk4, wk2, 5 )

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!------------------------------------------------------
! Add in last cross terms, plus z diffusion for u and v
!------------------------------------------------------
   do k=1,nlayer
     do j=1,ny
       do i=1,nx
         vdot(i,j,k,3) = vdot(i,j,k,3) + wk2(i,j,k)
         wdot(i,j,k,3) = wdot(i,j,k,3) + wk1(i,j,k)
       end do
     end do
   end do

!------------------
! Increment t7_time
!------------------
  call findtime(time_tmp2)
  t7_time = t7_time + (time_tmp2 - time_tmp)

!----------------------------------
! Z-derivatives for diffusion terms
!----------------------------------
  call d2bydz(u,wk1,3)
  call d2bydz(v,wk2,3)

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!-----------------------
! Add in diffusion terms
!-----------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        udot(i,j,k,3) = udot(i,j,k,3) + sck*wk1(i,j,k)
        vdot(i,j,k,3) = vdot(i,j,k,3) + sck*wk2(i,j,k)
      end do
    end do
  end do

!-------------------
! Increment t8_time
!-------------------
  call findtime(time_tmp2)
  t8_time = t8_time + (time_tmp2 - time_tmp)

!----------------------------
! Calculate 2nd x derivatives
!----------------------------
  call d2bydx(u,ex,wk1)
  call d2bydx(v,ey,wk2)
  call d2bydx(w,ez,wk3)

!-------------
! Measure time
!------------- 
  call findtime(time_tmp)

!--------------------
! Add derivative term
!--------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        udot(i,j,k,3) = udot(i,j,k,3) + sck*wk1(i,j,k)
        vdot(i,j,k,3) = vdot(i,j,k,3) + sck*wk2(i,j,k)
        wdot(i,j,k,3) = wdot(i,j,k,3) + sck*wk3(i,j,k)
      end do
    end do
  end do

!------------------------
! Increment t9_time term
!------------------------   
  call findtime(time_tmp2)
  t9_time = t9_time + (time_tmp2 - time_tmp)

!----------------------
! Calculate derivatives
!----------------------
   call d2bydy(u,ex,wk1)
   call d2bydy(v,ey,wk2)
   call d2bydy(w,ez,wk3)

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!----------------
! Add derivatives
!----------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        udot(i,j,k,3) = udot(i,j,k,3) + sck*wk1(i,j,k)
        vdot(i,j,k,3) = vdot(i,j,k,3) + sck*wk2(i,j,k)
        wdot(i,j,k,3) = wdot(i,j,k,3) + sck*wk3(i,j,k)
        ex(i,j,k) = 0.0d0
        ey(i,j,k) = 0.0d0
        ez(i,j,k) = 0.0d0
      end do
    end do
  end do

!-------------------
! Increment t10_time
!-------------------
  call findtime(time_tmp2)
  t10_time = t10_time + (time_tmp2 - time_tmp)

!-------------------------------
! Calculate z-derivatives of rho
!-------------------------------
  call d1upbydz(r,wk1,0)
  call d1bydz(r,wk2,0)

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!--------------------------------------------------------------------
! Boundary term for drho/dz is found by equating dw/dt to zero at top
! and bottom. 
! -------------------------------------------------------------------
  if (myrank .eq. 0) then
    do j=1,ny
      do i=1,nx
        wk1(i,j,1) = wdot(i,j,1,3)/t(i,j,1)
        wk2(i,j,1) = wdot(i,j,1,3)/t(i,j,1)
      end do
    end do
  end if
  if (myrank .eq. nproc-1) then
    do j=1,ny
      do i=1,nx
        wk1(i,j,nlayer) = wdot(i,j,nlayer,3)/t(i,j,nlayer)
        wk2(i,j,nlayer) = wdot(i,j,nlayer,3)/t(i,j,nlayer)
      end do
    end do
  end if

!------------------------------------
! Update rho and w with drho/dz terms
!------------------------------------
  do k=1,nlayer 
    do j=1,ny
      do i=1,nx
        rdot(i,j,k,3) = rdot(i,j,k,3) - w(i,j,k)*wk1(i,j,k)
        wdot(i,j,k,3) = wdot(i,j,k,3) - t(i,j,k)*wk2(i,j,k)
      end do
    end do
  end do

!--------------------------------------------------------
! Divide through by rho to finish off evolution equations
!--------------------------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        rr=1.0d0/(r(i,j,k)) 
        tdot(i,j,k,3) = tdot(i,j,k,3)*rr
        udot(i,j,k,3) = udot(i,j,k,3)*rr
        vdot(i,j,k,3) = vdot(i,j,k,3)*rr
        wdot(i,j,k,3) = wdot(i,j,k,3)*rr
      end do
    end do
  end do

!-------------------
! Increment t11_time
!------------------- 
  call findtime(time_tmp2)
  t11_time = t11_time + (time_tmp2 - time_tmp)

END SUBROUTINE hydro_deriv
