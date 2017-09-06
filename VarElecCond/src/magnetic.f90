!=============================================
! This subroutine calculates poloidal-toroidal
! potentials at the start of the code
!=============================================
SUBROUTINE mag_init()
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
  integer        :: i,j,k
  real(kind=dp)  :: seedfield,twopionx,twopiony

!---------------------------------------------------
! A good place to put in a dynamo field if necessary
!---------------------------------------------------
  if (dynamo) then
    if (myrank .eq. 0) then
      Print*,'Inserting seed field'
    end if
    seedfield = 0.003d0
    twopionx = 8.0d0*atan(1.0d0)/(nx*1.0d0)
    twopiony = 8.0d0*atan(1.0d0)/(ny*1.0d0)
    if (perfect) then
      do k=1,nlayer
        do j=1,ny
          do i=1,nx
            bx(i,j,k) = seedfield*cos(twopiony*j)
            by(i,j,k) = seedfield*cos(twopionx*i)
            bz(i,j,k) = 0.0d0
          end do
        end do
      end do
    else
      do k=1,nlayer
        do j=1,ny
          do i=1,nx
            bx(i,j,k) = 0.0d0
            by(i,j,k) = 0.0d0
            bz(i,j,k) = seedfield*cos(twopionx*i)*cos(twopiony*j)
          end do
        end do
      end do
    end if
  end if

!--------------------------------------
! Take the forward FFT of Bx, By and Bz
! Use ex, ey, ez as work arrays  
!--------------------------------------
  call forwardfft(bx,ex)
  call forwardfft(by,ey)
  call forwardfft(bz,ez)

!------------------------------------
! Dealias here - first in x direction
!------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,icutoffx
        ex((nx/2)+1+i,j,k) = 0.0d0
        ex((nx/2)+1-i,j,k) = 0.0d0
        ey((nx/2)+1+i,j,k) = 0.0d0
        ey((nx/2)+1-i,j,k) = 0.0d0
        ez((nx/2)+1+i,j,k) = 0.0d0
        ez((nx/2)+1-i,j,k) = 0.0d0
      end do
      ex((nx/2)+1,j,k) = 0.0d0
      ey((nx/2)+1,j,k) = 0.0d0
      ez((nx/2)+1,j,k) = 0.0d0
    end do
  end do

!-----------------------
! Dealias in y direction
!-----------------------
  do k=1,nlayer
    do i=1,nx
      ex(i,(ny/2)+1,k) = 0.0d0
      ey(i,(ny/2)+1,k) = 0.0d0
      ez(i,(ny/2)+1,k) = 0.0d0
      do j=1,icutoffy
        ex(i,(ny/2)+1+j,k) = 0.0d0
        ex(i,(ny/2)+1-j,k) = 0.0d0
        ey(i,(ny/2)+1+j,k) = 0.0d0
        ey(i,(ny/2)+1-j,k) = 0.0d0        
        ez(i,(ny/2)+1+j,k) = 0.0d0
        ez(i,(ny/2)+1-j,k) = 0.0d0     
      end do
    end do
  end do

!---------------------------------------------
! Use the FFT to work out the mean values
!---------------------------------------------
! Note: bzbar should be independent of z and t
!       Also need to be normalised by nx, ny
!---------------------------------------------
  do k=1, nlayer
    bxbar(k) = ex(1,1,k)
    bybar(k) = ey(1,1,k)
  end do
  bxbar = bxbar*xn
  bybar = bybar*xn
  bzbar = ez(1,1,1)*xn

!------------------------------------------------------
! Work out pol by inverting FFT(Bz) = FFT(-del_H^2 pol) 
! DC component set to zero
!------------------------------------------------------
  do k=1, nlayer
    pol(1,1,k) = 0.0d0
    do i=2,nx
      pol(i,1,k) = -1.0d0*ez(i,1,k)/akx2(i)
    end do
    do j=2,ny
      pol(1,j,k) = -1.0d0*ez(1,j,k)/aky2(j)
    end do
    do j=2,ny
      do i=2,nx
        pol(i,j,k) = -1.0d0*ez(i,j,k)/(akx2(i)+aky2(j))
      end do
    end do
  end do

!-------------------------
! Store FFT(dBy/dx) in wk3
!-------------------------
  do k=1,nlayer
    do j=1,ny
      wk3(1,j,k)=0.0d0
      wk3((nx/2)+1,j,k)=0.0d0
      do i=1,nx/2-1
        wk3(i+1,j,k)   = -ey(nx+1-i,j,k)*akx1(i+1)
        wk3(nx+1-i,j,k) =  ey(i+1,j,k)*akx1(nx+1-i)
      end do
    end do
  end do

!-------------------------
! Store FFT(dBx/dy) in wk4
!-------------------------
  do k=1,nlayer
    do j=1,ny/2-1
      do i=1,nx
        wk4(i,1,k)=0.0d0
        wk4(i,(ny/2)+1,k)=0.0d0
        wk4(i,j+1,k)   = -ex(i,ny+1-j,k)*aky1(j+1)  
        wk4(i,ny+1-j,k) =  ex(i,j+1,k)*aky1(ny+1-j) 
      end do
    end do
  end do

!-----------------------------------
! Work out wk3 = FFT(jz) = wk3 - wk4
!-----------------------------------
  do k=1, nlayer
    do j=1, ny
      do i=1, nx
        wk3(i,j,k)=wk3(i,j,k)-wk4(i,j,k)
      end do
    end do
  end do

!------------------------------------------------------
! Work out tor by inverting FFT(jz) = FFT(-del_H^2 tor) 
! DC component set to zero
!------------------------------------------------------
  do k=1, nlayer
    tor(1,1,k) = 0.0d0
    do i=2,nx
      tor(i,1,k) = -1.0d0*wk3(i,1,k)/akx2(i)
    end do
    do j=2,ny
      tor(1,j,k) = -1.0d0*wk3(1,j,k)/aky2(j)
    end do
    do j=2,ny
      do i=2,nx
        tor(i,j,k) = -1.0d0*wk3(i,j,k)/(akx2(i)+aky2(j))
      end do
    end do
  end do

!------------------------------
! Work out wk3=FFT(jz) from tor
!------------------------------
  do k=1, nlayer
    do j=1, ny
      do i=1, nx
        wk3(i,j,k)=-1.0d0*tor(i,j,k)*(akx2(i)+aky2(j)) 
      end do
    end do
  end do

!---------------------------------
! Work out dBz/dy in Fourier space
!---------------------------------
  do k=1,nlayer
    do j=1,ny/2-1
      do i=1,nx
        wk1(i,1,k)=0.0d0
        wk1(i,(ny/2)+1,k)=0.0d0
        wk1(i,j+1,k)   = -ez(i,ny+1-j,k)*aky1(j+1)  
        wk1(i,ny+1-j,k) =  ez(i,j+1,k)*aky1(ny+1-j) 
      end do
    end do
  end do

!---------------------------------
! Work out dBz/dx in Fourier space
!---------------------------------
  do k=1,nlayer
    do j=1,ny
      wk2(1,j,k)=0.0d0
      wk2((nx/2)+1,j,k)=0.0d0
      do i=1,nx/2-1
        wk2(i+1,j,k)   = -ez(nx+1-i,j,k)*akx1(i+1)
        wk2(nx+1-i,j,k) =  ez(i+1,j,k)*akx1(nx+1-i)
      end do
    end do
  end do  

!-------------------------------------------------
! Finished in Fourier space
! Invert wk1, wk2 and wk3 - store in ex, ey and ez
!-------------------------------------------------
  call inversefft(wk1,ex)
  call inversefft(wk2,ey)
  call inversefft(wk3,ez)

!---------------------------
! Work out dBy/dz and dBx/dz
!---------------------------
  call d1bydz(by,wk1,5)
  call d1bydz(bx,wk2,5)

!----------------------
! Calculate J^2 --> jsq
!----------------------
  jsq = 0.0d0
  do k=1, nlayer 
    do j=1, ny
      do i=1, nx
        jsq(i,j,k) = ez(i,j,k)*ez(i,j,k)                              &       
                   + ((ex(i,j,k)-wk1(i,j,k))*(ex(i,j,k)-wk1(i,j,k)))  &
                   + ((wk2(i,j,k)-ey(i,j,k))*(wk2(i,j,k)-ey(i,j,k)))
      end do
    end do
  end do  

!--------------------------
! Check boundary conditions
!--------------------------
  call mag_boundary()

!------------------
! Reset work arrays
!------------------
  ex = 0.0d0
  ey = 0.0d0
  ez = 0.0d0
  wk1 = 0.0d0
  wk2 = 0.0d0
  wk3 = 0.0d0
  wk4 = 0.0d0

END SUBROUTINE mag_init

!================================================
! This subroutine calculates the right-hand sides
! of the magnetic variables
!================================================
SUBROUTINE mag_deriv()
  use types
  use comms
  use params
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
  integer        :: i,j,k

!-------------------------------------
! Work out components of Ex, Ey and Ez
! Store in wk1, wk2, wk3
!-------------------------------------
  do k=1, nlayer
    do j=1, ny
      do i=1, nx
        wk1(i,j,k) = v(i,j,k)*bz(i,j,k)-w(i,j,k)*by(i,j,k)
        wk2(i,j,k) = w(i,j,k)*bx(i,j,k)-u(i,j,k)*bz(i,j,k)
        wk3(i,j,k) = u(i,j,k)*by(i,j,k)-v(i,j,k)*bx(i,j,k)
      end do
    end do
  end do

!-----------------------------------
! Take ex, ey, ez into Fourier space
!-----------------------------------
  call forwardfft(wk1,ex)
  call forwardfft(wk2,ey)
  call forwardfft(wk3,ez)

!-----------------------------------------------------
! Copy DC components of ex and ey to wk1_1D and wk2_1D
!-----------------------------------------------------
  do k=1, nlayer
    wk1_1D(k) = xn*ex(1,1,k)
    wk2_1D(k) = xn*ey(1,1,k)
  end do

!--------------------------------
! Find z derivatives of mean emfs
! Store in bxbardot and bybardot
!--------------------------------
  call d1bydz1D(wk2_1D,bxbardot(1,3),5)
  call d1bydz1D(wk1_1D,bybardot(1,3),5)

!-------------------------------------------
! Find second derivatives of bxbar and bybar
!-------------------------------------------
  call d2bydz1D(bxbar,wk1_1D,3)
  call d2bydz1D(bybar,wk2_1D,3)

!---------------------------
! Calculate mean derivatives
!---------------------------
  do k=1, nlayer
    bxbardot(k,3) = zck*wk1_1D(k) - bxbardot(k,3) 
    bybardot(k,3) = zck*wk2_1D(k) + bybardot(k,3)
  end do

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!------------------------------------------------
! Take x derivative of ex and ey in Fourier space
!------------------------------------------------
  do k=1,nlayer
    do j=1,ny
      wk1(1,j,k)=0.0d0
      wk1((nx/2)+1,j,k)=0.0d0
      wk2(1,j,k)=0.0d0
      wk2((nx/2)+1,j,k)=0.0d0
      do i=1,nx/2-1
        wk1(i+1,j,k)   = -ex(nx+1-i,j,k)*akx1(i+1)
        wk1(nx+1-i,j,k) =  ex(i+1,j,k)*akx1(nx+1-i)
        wk2(i+1,j,k)   = -ey(nx+1-i,j,k)*akx1(i+1)
        wk2(nx+1-i,j,k) =  ey(i+1,j,k)*akx1(nx+1-i)
      end do
    end do
  end do

!------------------------------------------------
! Take y derivative of ex and ey in Fourier space
!------------------------------------------------
  do k=1,nlayer
    do j=1,ny/2-1
      do i=1,nx
        wk3(i,1,k)=0.0d0
        wk3(i,(ny/2)+1,k)=0.0d0
        wk3(i,j+1,k)   = -ex(i,ny+1-j,k)*aky1(j+1)  
        wk3(i,ny+1-j,k) =  ex(i,j+1,k)*aky1(ny+1-j) 
        wk4(i,1,k)=0.0d0
        wk4(i,(ny/2)+1,k)=0.0d0
        wk4(i,j+1,k)   = -ey(i,ny+1-j,k)*aky1(j+1)  
        wk4(i,ny+1-j,k) =  ey(i,j+1,k)*aky1(ny+1-j) 
     end do
    end do
  end do

!--------------------------------------------------
! Combine these derivatives for evolution equations
! wk2 = dex/dy - dey/dx (in Fourier space)
! wk1 = dex/dx + dey/dy (in Fourier space)
!--------------------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        wk2(i,j,k) = wk3(i,j,k) - wk2(i,j,k)
        wk1(i,j,k) = wk1(i,j,k) + wk4(i,j,k)
      end do
    end do
  end do

!-------------------
! Increment t12_time
!------------------- 
  call findtime(time_tmp2)
  t12_time = t12_time + (time_tmp2 - time_tmp)


!------------------------------------------------------------------
! Work out wk3=d/dz(wk1) = d/dz(dex/dx + dey/dy) (in Fourier space)
!------------------------------------------------------------------
  call d1bydz(wk1,wk3,5)

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!-------------------------------------
! Take inverse(del_H^2) of wk3 and wk2
!-------------------------------------
  do k=1, nlayer
    wk2(1,1,k) = 0.0d0
    wk3(1,1,k) = 0.0d0
    do i=2,nx
      wk2(i,1,k) = wk2(i,1,k)/akx2(i)
      wk3(i,1,k) = wk3(i,1,k)/akx2(i)      
    end do
    do j=2,ny
      wk2(1,j,k) = wk2(1,j,k)/aky2(j)
      wk3(1,j,k) = wk3(1,j,k)/aky2(j)
    end do
    do j=2,ny
      do i=2,nx
        wk2(i,j,k) = wk2(i,j,k)/(akx2(i)+aky2(j))
        wk3(i,j,k) = wk3(i,j,k)/(akx2(i)+aky2(j))
      end do
    end do
  end do

!------------------------
! Increment t13_time term
!------------------------   
  call findtime(time_tmp2)
  t13_time = t13_time + (time_tmp2 - time_tmp)

!-----------------------------------------
! Work out wk1=d2pol/dz2 and wk4=d2tor/dz2
!-----------------------------------------
   call d2bydz(pol,wk1,3)
   call d2bydz(tor,wk4,3)

!-------------
! Measure time
!-------------
  call findtime(time_tmp)

!-------------------------------------------
! Store del_H^2 S and T in poldot and tordot
! Then complete update of poldot and tordot
!-------------------------------------------

  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        poldot(i,j,k,3) = zck*pol(i,j,k)*(akx2(i)+aky2(j))   &
                        + zck*wk1(i,j,k)+wk2(i,j,k)
        tordot(i,j,k,3) = zck*tor(i,j,k)*(akx2(i)+aky2(j))   &
                        + zck*wk4(i,j,k)-wk3(i,j,k)+ez(i,j,k)
      end do
    end do
  end do   

!------------------------
! Increment t14_time term
!------------------------   
  call findtime(time_tmp2)
  t14_time = t14_time + (time_tmp2 - time_tmp)

END SUBROUTINE mag_deriv

!====================
! Calculate B and J^2
!====================
SUBROUTINE bandj()
  use types
  use comms
  use params
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
  integer        :: i,j,k

!-------------
! Measure time
!------------- 
  call findtime(time_tmp)

!-------------------
! Set wk1 to dpol/dz
!-------------------
  call d1bydz(pol,wk1,5)

!---------------------------------------
! Set ex to d^2pol/dxdz in Fourier space
! Set wk3 to dtor/dx in Fourier space 
!---------------------------------------
  do k=1,nlayer
    do j=1,ny
      ex(1,j,k)=0.0d0
      ex((nx/2)+1,j,k)=0.0d0
      wk3(1,j,k)=0.0d0
      wk3((nx/2)+1,j,k)=0.0d0
      do i=1,nx/2-1
        ex(i+1,j,k)   = -wk1(nx+1-i,j,k)*akx1(i+1)
        ex(nx+1-i,j,k) =  wk1(i+1,j,k)*akx1(nx+1-i)
        wk3(i+1,j,k)   = -tor(nx+1-i,j,k)*akx1(i+1)
        wk3(nx+1-i,j,k) =  tor(i+1,j,k)*akx1(nx+1-i)
      end do
    end do
  end do  

!---------------------------------------
! Set ey to d^2pol/dydz in Fourier space
! Set wk2 to dtor/dy in Fourier space
!---------------------------------------
  do k=1,nlayer
    do j=1,ny/2-1
      do i=1,nx
        ey(i,1,k)=0.0d0
        ey(i,(ny/2)+1,k)=0.0d0
        ey(i,j+1,k)   = -wk1(i,ny+1-j,k)*aky1(j+1)  
        ey(i,ny+1-j,k) =  wk1(i,j+1,k)*aky1(ny+1-j) 
        wk2(i,1,k)=0.0d0
        wk2(i,(ny/2)+1,k)=0.0d0
        wk2(i,j+1,k)   = -tor(i,ny+1-j,k)*aky1(j+1)  
        wk2(i,ny+1-j,k) =  tor(i,j+1,k)*aky1(ny+1-j)
      end do
    end do
  end do

!---------------------------------------------------------
! Set ex=FFT(Bx), ey=FFT(By), ez=FFT(Bz)=FFT(-del_H^2 pol)
!---------------------------------------------------------
  do k=1,nlayer
    do j=1,ny
      do i=1,nx
        ex(i,j,k) = ex(i,j,k) + wk2(i,j,k)
        ey(i,j,k) = ey(i,j,k) - wk3(i,j,k)
        ez(i,j,k) = -1.0d0*pol(i,j,k)*(akx2(i)+aky2(j)) 
      end do
    end do
  end do

!-------------------------------------
! Evaluate Bx, By and Bz in real space
!-------------------------------------
  call inversefft(ex,bx)
  call inversefft(ey,by)
  call inversefft(ez,bz)

!------------------------------
! Add in mean parts of B
!------------------------------
! Work out wk3=FFT(jz) from tor
!------------------------------
  do k=1, nlayer
    do j=1, ny
      do i=1, nx
        bx(i,j,k) = bx(i,j,k)+bxbar(k)
        by(i,j,k) = by(i,j,k)+bybar(k)
        bz(i,j,k) = bz(i,j,k)+bzbar
        wk3(i,j,k)=-1.0d0*tor(i,j,k)*(akx2(i)+aky2(j)) 
      end do
    end do
  end do

!---------------------------------
! Work out dBz/dy in Fourier space
!---------------------------------
! Work out dBz/dx in Fourier space
!---------------------------------
  do k=1,nlayer
    do j=1,ny/2-1
      do i=1,nx
        wk1(i,1,k)=0.0d0
        wk1(i,(ny/2)+1,k)=0.0d0
        wk1(i,j+1,k)   = -ez(i,ny+1-j,k)*aky1(j+1)  
        wk1(i,ny+1-j,k) =  ez(i,j+1,k)*aky1(ny+1-j) 
      end do
    end do
    do j=1,ny
      wk2(1,j,k)=0.0d0
      wk2((nx/2)+1,j,k)=0.0d0
      do i=1,nx/2-1
        wk2(i+1,j,k)   = -ez(nx+1-i,j,k)*akx1(i+1)
        wk2(nx+1-i,j,k) =  ez(i+1,j,k)*akx1(nx+1-i)
      end do
    end do
  end do

!-------------------------------------------------
! Finished in Fourier space
! Invert wk1, wk2 and wk3 - store in ex, ey and ez
!-------------------------------------------------
  call inversefft(wk1,ex)
  call inversefft(wk2,ey)
  call inversefft(wk3,ez)

!---------------------------
! Work out dBy/dz and dBx/dz
!---------------------------
  call d1bydz(by,wk1,5)
  call d1bydz(bx,wk2,5)

!----------------------
! Calculate J^2 --> jsq
!----------------------
  jsq = 0.0d0
  do k=1, nlayer 
    do j=1, ny
      do i=1, nx
        jsq(i,j,k) = ez(i,j,k)*ez(i,j,k)                               &       
                   + ((ex(i,j,k)-wk1(i,j,k))*(ex(i,j,k)-wk1(i,j,k)))  &
                   + ((wk2(i,j,k)-ey(i,j,k))*(wk2(i,j,k)-ey(i,j,k)))
      end do
    end do
  end do  

!------------------
! Reset work arrays
!------------------
  ex = 0.0d0
  ey = 0.0d0
  ez = 0.0d0
  wk1 = 0.0d0
  wk2 = 0.0d0
  wk3 = 0.0d0
  wk4 = 0.0d0

!------------------------
! Increment t16_time term
!------------------------   
  call findtime(time_tmp2)
  t16_time = t16_time + (time_tmp2 - time_tmp)

END SUBROUTINE bandj

SUBROUTINE magdealias()
  use types
  use comms
  use params
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
  integer        :: i,j,k

!------------------------------------                                         
! Dealias here - first in x direction                                          
!------------------------------------                                          
  do k=1,nlayer
    do j=1,ny
      do i=1,icutoffx
        pol((nx/2)+1+i,j,k) = 0.0d0
        pol((nx/2)+1-i,j,k) = 0.0d0
        tor((nx/2)+1+i,j,k) = 0.0d0
        tor((nx/2)+1-i,j,k) = 0.0d0
      end do
      pol((nx/2)+1,j,k) = 0.0d0
      tor((nx/2)+1,j,k) = 0.0d0
    end do
  end do

!-----------------------
! Dealias in y direction                                                     
!-----------------------                                                      
  do k=1,nlayer
    do i=1,nx
      pol(i,(ny/2)+1,k) = 0.0d0
      tor(i,(ny/2)+1,k) = 0.0d0
      do j=1,icutoffy
        pol(i,(ny/2)+1+j,k) = 0.0d0
        pol(i,(ny/2)+1-j,k) = 0.0d0
        tor(i,(ny/2)+1+j,k) = 0.0d0
        tor(i,(ny/2)+1-j,k) = 0.0d0
      end do
    end do
  end do

END SUBROUTINE magdealias
