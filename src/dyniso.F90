!========================================================================
!
!  dyniso:     Fourier spectral, homogeneous, isotropic turbulence code 
!              with dynamic subgrid-scale model.
!
!              This version supports parallel execution on SGI platforms
!              and using the threaded version of FFTW 2.1.5
!   
!  Author:     S. Scott Collis
!
!  Written:    6-1-98
!
!  Copyright:  S. Scott Collis
!
!  9-18-98     Added constant CFL time steping (negative dt = cfl)
!  9-21-98     Added Rogollo statistics output in subroutine spectra
!              Added the RANF flag to link to the vectoral ranf function
!  9-22-98     Fixed restart problem, don't do FFT-1, FFT if possible
!  9-24-98     Added dynamic model
!  9-27-02     Ported to Linux + FFTW
!  9-27-02     Added OpenMP capability (there are problems with inline
!              format statement for OpenMP in the IFC)
!  2-09-12     Added threaded FFTW ability and some bug fixes
!========================================================================
! Dyniso:  DNS and LES of homogeneous, isotropic turbulence
! Copyright (C) 2002, 2012 S. Scott Collis 
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!    * Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!    * Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in
!      the documentation and/or other materials provided with the
!      distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!========================================================================

module intmod
      logical :: halt
end module intmod

module field
      integer :: nx, ny, nz, ndof, ndofles
      integer :: mx, my, mz
      integer :: px, py, pz, mpx
      integer :: nstep, istep, nout, scheme
      real :: nu, dt, t, cfl, velmax
      real :: lx, ly, lz
      real :: dx, dy, dz
      real, target, allocatable :: x(:), y(:), z(:)
      real, target, allocatable :: kx(:), ky(:), kz(:), kx2(:)
      real, target, allocatable :: pkx(:), pky(:), pkz(:)
      real, target, allocatable :: u(:,:,:,:), r(:,:,:,:), f(:,:,:,:)
      real, target, allocatable :: coef(:), coefp(:)
      integer :: iprint, ictype, dealias, stats

!.... Random number seed

#if   defined(CRAY)
      integer :: seed
#elif defined(RANF)
      integer(4) :: seed
#else
      integer :: seed
#endif

!.... Initial spectrum parameters

      integer :: ispec
      real :: kmax, v0

!.... LES parameters and storage

      integer :: LES   
      real, allocatable :: C(:,:,:), uh(:,:,:,:), kernel(:,:,:)
      real, allocatable :: Sijh(:,:,:,:), Lij(:,:,:,:), Mij(:,:,:,:)
      real, allocatable :: num(:,:,:), den(:,:,:)
      real :: Cs, delta_sq, ratio, filter_radius, filter_radius_sq, &
              alpha, alpha_sq, truncation_radius_sq
      integer :: filter_type, contraction, average, grid_filter

end module field

module timemod
#ifdef CRAY
      real, external :: second
      real :: cpu0, cpu, cpu2
#else
      real(4), external :: second
      real(4) :: cpu0, cpu, cpu2
#endif
      integer :: nthreads
end module timemod

module const
      real, parameter    :: zero   = 0.0
      real, parameter    :: pt25   = 0.25
      real, parameter    :: pt33   = 0.333333333333333333
      real, parameter    :: pt5    = 0.5
      real, parameter    :: one    = 1.0
      real, parameter    :: onept5 = 1.5
      real, parameter    :: two    = 2.0
      real, parameter    :: three  = 3.0
      real, parameter    :: pi     = 3.14159265358979323846
      complex, parameter :: czero  = (0.0,0.0)
      complex, parameter :: cone   = (1.0,0.0)
      complex, parameter :: iota   = (0.0,1.0)
end module const

program dyniso
!=============================================================================
!.... Simulate the decay of isotropic, homogeneous turbulence using a Fourier
!.... spectral method with anti-aliasing and the dynamic sub-grid scale model.
!
!.... Author:  S. Collis
!
!.... Date:    10-27-97
!
!============================================================================
      use field
      use const
      use timemod
      implicit none

!$    integer, external :: omp_get_num_threads
!$    integer, external :: omp_get_thread_num

      nthreads = 1

!$omp parallel
!$    if (omp_get_thread_num() == 0) then
!$      write(*,10) 'Running on ',omp_get_num_threads(), ' processors.'
!$    end if
!$    nthreads = omp_get_num_threads()
!$omp end parallel
  10  format(a,i4,a)

      call interupt()

      call setup()

      if (scheme.eq.1) then
        call euler(u, r)     ! Explicit Euler (debugging only)
      else if (scheme.eq.2) then
        call rk2(u, r)       ! RK2
      else if (scheme.eq.3) then
        call rk3(u, r)       ! RK3
      else
        call error('main$','Illegal value for scheme$')
      end if

      call post(u)

      stop
end program dyniso

subroutine setup
      use field
      use const
      use timemod
      implicit none
      integer :: i, j, k, idof, ier, i2

      namelist /in/ nx, ny, nz, Lx, Ly, Lz, scheme, nstep, nout, dt, nu, &
                    dealias, iprint, ictype, seed, kmax, v0, ispec, LES, &
                    Cs, stats

      namelist /dynamic/ average, contraction, filter_type, ratio, &
                         grid_filter

!.... start the clock

      cpu  = second()
      cpu0 = cpu

!.... set default input values

      nx = 32; ny = 32; nz = 32; nstep = 1; dt = 0.004; dealias = 1
      iprint = 0; Lx = -one; Ly = -one; Lz = -one; nu = 0.01189
      ndof = 3; ndofles = 8; ictype = 0; seed = 0; kmax = 4.756828460010884; 
      v0 = one; ispec = 1; nout = 50; LES = 1; Cs = 0.16124515
      scheme = 2; cfl=0.0; stats=1

!.... get namelist input

#ifdef VERBOSE
      write(*,in)
#endif
      read(*,in)
#ifdef VERBOSE
      write(*,in)
#endif

!.... open the output files

      if (ictype.eq.2) then
        open(11,file='cpu.dat',position='append',err=100)
        if (LES.eq.2) open(12,file='smag.dat',position='append',err=100)
        open(20,file='stat1.dat',position='append',err=100)
        open(21,file='stat2.dat',position='append',err=100)
        open(22,file='screen.dat',position='append',err=100)
      else
        open(11,file='cpu.dat',err=100)
        if (LES.eq.2) open(12,file='smag.dat',err=100)
        open(20,file='stat1.dat',err=100)
        open(21,file='stat2.dat',err=100)
        open(22,file='screen.dat',err=100)
      end if

!.... set time-step or CFL

      if (dt.lt.zero) then
        cfl = abs(dt)
      else
        cfl = zero
      end if

!.... sanity check

      if (mod(nx,2).ne.0.0) call error('setup$','Nx must be power of 2$')
      if (mod(ny,2).ne.0.0) call error('setup$','Ny must be power of 2$')
      if (mod(nz,2).ne.0.0) call error('setup$','Nz must be power of 2$')

!.... oddball wavenumbers

      mx = (nx+2)/2
      my = (ny+2)/2
      mz = (nz+2)/2

!.... Domain size logic
      
      if (Lx .lt. zero) then
        Lx = abs(Lx) * pi
      end if

      if (Ly .lt. zero) then
        Ly = abs(Ly) * Pi
      end if

      if (Lz .lt. zero) then
        Lz = abs(Lz) * pi
      end if

!.... Use 3/2 rule for anti-aliasing
      
      if (dealias.eq.1) then
        px = onept5 * nx
        py = onept5 * ny
        pz = onept5 * nz
      else
        px = nx
        py = ny
        pz = nz
      end if

      mpx = (px+2)/2
      truncation_radius_sq = 2.0/9.0

      write(*,10) nx,ny,nz
      write(*,20) px,py,pz
 10   format('Nx, Ny, Nz = (',i3,',',i3,',',i3,')')
 20   format('Px, Py, Pz = (',i3,',',i3,',',i3,')')

!.... make the mesh and initialize the wavenumbers

      allocate( x(nx), y(ny), z(nz), kx(mx), ky(ny), kz(nz), &
                kx2(2*mx), stat=ier )
      if (ier .ne. 0) call error('setup$','Insufficient Memory for mesh$')

      dx = Lx / real(nx)
      do i = 1, nx
        x(i) = (i-1)*dx
      end do
      i2 = 1
      do i = 1, mx
        kx(i) = two * pi / Lx * real(i-1)
        kx2(i2) = kx(i)
        i2 = i2 + 1
        kx2(i2) = kx(i)
        i2 = i2 + 1
      end do

      dy = Ly / real(ny)
      do j = 1, ny
        y(j) = (j-1)*dy
      end do
      ky(1) = zero
      do j = 2, (ny+2)/2
        ky(ny+2-j) = -two * pi / Ly * real(j-1)
        ky(j) = two * pi / Ly * real(j-1)
      end do

      dz = Lz / real(nz)
      do k = 1, nz
        z(k) = (k-1)*dz
      end do
      kz(1) = zero
      do k = 2, (nz+2)/2
        kz(nz+2-k) = -two * pi / Lz * real(k-1)
        kz(k) = two * pi / Lz * real(k-1)
      end do

!.... wave numbers for padded field

      allocate( pkx(mpx), pky(py), pkz(pz), stat=ier )
      if (ier .ne. 0) call error('setup$','Insufficient Memory for mesh$')

      do i = 1, mpx
        pkx(i) = two * pi / Lx * real(i-1)
      end do

      pky(1) = zero
      do j = 2, (py+2)/2
        pky(py+2-j) = -two * pi / Ly * real(j-1)
        pky(j) = two * pi / Ly * real(j-1)
      end do

      pkz(1) = zero
      do k = 2, (pz+2)/2
        pkz(pz+2-k) = -two * pi / Lz * real(k-1)
        pkz(k) = two * pi / Lz * real(k-1)
      end do

!.... initialize the FFT coefficients

#if   defined(CRAY)
      call mkwork( px, py, pz )
      allocate( coef(100 + 2*(nx + ny + nz)), &
               coefp(100 + 2*(px + py + pz)), stat=ier )
#elif defined(SGI_FFT3D) || defined(SGI_FFT)
      allocate( coef((nx+15)+2*(ny+15)+2*(nz+15)), &
               coefp((px+15)+2*(py+15)+2*(pz+15)), stat=ier )
#elif defined(INTEL_MKL)
      if (dealias .eq. 1) &
        call error('setup$','FFT must be power of 2 for MKL (no dealias)$')
      allocate(  coef(2*nx+4 + 3*ny + 3*nz ), &
                coefp(2*px+4 + 3*py + 3*pz ), stat=ier )
#elif defined(FFTW)
      allocate( coef(1), coefp(1), stat=ier )
      coef(1) = 1
      coefp(1) = 2
#endif
      if (ier .ne. 0) call error('setup$','Insufficient Memory for FFT$')
      call fft3di(nx, ny, nz, coef)
      call fft3di(px, py, pz, coefp)

!.... allocate memory for flow variables

      allocate( u(2*mx,ny,nz,ndof), r(2*mx,ny,nz,ndof), &
                f(2*mpx,py,pz,ndofles), stat=ier )
      if (ier .ne. 0) call error('setup$','Insufficient Memory for fields$')

!.... initialize data for first-touch memory allocation on O2k, otherwise
!.... its a good idea to initialize memory anyway!

      !$omp parallel do private(k,j,idof,i)
      do k = 1, nz
        do j = 1, ny
          do idof = 1,ndof
            do i = 1, 2*mx
              u(i,j,k,idof) = zero
            end do
          end do
        end do
      end do

      !$omp parallel do private(k,j,idof,i)
      do k = 1, nz
        do j = 1, ny
          do idof = 1,ndof
            do i = 1, 2*mx
              r(i,j,k,idof) = zero
            end do
          end do
        end do
      end do

      !$omp parallel do private(k,j,idof,i)
      do k = 1, pz
        do j = 1, py
          do idof = 1,ndofles
            do i = 1, 2*mpx
              f(i,j,k,idof) = zero
            end do
          end do
        end do
      end do

!.... set the initial condition

      t = zero
      istep = 0
      if (ictype.eq.0) then
        call taylor()
      else if (ictype.eq.1) then
        call turb(u,f)
      else if (ictype.eq.2) then
        call restart(u,f)
      else
        call error('setup$','Illegal value for ictype$')
      end if

!.... setup for LES

!     delta_sq = (dx**2 + dy**2 + dz**2)
      delta_sq = ((dx*dy*dz)**pt33)**2

      if (LES.gt.0) then
        allocate( C(2*mpx,py,pz) )
        if (LES.eq.1) then
          C = Cs
        else if (LES.eq.2) then
          average = 2; contraction = 1; filter_type = 0; ratio = pt5;
          grid_filter = 1
          write(*,dynamic)
          read(*,dynamic)
          write(*,dynamic)
          filter_radius = ratio * sqrt(2.0)/3.0
          filter_radius_sq = filter_radius**2
          alpha_sq = one / ratio**2
          alpha = sqrt( alpha_sq )
          allocate( uh(2*mpx,py,pz,ndof), Sijh(2*mpx,py,pz,6), &
                    Lij(2*mpx,py,pz,6), Mij(2*mpx,py,pz,6), &
                    num(2*mpx,py,pz), den(2*mpx,py,pz), kernel(mpx,py,pz) )
          if (filter_type .eq. 3) alpha_sq = 1.5 * alpha_sq
          call get_kernel
        else
          call error('setup$','Illegal value of LES$')
        end if
      else if (LES.lt.0) then
        call error('setup$','Illegal value of LES$')
      end if

      cpu2 = second()
      write(*,30) (cpu2 - cpu)/nthreads
 30   format('Initialize time   = ',1pe13.7)
      cpu = cpu2

      return
100   call error('setup$','Error opening output files$')
end subroutine setup

subroutine turb(ul, fl)

!.... Generate a real, incompressible, homogeneous, isotropic 
!.... field of turbulence using Alan Wray's method

      use field
      use const
      implicit none
      complex :: ul(mx,ny,nz,ndof), fl(mx,ny,nz,ndof)
      integer :: i, j, k, idof, im, jm, km, n, itmp, ijk
      real :: Rlambda, kmag, fact, avg, k12, theta1, theta2, phi
      complex :: alp, beta
      logical, external :: mask
      real, external :: E0, I_E0_dk, I_E0k2_dk  ! spectrum functions
#ifdef CRAY
      real(8), external :: ranf
#else
#ifdef RANF
      real(8), external :: ranf
#else
      real(8), external :: rand
#endif
#endif
      integer :: idum
      real    :: rdum
      if (seed.ne.0) then
#ifdef CRAY
        call ranset( seed )
#else
#ifdef RANF
        call ranset( seed )
#else
        call srand( seed )
#endif
#endif
      end if

!.... if kmax=2 then Rlambda = 1/nu

      if (ispec.eq.1) then
        Rlambda = I_E0_dk() * sqrt(20.0/3.0) / sqrt(two * I_E0k2_dk()) / nu
        write(*,"('R_lambda   = ',1pe13.6)") Rlambda
      end if

!.... output the input spectrum
      
      open(10,file="spec.in")
      do i = 1, nx/2
        write(10,"(2(1pe13.6,1x))") kx(i), E0(kx(i))
      end do
      close(10)

!.... currently running serial to avoid problems with random number generator

      ul = czero
      do k = 1, nz
        do j = 1, ny
          do i = 1, mx-1
#ifdef CRAY
            theta1 = ranf()
            theta2 = ranf()
            phi    = ranf()
#else
#ifdef RANF
            theta1 = ranf(i+1)
            theta2 = ranf(i+2)
            phi    = ranf(i+3)
#else
            theta1 = rand()
            theta2 = rand()
            phi    = rand()
!           call random_number( theta1 )
!           call random_number( theta2 )
!           call random_number(    phi )
#endif
#endif
            kmag   = sqrt( kx(i)**2 + ky(j)**2 + kz(k)**2 )
            k12    = sqrt( kx(i)**2 + ky(j)**2 )

            if ( kmag .eq. zero .or. mask(i,j,k) ) then
              fact = zero
            else
              fact = sqrt( E0( kmag ) / (two * pi) ) / kmag
            end if

!           write(*,"(8(1pe12.5,1x))") kx(i), ky(j), kz(k), kmag, &
!                                      E0(kmag), fact

!.... For random (no mean) helicity

            alp = fact * exp( two * iota * pi * theta1 ) * &
                    cos( two * pi * phi )
            beta = fact * exp( two * iota * pi * theta2 ) * &
                   sin( two * pi * phi )

!.... For maximum mean helicity

!           alp = fact * exp( two * iota * pi * theta1 ) / sqrt(two)
!           beta  = iota * alp

            if ( k12 .eq. zero ) then
              ul(i,j,k,1) = alp
              ul(i,j,k,2) = beta
              if ( k .eq. 1 ) ul(i,j,k,3) = zero 
            else
              ul(i,j,k,1) = ( alp * kmag * ky(j) + beta * kx(i) * kz(k) ) / &
                            ( kmag * k12 )
              ul(i,j,k,2) = ( beta * ky(j) * kz(k) - alp * kmag * kx(i) ) / &
                            ( kmag * k12 )
              if ( k .eq. 1 ) ul(i,j,k,3) = -beta * k12 / kmag
            end if
          end do
        end do
      end do

!.... satisfy conjugate symmetry

      do j = 2, my-1
        ul(1,j,1,3) = conjg( ul(1,ny+2-j,1,3) )
      end do

      do idof = 1, 2
        ul(1,1,1,idof) = zero
        do j = 2, my-1
          ul(1,j,1,idof) = conjg( ul(1,ny+2-j,1,idof) )
        end do
        do k = 2, mz-1
          ul(1,1,k,idof) = conjg( ul(1,1,nz+2-k,idof) )
          do j = 2, ny
            ul(1,j,k,idof) = conjg( ul(1,ny+2-j,nz+2-k,idof) )
          end do
        end do
      end do

!.... satisfy continuity

      do i = 2, mx-1
        ul(i,1,1,1) = zero
      end do

      !$omp parallel do private(j,i)
      do j = 2, ny
        do i = 1, mx-1
          ul(i,j,1,2) = -(kx(i)/ky(j))*ul(i,j,1,1)
        end do
      end do

      !$omp parallel do private(k,j,i)
      do k = 2, nz
        do j = 1, ny
          do i = 1, mx
            ul(i,j,k,3) = -( kx(i)*ul(i,j,k,1) + ky(j)*ul(i,j,k,2) ) / kz(k)
          end do
        end do
      end do

!.... zero the oddball wavenumber

      call oddball( ul )

!.... compute statistics

      if (stats.eq.1) call statistics( u, f )

      return
end subroutine turb

logical function mask( i, j, k )
      use field
      use const
      implicit none
      integer :: i, j, k
      if ( (kx(i)/nx)**2 + (ky(j)/ny)**2 + (kz(k)/nz)**2 .gt. &
           truncation_radius_sq ) then
        mask = .true.
      else
        mask = .false.
      end if
      return
end function mask

logical function pmask( i, j, k )
      use field
      use const
      implicit none
      integer :: i, j, k
      if ( (pkx(i)/nx)**2 + (pky(j)/ny)**2 + (pkz(k)/nz)**2 .gt. &
           truncation_radius_sq ) then
        pmask = .true.
      else
        pmask = .false.
      end if
      return
end function pmask

subroutine restart(ul, fl)
      use field
      use const
      implicit none
      complex :: ul(mx,ny,nz,ndof), fl(mx,ny,nz,ndof)
      integer :: inx, iny, inz
      integer :: i, j, k, idof
      write(*,"('Reading from restart.in ...')")
      open(10,file='restart.in',form='unformatted')
      read(10) istep, inx, iny, inz
      if (inx.ne.nx .or. iny.ne.ny .or. inz.ne.nz) then
        call error('restart$','Grid dimensions do not match$')
      end if
      read(10) t, ul
      close(10)
      if (stats.eq.1) call statistics( u, f )
      return
end subroutine restart

subroutine taylor
      use field
      use const
      implicit none
      integer :: i, j, k, idof

!.... Taylor-Green vortex:  See Orszag, Studies in Applied Mathematics, 50
!.....                      pp. 293-327

      !$omp parallel do private(k,j,i)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            u(i,j,k,1) =  cos(kx(2)*x(i)) * sin(ky(2)*y(j)) * cos(kz(2)*z(k))
            u(i,j,k,2) = -sin(kx(2)*x(i)) * cos(ky(2)*y(j)) * cos(kz(2)*z(k))
            u(i,j,k,3) = zero
          end do
        end do
      end do

!.... write out the initial condition

      call wgrid( 'out.xyz'//char(0), nx, ny, nz, x, y, z )
      call wdata( 'taylor.q'//char(0), ndof, nx, ny, nz, u, zero, zero, &
                  one/nu, t )

!.... Take the FFT of the initial field

      do idof = 1, ndof
        call fft3d(-1, nx, ny, nz, u(1,1,1,idof), mx, ny, coef)
      end do

!.... compute some statistics

      if (stats.eq.1) call statistics( u, f )

      return
end subroutine taylor

subroutine post(ul)
      use const
      use field
      use timemod
      implicit none
      complex :: ul(mx,ny,nz,ndof)
      integer :: idof, i, j, k
      
!.... save a restart file

      open(10,file='restart.out',form='unformatted')
      write(10) istep, nx, ny, nz
      write(10) t, ul
      close(10)
      
!.... output statistics for the final time-step

      call rhs3d(u, r, f, .false., uh, Sijh, Lij, Mij)
      call spectra(2, istep, u, f)

      write(*,10)
  10  format('Wrote restart.out file')

!.... inverse FFT

      do idof = 1, ndof
        call fft3d(1, nx, ny, nz, u(1,1,1,idof), mx, ny, coef)
      end do

!.... output the results in a Plot3d file
      call wgrid( 'out.xyz'//char(0), nx, ny, nz, x, y, z )
      call wdata( 'out.q'//char(0), ndof, nx, ny, nz, u, zero, zero, &
                  one/nu, t )
      write(*,20)
  20  format('Wrote Plot3d files')

!.... close the output files

      close(20); close(21); close(22); close(11)
      if (LES.eq.2) close(12)

!.... stop the clock

      cpu2 = second()
      write(*,30) (cpu2 - cpu)/nthreads
  30  format('Post process time = ',1pe13.7) 
      cpu = cpu2
      write(*,40) (second() - cpu0)/nthreads
  40  format('Total time        = ',1pe13.7)

      return
end subroutine post
   
subroutine euler(ul, rl)

!.... Explicit Euler time advancement

      use field
      use timemod
      use intmod
      implicit none
      complex :: ul(mx*ny*nz*ndof), rl(mx*ny*nz*ndof)

      real :: tn
      complex :: utmp, rtmp
      integer :: n, i, imax
      integer :: ier

      imax = mx*ny*nz*ndof

!.... Integrate in time using explicit Euler

      do n = 1, nstep

        istep = istep + 1

        tn = t
        call rhs3d(u, r, f, .true., uh, Sijh, Lij, Mij)
        if (stats.eq.1) call spectra(n, istep-1, u, f)

        if (cfl.ne.0.0) dt = cfl / velmax

        !$omp parallel do private(i)
        do i = 1, imax
          ul(i) = ul(i) + dt * rl(i)
        end do

        t = tn + dt
        if (stats.eq.1) call statistics( u, f )

        cpu2 = second()
        write(11,"(i5,2(1x,1pe10.3))") istep, t, (cpu2-cpu)/nthreads
        cpu = cpu2

        if (halt) exit
      end do

      return
end subroutine euler

subroutine rk2(ul, rl)

!.... Second-order Runge-Kutta time advancement

      use field
      use const
      use timemod
      use intmod
      implicit none
      complex :: ul(mx*ny*nz*ndof), rl(mx*ny*nz*ndof)

      real :: tn, fac
      integer :: n, i, imax
      integer :: ier

      imax = mx*ndof*ny*nz

!.... Integrate in time using explicit RK2

      do n = 1, nstep

        istep = istep + 1

        tn = t

        call rhs3d(u, r, f, .true., uh, Sijh, Lij, Mij)
        if (stats.eq.1) call spectra(n, istep-1, u, f)

        if (cfl.ne.0.0) dt = cfl / velmax

        fac = pt5 * dt
        !$omp parallel do private(i)
        do i = 1, imax
          rl(i) = ul(i) + fac * rl(i)
        end do

        t = tn + pt5 * dt
        call rhs3d(r, r, f, .false., uh, Sijh, Lij, Mij)
        
        !$omp parallel do private(i)
        do i = 1, imax
          ul(i) = ul(i) + dt * rl(i)
        end do

        t = tn + dt
        if (stats.eq.1) call statistics( u, f )

        cpu2 = second()
        write(11,"(i5,2(1x,1pe10.3))") istep, t, (cpu2-cpu)/nthreads
        call flush(11)
        cpu = cpu2

        if (halt) exit
      end do

      !open(10,file='tmp2.dat',form='unformatted')
      !write(10) istep, nx, ny, nz
      !write(10) t, ul
      !close(10)

      return
end subroutine rk2

subroutine rk3(ul, rl)
      use field
      use const
      use timemod
      use intmod
      implicit none
      complex :: ul(mx*ny*nz*ndof), rl(mx*ny*nz*ndof)

!.... Third-order Runge-Kutta coefficients are from Allen Wray, "Minimal
!.... Storage Time-Advancement Schemes" (I think its in JCP)

!.... RK3(1)

      real, parameter :: a1 = 0.2500000000000000, a2 = 0.5333333333333333333
      real, parameter :: b1 = 0.0000000000000000, b2 = 0.4166666666666666667
      real, parameter :: c1 = 0.7500000000000000

!.... RK3(2)

!     real, parameter :: a1 = 0.2500000000000000, a2 = 0.6666666666666666667 
!     real, parameter :: b1 = 0.1500000000000000, b2 = 0.4166666666666666667
!     real, parameter :: c1 = 0.6000000000000000

      real :: tn
      complex :: utmp, rtmp
      integer :: n, i, imax
      integer :: ier

      imax = mx*ny*nz*ndof

!.... Integrate in time using low-storage RK3

      do n = 1, nstep

        istep = istep + 1

        tn = t
        call rhs3d(u, r, f, .true., uh, Sijh, Lij, Mij)
        if (stats.eq.1) call spectra(n, istep-1, u, f)

        if (cfl.ne.0.0) dt = cfl / velmax

        !$omp parallel do private(i, utmp, rtmp)
        do i = 1, imax
          utmp = ul(i)
          rtmp = rl(i)
          ul(i) = utmp + a1 * dt * rtmp
          rl(i) = utmp + a2 * dt * rtmp
        end do

        t = tn + a2 * dt
        call rhs3d(r, r, f, .false., uh, Sijh, Lij, Mij)
        
        !$omp parallel do private(i, utmp, rtmp)
        do i = 1, imax
          utmp = ul(i)
          rtmp = rl(i)
!         ul(i) = utmp + b1 * dt * rtmp   ! Don't need for RK3(1)
          rl(i) = utmp + b2 * dt * rtmp
        end do

         t = tn + (a1 + b2) * dt
        call rhs3d(r, r, f, .false., uh, Sijh, Lij, Mij)

        !$omp parallel do private(i)
        do i = 1, imax
          ul(i) = ul(i) + c1 * dt * rl(i)
        end do

        t = tn + dt
        if (stats.eq.1) call statistics( u, f )

        cpu2 = second()
        write(11,"(i5,2(1x,1pe10.3))") istep, t, (cpu2-cpu)/nthreads
        cpu = cpu2

        if (halt) exit
      end do

      return
end subroutine rk3

subroutine rhs3d(ul, rl, fl, compute_max, uhl, Sijhl, Lijl, Mijl)
      use field
      use const
      implicit none
      complex :: ul(mx,ny,nz,ndof), rl(mx,ny,nz,ndof)
      complex :: fl(mpx,py,pz,ndofles), Mijl(mpx,py,pz,6)
      complex :: uhl(mpx,py,pz,ndof), Sijhl(mpx,py,pz,6), Lijl(mpx,py,pz,6)
      real    :: rtmp, nu_e, S, Sh, S33
      complex :: ctmp
      integer :: i, j, k, idof, jl, kl
      logical :: compute_max
      logical, external :: mask
      real    :: vtmp(2*mx)

      if(LES.gt.0) then

!.... apply the grid filter to the solution (Must do to better match T. Lund)

      if (LES.eq.2 .and. grid_filter.eq.1) then
        !$omp parallel do private(k,j,i)
        do k = 1, nz
          do j = 1, ny
            do i = 1, mx
              if (mask(i,j,k)) then
                ul(i,j,k,1) = zero
                ul(i,j,k,2) = zero
                ul(i,j,k,3) = zero
              end if
            end do
          end do
        end do
      end if

!.... copy the regular field to an 3/2 field

      call pad( ndof, nx, ny, nz, ul, px, py, pz, fl )

!.... compute the strain-rate tensor and put into fl

      !$omp parallel do private(k,j,i)
      do k = 1, pz
        do j = 1, py
          do i = 1, mpx
            fl(i,j,k,4) = iota * pt5*( pky(j)*fl(i,j,k,1) + pkx(i)*fl(i,j,k,2))
            fl(i,j,k,5) = iota * pt5*( pky(j)*fl(i,j,k,3) + pkz(k)*fl(i,j,k,2))
            fl(i,j,k,6) = iota * pt5*( pkz(k)*fl(i,j,k,1) + pkx(i)*fl(i,j,k,3))
            fl(i,j,k,7) = iota * pkx(i) * fl(i,j,k,1)
            fl(i,j,k,8) = iota * pky(j) * fl(i,j,k,2)
          end do
        end do
      end do

!.... Dynamic model

      if (LES.eq.2) then

!.... Apply the test filter to the velocity and strain-rate

        !$omp parallel do private(k,j,i)
        do k = 1, pz
          do j = 1, py
            do i = 1, mpx
              uhl(i,j,k,1)   = kernel(i,j,k) * fl(i,j,k,1)
              uhl(i,j,k,2)   = kernel(i,j,k) * fl(i,j,k,2)
              uhl(i,j,k,3)   = kernel(i,j,k) * fl(i,j,k,3)
              Sijhl(i,j,k,1) = kernel(i,j,k) * fl(i,j,k,7)
              Sijhl(i,j,k,2) = kernel(i,j,k) * fl(i,j,k,8)
              Sijhl(i,j,k,3) = kernel(i,j,k) * (-fl(i,j,k,7)-fl(i,j,k,8))
              Sijhl(i,j,k,4) = kernel(i,j,k) * fl(i,j,k,4)
              Sijhl(i,j,k,5) = kernel(i,j,k) * fl(i,j,k,5)
              Sijhl(i,j,k,6) = kernel(i,j,k) * fl(i,j,k,6)
            end do
          end do
        end do

!.... Inverse FFT of filtered velocity and strain-rate

        do idof = 1, ndof
          call fft3d(1, px, py, pz, uhl(1,1,1,idof), mpx, py, coefp)
        end do
        do idof = 1, 6
          call fft3d(1, px, py, pz, Sijhl(1,1,1,idof), mpx, py, coefp)
        end do

      end if

!.... Inverse FFT of 3/2 field

      do idof = 1, ndofles
        call fft3d(1, px, py, pz, fl(1,1,1,idof), mpx, py, coefp)
      end do

!.... Dynamic model

      if (LES.eq.2) then

!.... Form Lij, S_mag, S_mag_hat, Mij

      !$omp parallel do private(k,j,i,S33,S,Sh)
      do k = 1, pz
        do j = 1, py
          do i = 1, px
            Lij(i,j,k,1) = f(i,j,k,1)**2 - uh(i,j,k,1)**2
            Lij(i,j,k,2) = f(i,j,k,2)**2 - uh(i,j,k,2)**2
            Lij(i,j,k,3) = f(i,j,k,3)**2 - uh(i,j,k,3)**2
            Lij(i,j,k,4) = f(i,j,k,1)*f(i,j,k,2) - uh(i,j,k,1)*uh(i,j,k,2)
            Lij(i,j,k,5) = f(i,j,k,2)*f(i,j,k,3) - uh(i,j,k,2)*uh(i,j,k,3)
            Lij(i,j,k,6) = f(i,j,k,3)*f(i,j,k,1) - uh(i,j,k,3)*uh(i,j,k,1)

            S33 = -(f(i,j,k,7)+f(i,j,k,8))
            S = sqrt(two*(f(i,j,k,7)**2 + f(i,j,k,8)**2 + S33**2 + &
                     two*(f(i,j,k,4)**2 + f(i,j,k,5)**2 + f(i,j,k,6)**2)))
            Sh = sqrt(two*(Sijh(i,j,k,1)**2 + Sijh(i,j,k,2)**2 + &
                      Sijh(i,j,k,3)**2 + two*(Sijh(i,j,k,4)**2 + &
                      Sijh(i,j,k,5)**2 + Sijh(i,j,k,6)**2)))

            Mij(i,j,k,1) = alpha_sq * Sh * Sijh(i,j,k,1) - S * f(i,j,k,7)
            Mij(i,j,k,2) = alpha_sq * Sh * Sijh(i,j,k,2) - S * f(i,j,k,8)
            Mij(i,j,k,3) = alpha_sq * Sh * Sijh(i,j,k,3) - S * S33
            Mij(i,j,k,4) = alpha_sq * Sh * Sijh(i,j,k,4) - S * f(i,j,k,4)
            Mij(i,j,k,5) = alpha_sq * Sh * Sijh(i,j,k,5) - S * f(i,j,k,5)
            Mij(i,j,k,6) = alpha_sq * Sh * Sijh(i,j,k,6) - S * f(i,j,k,6)
          end do
        end do
      end do

!.... FFT of Lij and Mij

      do idof = 1, 6
        call fft3d(-1, px, py, pz, Lij(1,1,1,idof), mpx, py, coefp)
        call fft3d(-1, px, py, pz, Mij(1,1,1,idof), mpx, py, coefp)
      end do

!.... Filter Lij and Mij

      !$omp parallel do private(k,idof,j,i)
      do k = 1, pz
        do idof = 1, 6
          do j = 1, py
            do i = 1, mpx
              Lijl(i,j,k,idof) = kernel(i,j,k) * Lijl(i,j,k,idof)
              Mijl(i,j,k,idof) = kernel(i,j,k) * Mijl(i,j,k,idof)
            end do
          end do
        end do
      end do

!.... Inverse FFT of Lij and Mij

      do idof = 1, 6
        call fft3d(1, px, py, pz, Lij(1,1,1,idof), mpx, py, coefp)
        call fft3d(1, px, py, pz, Mij(1,1,1,idof), mpx, py, coefp)
      end do

!.... Perform contraction

      if (contraction .eq. 0) then  ! Germano et al.
        !$omp parallel do private(k,j,i,S33)
        do k = 1, pz
          do j = 1, py
            do i = 1, px
              S33 = -(f(i,j,k,7)+f(i,j,k,8))
              num(i,j,k) = Lij(i,j,k,1)*f(i,j,k,7) + Lij(i,j,k,2)*f(i,j,k,8) &
                         + Lij(i,j,k,3)*S33 + two * (Lij(i,j,k,4)*f(i,j,k,4) &
                         + Lij(i,j,k,5)*f(i,j,k,5) + Lij(i,j,k,6)*f(i,j,k,6) )
              den(i,j,k) = Mij(i,j,k,1)*f(i,j,k,7) + Mij(i,j,k,2)*f(i,j,k,8) &
                         + Mij(i,j,k,3)*S33 + two * (Mij(i,j,k,4)*f(i,j,k,4) &
                         + Mij(i,j,k,5)*f(i,j,k,5) + Mij(i,j,k,6)*f(i,j,k,6) )
            end do
          end do
        end do
      else if ( contraction .eq. 1) then  ! Lilly 
        !$omp parallel do private(k,j,i)
        do k = 1, pz
          do j = 1, py
            do i = 1, px
              num(i,j,k) = Lij(i,j,k,1) * Mij(i,j,k,1) + &
                           Lij(i,j,k,2) * Mij(i,j,k,2) + &
                           Lij(i,j,k,3) * Mij(i,j,k,3) + two * ( &
                           Lij(i,j,k,4) * Mij(i,j,k,4) + &
                           Lij(i,j,k,5) * Mij(i,j,k,5) + &
                           Lij(i,j,k,6) * Mij(i,j,k,6) )
              den(i,j,k) = Mij(i,j,k,1) * Mij(i,j,k,1) + &
                           Mij(i,j,k,2) * Mij(i,j,k,2) + &
                           Mij(i,j,k,3) * Mij(i,j,k,3) + two * ( &
                           Mij(i,j,k,4) * Mij(i,j,k,4) + &
                           Mij(i,j,k,5) * Mij(i,j,k,5) + &
                           Mij(i,j,k,6) * Mij(i,j,k,6) )
            end do
          end do
        end do
      else
        call error("rhs3d$","Illegal value of contraction$")
      end if

!.... compute C (this needs OpenMP directives)

      if (average.eq.0) then
        C = -pt5 / delta_sq * num / den
      else if (average.eq.1) then
        do k = 1, pz
          C(:,:,k) = -pt5 / delta_sq * sum( num(:,:,k) ) / sum( den(:,:,k) )
        end do
      else if (average.eq.2) then
        C = -pt5 / delta_sq * sum(num)/sum(den)
        write(12,"(2(1pe13.6,1x))") t, C(1,1,1)
        call flush(12)
      else
        call error("rhs3d$","Illegal value of average$")
      end if

      end if  ! LES.eq.2

!.... find max velocity for CFL constraint

      if (compute_max) then
        velmax = zero
        !$omp parallel do private(k,j,i), reduction(max:velmax)
        do k = 1, pz
          do j = 1, py
            do i = 1, px
              velmax = max( velmax, nx*abs(f(i,j,k,1)) + ny*abs(f(i,j,k,2)) + &
                            nz*abs(f(i,j,k,3)) )
            end do
          end do
        end do
      end if
       
!.... Form nonlinear products (in physical space)

      s = zero
      !$omp parallel do private(k,j,i,S33,S,nu_e)
      do k = 1, pz
        do j = 1, py
          do i = 1, px
            S33 = -(f(i,j,k,7)+f(i,j,k,8))
            S = sqrt( two*(f(i,j,k,7)**2 + f(i,j,k,8)**2 + s33**2 + &
                      two*(f(i,j,k,4)**2 + f(i,j,k,5)**2 + f(i,j,k,6)**2)) )
            nu_e = C(i,j,k) * delta_sq * S
            f(i,j,k,4) = f(i,j,k,1) * f(i,j,k,2) - two * nu_e * f(i,j,k,4)
            f(i,j,k,5) = f(i,j,k,2) * f(i,j,k,3) - two * nu_e * f(i,j,k,5)
            f(i,j,k,6) = f(i,j,k,1) * f(i,j,k,3) - two * nu_e * f(i,j,k,6)
            f(i,j,k,1) = f(i,j,k,1) * f(i,j,k,1) - two * nu_e * f(i,j,k,7)     
            f(i,j,k,2) = f(i,j,k,2) * f(i,j,k,2) - two * nu_e * f(i,j,k,8)
            f(i,j,k,3) = f(i,j,k,3) * f(i,j,k,3) - two * nu_e * S33
          end do
        end do
      end do

!.... FFT of 3/2 field

      do idof = 1, 2*ndof
        call fft3d(-1, px, py, pz, fl(1,1,1,idof), mpx, py, coefp)
      end do

!.... convert from anti-aliased 3/2 field to regular field (in place)

      call unpad( 2*ndof, nx, ny, nz, px, py, pz, fl )
     
      else   ! DNS

!.... copy the regular field to an 3/2 field

      call pad( ndof, nx, ny, nz, ul, px, py, pz, fl )

!.... Inverse FFT of 3/2 field

      do idof = 1, ndof
        call fft3d(1, px, py, pz, fl(1,1,1,idof), mpx, py, coefp)
      end do

!.... Compute maximum velocity for CFL control (should be in parallel!)

      if (compute_max) then
        velmax = zero
        !$omp parallel do private(k,j,i), reduction(max:velmax)
        do k = 1, pz
          do j = 1, py
            do i = 1, px
              velmax = max( velmax, nx*abs(f(i,j,k,1)) + ny*abs(f(i,j,k,2)) + &
                            nz*abs(f(i,j,k,3)) )
            end do
          end do
        end do
      end if

!.... Form nonlinear products (in physical space)

      !$omp parallel do private(k,j,i), shared(f)
      do k = 1, pz
        do j = 1, py
#ifdef VECTORIZE
          f(:,j,k,4) = f(:,j,k,1) * f(:,j,k,2)   !  uv
          f(:,j,k,5) = f(:,j,k,2) * f(:,j,k,3)   !  vw
          f(:,j,k,6) = f(:,j,k,1) * f(:,j,k,3)   !  uw
#else
          do i = 1, px
            f(i,j,k,4) = f(i,j,k,1) * f(i,j,k,2)   !  uv
            f(i,j,k,5) = f(i,j,k,2) * f(i,j,k,3)   !  vw
            f(i,j,k,6) = f(i,j,k,1) * f(i,j,k,3)   !  uw
          end do
#endif
        end do
      end do

      !$omp parallel do private(k,j,i), shared(f)
      do k = 1, pz
        do j = 1, py
          do i = 1, px
            f(i,j,k,1) = f(i,j,k,1) ** 2           !  u^2
            f(i,j,k,2) = f(i,j,k,2) ** 2           !  v^2
            f(i,j,k,3) = f(i,j,k,3) ** 2           !  w^2
          end do
        end do
      end do

!.... FFT of 3/2 field

      do idof = 1, 2*ndof
        call fft3d(-1, px, py, pz, fl(1,1,1,idof), mpx, py, coefp)
      end do

!.... convert from anti-aliased 3/2 field to regular field (in place)

      call unpad( 2*ndof, nx, ny, nz, px, py, pz, fl )

      end if  ! DNS

!.... Form the convection term by taking the divergence of f

      !$omp parallel do private(k,j,i)
      do k = 1, nz
        do j = 1, ny
          do i = 1, mx
            fl(i,j,k,1) = -iota * ( kx(i) * fl(i,j,k,1) + ky(j) * &
                          fl(i,j,k,4) + kz(k) * fl(i,j,k,6) )
            fl(i,j,k,2) = -iota * ( kx(i) * fl(i,j,k,4) + ky(j) * &
                          fl(i,j,k,2) + kz(k) * fl(i,j,k,5) )
            fl(i,j,k,3) = -iota * ( kx(i) * fl(i,j,k,6) + ky(j) * &
                          fl(i,j,k,5) + kz(k) * fl(i,j,k,3) )
          end do
        end do
      end do

!.... Form the viscous term

      !$omp parallel do private(k,j,i,rtmp,vtmp)
      do k = 1, nz
        do j = 1, ny
#ifdef VECTORIZE
          vtmp(:) = -nu * ( kx2(:)**2 + ky(j)**2 + kz(k)**2 )
          r(:,j,k,1) = vtmp(:) * u(:,j,k,1)
          r(:,j,k,2) = vtmp(:) * u(:,j,k,2)
          r(:,j,k,3) = vtmp(:) * u(:,j,k,3)
#else
          do i = 1, mx
            rtmp = -nu * ( kx(i)**2 + ky(j)**2 + kz(k)**2 )
            rl(i,j,k,1) = rtmp * ul(i,j,k,1)
            rl(i,j,k,2) = rtmp * ul(i,j,k,2)
            rl(i,j,k,3) = rtmp * ul(i,j,k,3)
          end do
#endif
        end do
      end do

!.... Form the pressure gradient term (be careful about the mean pressure)

      !$omp parallel do private(k,j,i,ctmp)
      do k = 2, nz
        do j = 1, ny
          do i = 1, mx
            ctmp = ( kx(i) * fl(i,j,k,1) + ky(j) * fl(i,j,k,2) + &
                     kz(k) * fl(i,j,k,3) ) / &
                   ( kx(i)**2 + ky(j)**2 + kz(k)**2 )
            fl(i,j,k,1) = fl(i,j,k,1) - kx(i) * ctmp
            fl(i,j,k,2) = fl(i,j,k,2) - ky(j) * ctmp
            fl(i,j,k,3) = fl(i,j,k,3) - kz(k) * ctmp
          end do
        end do
      end do

      k = 1
      !$omp parallel do private(j,i,ctmp)
      do j = 2, ny
        do i = 1, mx
          ctmp = ( kx(i) * fl(i,j,k,1) + ky(j) * fl(i,j,k,2) ) / &
                 ( kx(i)**2 + ky(j)**2 )
          fl(i,j,k,1) = fl(i,j,k,1) - kx(i) * ctmp
          fl(i,j,k,2) = fl(i,j,k,2) - ky(j) * ctmp
        end do
      end do

      k = 1; j = 1
#ifdef VECTORIZE
      f(:,j,k,1) = zero
#else
      !$omp parallel do private(i)
      do i = 1, mx
        fl(i,j,k,1) = zero
      end do
#endif

!.... Include the nonlinear term

      !$omp parallel do private(j,i)
      do k = 1, nz
        do j = 1, ny
#ifdef VECTORIZE
          r(:,j,k,1) = r(:,j,k,1) + f(:,j,k,1)
          r(:,j,k,2) = r(:,j,k,2) + f(:,j,k,2)
          r(:,j,k,3) = r(:,j,k,3) + f(:,j,k,3)
#else
          do i = 1, 2*mx
            r(i,j,k,1) = r(i,j,k,1) + f(i,j,k,1)
            r(i,j,k,2) = r(i,j,k,2) + f(i,j,k,2)
            r(i,j,k,3) = r(i,j,k,3) + f(i,j,k,3)
          end do
#endif
        end do
      end do

!     if (iprint .eq. 1) call cprint( ndof, nx, ny, nz, rl )

!.... make sure to zero the odd-ball wavenumbers

      call oddball( rl )

      return
end subroutine rhs3d

subroutine get_kernel
      use field
      use const
      implicit none
      integer :: i, j, k
      real, parameter :: eps = 1.0e-10    ! from Lund's code
      
      if (filter_type .eq. 0) then        ! spherical cutoff

        do k = 1, pz
          do j = 1, py
            do i = 1, mpx
              if ( (pkx(i)/nx)**2 + (pky(j)/ny)**2 + (pkz(k)/nz)**2 &
                   .le. filter_radius_sq ) then
                kernel(i,j,k) = one
              else
                kernel(i,j,k) = zero
              end if
            end do
          end do
        end do

      else if (filter_type .eq. 1) then   ! exact top-hat

        do k = 1, pz
          do j = 1, py
            do i = 1, mpx
              kernel(i,j,k) =  sin( pi * ( pkx(i)/nx + eps ) * alpha ) * &
                               sin( pi * ( pky(j)/ny + eps ) * alpha ) * &
                               sin( pi * ( pkz(k)/nz + eps ) * alpha ) / &
                               ( (pi*alpha)**3 * ( pkx(i)/nx + eps ) *   &
                               ( pky(j)/ny + eps ) * ( pkz(k)/nz + eps ) )
            end do
          end do
        end do

      else if ( filter_type .eq. 2) then  ! Gaussian

        do k = 1, pz
          do j = 1, py
            do i = 1, mpx
              kernel(i,j,k) =  exp( -pi**2/6.0 * ( (pkx(i)/nx)**2 + &
                               (pky(j)/ny)**2 + (pkz(k)/nz)**2 ) * alpha_sq )
            end do
          end do
        end do

      else if ( filter_type .eq. 3) then  ! Trapezoidal

        do k = 1, pz
          do j = 1, py
            do i = 1, mpx
              kernel(i,j,k) = ( 1 + cos( pi * ( pkx(i)/nx ) * alpha ) ) * &
                              ( 1 + cos( pi * ( pky(j)/ny ) * alpha ) ) * &
                              ( 1 + cos( pi * ( pkz(k)/nz ) * alpha ) ) * 0.125
            end do
          end do
        end do

      else if ( filter_type .eq. 4) then  ! Simpsons

        do k = 1, pz
          do j = 1, py
            do i = 1, mpx
              kernel(i,j,k) = ( 2 + cos( pi * ( pkx(i)/nx ) * alpha ) ) * &
                              ( 2 + cos( pi * ( pky(j)/ny ) * alpha ) ) * &
                              ( 2 + cos( pi * ( pkz(k)/nz ) * alpha ) ) / 27.0
            end do
          end do
        end do

      else
        call error('get_kernel$','Illegal value for filter_type$')
      end if

      return
end subroutine get_kernel

subroutine statistics( ul, fl )
      use field
      use const
      implicit none
      complex :: ul(mx,ny,nz,ndof), fl(mpx,py,pz,ndof)
      integer :: i, j, k, idof
      integer :: k0, k1, k2, ke, kb, ierr, temp, pid
      integer, external :: kill
      real :: urms1, urms2, urms3, vrms, E, kmag, L, lambda, R_L, R_lambda, D
      real :: savg, sx, sx1, sx2, sy, sy1, sy2, sz, sz1, sz2, enstrophy, div
      real :: eps, ups, eta, q, tke, Leps, Teps
      character*80 :: base='spec', fname

!.... compute the rms velocity

      urms1 = zero; urms2 = zero; urms3 = zero
      !$omp parallel do private(k,j,i), reduction(+: urms1,urms2,urms3)
      do k = 1, nz
        do j = 1, ny
          i = 1
          urms1 = urms1 + ul(i,j,k,1) * conjg(ul(i,j,k,1))
          urms2 = urms2 + ul(i,j,k,2) * conjg(ul(i,j,k,2))
          urms3 = urms3 + ul(i,j,k,3) * conjg(ul(i,j,k,3))
          do i = 2, mx
            urms1 = urms1 + two * ul(i,j,k,1) * conjg(ul(i,j,k,1))
            urms2 = urms2 + two * ul(i,j,k,2) * conjg(ul(i,j,k,2))
            urms3 = urms3 + two * ul(i,j,k,3) * conjg(ul(i,j,k,3))
          end do
        end do
      end do
      q    = sqrt( urms1 + urms2 + urms3 )
      tke  = pt5 * q**2
      vrms = sqrt( pt33 * q**2 )
      urms1 = sqrt( urms1 )
      urms2 = sqrt( urms2 )
      urms3 = sqrt( urms3 )

!.... compute the integral scale and Taylor microscale, defined in Orszag &
!.... Patterson "Numerical Simulation of Turbulence," Springer, 1972.

      L = zero
      lambda = zero
      !$omp parallel do private(k,j,i,kmag), reduction(+: L,lambda)
      do k = 1, nz
        do j = 1, ny
          i = 1
          kmag = sqrt( kx(i)**2 + ky(j)**2 + kz(k)**2 )
          if (kmag.ne.zero) then
            L = L + one/kmag * ( ul(i,j,k,1) * conjg(ul(i,j,k,1)) + &
                                 ul(i,j,k,2) * conjg(ul(i,j,k,2)) + &
                                 ul(i,j,k,3) * conjg(ul(i,j,k,3)) )
            lambda = lambda +  kmag**2 * ( &
                                 ul(i,j,k,1) * conjg(ul(i,j,k,1)) + &
                                 ul(i,j,k,2) * conjg(ul(i,j,k,2)) + &
                                 ul(i,j,k,3) * conjg(ul(i,j,k,3)) )
          end if
          do i = 2, mx
            kmag = sqrt( kx(i)**2 + ky(j)**2 + kz(k)**2 )
            L = L + two/kmag * ( ul(i,j,k,1) * conjg(ul(i,j,k,1)) + &
                                 ul(i,j,k,2) * conjg(ul(i,j,k,2)) + &
                                 ul(i,j,k,3) * conjg(ul(i,j,k,3)) )
            lambda = lambda +  two*kmag**2 * ( &
                                 ul(i,j,k,1) * conjg(ul(i,j,k,1)) + &
                                 ul(i,j,k,2) * conjg(ul(i,j,k,2)) + &
                                 ul(i,j,k,3) * conjg(ul(i,j,k,3)) )
          end do
        end do
      end do
      L = 0.25 * pi / vrms * L
      lambda = sqrt( 15.0 * vrms**2 / lambda )

!.... compute the band averaged energy spectrum every nout time steps

      if (mod(istep,nout).eq.0 .and. .false.) then
        call makename(base,istep,fname)
        open(10,file=fname)
        write(10,"('# t = ', 1pe13.6)") t
        k1 = 0
        k2 = sqrt( kx(mx)**2 + ky(my)**2 + kz(mz)**2 )
        do k0 = k1+1, k2-1, 2
          kb = k0 - 1
          ke = k0 + 1
          E = zero
          D = zero
          !$omp parallel do private(k,j,i,kmag), reduction(+: E, D)
          do k = 1, nz
            do j = 1, ny
              i = 1
              kmag = sqrt( kx(i)**2 + ky(j)**2 + kz(k)**2 )
              if (kmag .ge. kb .and. kmag .lt. ke) then
                E = E + pt5 * ( ul(i,j,k,1) * conjg(ul(i,j,k,1)) + &
                                ul(i,j,k,2) * conjg(ul(i,j,k,2)) + &
                                ul(i,j,k,3) * conjg(ul(i,j,k,3)) )
                D = D + nu * kmag**2 * ( &
                                ul(i,j,k,1) * conjg(ul(i,j,k,1)) + &
                                ul(i,j,k,2) * conjg(ul(i,j,k,2)) + &
                                ul(i,j,k,3) * conjg(ul(i,j,k,3)) )
              end if
              do i = 2, mx
                kmag = sqrt( kx(i)**2 + ky(j)**2 + kz(k)**2 )
                if (kmag .ge. kb .and. kmag .lt. ke) then
                  E = E +       ( ul(i,j,k,1) * conjg(ul(i,j,k,1)) + &
                                  ul(i,j,k,2) * conjg(ul(i,j,k,2)) + &
                                  ul(i,j,k,3) * conjg(ul(i,j,k,3)) )
                  D = D + two * nu * kmag**2 * ( &
                                  ul(i,j,k,1) * conjg(ul(i,j,k,1)) + &
                                  ul(i,j,k,2) * conjg(ul(i,j,k,2)) + &
                                  ul(i,j,k,3) * conjg(ul(i,j,k,3)) )
                end if
              end do
            end do
          end do
          E = pt5 * E
          D = pt5 * D
          write(10,"(i4,4(1x,1pe13.6))") k0, E, D
        end do
        close(10)
      end if

!.... form the velocity-derivatives and the vorticity components

      !$omp parallel do private(k,j,i)
      do k = 1, nz
        do j = 1, ny
          do i = 1, mx
            fl(i,j,k,1) = iota*( ky(j)*ul(i,j,k,3) - kz(k)*ul(i,j,k,2) )
            fl(i,j,k,2) = iota*( kz(k)*ul(i,j,k,1) - kx(i)*ul(i,j,k,3) )
            fl(i,j,k,3) = iota*( kx(i)*ul(i,j,k,2) - ky(j)*ul(i,j,k,1) )
            fl(i,j,k,4) = iota * kx(i) * ul(i,j,k,1)
            fl(i,j,k,5) = iota * ky(j) * ul(i,j,k,2)
            fl(i,j,k,6) = iota * kz(k) * ul(i,j,k,3)
          end do
        end do
      end do

      sx = zero; sy = zero; sz = zero; enstrophy = zero; div = zero

      if (.false.) then          ! need to do this in wave space

!.... Convert to physical space (This is expensive, I should be able to compute
!.... most statistics in wave-space)

      do idof = 1, 2*ndof
        call fft3d(1, nx, ny, nz, fl(1,1,1,idof), mpx, py, coef)
      end do

!.... compute the average skewness, enstrophy, and RMS divergence

!.... Note that pointwise averaging in space is entirely equivalent to an 
!.... integral average for discrete Fourier series

      sx1 = zero
      sx2 = zero
      sy1 = zero
      sy2 = zero
      sz1 = zero
      sz2 = zero
      !$omp  parallel do private(k,j,i), &
      !$omp& reduction(+: enstrophy,div,sx1,sx2,sy1,sy2,sz1,sz2)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            sx1 = sx1 + f(i,j,k,4)**3
            sx2 = sx2 + f(i,j,k,4)**2
            sy1 = sy1 + f(i,j,k,5)**3
            sy2 = sy2 + f(i,j,k,5)**2
            sz1 = sz1 + f(i,j,k,6)**3
            sz2 = sz2 + f(i,j,k,6)**2
            enstrophy = enstrophy + f(i,j,k,1)**2 + f(i,j,k,2)**2 + &
                                    f(i,j,k,3)**2
            div = div + (f(i,j,k,4) + f(i,j,k,5) + f(i,j,k,6))**2
          end do
        end do
      end do

      sx1 = sx1 / real(nx*ny*nz)           ! skewness components
      sx2 = sx2 / real(nx*ny*nz)
      if (sx2.ne.zero) then
        sx = -sx1 / sx2**onept5
      else
        sx = zero
      end if

      sy1 = sy1 / real(nx*ny*nz)
      sy2 = sy2 / real(nx*ny*nz)
      if (sy2.ne.zero) then
        sy = -sy1 / sy2**onept5
      else
        sy = zero
      end if

      sz1 = sz1 / real(nx*ny*nz)
      sz2 = sz2 / real(nx*ny*nz)
      if (sz2.ne.zero) then
        sz = -sz1 / sz2**onept5
      else
        sz = zero
      end if

      end if

!.... since the flow is isotropic, might as well average in all directions

      savg = pt33 * ( sx + sy + sz )

      enstrophy = enstrophy / real(nx*ny*nz)
      div = sqrt( div / real(nx*ny*nz) )

      R_L = L * vrms / nu                  ! Integral scale Reynolds number
      R_lambda = lambda * vrms / nu        ! Talor microscale Reynolds number

      write(20,10) t, vrms, urms1, urms2, urms3, &
                   L, lambda, R_L, R_lambda, sx, enstrophy
  10  format(12(1pe13.6,1x))
      call flush(20)

      eps = 5.0 * nu * q**2 / lambda**2    ! turbulence dissipation rate
      ups = (nu * eps)**pt25               ! small-eddy velocity-scale
      eta = (nu**3 / eps)**pt25            ! small-eddy length-scale

      Leps = q**3 / eps                    ! large-eddy length-scale
      Teps = q**2 / eps                    ! large-eddy time-scale

      write(21,10) t, savg, sx, sy, sz, eps, ups, eta, Leps, Teps
      call flush(21)

!.... Interface to the GUI if available

#ifdef GUI
      call sendtoui(istep, t, eps, vrms, urms1, urms2, urms3, ups, L, &
           lambda, eta, Leps, Teps, R_L, R_lambda, enstrophy, savg, sx, sy, sz)
      open(unit=80, file='pid', status='old', action='read', iostat=ierr)
      if (ierr.eq.0) then
        read (80,*) pid
        close(80)
        temp = kill(pid,2)
      end if 
#endif

      return
end subroutine statistics

!=============================================================================
subroutine spectra(n, il, ul, du)
!=============================================================================
!
!.... Computes the 3d and 1d spectra and statistics
!
!.... du is the time derivative with nonlinear term and pressure only
!
!=============================================================================
      use field
      use const
      implicit none
      complex :: ul(mx,ny,nz,ndof), du(mpx,py,pz,ndof)
      integer :: n, il, i, j, k, idof

      logical, external :: mask
      integer :: nshells, shell(nx/2)
      real :: fact(mx), kk(mx), uu(mx), ww(mx), w2(mx), uw(mx)
      real :: udu(mx), wdw(mx)
      real :: sum_ek, sum_dk, sum_d2k, sum_hk, sum_tek, sum_tdk, de
      real :: skewness, R_lambda
      real, allocatable :: sample(:), ek(:), dk(:), d2k(:), hk(:)
      real, allocatable :: tek(:), tdk(:)

      real :: mean(3), uiuj(nx/2), usum(nx/2), xsum
      real :: xrij(nx/2+1,3), yrij(ny/2+1,3), zrij(nz/2+1,3)
      character*80 :: base='spectra', fname
!=============================================================================

      nshells = max(nx,ny,nz)
      allocate( sample(nshells), ek(nshells), dk(nshells), &
                d2k(nshells), hk(nshells), tek(nshells), tdk(nshells) )

      xrij = zero; yrij = zero; zrij = zero
      ek=zero; dk=zero; hk=zero; d2k=zero; sample=zero; tek=zero; tdk=zero

!.... Three-dimensional spectral

     do k = 1, nz
        do j = 1, ny

!         if (mask(1,j,k)) then
!           fact(1) = zero
!         else
!           fact(1) = one
!         end if
!         do i = 2, mx
!           if (mask(i,j,k)) then
!             fact(i) = zero
!           else
!             fact(i) = two
!           end if
!         end do

          fact = two
          fact(1) = one

          do i = 1, nx/2
            kk(i) = kx(i)**2 + ky(j)**2 + kz(k)**2
            shell(i) = 1 + int(sqrt(kk(i)) + pt5)
!           write(*,"(i4,1x,i4)") i, shell(i)
            uu(i) = fact(i) * ( abs(ul(i,j,k,1))**2 + abs(ul(i,j,k,2))**2 + &
                                abs(ul(i,j,k,3))**2 )
            ww(i) = kk(i) * uu(i)
            w2(i) = kk(i) * ww(i)
            uw(i) = fact(i)*two*aimag(kx(i)*ul(i,j,k,2)*conjg(ul(i,j,k,3)) + &
                                      ky(j)*ul(i,j,k,3)*conjg(ul(i,j,k,1)) + &
                                      kz(k)*ul(i,j,k,1)*conjg(ul(i,j,k,2)))
            udu(i) = fact(i) * real( ul(i,j,k,1) * conjg(du(i,j,k,1)) + &
                                     ul(i,j,k,2) * conjg(du(i,j,k,2)) + &
                                     ul(i,j,k,3) * conjg(du(i,j,k,3)) )
            wdw(i) = kk(i) * udu(i)
          end do
          do i = 1, nx/2
            sample(shell(i)) = sample(shell(i)) + fact(i)  ! shell sample
            ek(shell(i))     = ek(shell(i))     + uu(i)    ! 2 * energy sum
            dk(shell(i))     = dk(shell(i))     + ww(i)    ! enstrophy sum
            d2k(shell(i))    = d2k(shell(i))    + w2(i)    ! palenstrophy sum
            hk(shell(i))     = hk(shell(i))     + uw(i)    ! helicity sum
            tek(shell(i))    = tek(shell(i))    + udu(i)   ! energy transfer
            tdk(shell(i))    = tdk(shell(i))    + wdw(i)   ! enstrophy transfer
          end do
        end do

!.... One-dimensional spectra

        fact = two
        fact(1) = one
        do idof = 1, ndof
          usum = zero
          do j = 1, ny
            do i = 1, nx/2
              uiuj(i) = abs( ul(i,j,k,idof) )**2
              uu(i) = fact(i) * uiuj(i)
              usum(i) = usum(i) + uiuj(i)
            end do
            xsum = sum( uiuj(2:nx/2) )
            if(j.eq.1) then
              yrij(j,idof) = yrij(j,idof) + uiuj(1) + 2 * xsum
            else if (j .lt. ny/2+1) then
              yrij(j,idof) = yrij(j,idof) + uiuj(1) + xsum
            else if (j .gt. ny/2+1) then
              yrij(ny+2-j,idof) = yrij(ny+2-j,idof) + xsum
            end if
          end do
          
          xrij(1:nx/2,idof) = xrij(1:nx/2,idof) + usum(1:nx/2)
          
          xsum = sum( usum(2:nx/2) )
          
          if (k.eq.1) then
            zrij(k,idof) =                usum(1) + 2 * xsum
          else if (k.lt.nz/2+1) then
            zrij(i,idof) = zrij(i,idof) + usum(1) + xsum
          else if (k.gt.nx/2+1) then
            zrij(nz+2-k,idof) = zrij(nz+2-k,idof) + xsum
          endif
        end do

      end do

!.... compute final statistics

      do idof = 1, 3
        mean(idof) = xrij(1,idof) + 2 * sum( xrij(2:nx/2,idof) )
      end do

      sum_ek  = sum( ek )
      sum_dk  = sum( dk )
      sum_d2k = sum( d2k )
      sum_hk  = sum( hk )
      sum_tek = sum( tek )
      sum_tdk = sum( tdk )

!.... write out to screen standard quantities 

      if (n.eq.1) then
        write(*,"('Step    time      sum_ek      sum_dk    sum_hk    ', &
                 &'skewness    R_lambda  power_in    sum_tke    ', &
                 &'sum_tdk      de')")
      endif

      skewness = -6.0*sqrt(15.0)/7.0 * sum_tdk / sum_dk**1.5
      if (nu.ne.0) then
        de = 1.0/(1.0-(-2.0*nu)*(sum_tdk-nu*sum_d2k)*0.5*sum_ek/ &
             (-nu*sum_dk)**2)
        R_lambda = sqrt(20.0/(3.0*(2.0*nu*sum_dk/2.0)*nu))*(sum_ek/2.0)
      else
        de = zero
        R_lambda = zero
      end if

      write(*,"(i4,10(1pe11.3))") il, t, sum_ek, sum_dk, sum_hk, &
                                  skewness, R_lambda, zero, sum_tek, &
                                  sum_tdk, de

      write(22,"(i4,10(1x,1pe20.12))") il, t, sum_ek, sum_dk, sum_hk, &
                                       skewness, R_lambda, zero, sum_tek, &
                                       sum_tdk, de
      call flush(22)

!.... output 3d spectra

      if (mod(il,nout).eq.0) then
        call makename(base,il,fname)
        open(10,file=fname)

!       write(10,"('# ',i4,10(1x,1pe13.6))") il, t, mean
!       do i = 1, nx/2
!         write(10,"(i4,10(1x,1pe13.6))") i, xrij(i,:), yrij(i,:), zrij(i,:)
!       end do
!       close(10)

        write(10,"('# ',i4,10(1x,1pe13.6))") il, t, mean
        do i = 1, nshells
          if (sample(i).ne.zero) then
            write(10,"(i4,10(1x,1pe13.6))") i-1, 4.0*pi*((i-1)**2)/ &
                                            sample(i)*(0.5*ek(i)),  &
                                            4.0*pi*((i-1)**2)/ &
                                            sample(i)*dk(i)*nu
          end if
        end do
        close(10)

      endif

      deallocate( sample, ek, dk, d2k, hk, tek, tdk )

      return
end subroutine spectra

subroutine wgrid( name, nx, ny, nz, x, y, z )
      implicit none
      character*80 :: name
      integer :: nx, ny, nz
      real :: x(nx), y(ny), z(nz)
      integer :: i, j, k, idof

      open(10,file=name,form='unformatted')
      write(10) int(nx,4), int(ny,4), int(nz,4)
      write(10) &
        (((real(x(i),4), i=1,nx), j=1,ny), k=1,nz), &
        (((real(y(j),4), i=1,nx), j=1,ny), k=1,nz), &
        (((real(z(k),4), i=1,nx), j=1,ny), k=1,nz)
      close(10)

      return
end subroutine wgrid

subroutine wdata( name, ndof, nx, ny, nz, u, fsmach, alpha, re, time)
      use const
      implicit none
      character*80 :: name
      integer :: ndof, nx, ny, nz
      real :: u( 2*((nx+2)/2), ny, nz, ndof ), fsmach, alpha, re, time
      integer :: i, j, k, idof

!.... for an incompressible field, rho and T are 1.0

      open(10,file=name,form='unformatted')
      write(10) int(nx,4), int(ny,4), int(nz,4)
      write(10) real(fsmach,4), real(alpha,4), real(re,4), real(time,4)
      write(10) &
        ((( real(one,4)       , i=1,nx), j=1,ny), k=1,nz ), &
        ((( real(u(i,j,k,1),4), i=1,nx), j=1,ny), k=1,nz ), &
        ((( real(u(i,j,k,2),4), i=1,nx), j=1,ny), k=1,nz ), &
        ((( real(u(i,j,k,3),4), i=1,nx), j=1,ny), k=1,nz ), &
        ((( real(one,4)       , i=1,nx), j=1,ny), k=1,nz )
      close(10)

      return
end subroutine wdata

subroutine error(name,msg)
      implicit none

      integer loc
      character*80 name, msg

      loc = index(name,'$')-1
      write(*,"(/,'*****************************************************')")
      write(*,"('Error in --> ',a)") name(1:loc)
      loc = index(msg,'$')-1
      write(*,"('-----------> ',a)") msg(1:loc)
      write(*,"('*****************************************************',/)")

      stop
end subroutine error

#ifdef CRAY

!.... USE CRAY FFT routines

subroutine fft3di(nx, ny, nz, coef)
      implicit none
      integer :: nx, ny, nz
      real :: coef(100 + 2*(nx + ny + nz)), dum
      call scfft3d(0, nx, ny, nz, 0.0, dum, 1, 1, dum, 1, 1, coef, dum, 0)
      return
end subroutine fft3di

module ffttmp
      integer :: isys = 1
      real, allocatable :: work(:)
end module ffttmp

subroutine mkwork( nx, ny, nz )
      use ffttmp
      implicit none
      integer :: nx, ny, nz
      integer, parameter :: NCPUS=4

      if (isys.eq.0) then
        allocate( work(512 * max(nx,ny,nz)) )
      else
        allocate( work(4 * NCPUS * max(nx*ny,ny*nx,nz*nx)) )
      end if
end subroutine mkwork

subroutine fft3d(sign, nx, ny, nz, u, ldx, ldy, coef)
      use const
      use ffttmp
      implicit none
      integer :: sign, nx, ny, nz, ldx, ldy
      real :: u(2*ldx, ldy, nz), coef(100 + 2*(nx + ny + nz)), alpha
      integer :: i, j, k, il
      if (sign.eq.-1) then
        alpha = one/real(nx*ny*nz)        
        call scfft3d(sign, nx, ny, nz, alpha, u, 2*ldx, ldy, u, ldx, ldy, &
                     coef, work, isys)
      else
        alpha = one
        call csfft3d(sign, nx, ny, nz, alpha, u, ldx, ldy, u, 2*ldx, ldy, &
                     coef, work, isys)
      end if
      return
end subroutine fft3d

#elif defined(SGI_FFT3D)

!.... Use SGI three-dimensional FFT routines (slow?)

subroutine fft3di(nx, ny, nz, coef)
      implicit none
      integer :: nx, ny, nz
      real :: coef((nx+15)+2*(ny+15)+2*(nz+15))
#ifdef R8
      call dzfft3dui( nx, ny, nz, coef )
#else
      call scfft3dui( nx, ny, nz, coef )
#endif
      return
end subroutine fft3di

subroutine fft3d(sign, nx, ny, nz, u, ldx, ldy, coef)
      use const
      implicit none
      integer :: sign, nx, ny, nz, ldx, ldy
      real :: u(2*ldx, ldy, nz), coef( (nx+15)+2*(ny+15)+2*(nz+15) ), alpha
      integer :: i, j, k, il
#ifdef R8
      if (sign.eq.-1) then
        call dzfft3du(sign, nx, ny, nz, u, 2*ldx, ldy, coef)
        alpha = one/real(nx*ny*nz)
        call dscal3d(nx, ny, nz, alpha, u, 2*ldx, ldy)
      else
        call zdfft3du(sign, nx, ny, nz, u, 2*ldx, ldy, coef)
      end if
#else
      if (sign.eq.-1) then
        call scfft3du(sign, nx, ny, nz, u, 2*ldx, ldy, coef)
        alpha = one/real(nx*ny*nz)
        call sscal3d(nx, ny, nz, alpha, u, 2*ldx, ldy)
      else
        call csfft3du(sign, nx, ny, nz, u, 2*ldx, ldy, coef)
      end if
#endif
      return
end subroutine fft3d

#elif defined(SGI_FFT)

!.... Use SGI one-dimensional FFT routines

subroutine fft3di(nx, ny, nz, coef)
      implicit none
      integer :: nx, ny, nz
      real :: coef((nx+15)+2*(ny+15)+2*(nz+15))
#ifdef R8
      call dzfft1dui( nx, coef(1) )
      call zfft1di( ny, coef(nx+15+1) )
      call zfft1di( nz, coef((nx+15)+2*(ny+15)+1) )
#else
      call scfft1dui( nx, coef(1) )
      call cfft1di( ny, coef(nx+15+1) )
      call cfft1di( nz, coef((nx+15)+2*(ny+15)+1) )
#endif
      return
end subroutine fft3di

subroutine fft3d(sign, nx, ny, nz, u, ldx, ldy, coef)
      use const
      implicit none
      integer :: sign, nx, ny, nz, ldx, ldy, mx
      real :: u(2*ldx, ldy, nz), coef( (nx+15)+2*(ny+15)+2*(nz+15) )
      integer :: i, j, k, il
      mx = (nx+2)/2

!.... sign.eq.-1 is the forward transform, sign.eq.1 is the inverse transform

#ifdef R8
      if (sign .eq. -1) then
#ifdef DEBUG
        do i = 1, nx
          do j = 1, ny
            do k = 1, nz
              write(*,10) i, j, k, u(i,j,k)
  10          format(3(i4,1x),1pe13.6)
            end do
          end do
        end do
#endif
        !$omp parallel do private(k,j), shared(sign, nx)
        do k = 1, nz
          do j = 1, ny
            call dzfft1du( sign, nx, u(1,j,k), 1, coef(1) )
            call dscal1d( nx, (one/real(nx)), u(1,j,k), 1 )
          end do
        end do
        !$omp parallel do private(k,i,il), shared(sign, ny, ldx)
        do k = 1, nz
          do i = 1, mx
            il = 1 + (i-1)*2
            call zfft1d( sign, ny, u(il,1,k), ldx, coef(nx+15+1) )
            call zscal1d( ny, (one/real(ny)), u(il,1,k), ldx )
          end do
        end do
        !$omp parallel do private(j,i,il), shared(sign, nz, ldx)
        do j = 1, ny
          do i = 1, mx
            il = 1 + (i-1)*2
            call zfft1d( sign, nz, u(il,j,1), ldx*ldy, &
                         coef((nx+15)+2*(ny+15)+1) )
            call zscal1d( nz, (one/real(nz)), u(il,j,1), ldx*ldy )
          end do
        end do
      else if (sign .eq. 1) then
        !$omp parallel do private(j,i,il), shared(sign, nz, ldx, ldy)
        do j = 1, ny
          do i = 1, mx
            il = 1 + (i-1)*2
            call zfft1d( sign, nz, u(il,j,1), ldx*ldy, &
                         coef((nx+15)+2*(ny+15)+1) )
          end do
        end do
        !$omp parallel do private(k,i,il), shared(sign, ny, ldx)
        do k = 1, nz
          do i = 1, mx
            il = 1 + (i-1)*2
            call zfft1d( sign, ny, u(il,1,k), ldx, coef(nx+15+1) )
          end do
        end do
        !$omp parallel do private(k,j), shared(sign, nx)
        do k = 1, nz
          do j = 1, ny
            call zdfft1du( sign, nx, u(1,j,k), 1, coef(1) )
          end do
        end do
      else
        call error('fft3d$','Illegal value of sign$')
      end if
#else
      if (sign .eq. -1) then
        !$omp parallel do private(k,j), shared(sign, nx)
        do k = 1, nz
          do j = 1, ny
            call scfft1du( sign, nx, u(1,j,k), 1, coef(1) )
            call sscal1d( nx, (one/real(nx)), u(1,j,k), 1 )
          end do
        end do
        !$omp parallel do private(k,i,il), shared(sign, ny, ldx)
        do k = 1, nz
          do i = 1, mx
            il = 1 + (i-1)*2
            call cfft1d( sign, ny, u(il,1,k), ldx, coef(nx+15+1) )
            call cscal1d( ny, (one/real(ny)), u(il,1,k), ldx )
          end do
        end do
        !$omp parallel do private(j,i,il), shared(sign, nz, ldx)
        do j = 1, ny
          do i = 1, mx
            il = 1 + (i-1)*2
            call cfft1d( sign, nz, u(il,j,1), ldx*ldy, &
                         coef((nx+15)+2*(ny+15)+1) )
            call cscal1d( nz, (one/real(nz)), u(il,j,1), ldx*ldy )
          end do
        end do
      else if (sign .eq. 1) then
        !$omp parallel do private(j,i,il), shared(sign, nz, ldx, ldy)
        do j = 1, ny
          do i = 1, mx
            il = 1 + (i-1)*2
            call cfft1d( sign, nz, u(il,j,1), ldx*ldy, &
                         coef((nx+15)+2*(ny+15)+1) )
          end do
        end do
        !$omp parallel do private(k,i,il), shared(sign, ny, ldx)
        do k = 1, nz
          do i = 1, mx
            il = 1 + (i-1)*2
            call cfft1d( sign, ny, u(il,1,k), ldx, coef(nx+15+1) )
          end do
        end do
        !$omp parallel do private(k,j), shared(sign, nx)
        do k = 1, nz
          do j = 1, ny
            call csfft1du( sign, nx, u(1,j,k), 1, coef(1) )
          end do
        end do
      else
        call error('fft3d$','Illegal value of sign$')
      end if
#endif
      return
end subroutine fft3d

#elif defined(INTEL_MKL)

!.... Use Intel MKL FFT routines (only power of 2 -- dealiasing won't work)

subroutine fft3di(nx, ny, nz, coef)
      implicit none
      integer :: nx, ny, nz
      real :: coef( 2*nx+4 + 3*ny + 3*nz )
      real :: tmp(1)
#ifdef R8
      call dzfft1d( tmp, nx, 0, coef(1) )
      call zfft1d ( tmp, ny, 0, coef(2*nx+4 + 1) )
      call zfft1d ( tmp, nz, 0, coef(2*nx+4 + 3*ny + 1)
#else
      integer isign
      isign = 0
      write(*,*) nx, ny, nz
      call scfft1d( coef(1), nx, isign, coef(1) )
      call cfft1d ( coef(1), ny, isign, coef(2*nx+4 + 1) )
      call cfft1d ( coef(1), nz, isign, coef(2*nx+4 + 3*ny + 1) )
#endif
      return
end subroutine fft3di

subroutine fft3d(sign, nx, ny, nz, u, ldx, ldy, coef)
      use const
      implicit none
      integer :: sign, nx, ny, nz, ldx, ldy
      real :: u(2*ldx, ldy, nz), coef( (nx+15)+2*(ny+15)+2*(nz+15) ), alpha
      integer :: i, j, k, il

!.... skeleton, need to implement this...

      return
end subroutine fft3d

#elif defined(FFTW)

!.... Use FFTW routines (need to update to threaded versions)

module ffttmp
      integer FFTW_FORWARD,FFTW_BACKWARD
      parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

      integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
      parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

      integer FFTW_ESTIMATE,FFTW_MEASURE
      parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

      integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
      parameter (FFTW_OUT_OF_PLACE=0)
      parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

      integer FFTW_THREADSAFE
      parameter (FFTW_THREADSAFE=128)

      ! Constants for the MPI wrappers:
      integer FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
      integer FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
      parameter(FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
      parameter(FFTW_SCRAMBLED_INPUT=8192)
      parameter(FFTW_SCRAMBLED_OUTPUT=16384)

      ! make sure that the pointer is large enough
      integer(8) :: fplan(2), bplan(2)
end module ffttmp

!> Initialize FFTW for 3d FFT
subroutine fft3di(nx, ny, nz, rplan)
      use ffttmp
      implicit none
      integer :: nx, ny, nz
      real :: rplan
      integer :: iplan
      iplan = rplan
      call rfftw3d_f77_create_plan( fplan(iplan), nx, ny, nz, &
        FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE + FFTW_IN_PLACE )
      call rfftw3d_f77_create_plan( bplan(iplan), nx, ny, nz, &
        FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE + FFTW_IN_PLACE )
      return
end subroutine fft3di

!> \note If you are using OpenMP, then make sure to build FFTW with 
!! --enable-threads --with-openmp
!! \note You can also using --enable-float to make a single precision version
subroutine fft3d(sign, nx, ny, nz, u, ldx, ldy, rplan)
      use const
      use ffttmp
      use timemod
      implicit none
      integer :: sign, nx, ny, nz, ldx, ldy
      real :: u(2*ldx, ldy, nz), rplan, alpha
      integer :: i, j, k, il, iplan
      iplan = rplan
      if (sign.eq.-1) then
#ifdef DEBUG
        do i = 1, nx
          do j = 1, ny
            do k = 1, nz
              write(*,10) i, j, k, u(i,j,k)
  10          format(3(i4,1x),1pe13.6)
            end do
          end do
        end do
#endif
#ifndef _OPENMP 
        call rfftwnd_f77_one_real_to_complex( fplan(iplan), u, u )
#else
        call rfftwnd_f77_threads_one_real_to_complex(nthreads, &
        fplan(iplan), u, u )
#endif
        alpha = one/real(nx*ny*nz)
        u = u * alpha
      else
#ifndef _OPENMP 
        call rfftwnd_f77_one_complex_to_real( bplan(iplan), u, u )
#else
        call rfftwnd_f77_threads_one_complex_to_real(nthreads, &
        bplan(iplan), u, u )
#endif
      end if
      return
end subroutine fft3d

#endif

subroutine pad( ndof, nx, ny, nz, u, px, py, pz, up)

!.... Takes a field u and pads it into up 
!.... (u and up cannot share the same memory)

      use const
      implicit none
      integer :: ndof, nx, ny, nz, px, py, pz
      complex :: u((nx+2)/2,ny,nz,ndof), up((px+2)/2,py,pz,ndof)
      integer :: i, j, k, idof, mx, mpx, jl, kl

      mx  = (nx+2)/2
      mpx = (px+2)/2

      do idof = 1, ndof

        !$omp parallel do private(i,j,k,jl,kl)
        do k = 1, (nz+2)/2
          do j = 1, (ny+2)/2
            do i = 1, mx
              up(i,j,k,idof) = u(i,j,k,idof)
            end do
            do i = mx+1,mpx
              up(i,j,k,idof) = czero
            end do
          end do
          do j = ny, (ny+2)/2+1, -1
            jl = py - (ny-j)
            do i = 1, mx
              up(i,jl,k,idof) = u(i,j,k,idof)
            end do
            do i = mx+1, mpx
              up(i,jl,k,idof) = czero
            end do
          end do
          do j = (ny+2)/2+1, py-ny+(ny+2)/2
            do i = 1, mpx
              up(i,j,k,idof) = czero
            end do
          end do
        end do

        !$omp parallel do private(i,j,k,jl,kl)
        do k = nz, (nz+2)/2+1, -1
          kl = pz - (nz-k) 
          do j = 1, (ny+2)/2
            do i = 1, mx
              up(i,j,kl,idof) = u(i,j,k,idof)
            end do
            do i = mx+1,mpx
              up(i,j,kl,idof) = czero
            end do
          end do
          do j = ny, (ny+2)/2+1, -1
            jl = py - (ny-j)
            do i = 1, mx
              up(i,jl,kl,idof) = u(i,j,k,idof)
            end do
            do i = mx+1, mpx
              up(i,jl,kl,idof) = czero
            end do
          end do
          do j = (ny+2)/2+1, py-ny+(ny+2)/2
            do i = 1, mpx
              up(i,j,kl,idof) = czero
            end do
          end do
        end do

        !$omp parallel do private(i,j,k,jl,kl)
        do k = (nz+2)/2+1, pz-nz+(nz+2)/2
          do j = 1, py
            do i = 1, mpx
              up(i,j,k,idof) = czero
            end do
          end do
        end do

      end do

      return
end subroutine pad

subroutine unpad( ndof, nx, ny, nz, px, py, pz, up)

!.... This unpads the field up in place

      use const
      implicit none
      integer :: ndof, nx, ny, nz, px, py, pz
      complex :: up((px+2)/2,py,pz,ndof)
      integer :: i, j, k, idof, mx, mpx, jl, kl

      mx  = (nx+2)/2
      mpx = (px+2)/2

      do idof = 1, ndof

        !$omp parallel do private(i,j,k,jl,kl)
        do k = 1, (nz+2)/2
          do j = 1, (ny+2)/2
            do i = 1, mx
              up(i,j,k,idof) = up(i,j,k,idof)
            end do
            do i = mx+1,mpx
              up(i,j,k,idof) = czero
            end do
          end do
          do j = (ny+2)/2+1, ny
            jl = py - (ny-j)
            do i = 1, mx
              up(i,j,k,idof) = up(i,jl,k,idof)
            end do
            do i = mx+1, mpx
              up(i,j,k,idof) = czero
            end do
          end do
          do j = ny+1, py
            do i = 1, mpx
              up(i,j,k,idof) = czero
            end do
          end do
        end do

        !$omp parallel do private(i,j,k,jl,kl)
        do k = (nz+2)/2+1, nz
          kl = pz - (nz-k)
          do j = 1, (ny+2)/2
            do i = 1, mx
              up(i,j,k,idof) = up(i,j,kl,idof)
            end do
            do i = mx+1,mpx
              up(i,j,k,idof) = czero
            end do
          end do
          do j = (ny+2)/2+1, ny
            jl = py - (ny-j)
            do i = 1, mx
              up(i,j,k,idof) = up(i,jl,kl,idof)
            end do
            do i = mx+1, mpx
              up(i,j,k,idof) = czero
            end do
          end do
          do j = ny+1, py
            do i = 1, mpx
              up(i,j,kl,idof) = czero
            end do
          end do
        end do

        !$omp parallel do private(i,j,k,jl,kl)
        do k = nz+1, pz
          do j = 1, py
            do i = 1, mpx
              up(i,j,k,idof) = czero
            end do
          end do
        end do

      end do

      return
end subroutine unpad

! Note:  to list available defines use 
!   `echo "" | gfortran -dM -E -fopenmp - | sort`
subroutine interupt
      use intmod
      implicit none
      integer(4) :: flag=-1, i
      integer(4), parameter :: SIGCONT=25, SIGINT=2, SIGKILL=9, SIGTERM=15
#if __GNUC__ || __GFORTRAN__
      external :: handler
      intrinsic signal
      halt = .false.
      i = signal( SIGINT,  handler )
#else
      integer(4), external :: signal, handler
      halt = .false.
      i = signal( SIGCONT, handler, flag )
#endif
      return
end subroutine interupt

subroutine handler
      use intmod
      implicit none
      if (halt) call exit(0)
      write(*,"('Interupt:  attempting to finish current time-step...')")
      write(*,"('  An additional Ctrl-C will terminate immediately')")
      halt = .true.
      return
end subroutine handler

subroutine cprint(ndof, nx, ny, nz, u)
      use const
      implicit none
      integer :: nx, ny, nz, ndof
      complex :: u( (nx+2)/2, ny, nz, ndof ), utmp
      integer :: i, j, k, idof

      write(*,*)
      do idof = 1, ndof
        do k = 1, nz
          do j = 1, ny
            do i = 1, (nx+2)/2
              if (abs(u(i,j,k,idof)) .gt. 1.0e-12) then
                utmp = u(i,j,k,idof)
              else
                utmp = czero
              end if
              write(*,"(4(1x,i5),6(1x,1pe11.4))") i, j, k, idof, utmp
            end do
          end do
        end do
      end do

      return
end subroutine cprint

subroutine rprint(ndof, nx, ny, nz, u)
      use const
      implicit none
      integer :: nx, ny, nz, ndof
      real :: u(2*((nx+2)/2), ny, nz, ndof ), utmp
      integer :: i, j, k, idof

      write(*,*)
      do idof = 1, ndof
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if (abs(u(i,j,k,idof)) .gt. 1.0e-12) then
                utmp = u(i,j,k,idof)
              else
                utmp = zero
              end if
              write(*,"(4(1x,i5),6(1x,1pe11.4))") i, j, k, idof, utmp
            end do
          end do
        end do
      end do

      return
end subroutine rprint

subroutine oddball( rl )
      use field
      use const
      implicit none
      complex :: rl(mx,ny,nz,ndof)
      integer :: i, j, k, idof

      if (mod(nx,2) .eq. 0) then
        !$omp parallel do private(idof,k,j)
        do idof = 1, ndof
          do k = 1, nz
            do j = 1, ny
              rl(mx,j,k,idof) = czero
            end do
          end do
        end do
      end if

      if (mod(ny,2) .eq. 0) then
        !$omp parallel do private(idof,k,i)
        do idof = 1, ndof
          do k = 1, nz
            do i = 1, mx
              rl(i,my,k,idof) = czero
            end do
          end do
        end do
      end if
        
      if (mod(nz,2) .eq. 0) then
        !$omp parallel do private(idof,j,i)
        do idof = 1, ndof
          do j = 1, ny
            do i = 1, mx
              rl(i,j,mz,idof) = czero
            end do
          end do
        end do
      end if

      return
end subroutine oddball

subroutine compute_q( ul, fl )
      use field
      use const
      implicit none
      complex :: ul(mx,ny,nz,ndof), fl(mpx,py,pz,ndof)
      integer :: i, j, k, idof
      real :: ek
      
!.... compute the rms velocity

      ek = zero
      do k = 1, nz
        do j = 1, ny
          i = 1
          ek = ek + ( abs(ul(i,j,k,1))**2 + abs(ul(i,j,k,2))**2 + &
                      abs(ul(i,j,k,3))**2 )
          do i = 2, mx-1
            ek = ek + two * ( abs(ul(i,j,k,1))**2 + abs(ul(i,j,k,2))**2 + &
                              abs(ul(i,j,k,3))**2 )
          end do
        end do
      end do

      write(*,*) t, ek

      return
      end

subroutine compute_div( ul )
      use field
      use const
      implicit none
      complex :: ul(mx,ny,nz,ndof)
      integer :: i, j, k, idof
      complex :: div, ctmp
      
      div = czero
      do k = 1, nz
        do j = 1, ny
          do i = 1, mx
            ctmp = kx(i)*ul(i,j,k,1) + ky(j)*ul(i,j,k,2) + kz(k)*ul(i,j,k,3)
            if (abs(ctmp).gt.abs(div)) div = ctmp
          end do
        end do
      end do
      write(*,*) 'Maximum divergence in initial field ',div

      return
end subroutine compute_div

real function E0(k)

!.... Compute a model turbulence spectrum (same as used by A. Wray)

      use field
      use const
      implicit none
      real :: k, arg
      integer :: i

!.... Comte-Bellot & Corrsin spectra parameters, the value of v_0 is used
!.... to match to the initial eneergy of the experimental spectra

      real :: a(0:6), c1, c2
!     real :: l = 50.8, m = 5.08, t0=42.0, u_inf=1000.0  ! matches JFM

      real :: l = 55.0, m = 5.08, t0=42.0, u_inf=1000.0  ! from Germiso
!.... v0 = 27.1893

      if (ispec.eq.0) then
        if (k .ge. 4.0 .and. k .lt. 6) then
          E0 = 0.75
        else
          E0 = zero
        end if
      else if (ispec.eq.1) then
        E0 = 16.0 * sqrt(two/pi) * v0**2 * kmax**(-5) * k**4 * &
             exp(-two*(k/kmax)**2)
      else if (ispec.eq.2) then
        a(0) =  0.56102E+01;   a(1) = -0.11236E+01;   a(2) = -0.30961E+00
        a(3) =  0.33172E+00;   a(4) = -0.10959E+00;   a(5) = -0.22320E-01   
        a(6) =  0.66575E-02  
        
        c1 = 2.0 * pi / l
        c2 = c1 / (v0**2)

        if (k.eq.zero) then
          E0 = zero
        else
          arg = zero
          do i = 0, 6
            arg = arg + a(i)*log(c1*k)**i
          end do
          E0 = c2 * exp(arg)
        end if
      else if (ispec.eq.3) then
        c1 = one / ( one - int( sqrt(two)/3.0 * nx )**(-2.0/3.0) )
        if (k.eq.zero) then
          E0 = zero
        else
          E0 = c1 * k**(-5.0/3.0)
        end if
      else
        call error("E0$","Illegal value of ispec$")
      end if

      return
end function E0

real function I_E0_dk()

!.... Compute a model turbulence spectrum (same as used by A. Wray)

      use field
      use const
      implicit none
      real :: k

      I_E0_dk = 1.5
      return

end function I_E0_dk

real function I_E0k2_dk()

!.... Compute a model turbulence spectrum (same as used by A. Wray)

      use field
      use const
      implicit none
      real :: k
      real, external :: I_E0_dk

      I_E0k2_dk = kmax**2 * 5.0/4.0 * I_E0_dk()
      return

end function I_E0k2_dk

subroutine makename(base,iver,fname)

!.... put a version number on a filename

      character*80 base, fname

      length = index(base,' ')
      fname = base
      if (iver .lt. 10) then
        write(fname(length:80),"('.',i1)") iver
      else if (iver .lt. 100) then
        write(fname(length:80),"('.',i2)") iver
      else if (iver .lt. 1000) then
        write(fname(length:80),"('.',i3)") iver
      else if (iver .lt. 10000) then
        write(fname(length:80),"('.',i4)") iver
      else if (iver .lt. 100000) then
        write(fname(length:80),"('.',i5)") iver
      else if (iver .lt. 1000000) then
        write(fname(length:80),"('.',i6)") iver
      else
        call error('makename$','version number too large$')
      end if

      return
end subroutine makename
