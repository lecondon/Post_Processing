subroutine sieve(is_prime, n_max)
! =====================================================
! Uses the sieve of Eratosthenes to compute a logical
! array of size n_max, where .true. in element i
! indicates that i is a prime.
! =====================================================
    integer, intent(in)   :: n_max
    logical, intent(out)  :: is_prime(n_max)
    integer :: i
    is_prime = .true.
    is_prime(1) = .false.
    do i = 2, int(sqrt(real(n_max)))
        if (is_prime (i)) is_prime (i * i : n_max : i) = .false.
    end do
    return
end subroutine

subroutine logical_to_integer(prime_numbers, is_prime, num_primes, n)
! =====================================================
! Translates the logical array from sieve to an array
! of size num_primes of prime numbers.
! =====================================================
    integer                 :: i, j=0
    integer, intent(in)     :: n
    logical, intent(in)     :: is_prime(n)
    integer, intent(in)     :: num_primes
    integer, intent(out)    :: prime_numbers(num_primes)
    do i = 1, size(is_prime)
        if (is_prime(i)) then
            j = j + 1
            prime_numbers(j) = i
        end if
    end do
end subroutine

subroutine clmdaily(hr_st, hr_ed,nx, ny, nz, nclm, leap, clmdr,dname,run)
! =====================================================
! Calculated Daily Average from PARFLOW CLM OUT
! Seeing how this works
! =====================================================
    ! Input variables
    integer, intent(in)     :: leap
    integer, intent(in)     :: hr_st, hr_ed
    integer, intent(in)     :: nx, ny, nz, nclm
    character*200, intent(in)   :: clmdr, dname, run

    real*8  :: dx2, dy2, dz2
    integer  :: i, j, k, nni, nnj, nnk, ix, iy, iz, ns, rx, ry, rz
    integer  :: nnx, nny, nnz    
    integer  :: ijk, namelength, xtent,ytent,ztent, cnt
    integer :: t,counter, day, month, hour, c1, dayofyear, out
    real :: t1, t2, t3
    ! output variables
    real*8, dimension(nx,ny,nz)     :: press, sat
    real*8, dimension(nx,ny,nclm)   :: CLM
    ! daily variables
    real*8, dimension(nx,ny) :: flow, sm, WTd, Storage
    real*8, dimension(nx,ny) :: SWE, Tsoil, trans_veg, LH, SH
    real*8, dimension(nx,ny) :: Tgrnd, Irr, evap_grnd, LWout, evap_tot
    real*8, dimension(nx,ny) :: evap_soi, netrad, evap_veg
    ! hourly variables
    real*8, dimension(nx,ny) :: SWEh, Tsoilh, trans_vegh, LHh, SHh, Tgrndh
    real*8, dimension(nx,ny) :: evap_grndh, LWouth, evap_toth, evap_soih
    real*8, dimension(nx,ny) :: netradh, evap_vegh, Irrh
    ! monthly variables
    real*8, dimension(nx,ny) :: SWEm, Tsoilm, trans_vegm, LHm, SHm, Tgrndm
    real*8, dimension(nx,ny) :: evap_grndm, LWoutm, evap_totm, evap_soim
    real*8, dimension(nx,ny) :: netradm, evap_vegm, Irrm
    ! yearly variables
    real*8, dimension(nx,ny) :: SWEy, Tsoily, trans_vegy, LHy, SHy, Tgrndy
    real*8, dimension(nx,ny) :: evap_grndy, LWouty, evap_toty, evap_soiy
    real*8, dimension(nx,ny) :: netrady, evap_vegy, Irry

    integer*4, dimension(12) :: days
    character*200 :: pname, fname, filenum
    character( len= 3), dimension( 12) ::  monthname = [ 'OCT', 'NOV', 'DEC', &
    'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN','JUL','AUG', 'SEP']

    ! set month days based on is_leap
    if (leap.eq.1) then
        days = [ 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30]
    else
        days = [ 31, 30, 31, 31, 29, 31, 30, 31, 30, 31, 31, 30]
    end if


    !Process Hourly
    cnt = 1
    do out = hr_st,hr_ed

    ! clear out hourly data
    SWEh = 0.0d0
    Tsoilh = 0.0d0
    trans_vegh = 0.0d0
    LHh = 0.0d0
    SHh = 0.0d0
    Tgrndh = 0.0d0
    evap_grndh = 0.0d0
    LWouth = 0.0d0
    evap_toth = 0.0d0
    evap_soih = 0.0d0
    netradh = 0.0d0
    evap_vegh = 0.0d0
    Irrh = 0.0d0

    pname=trim(adjustl(clmdr))//trim(adjustl(run))//'.out.clm_output'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.C.pfb'
    print*,trim(adjustl(clmdr))
    print*,'pname is'
    print*,pname
    print*,'fname is'
    print*,fname
    call pfb_read(CLM,fname,nx,ny,nz)

    ! assign fluxes for hourly CLM outputs
    SWEh(:,:) = SWEh(:,:) + CLM(:,:,11)
    Tsoilh(:,:) = Tsoilh(:,:) + CLM(:,:,14)
    trans_vegh(:,:) = trans_vegh(:,:) + CLM(:,:,9)
    LHh(:,:) = LHh(:,:) + CLM(:,:,1)
    SHh(:,:) = SHh(:,:) + CLM(:,:,3)
    Tgrndh(:,:) = Tgrndh(:,:) + CLM(:,:,12)
    evap_grndh(:,:) = evap_grndh(:,:) + CLM(:,:,6)
    LWouth(:,:) = LWouth(:,:) + CLM(:,:,2)
    evap_toth(:,:) = evap_toth(:,:) + CLM(:,:,5)
    evap_soih(:,:) = evap_soih(:,:) + CLM(:,:,7)
    netradh(:,:) = netradh(:,:) + CLM(:,:,4)
    evap_vegh(:,:) = evap_vegh(:,:) + CLM(:,:,8)
    Irrh(:,:) = Irrh(:,:) + CLM(:,:,12)

    pname = 'SWE.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) SWEh
    close(10)

    pname = 'Tsoil.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Tsoilh
    close(10)

    pname = 'trans_veg.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) trans_vegh
    close(10)

    pname = 'LH.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) LHh
    close(10)

    pname = 'SH.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) SHh
    close(10)

    pname = 'Tgrnd.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Tgrndh
    close(10)

    pname = 'evap_grnd.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_grndh
    close(10)

    pname = 'LWout.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) LWouth
    close(10)

    pname = 'evap_tot.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_toth
    close(10)

    pname = 'evap_soi.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_soih
    close(10)

    pname = 'netrad.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) netradh
    close(10)

    pname = 'evap_veg.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_vegh
    close(10)

    pname = 'Irr.hourly'
    write(filenum,'(i5.5)') out
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Irrh
    close(10)

    print*, cnt
    cnt = cnt+1

    end do

    ! set running counter for PF file number
    counter = 1
    dayofyear = 1

    ! clear out yearly avgs
    SWEy = 0.0d0
    Tsoily =  0.0d0
    trans_vegy =  0.0d0
    LHy =  0.0d0
    SHy = 0.0d0
    Tgrndy = 0.0d0
    evap_grndy = 0.0d0
    LWouty = 0.0d0
    evap_toty = 0.0d0
    evap_soiy = 0.0d0
    netrady = 0.0d0
    evap_vegy = 0.0d0
    Irry = 0.0d0

    do month = 1, 12

    ! clear out monthly avgs
    SWEm = 0.0d0
    Tsoilm =  0.0d0
    trans_vegm =  0.0d0
    LHm =  0.0d0
    SHm = 0.0d0
    Tgrndm = 0.0d0
    evap_grndm = 0.0d0
    LWoutm = 0.0d0
    evap_totm = 0.0d0
    evap_soim = 0.0d0
    netradm = 0.0d0
    evap_vegm = 0.0d0
    Irrm = 0.0d0

    ! loop over month
    do day = 1, days(month)

    ! zero out daily avg

    ! clear out daily avgs
    SWE = 0.0d0
    Tsoil = 0.0d0
    trans_veg = 0.0d0
    LH = 0.0d0
    SH = 0.0d0
    Tgrnd = 0.0d0
    evap_grnd = 0.0d0
    LWout = 0.0d0
    evap_tot = 0.0d0
    evap_soi = 0.0d0
    netrad = 0.0d0
    evap_veg = 0.0d0
    Irr = 0.0d0

    print*,'is this working yet?'

    CALL CPU_TIME(t1)
    ! loop over day
    do hour = 1, 24

    pname=trim(adjustl(clmdr))//trim(adjustl(run))//'.out.clm_output'
    write(filenum,'(i5.5)') counter
    fname=trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.C.pfb'
    print*,'fname is'
    print*,fname
    call pfb_read(CLM,fname,nx,ny,nz)

    print*,'looped over days'

    ! assign fluxes and states from single CLM file
    SWE(:,:) = SWE(:,:) + CLM(:,:,11) / 24.0d0                      !avg
    Tsoil(:,:) = Tsoil(:,:) + CLM(:,:,14) / 24.0d0                  !avg
    ! should I multiply these two fluxes by one hour (3600 s)?
    trans_veg(:,:) = trans_veg(:,:) + CLM(:,:,9)
    LH(:,:) = LH(:,:) + CLM(:,:,1) / 24.0d0
    SH(:,:) = SH(:,:) + CLM(:,:,3) / 24.0d0
    Tgrnd(:,:) = Tgrnd(:,:) + CLM(:,:,12) / 24.0d0
    evap_grnd(:,:) = evap_grnd(:,:) + CLM(:,:,6)
    LWout(:,:) = LWout(:,:) + CLM(:,:,2) / 24.0d0
    evap_tot(:,:) = evap_tot(:,:) + CLM(:,:,5)
    evap_soi(:,:) = evap_soi(:,:) + CLM(:,:,7)
    netrad(:,:) = netrad(:,:) + CLM(:,:,4) / 24.0d0
    evap_veg(:,:) = evap_veg(:,:) + CLM(:,:,8)
    Irr(:,:) = Irr(:,:) + CLM(:,:,12)

    print*,'Fluxes Assigned'

    counter = counter + 1

    end do !hr

    !  write daily averages
    pname = 'SWE.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) SWE
    close(10)

    pname = 'Tsoil.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Tsoil
    close(10)

    pname = 'trans_veg.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) trans_veg
    close(10)

    pname = 'LH.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) LH
    close(10)

    pname = 'SH.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) SH
    close(10)

    pname = 'Tgrnd.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Tgrnd
    close(10)

    pname = 'evap_grnd.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_grnd
    close(10)

    pname = 'LWout.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) LWout
    close(10)

    pname = 'evap_tot.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_tot
    close(10)

    pname = 'evap_soi.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_soi
    close(10)

    pname = 'netrad.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) netrad
    close(10)

    pname = 'evap_veg.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_veg
    close(10)

    pname = 'Irr.daily'
    write(filenum,'(i3.3)') dayofyear
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Irr
    close(10)
 
    print*,'daily averages written'
 
    dayofyear = dayofyear + 1
 
    ! compute monthly averages
    SWEm = SWEm + SWE / float(days(month))
    Tsoil = Tsoilm + Tsoil / float(days(month))
    trans_vegm = trans_vegm + trans_veg
    LHm = LHm + LH / float(days(month))
    SHm = SHm + SH / float(days(month))
    Tgrndm = Tgrndm + Tgrnd / float(days(month))
    evap_grndm = evap_grndm + evap_grnd
    LWoutm = LWoutm + LWout / float(days(month))
    evap_totm = evap_totm + evap_tot
    evap_soim = evap_soim + evap_soi
    netradm = netradm + netrad / float(days(month))
    evap_vegm = evap_vegm + evap_veg
    Irrm = Irrm + Irr
 
    CALL CPU_TIME(t3)
    print*, counter-1, dayofyear-1 , day, monthname(month), t3-t1

    end do !day

    ! write monthly averages
    pname = 'SWE.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) SWEm
    close(10)

    pname = 'Tsoil.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Tsoilm
    close(10)

    pname = 'trans_veg.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) trans_vegm
    close(10)

    pname = 'LH.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) LHm
    close(10)

    pname = 'SH.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) SHm
    close(10)


    pname = 'Tgrnd.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Tgrndm
    close(10)

    pname = 'evap_grnd.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_grndm
    close(10)

    pname = 'LWout.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) LWoutm
    close(10)

    pname = 'evap_tot.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_totm
    close(10)

    pname = 'evap_soi.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_soim
    close(10)

    pname = 'netrad.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) netradm
    close(10)

    pname = 'evap_veg.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_vegm
    close(10)

    pname = 'Irr.monthly'
    write(filenum,'(i2.2)') month
    fname=trim(adjustl(dname))//'/'//trim(adjustl(pname))//'.'//trim(adjustl(filenum))//'.bin'
    open (10,file=trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Irrm
    close(10)

    print*,'monthly averages written'

    ! compute yearly averages
    SWEy = SWEy + SWEm*(float(days(month)) / 365.d0)
    Tsoily = Tsoily +Tsoilm*(float(days(month)) / 365.d0)
    trans_vegy = trans_vegy + trans_vegm
    LHy = LHy + LHm*(float(days(month)) / 365.d0)
    SHy = SHy + SHm *(float(days(month))/ 365.d0)
    Tgrndy = Tgrndy + Tgrndm*(float(days(month)) / 365.d0)
    evap_grndy = evap_grndy + evap_grndm
    LWouty = LWouty + LWoutm*(float(days(month))/365.d0)
    evap_toty = evap_toty + evap_totm
    evap_soiy = evap_soiy + evap_totm
    netrady = netrady + netradm*(float(days(month))/365.d0)
    evap_vegy = evap_vegy + evap_vegm
    Irry = Irry + Irrm

    end do !mo

    print*,'yearly averages computed'

    ! write annual avg
    fname = 'SWE.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) SWEy
    close(10)

    fname = 'trans_veg.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) trans_vegy
    close(10)

    fname = 'LH.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) LHy
    close(10)

    fname = 'SH.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) SHy
    close(10)

    fname = 'Tgrnd.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Tgrndy
    close(10)

    fname = 'evap_grnd.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_grndy
    close(10)

    fname = 'LWout.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) LWouty
    close(10)

    fname = 'evap_tot.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_toty
    close(10)

    fname = 'evap_soi.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_soiy
    close(10)

    fname = 'netrad.yearly.bin'
    open(10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) netrady
    close(10)

    fname = 'evap_veg.yearly.bin'
    open(10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) evap_vegy
    close(10)

    fname = 'Irr.yearly.bin'
    open (10,file=trim(adjustl(dname))//'/'//trim(adjustl(fname)),form='unformatted',access='stream')
    write(10) nx, ny, 1
    write(10) Irry
    close(10)

end subroutine

