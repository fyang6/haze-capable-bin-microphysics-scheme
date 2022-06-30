subroutine write_fields3D_micro

  use vars
  use microphysics

  implicit none

  character *80 filename, long_name
  character *8 name
  character *10 timechar
  character *4 rankchar
  character *6 filetype
  character *10 units

  integer :: i,j,k,nfields, nfields1
  real(4) tmp(nx,ny,nzm)

  character *3 binumber
  integer m
! fields 
  nfields = 8 + naerosol + 2*ncloud
  nfields1 = 0

  if (masterproc) then

    write(rankchar,'(i4)') nsubdomains
    write(timechar,'(i10)') nstep
    do k=1,11-lenstr(timechar)-1
      timechar(k:k) = '0'
    end do

    if (RUN3D) then
      if (save3Dbin) then
        filetype = '.bin3D'
      else
        filetype = '.com3D'
      end if
      filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_micro_'// &
            rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype
      open(46,file=filename,status='unknown',form='unformatted')
    else
      print*,"Not 3D simulation [ToDo]"
      stop
    end if

    if (save3Dbin) then
      write(46) nx, ny, nzm, naerosol, ncloud, &
                nsubdomains, nsubdomains_x, nsubdomains_y, nfields
      do k=1, nzm
        write(46) z(k)
      end do
      do k=1,nzm
        write(46) pres(k)
      end do
      do k=1,nzm
        write(46) rho(k)
      end do
      write(46) dx
      write(46) dy
      write(46) nstep*dt/(3600.*24.) + day0
      do k=1,naerosol
        write(46) raerosol(k)
      end do
      do k=1,ncloud
        write(46) rcloud(k)
      end do
    else
      print*,"Not 3Dbin in micro [ToDo]"
      stop
    end if

  end if

  nfields1 = nfields1 + 1
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        tmp(i,j,k) = qt(i,j,k)*1000.*rho(k)
      end do
    end do
  end do
  name='QT'
  long_name = 'Total Water Vapor Conent'
  units = 'g/m^3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                  save3Dbin,dompi,rank,nsubdomains)

  nfields1 = nfields1 + 1
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        tmp(i,j,k) = ssatw(i,j,k)
      end do
    end do
  end do
  name='SSW'
  long_name = 'Saturation Ratio (water)'
  units = '%'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                  save3Dbin,dompi,rank,nsubdomains)

  nfields1 = nfields1 + 1
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        !tmp(i,j,k) = qcl(i,j,k)*1000.*rho(k)
        tmp(i,j,k) = qcl(i,j,k)*1000.
      end do
    end do
  end do
  name='QC'
  long_name = 'Cloud Water Content'
  units = 'g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                  save3Dbin,dompi,rank,nsubdomains)

  nfields1 = nfields1 + 1
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        !tmp(i,j,k) = qpl(i,j,k)*1000.*rho(k)
        tmp(i,j,k) = qpl(i,j,k)*1000.
      end do
    end do
  end do
  print*,rank,tmp(1,1,17)
  name='QR'
  long_name = 'Rain Water Content'
  units = 'g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                  save3Dbin,dompi,rank,nsubdomains)

  nfields1 = nfields1 + 1
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        tmp(i,j,k) = qna(i,j,k)
      end do
    end do
  end do
  name='NA'
  long_name = 'Aerosol Number Concentration'
  units = '#/mg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                  save3Dbin,dompi,rank,nsubdomains)

nfields1 = nfields1 + 1
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        tmp(i,j,k) = qnh(i,j,k)
      end do
    end do
  end do
  name='NH'
  long_name = 'Haze Number Concentration'
  units = '#/mg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                  save3Dbin,dompi,rank,nsubdomains)

  nfields1 = nfields1 + 1
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        tmp(i,j,k) = qnc(i,j,k)
      end do
    end do
  end do
  name='NC'
  long_name = 'Cloud Droplet Number Concentration'
  units = '#/mg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                  save3Dbin,dompi,rank,nsubdomains)

  nfields1 = nfields1 + 1
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        tmp(i,j,k) = qnr(i,j,k)
      end do
    end do
  end do
  name='NR'
  long_name = 'Rain Drop Number Concentration'
  units = '#/mg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                  save3Dbin,dompi,rank,nsubdomains)

  do m = 1,naerosol
    write(binumber,'(i3)') m
    do k=1,3-lenstr(binumber)
      binumber(k:k) = '0'
    end do

    nfields1 = nfields1 + 1
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,k) = Naer(i,j,k,m)*1.0e-6
        end do
      end do
    end do
    name='FNA'//binumber(1:3)
    long_name = 'Aerosol Number Concentration in Bin'//binumber(1:3)
    units = '#/mg/bin'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                    save3Dbin,dompi,rank,nsubdomains)
  end do

  do m = 1,ncloud
    write(binumber,'(i3)') m
    do k=1,3-lenstr(binumber)
      binumber(k:k) = '0'
    end do

    nfields1 = nfields1 + 1
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,k) = Ncld(i,j,k,m)*1.0e-6
        end do
      end do
    end do
    name='FNC'//binumber(1:3)
    long_name = 'Cloud Droplet Number Concentration in Bin'//binumber(1:3)
    units = '#/mg/bin'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                    save3Dbin,dompi,rank,nsubdomains)
  end do

  do m = 1,ncloud
    write(binumber,'(i3)') m
    do k=1,3-lenstr(binumber)
      binumber(k:k) = '0'
    end do

    nfields1 = nfields1 + 1
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,k) = Mcld(i,j,k,m)
        end do
      end do
    end do
    name='FMC'//binumber(1:3)
    long_name = 'Cloud Mass Concentration in Bin'//binumber(1:3)
    units = 'g/kg/bin'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                    save3Dbin,dompi,rank,nsubdomains)
  end do

  call task_barrier()

  if (nfields.ne.nfields1) then
    if(masterproc) print*,"write_fields3D_micro error: nfields"
    call task_abort()
  end if
  if(masterproc) then
    close(46)
    if(RUN3D.or.save3Dsep) then
      if(dogzip3D) call systemf('gzip -f'//filename)
      print*, 'Writting 3D data. file:'//filename
    else
      print*, 'Appending 3D data. file:'//filename
    end if
  end if
end subroutine write_fields3D_micro
