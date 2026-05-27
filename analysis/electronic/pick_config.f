      character*80 ofile, odirctry(20), native
      parameter(namax=1000)
      dimension rba(3,3)
      dimension nhk(namax)
      dimension is(namax)
      real x(3,namax), v(3,namax)
      character dummy


c------------------------------------------------
      noffil = 1
      odirctry(1) = 'data1'
c      odirctry(1) = '.'

      write(*,*) 'Step # ?'
      read(*,*) npick
      
c------------------------------------------------
c    get supercell edges
      call celledg( noffil, ofile, odirctry, rba )
c------------------------------------------------
      audang =  0.529177249d0

 
c--- store zero for atomic charge
      nfile = 0
  100 continue
      nfile = nfile + 1
c--- allocate I/O files
      native = '/qm_ion.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 1, err=998, file=ofile, status='old', action='read' )
      write(*,*) 'open : ', ofile

      native = '/md_vel.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 3, err=997, file=ofile, status='old', action='read' )
      write(*,*) 'open : ', ofile
      ifvel = 0
      go to 996
  997 continue
      ifvel = 1
  996 continue


      read(1,'(a1)') dummy
      if( ifvel.eq.0 ) read(3,'(a1)') dummy
c--- read atomic coordinates
   1  read(1,'(100i7)',end=9999) nstep,
     &               ntype, ( nhk(it), it = 1, ntype )
      ntot = 0
      do it = 1, ntype
         ntot = ntot + nhk(it)
      enddo
      read(1,*) atconfx
c      atconfx = 1.0
      read(1,'(9f8.5)')
     &         ( ( x(i,ia), i = 1, 3 ), ia = 1, ntot )
c      if( nstep.lt.nstop ) go to 1
      
      ii = 0
      do it = 1, ntype
      do i = 1, nhk(it)
         ii = ii + 1
         is(ii) = it
      do ix = 1, 3
         x(ix,ii) = x(ix,ii)*atconfx
      enddo
      enddo
      enddo

c--- read atomic velocities
      if( ifvel.eq.0 ) then
      read(3,'(11i7)',end=9999) nstep3, ntot
      read(3,*) atconfx
      read(3,'(9f8.5)')
     &         ( ( v(i,ia), i = 1, 3 ), ia = 1, ntot )
      else
        atconfx = 0.0
      endif
      if( nstep.lt.npick ) go to 1
c      write(*,*) nstep

      ii = 0
      do it = 1, ntype
      do i = 1, nhk(it)
         ii = ii + 1
c         is(ii) = it
         sdr1 = v(1,ii)*atconfx
         sdr2 = v(2,ii)*atconfx
         sdr3 = v(3,ii)*atconfx
         SX  = rba(1,1)*SDR1+rba(1,2)*SDR2+rba(1,3)*SDR3
         SY  = rba(2,1)*SDR1+rba(2,2)*SDR2+rba(2,3)*SDR3
         SZ  = rba(3,1)*SDR1+rba(3,2)*SDR2+rba(3,3)*SDR3
         v(1,ii) = sx/audang
         v(2,ii) = sy/audang
         v(3,ii) = sz/audang
      enddo
      enddo

c--- create image of atoms in x & y & z directions
c    Note: Atoms resides in the range of
c             ext_x1 < x < ext_x2
c             ext_y1 < y < ext_y2
c             ext_z1 < z < ext_z2
c          are treated as the real atom.
c          Atoms outside this range are image atoms.

      write(*,*) nstep
c--- coordinates
      open(10,file='Config.dat',status='unknown')
      write(10,'(i7)') ntot
      do i = 1, ntot
         write(10,'(i2,3f9.6)') is(i),(x(j,i),j=1,3) 
      end do
c--- coordinates
      if( ifvel.eq.0 ) then
      open(11,file='Veloc.dat',status='unknown')
      write(11,'(i7)') ntot
      vmax = 0.d0
      do i = 1, ntot
         vmax = max( vmax, abs(v(1,i)), abs(v(2,i)), abs(v(3,i)) )
      end do
      write(11,*) vmax
      do i = 1, ntot
         write(11,'(i2,3f10.6)') is(i),(v(j,i)/vmax,j=1,3) 
      end do
      end if
      stop

 9999 continue
      close(1)
c      close(2)
      close(3)
      if( nfile.lt.noffil ) go to 100
      stop
  998 continue
      write(*,*) 'cannot find : ', ofile
      end




      subroutine celledg( noffil, ofile, odirctry, rba )
c------------------------------------------------
c    get supercell edges
c------------------------------------------------
c      implicit real*8(a-h, o-z)
      parameter( audang =  0.529177249d0 )
      character*80 ofile, odirctry(20), native
      character*1 dummy
      dimension rba(3,3)

      count = 0.d0
      ava = 0.d0
      avb = 0.d0
      avc = 0.d0
      avbc = 0.d0
      avca = 0.d0
      avab = 0.d0
      nfile = 1
   10 continue
      native = '/md_box.d'
      call getfname( odirctry(nfile), ofile, native )
      write(*,*) 'open : ', ofile
      open( 1, err=998, file=ofile, status='old' )
 1000 read(1,'(a1)') dummy
      if( dummy.eq.'#' ) go to 1000
      backspace 1
c----------------------------------------------------------
    1 read(1,*,end=999) ncstp,
     &           dltca, dltcb, dltcc, angalf, angbet, anggam
      count = count + 1.d0
      ava = ava + dltca
      avb = avb + dltcb
      avc = avc + dltcc
      avbc = avbc + angalf
      avca = avca + angbet
      avab = avab + anggam

 999  continue
      close(1)
      nfile = nfile + 1
      if( nfile.le.noffil ) go to 10
c----------------------------------------------------------
      dltca = ava/count * audang
      dltcb = avb/count * audang
      dltcc = avc/count * audang
      avbc = avbc/count
      avca = avca/count
      avab = avab/count
      angcon = acos(-1.d0)/180.d0
      angalf = avbc * angcon
      angbet = avca * angcon
      anggam = avab * angcon
C                           --- unit vectors parallel to cell vectors ---
                                   E1X = 1.0
                                   E1Y = 0.0
                                   E1Z = 0.0
                                   E2X = COS(anggam)
                                   E2Y = SIN(anggam)
                                   E2Z = 0.0
                                   E3X = COS(angbet)
                                   E3Y = COS(angalf) - E3X*E2X
                                   E3Y = E3Y/E2Y
                                   E3Z = 1.0 - E3X*E3X - E3Y*E3Y
                                   E3Z = SQRT(E3Z)

                                   rba(1,1) = dltca*E1X
                                   rba(2,1) = dltca*E1Y
                                   rba(3,1) = dltca*E1Z
                                   rba(1,2) = dltcb*E2X
                                   rba(2,2) = dltcb*E2Y
                                   rba(3,2) = dltcb*E2Z
                                   rba(1,3) = dltcc*E3X
                                   rba(2,3) = dltcc*E3Y
                                   rba(3,3) = dltcc*E3Z

      write(*,*) 'cell edges in [A]'
      write(*,'(3f10.5)') rba(1,1), rba(2,1), rba(3,1)
      write(*,'(3f10.5)') rba(1,2), rba(2,2), rba(3,2)
      write(*,'(3f10.5)') rba(1,3), rba(2,3), rba(3,3)

      return
 998  continue
      write(*,*) 'missing fort.47'
      stop
      end




      subroutine getfname( odirctry, ofile, native )
      character*80 ofile, odirctry, native

      do i = 1, 80
         if( odirctry(i:i).ne.' ' ) then
             io1 = i
             go to 1
         endif
      enddo
    1 continue
      do i = 80, 1, -1
         if( odirctry(i:i).ne.' ' ) then
             io2 = i
             go to 2
         endif
      enddo
    2 continue
      do i = 1, 80
         if( native(i:i).ne.' ' ) then
             in1 = i
             go to 3
         endif
      enddo
    3 continue
      do i = 80, 1, -1
         if( native(i:i).ne.' ' ) then
             in2 = i
             go to 4
         endif
      enddo
    4 continue
      ofile = ( odirctry(io1:io2)//native(in1:in2) )
c      write(*,*) '1:', ofile
c      write(*,*) '2:', odirctry(io1:io2)
c      write(*,*) '3:', native(in1:in2)

      return
      end




      BLOCK DATA SETATM
c-----------------------------------------------------------------------
c     dmassn : mass number
c-----------------------------------------------------------------------
      implicit real*8(a-h, o-z)
      common/dmassn/ dmassn(103)

      data ( dmassn(i), i = 1, 56 ) /
     &    1.0079400d0,   4.0026020d0,   6.9410000d0,   9.0121820d0,
     &   10.8110000d0,  12.0100000d0,  14.0067400d0,  15.9994000d0,
     &   18.9984032d0,  20.1797000d0,  22.9897680d0,  24.3050000d0,
     &   26.9815390d0,  28.0855000d0,  30.9737620d0,  32.0660000d0,
     &   35.4527000d0,  39.9480000d0,  39.0983000d0,  40.0780000d0,
     &   44.9559100d0,  47.8800000d0,  50.9415000d0,  51.9961000d0,
     &   54.9380500d0,  55.8470000d0,  58.9332000d0,  58.6900000d0,
     &   63.5460000d0,  65.3900000d0,  69.7230000d0,  72.6100000d0,
     &   74.9215900d0,  78.9600000d0,  79.9040000d0,  83.8000000d0,
     &   85.4678000d0,  87.6200000d0,  88.9058500d0,  91.2240000d0,
     &   92.9063800d0,  95.9400000d0,  98.9100000d0, 101.0700000d0,
     &  102.9055000d0, 106.4200000d0, 107.8682000d0, 112.4110000d0,
     &  114.8200000d0, 118.7100000d0, 121.7500000d0, 127.6000000d0,
     &  126.9044700d0, 131.2900000d0, 132.9054300d0, 137.3270000d0 /
      data ( dmassn(i), i = 57, 103 ) /
     &  138.9055000d0, 140.1150000d0, 140.9076500d0, 144.2400000d0,
     &  145.0000000d0, 150.0600000d0, 151.9650000d0, 157.2500000d0,
     &  158.9253400d0, 162.5000000d0, 164.9303200d0, 167.2600000d0,
     &  168.9342100d0, 173.0400000d0, 174.9670000d0, 178.4900000d0,
     &  180.9479000d0, 183.8500000d0, 186.2070000d0, 190.2000000d0,
     &  192.2200000d0, 195.0800000d0, 196.9665400d0, 200.5900000d0,
     &  204.3833000d0, 207.2000000d0, 208.9803700d0, 210.0000000d0,
     &  210.0000000d0, 222.0000000d0, 223.0000000d0, 226.0000000d0,
     &  227.0000000d0, 232.0381000d0, 231.0358800d0, 238.0289000d0,
     &  237.0500000d0, 244.0000000d0, 243.0000000d0, 247.0000000d0,
     &  247.0000000d0, 251.0000000d0, 254.0000000d0, 257.0000000d0,
     &  256.0000000d0, 254.0000000d0, 257.0000000d0 /

      END

