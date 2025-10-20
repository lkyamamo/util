
!---------------------------------------------------------------
!                                                     2017/11/28
!   
!   This program is based on dosg-region2.f90.
!
!   You can choose whether or not DOS is multipled 
!   by the ratio of atoms 
!   => lmulti = .true. : DOS is multipled by the ratio of atoms
!   => lmulti = .false. : DOS is NOT multipled by the ratio of atoms
!
!   Note:
!
!                                                 Yuichi Uchida
!---------------------------------------------------------------
!---------------------------------------------------------------
!     Electronic density of states
!---------------------------------------------------------------
      parameter( mxtype = 30 )
      parameter( ntmax = 2000 )
      parameter( lmax = 3 )
      parameter( nbandx = 2000 )
      parameter( mshmx = 1000 )
      parameter( mxregion = 20 )
      implicit real*8(a-h, o-z)
      dimension dos(-mshmx:mshmx)
      dimension pdos(-mshmx:mshmx,0:mxtype,mxregion)
      dimension pdoshyb(-mshmx:mshmx,0:mxtype,0:mxtype,mxregion)
      dimension pdosl(-mshmx:mshmx,0:lmax,0:mxtype,mxregion)
      dimension weight(nbandx,0:mxtype,mxregion)
      dimension weightl(nbandx,0:lmax,0:mxtype,mxregion)
      integer :: weightla(nbandx,10,0:ntmax)
      integer :: nofatom(mxregion)
      integer :: id_atom(ntmax,mxregion)
      integer :: is(0:ntmax)
      dimension natom(0:mxtype)
      dimension natom_region(mxtype,mxregion)
      dimension nhk1(mxtype)
      dimension nhk2(mxtype)
      dimension egvl(nbandx)
      dimension nmulao(mxtype), lmulao(10,mxtype)
      dimension wgtex(10,mxtype)
      dimension gaussian(-mshmx:mshmx)
      character*80 ofile, odirctry(20), native, ofbase
      character dummy
      character num(0:9)
      logical  lcpdos, lgaussfilter
      data num / '0','1','2','3','4','5','6','7','8','9' /
!--------------------------------------------------------
      parameter( mxatom = 200 )
      integer :: isort(1:mxatom)
      logical :: lmulti
      
!--------------------------------------------------------

!------------------------------------------------
!--- No. of files to open
      nopen = 1
!      odirctry(1) = 'data9'
      odirctry(1) = '../data/'

      lmulti = .true.

      nini  = 0
      nskip = 1

      dele = 0.005d0
!      ----- in [Ryd.]

      deleh = dele * 0.5

      lgaussfilter = .true.
      wgauss = 0.01d0
!      ----- in [Ryd.]  Gaussian broadening
!------------------------------------------------
      ifop = 1
      nfile = ifop
      native = '/qm_ion.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 14, err=998, file=ofile, status='old' )
      write(*,*) 'open : ', ofile
      read(14,'(a1)') dummy
      read(14,'(90i7)') ncstp, ntype, (natom(it), it = 1, ntype)
      close(14)
      ntot = 0
      do it = 1, ntype
         ntot = ntot + natom(it)
      enddo
       

      nregions = 4    !---No. of regions

      nofatom(1) = 2   !---No. of atoms in region1
      nofatom(2) = 2  !---No. of atoms in region2
      nofatom(3) = 2  !---No. of atoms in region3

      nleft = 0
      do i = 1, nregions-1
        nleft = nleft + nofatom(i) 
      end do
      
      nofatom(nregions) = ntot - nleft !---No. of atoms in last region
    
!-----specify id of atoms (start)--------!

      !Note: write id of each region in ascending order  

      id_atom(1:nofatom(1),1) = (/97,103/)  !---id of atoms in region1
      id_atom(1:nofatom(2),2) = (/78,79/)  !---id of atoms in region2
      id_atom(1:nofatom(3),3) = (/27,28/)  !---id of atoms in region3

      ! gather all specified atoms
      nid = 0
      do nr = 1, nregions-1
         do id = 1, nofatom(nr)
            nid = nid + 1
            isort(nid) = id_atom(id,nr)
         end do
      end do
      ! sort in ascending order
      do i = 1, nid-1
         do j = i, nid
            if(isort(i) > isort(j)) then
               k = isort(i)
               isort(i) = isort(j)
               isort(j) = k
            end if
         end do
      end do
      ! make id_atom of last region (the others)
      ii = 0
      j = 1
      do i = 1, ntot
         if( i /= isort(j) ) then
            ii = ii + 1
            id_atom(ii,nregions) = i
         else
            j = j + 1
         end if
      end do

      ! check specified atoms
      do ir = 1, nregions
         do id = 1, nofatom(ir)
         write( unit=6, fmt="(a,i3,a,i5,a,i5)", iostat=ios ) &
          &  " ir=",ir,"  id=",id, "  id_atom(id,ir)=", id_atom(id,ir)
         end do
      end do

!-----specify id of atoms (end)--------!          

!------------------------------------------------


      !---set Gaussian 
      if( lgaussfilter ) then
      ngauss = wgauss/dele * 4
      w = 0.d0
      do i = -ngauss, ngauss
         e = dele * dble(i)
         egs = e/wgauss
         gus = exp(-egs**2)
         gaussian(i) = gus
         w = w + gus
      end do 
      do i = -ngauss, ngauss
         gaussian(i) = gaussian(i)/w
      end do
      end if



      do i = -mshmx, mshmx
         dos(i) = 0.0
      enddo
      do ir = 1, mxregion
      do it = 0, mxtype
      do i = -mshmx, mshmx
         pdos(i,it,ir) = 0.0
      enddo
      enddo
      do it = 0, mxtype
      do l =  0, lmax
      do i = -mshmx, mshmx
         pdosl(i,l,it,ir) = 0.0
      enddo
      enddo
      enddo
      do it = 0, mxtype
      do jt = 0, mxtype
      do i = -mshmx, mshmx
         pdoshyb(i,jt,it,ir) = 0.0
      enddo
      enddo
      enddo
      enddo

      ii1 = 0
      count = 0.0
      av_efermi=0.0
      pcount = 0.0
      ncstp1 = 0
      ifop = 1
!------------------------------------------------
    2 continue
      nfile = ifop

      native = '/qm_eig.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 1, err=998, file=ofile, status='old' )
      write(*,*) 'open : ', ofile

      native = '/qm_fer.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 2, err=998, file=ofile, status='old' )
      write(*,*) 'open : ', ofile

!      if( lcpdos ) then
!         native = '/qm_pds.d'
!         call getfname( odirctry(nfile), ofile, native )
!         open( 3, err=998, file=ofile, status='old' )
!         write(*,*) 'open : ', ofile

         native = '/qm_pda.d'
         call getfname( odirctry(nfile), ofile, native )
         open( 7, err=998, file=ofile, status='old' )
         write(*,*) 'open : ', ofile
!      end if

      native = '/qm_ion.d'
      call getfname( odirctry(nfile), ofile, native )
      open( 4, err=998, file=ofile, status='old' )
      write(*,*) 'open : ', ofile

      go to 887
  998 continue
      write(*,*) ofile,'does not exist.'
      stop
  887 continue
      read(1,'(a1)') dummy
      read(2,'(a1)') dummy
!      if( lcpdos ) then
!          read(3,'(a1)') dummy
          read(7,'(a1)') dummy
!      end if

!     -----get No. of atoms
      read(4,'(a1)') dummy
      read(4,'(90i7)',end=999) ncstp, ntype, (natom(it), it = 1, ntype)
      ntot = 0
      do it = 1, ntype
         ntot = ntot + natom(it)
      enddo
      close(4)
      write(*,*) 'total No. of atoms:', ntot


    3 continue
!      if( lcpdos ) then
!--- lcpdos-start
!      read(3,*,end=999) ncstp1, ntype, nband, lmx
!      if( lmx < 0 ) then
!          do it = 1, ntype
!             read(3,*) nmulao(it), (lmulao(l,it),l=1,nmulao(it))
!          enddo
!          lmx = 0
!          do it = 1, ntype
!          do l = 1, nmulao(it)
!             lmx = max( lmx, lmulao(l,it) )
!          end do
!          end do
!
!          do ib = 1, nband
!             read(3,*,end=999) idummy, 
!     &           ((wgtex(l,it), l=1,nmulao(it)),it=1,ntype)
!             do it = 1, ntype
!             do l = 0, lmax
!                weightl(ib,l,it) = 0.0
!             end do
!             end do
!             do it = 1, ntype
!             do l = 1, nmulao(it)
!                weightl(ib,lmulao(l,it),it)
!     &        = weightl(ib,lmulao(l,it),it) + wgtex(l,it)
!             end do
!             end do
!          end do
!        else
!          do ib = 1, nband
!             read(3,*,end=999) idummy, 
!     &           ((weightl(ib,l,it), l=0,lmx),it=1,ntype)
!          end do
!      end if

!---
      read(7,'(90i7)',end=999) ncstp1, ntype, (natom(it),it = 1, ntype), &
     & nband
      ntot = 0
      do it = 1, ntype
         ntot = ntot + natom(it)
      enddo
      nhk1(1) = 1
      nhk2(1) = nhk1(1) + natom(1) - 1
      do it = 2, ntype
         nhk1(it) = nhk2(it-1) + 1
         nhk2(it) = nhk1(it) + natom(it) - 1
      end do
      do it = 1, ntype
      do i = nhk1(it), nhk2(it)
         is(i) = it
      end do
      end do

          do it = 1, ntype
             read(7,*) nmulao(it), (lmulao(l,it),l=1,nmulao(it))
          enddo
          lmx = 0
          do it = 1, ntype
          do l = 1, nmulao(it)
             lmx = max( lmx, lmulao(l,it) )
          end do
          end do

      do ib = 1, nband
         do i = 1, ntot
            it = is(i)
            do l = 1, nmulao(it)
               read(7,*) weightla(ib,l,i)
            end do
         end do
      end do


      do ir = 1, nregions

         !---calculate the number of atoms in each region
         do it = 1, ntype
            natom_region(it,ir) = 0
         end do
         do n = 1, nofatom(ir)
            i = id_atom(n,ir)
            it = is(i)
            natom_region(it,ir) = natom_region(it,ir) + 1
         end do

         !---calculate the weights in each region
          do ib = 1, nband
             do it = 1, ntype
             do l = 0, lmax
                weightl(ib,l,it,ir) = 0.0
             end do
             end do
             do n = 1, nofatom(ir)
                i = id_atom(n,ir)
                it = is(i)
             do l = 1, nmulao(it)
                weightl(ib,lmulao(l,it),it,ir) &
     &        = weightl(ib,lmulao(l,it),it,ir) + weightla(ib,l,i)
             end do
             end do
          end do
      end do

      !---normalization for each band
      do ib = 1, nband

         do ir = 1, nregions
         do it = 1, ntype
            weight(ib,it,ir) = 0.0
            do l = 0, lmx
               weight(ib,it,ir) = weight(ib,it,ir) + weightl(ib,l,it,ir)
            end do
         end do
         end do

         sunp = 0.0
         do ir = 1, nregions
         do it = 1, ntype
            sunp = sunp + weight(ib,it,ir)
         enddo
         enddo
         do ir = 1, nregions
         do it = 1, ntype
            weight(ib,it,ir) = weight(ib,it,ir)/sunp
            do l = 0, lmx
               weightl(ib,l,it,ir) = weightl(ib,l,it,ir)/sunp
            end do
         end do
         end do

      end do
!--- lcpdos-end
!      end if

!      iterrec = -1
    1 continue
      read(1,*,end=999) ncstp, idummy, nband
      do ib = 1, nband
         read(1,*) ii, egvl(ib)
      enddo
      read(2,*,end=999) ncstp2, idummy, efermi
!      if( iterrec /= idummy ) then
!          iterrec = idummy
!          go to 1
!      end if
!      write(*,*) iterrec


      if( mod(ncstp,nskip).eq.0 .and. ncstp.ge.nini ) then

          count = count + 1.0
!------------------------------------------!
          av_efermi=av_efermi+efermi
!------------------------------------------!

          do ib=1, nband

             egvlff = egvl(ib)!-efermi
             idos = ( egvlff + deleh )/dele
             if( egvlff + deleh.lt.0.0 ) idos = idos - 1

             if( .not.lgaussfilter ) then
                 dos(idos) = dos(idos) + 1.0
               else
                 do i = -ngauss, ngauss
                    idosi = idos + i
                    if( idosi >= -mshmx .and. idosi <= mshmx ) then
                        dos(idosi) = dos(idosi) + gaussian(i)
                    end if
                 end do
             end if

          enddo

!          if( lcpdos ) then

              pcount = pcount + 1.0
              do ib=1, nband

                 egvlff = egvl(ib)!-efermi
                 idos = ( egvlff + deleh )/dele
                 if( egvlff + deleh.lt.0.0 ) idos = idos - 1

             if( .not.lgaussfilter ) then
                 do ir = 1, nregions
                 do it = 1, ntype
                    pdos(idos,it,ir) = pdos(idos,it,ir) + weight(ib,it,ir)
                    do l = 0, lmx
                       pdosl(idos,l,it,ir) = pdosl(idos,l,it,ir) &
     &                                     + weightl(ib,l,it,ir)
                    end do
                 end do
                 end do
               else
                 do i = -ngauss, ngauss
                    idosi = idos + i
                    if( idosi >= -mshmx .and. idosi <= mshmx ) then
                        do ir = 1, nregions
                        do it = 1, ntype
                           pdos(idosi,it,ir) = pdos(idosi,it,ir) &
     &                                    + weight(ib,it,ir) * gaussian(i)
                           do l = 0, lmx
                              pdosl(idosi,l,it,ir) = pdosl(idosi,l,it,ir) &
     &                                  + weightl(ib,l,it,ir) * gaussian(i)
                           end do
                        end do
                        end do
                    end if
                 end do
             end if

!---           decomposed by hybridization
             if( .not.lgaussfilter ) then
                 do ir = 1, nregions
                 do it = 1, ntype
                 do jt = 1, ntype
                    pdoshyb(idos,jt,it,ir) = pdoshyb(idos,jt,it,ir) &
     &                                  + weight(ib,jt,ir) * weight(ib,it,ir)
                 enddo
                 enddo
                 enddo
               else
                 do i = -ngauss, ngauss
                    idosi = idos + i
                    if( idosi >= -mshmx .and. idosi <= mshmx ) then
                        do ir = 1, nregions
                        do it = 1, ntype
                        do jt = 1, ntype
                           pdoshyb(idosi,jt,it,ir) = pdoshyb(idosi,jt,it,ir) &
     &                   + weight(ib,jt,ir) * weight(ib,it,ir) * gaussian(i)
                        end do
                        end do
                        end do
                    end if
                 end do
             end if

              enddo

          endif

!      endif
      if( ncstp.ge.ncstp1 ) then
          go to 3
        else
          go to 1
      end if
  999 continue
      close(1)
      close(2)
      close(3)
      close(7)
      ifop = ifop + 1
      if( ifop.le.nopen ) then
          go to 2
      endif
!-------------------------------------------------------
      do i = -mshmx, mshmx
         if( dos(i).gt.1e-5) then
             ii1 = i
             go to 30
         endif
      enddo
  30  continue
      if( ii1.gt.-mshmx ) ii1 = ii1 - 1
      do i = mshmx, -mshmx, -1
         if( dos(i).gt.1e-5 ) then
             ii2 = i
             go to 40
         endif
      enddo
  40  continue
      write(*,*) count, pcount
      if( count.lt.0.5 ) stop
!--------------------------------------!
     av_efermi=av_efermi/count*13.6058
!--------------------------------------!

      dele = dele*13.6058
!         [Ryd.] --> [eV]

      !---total DOS
      count = count*real(ntot)*dele*0.5
      open(10, file='DOS.dat')
      do i = ii1, ii2
!         write(10,'(100e15.7)') real(i)*dele, dos(i)/count
         write(10,'(100e15.7)') real(i)*dele-av_efermi, dos(i)/count
      enddo
      endfile(10)
      close(10)

      if( pcount.gt.0.5 ) then
          pcount = pcount*dele*0.5
          do ir = 1, nregions

             !--- total DOS in each region
             if( ir < 10 ) then
                 ofbase = 'DOS'//num(0)//num(ir)//'.dat'
               else
                 ir10 = ir/10
                 ir1  = mod(ir,10)
                 ofbase = 'DOS'//num(ir10)//num(ir1)//'.dat'
             end if

             do i = ii1, ii2
                dos(i) = 0.0
             end do
             do it = 1, ntype
                do i = ii1, ii2
                   dos(i) = dos(i) + pdos(i,it,ir)
                end do
             end do

          open(10, file=ofbase(1:len_trim(ofbase)) )
          do i = ii1, ii2
             if ( lmulti ) then
                write(10,'(100e15.7)') real(i)*dele-av_efermi, &
     &            dos(i)/pcount/real(ntot) !--multiply ratio
             else 
                write(10,'(100e15.7)') real(i)*dele-av_efermi, &
     &            dos(i)/pcount/real(nofatom(ir))
             end if
          enddo
          endfile(10)
          close(10)


          do it = 1, ntype
          if( natom_region(it,ir) > 0 ) then

             if( ir < 10 ) then
                 ofbase = 'DOS'//num(0)//num(ir)//'.dat'
               else
                 ir10 = ir/10
                 ir1  = mod(ir,10)
                 ofbase = 'DOS'//num(ir10)//num(ir1)//'.dat'
             end if
             if( it < 10 ) then
                 ofbase = ofbase(1:len_trim(ofbase))//'.type='//num(it)
               else
                 it10 = it/10
                 it1  = mod(it,10)
                 ofbase = ofbase(1:len_trim(ofbase))//'.type='//num(it10)//num(it1)
             end if

          open(10, file=ofbase(1:len_trim(ofbase)) )
          do i = ii1, ii2
             if ( lmulti ) then
                write(10,'(100e15.7)') real(i)*dele-av_efermi, &
     &            pdos(i,it,ir)/pcount/real(natom(it))!--multiply ratio
             else
                write(10,'(100e15.7)') real(i)*dele-av_efermi, &
     &            pdos(i,it,ir)/pcount/real(natom_region(it,ir))
             end if
          enddo
          endfile(10)
          close(10)

          do l = 0, lmx

          open(10, file=ofbase(1:len_trim(ofbase))//'.l='//num(l) )
          do i = ii1, ii2
             if ( lmulti ) then
                write(10,'(100e15.7)') real(i)*dele-av_efermi, &
     &            pdosl(i,l,it,ir)/pcount/real(natom(it))!-multiply ratio
             else
                write(10,'(100e15.7)') real(i)*dele-av_efermi, &
     &            pdosl(i,l,it,ir)/pcount/real(natom_region(it,ir))
             end if
          enddo
          endfile(10)
          close(10)

          end do

!---           decomposed by hybridization
          do jt = 1, ntype
          if( natom_region(jt,ir) > 0 ) then
          if( jt < 10 ) then
              open(10, file=ofbase(1:len_trim(ofbase))//'.type='//num(jt) )
            else
              jt10 = jt/10
              jt1  = mod(jt,10)
              open(10, file=ofbase(1:len_trim(ofbase)) &
     &                    //'.type='//num(jt10)//num(jt1) )
          end if
          do i = ii1, ii2
             if ( lmulti ) then
                write(10,'(100e15.7)') real(i)*dele-av_efermi, &
     &            pdoshyb(i,jt,it,ir)/pcount/real(natom(it))
             else
                write(10,'(100e15.7)') real(i)*dele-av_efermi, &
     &            pdoshyb(i,jt,it,ir)/pcount/real(natom_region(it,ir))
             end if
          enddo
          endfile(10)
          close(10)
          end if
          enddo

          end if
          end do
          end do
      endif
      write(*,*) "*********************************************"
      if (lmulti) then
      write(*,*) "lmulti = .true."
      write(*,*) "DOS is multipled by the ratio of atoms"
      else
      write(*,*) "lmulti = .false."
      write(*,*) "DOS is NOT multipled by the ratio of atoms"
      end if
      write(*,*) "*********************************************"
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
!      write(*,*) '1:', ofile
!      write(*,*) '2:', odirctry(io1:io2)
!      write(*,*) '3:', native(in1:in2)

      return
      end

