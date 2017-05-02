       program gofr_tutorial
       implicit none
       integer, parameter :: dim=3,dbl=8
       real(kind=dbl), parameter :: au_convert=.529177_dbl
       integer, dimension(:,:,:), allocatable :: num
       real(kind=dbl), dimension(:,:,:), allocatable :: gofr
       real(kind=dbl), dimension(:,:), allocatable :: r,s
       real(kind=dbl), dimension(dim) :: sij,rijvec,si,temp
       real(kind=dbl), dimension(dim,dim) :: box,boxinv
       integer :: nframe,natoms,ntypes,&
                  i,j,iframe,iatom,jatom,idim,maxbin,ibin,&
                  itype,jtype,jnum,inum,jdim
       integer, dimension(:), allocatable :: atomtypes,atomcount
       real(kind=dbl) :: delta,volume,rijsq,maxdist,rcutsq,rij,rupper,&
                   rmin, rlower,shellvolume,pi,rho,idealnum,dist,factor
       character(len=50) :: filename
       character(len=1), dimension(:), allocatable :: name
!
! open gofr input file that gives molecular information
!
       call getarg(1,filename)
       open(unit=1,file=filename)
       read(1,*) ntypes
!      write(6,*) ntypes
       allocate(atomtypes(ntypes))
       allocate(name(ntypes))
       allocate(atomcount(ntypes))
!      write(6,*) ntypes
!
! now open the coordinate file
!
       call getarg(2,filename)
       open(unit=2,file=filename)
       read(2,*) nframe
       write(6,*) nframe
       read(2,*) natoms
       write(6,*) natoms
       atomcount=0
       call getarg(3,filename)
       open(unit=3,file=filename)
 !
 ! read in each atom, its type and name
 ! at the same time, calculate the total number of
 ! atoms that are of the same type
 !
       do i=1,ntypes
       	read(1,*) atomtypes(i), atomcount(i), name(i)
        write(6,*) atomtypes(i), atomcount(i), name(i)
       enddo
       allocate(r(dim,natoms))
       allocate(s(dim,natoms))
       read(1,*) delta
       read(1,*) maxdist
       maxbin=int(maxdist/delta) + 1
       allocate(num(ntypes,ntypes,maxbin))
       allocate(gofr(ntypes,ntypes,maxbin))
       gofr=0._dbl
!
! danger, all of the results depend on the following line !!!
!
       pi=3.14159_dbl
!
! 
! read box matrix, assumed to have all cell angles 90
!
! 
! start loop over number of frames
!
ifl:   do iframe=1,nframe
       if(mod(iframe,10).eq.0) write(6,*) iframe
!
! read all the atomic positions 
!
        r=0._dbl
       box=0._dbl
       read(2,*)
       read(3,*)
       do idim=1,3
         read(3,*)  (box(idim,jdim),jdim=1,dim)
       enddo
	do iatom=1,natoms
	 read(2,*) r(:,iatom)
	enddo
!
! all group codes units are in atomic units, so convert
!
! NO CONVERSION for PW-RELAX
        r=r*au_convert
 	box=box*au_convert
!
! assume cell angles are 90 and get the square of half of the smallest
!
        do idim=1,dim
	 temp(idim)=box(idim,idim)
	enddo
	rcutsq=minval(temp)
	rcutsq=rcutsq*.5_dbl
	rcutsq=rcutsq*rcutsq
        rmin=0.66_dbl
!
! now scale all coordinates
!
	call gethinv(box,boxinv,volume)
	do iatom=1,natoms
	 s(:,iatom)=matmul(boxinv,r(:,iatom))
	enddo
!
! now loop over all molecules (intermolecular gofr)
!
        num=0
iatoml:	 do iatom=1,natoms
!
! iatom is the index of the current atom
!
           if(iatom.le.atomcount(1)) then
              i=1
           else 
              i=2
           endif
!         write(6,*) iatom,i
	  si=s(:,iatom)
jatoml:	   do jatom=1,natoms
            if(jatom.le.atomcount(1)) then
              j=1
            else 
              j=2
            endif
!
! jatom is the index of the j'th atom
!
!
! check to see the iatom is not the same as jatom
!
diffl:      if(jatom/=iatom) then
	     sij=s(:,jatom) - si
!
! minimum image 
!
	     sij=sij-dnint(sij)
	     rijvec=matmul(box,sij)
	     rijsq=dot_product(rijvec,rijvec)
cutmin:	     if(rijsq>rmin) then
cutl:	     if(rijsq<rcutsq) then
	      rij=sqrt(rijsq)
              if(atomtypes(i).eq.atomtypes(j).and.rij.lt.1.1) then
                write(6,*) r(:,iatom)
                write(6,*) r(:,jatom)
                write(6,*) iframe,iatom,jatom,rij
                stop
              endif
	      ibin=int(rij/delta+.5_dbl) + 1
	      if(ibin<=maxbin) &
	       num(atomtypes(i),atomtypes(j),ibin)=  &
	            num(atomtypes(i),atomtypes(j),ibin) + 1
             endif cutl
             endif cutmin
	    endif diffl
	   enddo jatoml
	 enddo iatoml
!
! now that this frame is over, get the running gofr
!
        factor=(4._dbl/3._dbl)*pi
itypel: do itype=1,ntypes
!
! get the number of atoms of type itype
!
	 inum=atomcount(itype)
jtypel:	 do jtype=1,ntypes
!
! get the number of atoms of type jtype
! if the types are the same, then the density should be
! calculated with one less atom
!
          if(itype==jtype) then
           jnum=atomcount(jtype) -1
	  else
	   jnum=atomcount(jtype)
	  endif
	  rho=real(jnum)/volume
ibinl:	  do ibin=1,maxbin
	   dist=delta*(real(ibin-1)) 
	   rupper=(dist + delta)**3
	   rlower=dist**3
	   shellvolume=factor*(rupper-rlower)
	   idealnum=shellvolume*rho
!
! divide the total number of atoms num by the ideal gas
! number. Then divide by inum since we want the average
! per particle of type itype
!
	   gofr(itype,jtype,ibin)=gofr(itype,jtype,ibin) + &
	    num(itype,jtype,ibin)/(idealnum * inum)
	  enddo ibinl
	 enddo jtypel
	enddo itypel
       enddo ifl
!
! frame loop is over. Now, divide by the number of frames
! to get the average.
!
       do itype=1,ntypes
        do jtype=itype,ntypes
!         write(6,*) itype,jtype,name(itype),name(jtype)
!
! construct the output filename based on the atom names
!
	 filename=name(itype)//'-'//name(jtype)//'.gofr'
	 open(unit=10,file=filename)
	 do ibin=1,maxbin
	  dist=delta*(real(ibin-1)) + .5_dbl*delta
	  write(10,*) dist,gofr(itype,jtype,ibin)/real(nframe)
	 enddo
	 close(10)
	enddo
       enddo
       end

!
! routine to get the inverse and determinate of a 3x3 matrix
!
      subroutine gethinv(hmat,hmati,deth)
      implicit none
      integer, parameter :: dim=3,dbl=8
      real(kind=dbl), dimension(dim,dim) :: hmat,hmati
      real(kind=dbl) :: deth
      deth = HMAT(1,1)*(HMAT(2,2)*HMAT(3,3)-HMAT(2,3)*HMAT(3,2)) &
           + HMAT(1,2)*(HMAT(2,3)*HMAT(3,1)-HMAT(2,1)*HMAT(3,3)) &
           + HMAT(1,3)*(HMAT(2,1)*HMAT(3,2)-HMAT(2,2)*HMAT(3,1))
       HMATI(1,1) = (HMAT(2,2)*HMAT(3,3)-HMAT(2,3)*HMAT(3,2))/deth
       HMATI(2,2) = (HMAT(1,1)*HMAT(3,3)-HMAT(1,3)*HMAT(3,1))/deth
       HMATI(3,3) = (HMAT(1,1)*HMAT(2,2)-HMAT(1,2)*HMAT(2,1))/deth
       HMATI(1,2) = (HMAT(1,3)*HMAT(3,2)-HMAT(1,2)*HMAT(3,3))/deth
       HMATI(2,1) = (HMAT(3,1)*HMAT(2,3)-HMAT(2,1)*HMAT(3,3))/deth
       HMATI(1,3) = (HMAT(1,2)*HMAT(2,3)-HMAT(1,3)*HMAT(2,2))/deth
       HMATI(3,1) = (HMAT(2,1)*HMAT(3,2)-HMAT(3,1)*HMAT(2,2))/deth
       HMATI(2,3) = (HMAT(1,3)*HMAT(2,1)-HMAT(2,3)*HMAT(1,1))/deth
       HMATI(3,2) = (HMAT(3,1)*HMAT(1,2)-HMAT(3,2)*HMAT(1,1))/deth
       RETURN
       END
