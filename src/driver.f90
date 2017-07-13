!>
!> driver for ADT based distance calculation
!>
program driver
  use DistCalc, only : initDistanceCalc,MinDist
  use cartwrite
  implicit none
  character*40 :: facetfile,dreqstr
  real*8 :: dreq
  integer :: nnode,ntri,nquad
  integer :: i,j,k,l,n,cindx
  integer :: jmax,kmax,lmax
  !
  real*8 :: xmin(3),xmax(3),dx(3),x1(3),v1(3)
  real*8, allocatable :: xsurf(:,:),x(:,:,:,:),dist(:,:,:)
  integer, allocatable :: tri(:,:)
  integer, allocatable :: quad(:,:)
  real*8 :: CP(3),lenx,leny,lenz
  !
  if (iargc() < 2) then
     write(6,*) "Usage: adtdist <facetfile> <dreq>"
     stop
  endif
  !
  call getarg(1,facetfile)
  call getarg(2,dreqstr)
  read(dreqstr,*) dreq
  !
  open(unit=1,file=facetfile,form='formatted')
  read(1,*) nnode
  allocate(xsurf(3,nnode))
  read(1,*) xsurf
  read(1,*) nquad
  if (nquad >0) then
    allocate(quad(4,nquad))
    read(1,*) quad
  endif
  read(1,*) ntri
  if (ntri > 0) then
   allocate(tri(3,ntri))
   read(1,*) tri
  endif
  !
  xmin=1E15
  xmax=-1E15
  do i=1,nnode
     do j=1,3
        xmin(j)=min(xmin(j),xsurf(j,i))
        xmax(j)=max(xmax(j),xsurf(j,i))
     enddo
  enddo
  !
  ! notional structured grid to compute distance
  ! on (hard code size now
  !
  jmax=64
  kmax=64
  lmax=64
  xmin=xmin-dreq
  xmax=xmax+dreq
  dx=(xmax-xmin)/(jmax+1)
  allocate(x(3,jmax,kmax,lmax))
  allocate(dist(jmax,kmax,lmax))
  !
  do l=1,lmax
     do k=1,kmax
        do j=1,jmax
           x(1,j,k,l)=(j-1)*dx(1)+xmin(1)
           x(2,j,k,l)=(k-1)*dx(2)+xmin(2)
           x(3,j,k,l)=(l-1)*dx(3)+xmin(3)
        enddo
     enddo
  enddo

   lenx=(jmax-1)*dx(1)
   leny=(kmax-1)*dx(2)
   lenz=(lmax-1)*dx(3)

   print *,lenx,leny,lenz,jmax,kmax,lmax,dx
  !
  call initDistanceCalc(xsurf,tri,quad,nnode,ntri,nquad)
  !
  ! compute the distances for each point one by one now
  !
     do l=1,lmax
        write(6,*) dble(l)/lmax*100,'% complete'
        do k=1,kmax
           do j=1,jmax           
              call MinDist(x(:,j,k,l),dist(j,k,l),CP)
           enddo
        enddo
     end do
     !
     open(unit=1,file='x.p3d',form='unformatted')
     write(1) jmax,kmax,lmax
     write(1) ((((x(n,j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),n=1,3)
     close(1)
     !
     open(unit=1,file='q.p3d',form='unformatted')
     write(1) jmax,kmax,lmax
     write(1) 0.5d0,0d0,100d0,0d0
     write(1) ((((dist(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),n=1,5)
     close(1)
   !


   call Writevtrfile(dist,1,"distance",&
    jmax,kmax,lmax,dx(1),dx(2),dx(3),xmin(1),xmin(2),xmin(3))

   deallocate(xsurf,x,dist)
   if (allocated(tri)) deallocate(tri)
   if (allocated(quad)) deallocate(quad)
   !
 end program driver
