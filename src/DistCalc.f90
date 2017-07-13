module DistCalc
  integer, allocatable :: bface(:,:)
  integer, allocatable :: adtIntegers(:,:)
  real*8, allocatable :: adtExtents(:,:)
  real*8, allocatable :: adtReals(:,:)
  real*8, allocatable :: coord(:,:)
  real*8, allocatable :: bfaceinfo(:,:)
  integer :: nelem
  integer :: initdist=0
  real*8  :: dreqinput

  contains
    subroutine initDistanceCalc(xsurf,tri,quad,nnode,ntri,nquad,dr0)
      implicit none
      integer, intent(in) :: nnode !< number of nodes in tesselation
      integer, intent(in) :: ntri  !< number of triangles
      integer, intent(in) :: nquad !< number of quads
      real*8, intent(in) :: xsurf(3,nnode) !< grid coords
      integer, intent(in) :: tri(3,ntri)   !< triangle connectivity
      integer, intent(in) :: quad(4,nquad) !< quadrilatera connectivity 
      real*8, optional, intent(in) :: dr0  !< optional bound to search in
      integer :: i,j,n
      !
      dreqinput=1E15
      initdist=1
      nelem=ntri+2*nquad
      allocate(bface(3,nelem))
      !
      do i=1,ntri
         bface(:,i)=tri(:,i)
      enddo
      !
      j=ntri
      do i=1,nquad
         j=j+1
         bface(:,j)=[quad(1,i),quad(2,i),quad(3,i)]
         j=j+1
         bface(:,j)=[quad(1,i),quad(3,i),quad(4,i)]
      enddo
      !
      ! initialize the ADT
      !
      allocate(bfaceinfo(22,nelem))
      call tribfdata(xsurf,bface,bfaceinfo,nnode,nelem)
      allocate(coord(6,nelem),adtIntegers(4,nelem),adtReals(6,nelem),adtExtents(6,2))
      do i=1,nelem
         coord(1:3,i)=1E15
         coord(4:6,i)=-1E15
         do j=1,3
            do n=1,3
               coord(n,i)=min(coord(n,i),xsurf(n,bface(j,i)))
               coord(n+3,i)=max(coord(n+3,i),xsurf(n,bface(j,i)))
            enddo
         enddo
      enddo
      write(6,*) 'nelem=',nelem
      call buildADT(coord,adtIntegers,adtReals,adtExtents,nelem,6)
      if (present(dr0)) then
         dreqinput=dr0*dr0*1.1d0*1.1d0
      end if
    end subroutine initDistanceCalc
    !
    subroutine MinDist(xp,dd,CP)
      real*8, intent(in) :: xp(3)  !< query point
      real*8, intent(out) :: dd    !< signed distance
      real*8, intent(out) :: CP(3) !< closest point on the surface
      integer :: cindx
      dd=dreqinput
      call searchADT_mindist(dd,CP,adtIntegers,adtReals,adtExtents,bfaceinfo,&
                 xp,cindx,nelem,6)
    end subroutine MinDist
    !
    subroutine DistCleanup
      if (initdist==1) then 
       deallocate(bface,bfaceinfo,coord,adtIntegers,adtReals,adtExtents)
      endif
      initdist=0 
    end subroutine DistCleanup
    !
  end module DistCalc
      
      
    
    
