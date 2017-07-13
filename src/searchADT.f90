!=======================================================================
subroutine searchADT_mindist(dmin,CPval,adtIntegers,adtReals,adtExtents,coord,&
     xsearch,cellIndex,nelem,ndim)
!
! subroutine arguments
! 
implicit none
!
integer, intent(in) :: nelem,ndim
real*8, intent(in) :: coord(22,nelem)
!
real*8, intent(inout) :: dmin
integer, intent(in) :: adtIntegers(4,nelem)
real*8, intent(in) :: adtReals(ndim,nelem)
real*8, intent(in) :: adtExtents(ndim,2)
!
real*8, intent(in) :: xsearch(ndim/2)
integer, intent(inout) :: cellIndex
real*8, intent(out) :: CPVal(3)
!
! local variables
!
integer :: rootNode
integer :: i
character*(80) :: integer_string,logfile,integer_string1
real*8 :: ddist,dsgn,dmag
integer :: itrack
integer, save :: ivol=0
real*8 :: TOL=1E-6
!
! begin
!
! first create search space
!
cellIndex=0
rootNode=1
!
! search for intersections with elements now
!
ivol=ivol+1
dsgn=1d0
dmag=0d0
!
call searchIntersections_mindist(dmin,dsgn,dmag,CPval,cellIndex,adtIntegers,adtReals,coord,&
        0,rootNode,xsearch,nelem,ndim)
dmin=dsgn*sqrt(dmin)
!
return
end subroutine searchADT_mindist
!
!=======================================================================
recursive subroutine searchIntersections_mindist(dmin,dsgn,dmag,CPval,cellIndex,adtIntegers,adtReals,&
     coord,level,node,xsearch,nelem,ndim) !LC0,UC0,LC,UC,nelem,ndim)
!
implicit none
!
real*8, intent(inout) :: dmin,dsgn,dmag
real*8, intent(inout) :: CPval(3)
integer, intent(in) :: nelem,ndim,level,node
integer, intent(in) :: adtIntegers(4,nelem)
real*8, intent(in) :: adtReals(ndim,nelem)
real*8, intent(in) :: coord(22,nelem)
real*8, intent(in) :: xsearch(ndim/2)
!
integer, intent(inout) :: cellIndex
!
! local variables
!
logical :: checkRegionOverlap,IsElementInRegion
logical :: Iin,Iover
integer :: d,nodeChild,dimcut
real*8  :: element(22),coordmid
real*8  :: LCC(ndim),UCC(ndim)
real*8  :: TOL=1E-6
real*8  :: dd,boxDist2
real*8  :: CP(3),s,t
real*8  :: dd2,dtest,vmag,emag
!
! begin
!
element=coord(:,adtIntegers(1,node))
! distance to the surface of the circumsphere of the element
dd=sqrt(dot_product(xsearch-element(10:12),xsearch-element(10:12)))-element(13)
if (dd < sqrt(dmin) ) then
   call ComputePointTriangleDist(element(14:16),element(17:19),element(20:22),&
        element(1:9),xsearch,dd2,CP,s,t)
   !
   ! determine sign of distance based on best visible
   ! element
   !
   if (abs(dd2-dmin) < 1e-12) then
      emag=sqrt(dot_product(element(7:9),element(7:9)))
      vmag=sqrt(dot_product(xsearch-CP,xsearch-CP))
      dtest=abs(dot_product(xsearch-CP,element(7:9))/vmag/emag)
      if (dtest > dmag) then
         dsgn=sign(1d0,dot_product(xsearch-CP,element(7:9)))
         dmag=dtest
      endif
   endif
   if (dd2 < dmin-1e-12) then
      dmin=dd2
      cellIndex=adtIntegers(1,node)
      CPval=CP
      vmag=sqrt(dot_product(xsearch-CP,xsearch-CP))
      emag=sqrt(dot_product(element(7:9),element(7:9)))
      dmag=abs(dot_product(xsearch-CP,element(7:9))/vmag/emag)
      dsgn=sign(1d0,dot_product(xsearch-CP,element(7:9)))
   endif
endif
!
do d=2,3
   if (adtIntegers(d,node) > 0) then
      nodeChild=adtIntegers(4,adtIntegers(d,node))
      if (boxDist2(xsearch,adtReals(:,nodeChild)) < dmin ) then
         call searchIntersections_mindist(dmin,dsgn,dmag,CPval,cellIndex,adtIntegers,adtReals,coord,&
              level+1,nodeChild,xsearch,nelem,ndim)
      endif
   endif
enddo
!
return
end subroutine searchIntersections_mindist
!
function boxDist2(xp,box)
  real*8, intent(in) :: xp(3)
  real*8, intent(in) :: box(6)
  real*8 :: boxDist2
  real*8 :: dx(3)
  integer :: j

  do j=1,3
     dx(j)=maxval((/box(j)-xp(j),0d0,xp(j)-box(j+3)/))
  enddo
  boxDist2=dx(1)**2+dx(2)**2+dx(3)**2
end function boxDist2
