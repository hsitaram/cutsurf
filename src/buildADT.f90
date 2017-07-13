!=========================================================================
subroutine buildADT(coord,adtIntegers,adtReals,adtExtents,nelem,ndim)
implicit none
!
!  Subroutine : buildADT
!
!  Builds Alternating Digital Tree for a given
!  set of elemets
!
!  e.g : an element in 3-D is a bounding box 
!        ndim = 6 
!        coord = [lower corner coord;upper corner coord] per element
!
!
!  Arguments:
!
!  inputs:
!
!  nelem             - number of elements
!  ndim              - element dimension
!  coord(ndim,nelem) - element coordinates
!  
!  outputs:
!
!  adtIntegers(4,nelem) - [elementIndex in ADT, left child, right child,
!                          tree position of each element] 
!  adtReals(nelem)      - median coordinate 
!  adtExtents(ndim,2)   - extents of search space
!
!
!  subroutine arguments
!
integer, intent(in) :: nelem,ndim
real*8, intent(in) :: coord(ndim,nelem)
!
integer, intent(inout) :: adtIntegers(4,nelem)
real*8, intent(inout) :: adtReals(ndim,nelem)
real*8, intent(inout) :: adtExtents(ndim,2)
!
! local variables 
!
integer :: i,nd, level
integer, allocatable :: ElementsAvailable(:)
real*8, allocatable :: adtWork(:)
!real*8, allocatable :: adtRealsTemp(:)
real*8 :: margin, dcoord
integer :: adtcount,side,parent,nav
!
! begin
!
!
! Allocate memory
!
!do i=1,nelem
! write(6,100) coord(:,i)
!enddo
100 format(6(1X,F10.4))
allocate(ElementsAvailable(nelem))
allocate(adtWork(nelem))
!allocate(adtRealsTemp(nelem))
!
! Determine extent of elements
!
nd=ndim/2
do i=1,ndim
   adtExtents(i,1)=minval(coord(i,:))
   adtExtents(i,2)=maxval(coord(i,:))
enddo
!
! Determine extent of ADT (add small margin 
! to avoid boundary problems)
!
margin=0.0
do i=1,ndim
   dcoord=margin*(adtExtents(i,2)-adtExtents(i,1))
   adtExtents(i,1)=adtExtents(i,1)-dcoord
   adtExtents(i,2)=adtExtents(i,2)+dcoord
enddo
!
!do i=1,nd
!   adtExtents(nd+i,1)=adtExtents(i,1)
!   adtExtents(nd+i,2)=adtExtents(i,2)
!enddo
!
! Built ADT by recursion
!
do i=1,nelem
   ElementsAvailable(i)=i
enddo
!
adtcount=0
side=0
parent=0
level=0
nav=nelem
!
call BuildADTrecursion(coord,adtReals,adtWork, &
     adtIntegers,ElementsAvailable, &
     adtcount,side,parent,level,ndim,nelem,nav)
!
!write(6,199) adtReals
!199 format('adtReals=',8(1X,F10.4))
!write(6,*) '--------------------------'
!write(6,*) 'adtCount=',adtCount
!do i=1,nelem
!   write(6,*) adtIntegers(1:3,i)
!enddo
!write(6,*) '--------------------------'
!
! create inverse map
!
do i=1,nelem
   adtIntegers(4,adtIntegers(1,i))=i
!   adtReals(i)=adtRealsTemp(i)
enddo
!
! Release memory 
!
deallocate(ElementsAvailable)
deallocate(adtWork)

return
end subroutine buildADT

!===========================================================
recursive subroutine BuildADTrecursion(coord,adtReals,adtWork, &
     adtIntegers,ElementsAvailable, &
     adtcount,side,parent,level,ndim,nelem,nav)
!
implicit none
!
! Subroutine: BuildADTrecursion
!
integer, intent(in) :: ndim,nelem,nav
integer, intent(in) :: side,parent,level
real*8,  intent(in) :: coord(ndim,nelem)
!
integer, intent(inout) :: adtCount
real*8,  intent(inout) :: adtReals(ndim,nelem)
real*8,  intent(inout) :: adtWork(nelem)
integer, intent(inout) :: adtIntegers(4,nelem),elementsAvailable(nav)
!
!
! Local variables
!
integer :: i,j,dimcut,nleft,parentToChild,nd
real*8  :: coordmid
!
!
!
nd=ndim/2
!
!write(6,*) 'elements=',elementsAvailable
!
if (nav > 1) then

   dimcut=1+mod(level,ndim)
   ! Coordinates along dimension dimcut:
   do i=1,nav
      adtWork(i)=coord(dimcut,ElementsAvailable(i))
   enddo
   ! Re-order ElementsAvailable with nleft 
   ! elements to the left of median of adtWork
   call median(ElementsAvailable,adtWork,nav,coordmid)

   ! Add median element to ADT (with no children)
   nleft=(nav+1)/2
   adtcount=adtcount+1
   adtIntegers(1,adtcount)=ElementsAvailable(nleft)
   adtIntegers(2,adtcount)=0
   adtIntegers(3,adtcount)=0
   !
   ! find minimum and maximum bounds of the elements
   ! contained in this leaf
   !
   adtReals(1:nd,adtCount)=1E15
   adtReals(nd+1:ndim,adtCount)=-1E15
   !
   do i=1,nav
      do j=1,nd
         adtReals(j,adtCount)=min(adtReals(j,adtCount),coord(j,ElementsAvailable(i)))
         adtReals(j+nd,adtCount)=max(adtReals(j+nd,adtCount),coord(j+nd,ElementsAvailable(i)))
      enddo
   enddo
   !
   ! Specify that new element is child of a parent (unless root)
   !
   if (side > 1) then
      adtIntegers(side,parent)=ElementsAvailable(nleft)
   endif

   parentToChild=adtcount

   ! Continue building left side of tree
    if (nleft > 1) &
    call BuildADTrecursion(coord,adtReals,adtWork,&
         adtIntegers,ElementsAvailable(1), &
         adtcount,2,parentToChild,level+1,ndim,nelem,nleft-1)

   ! Continue building right side of tree
   call BuildADTrecursion(coord,adtReals,adtWork, &
        adtIntegers,ElementsAvailable(nleft+1), &
        adtcount,3,parentToChild,level+1,ndim,nelem,nav-nleft)

elseif (nav==1) then

   ! Region can not be divided: add point with no child
   adtcount=adtcount+1
   adtIntegers(1,adtcount)=ElementsAvailable(1)
   adtIntegers(2,adtcount)=0
   adtIntegers(3,adtcount)=0
   do j=1,ndim
      adtReals(j,adtCount)=coord(j,ElementsAvailable(1))
   enddo
   if (side > 1) then
      adtIntegers(side,parent)=ElementsAvailable(1)
   endif
endif
!write(6,*) 'adtCount=',adtCount,adtIntegers(1,adtCount)
!write(6,199) adtReals
!199 format('adtReals=',8(1X,F10.4))

return
end subroutine BuildADTrecursion
!
!======================================================================
!
SUBROUTINE median(ix, x, n, xmed)

! Find the median of X(1), ... , X(N), using as much of the quicksort
! algorithm as is needed to isolate it.
! N.B. On exit, the array X is partially ordered.

!     Latest revision - 26 November 1996
IMPLICIT NONE

INTEGER, INTENT(IN)                :: n
INTEGER, intent(IN OUT)            :: ix(n)
REAL*8, INTENT(IN OUT)             :: x(n)
REAL*8, INTENT(OUT)                  :: xmed

! Local variables

REAL*8    :: temp, xhi, xlo, xmax, xmin
LOGICAL :: odd
INTEGER :: hi, lo, nby2, nby2p1, mid, i, j, k,itemp

nby2 = n / 2
nby2p1 = nby2 + 1
odd = .true.

!     HI & LO are position limits encompassing the median.

IF (n == 2 * nby2) odd = .false.
lo = 1
hi = n
IF (n < 3) THEN
  IF (n < 1) THEN
    xmed = 0.0
    RETURN
  END IF
  xmed = x(1)
  IF (n == 1) RETURN
  xmed = 0.5*(xmed + x(2))
   if (x(2) < x(1)) then
    temp=x(1)
    x(1)=x(2)
    x(2)=temp
    itemp=ix(1)
    ix(1)=ix(2)
    ix(2)=itemp
   endif
  RETURN
END IF

!     Find median of 1st, middle & last values.

10 mid = (lo + hi)/2
xmed = x(mid)
xlo = x(lo)
xhi = x(hi)
IF (xhi < xlo) THEN          ! Swap xhi & xlo
  temp = xhi
  xhi = xlo
  xlo = temp
END IF
IF (xmed > xhi) THEN
  xmed = xhi
ELSE IF (xmed < xlo) THEN
  xmed = xlo
END IF

! The basic quicksort algorithm to move all values <= the sort key (XMED)
! to the left-hand end, and all higher values to the other end.

i = lo
j = hi
50 DO
  IF (x(i) >= xmed) EXIT
  i = i + 1
END DO
DO
  IF (x(j) <= xmed) EXIT
  j = j - 1
END DO
IF (i < j) THEN
  temp = x(i)
  x(i) = x(j)
  x(j) = temp

  itemp=ix(i)
  ix(i)=ix(j)
  ix(j)=itemp

  i = i + 1
  j = j - 1

!     Decide which half the median is in.

  IF (i <= j) GO TO 50
END IF

IF (.NOT. odd) THEN
  IF (j == nby2 .AND. i == nby2p1) GO TO 130
  IF (j < nby2) lo = i
  IF (i > nby2p1) hi = j
  IF (i /= j) GO TO 100
  IF (i == nby2) lo = nby2
  IF (j == nby2p1) hi = nby2p1
ELSE
  IF (j < nby2p1) lo = i
  IF (i > nby2p1) hi = j
  IF (i /= j) GO TO 100

! Test whether median has been isolated.

  IF (i == nby2p1) RETURN
END IF
100 IF (lo < hi - 1) GO TO 10

IF (.NOT. odd) THEN
  if (x(nby2p1) < x(nby2)) then
   temp=x(nby2p1)
   x(nby2p1)=x(nby2)
   x(nby2)=temp
   itemp=ix(nby2p1)
   ix(nby2p1)=ix(nby2)
   ix(nby2)=itemp 
  endif
  xmed = 0.5*(x(nby2) + x(nby2p1))
  RETURN
END IF
temp = x(lo)
itemp = ix(lo)
IF (temp > x(hi)) THEN
  x(lo) = x(hi)
  x(hi) = temp
  ix(lo)=ix(hi)
  ix(hi)=itemp
END IF
xmed = x(nby2p1)
RETURN

! Special case, N even, J = N/2 & I = J + 1, so the median is
! between the two halves of the series.   Find max. of the first
! half & min. of the second half, then average.

130 xmax = x(1)
DO k = lo, j
  xmax = MAX(xmax, x(k))
END DO
xmin = x(n)
DO k = i, hi
  xmin = MIN(xmin, x(k))
END DO
xmed = 0.5*(xmin + xmax)

RETURN
END SUBROUTINE median
