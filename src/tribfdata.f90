!==================================================================================
subroutine tribfdata(xbgeom,bface,bfaceinfo,nnode,nbface)
implicit none
!
! Local variables
!
integer, intent(in) :: nbface,nnode
real*8, intent(in) :: xbgeom(3,nnode)
integer, intent(in) :: bface(3,nbface)
real*8, intent(out) :: bfaceinfo(22,nbface)
!
integer :: ibface,BF(3)
real*8 :: a,b,c,det,detinv,denom,a2,b2,c2,fA,fB,fC
real*8 :: AA(3),BB(3),CC(3),E0(3),E1(3),n(3)
!
do ibface=1,nbface
   BF=bface(:,ibface)
   AA=xbgeom(:,BF(1))
   BB=xbgeom(:,BF(2))
   CC=xbgeom(:,BF(3))
   E0=CC-BB
   E1=AA-BB
   call DotProd(a,E0,E0)
   call DotProd(b,E0,E1)
   call DotProd(c,E1,E1)
   det=a*c-b*b
   detinv=1.0/det
   denom=a-2*b+c
   bfaceinfo(1,ibface)=a
   bfaceinfo(2,ibface)=b
   bfaceinfo(3,ibface)=c
   bfaceinfo(4,ibface)=det
   bfaceinfo(5,ibface)=detinv
   bfaceinfo(6,ibface)=denom
   call CrossProd(n,E0,E1)
   bfaceinfo(7,ibface)=n(1)
   bfaceinfo(8,ibface)=n(2)
   bfaceinfo(9,ibface)=n(3)
   a2=a
   c2=c
   b2=denom
   fA=a2*(b2+c2-a2)
   fB=b2*(a2+c2-b2)
   fC=c2*(a2+b2-c2)
   denom=2*(a2*b2+a2*c2+b2*c2)-(a2*a2+b2*b2+c2*c2)
   bfaceinfo(10:12,ibface)=(fA*AA+fB*BB+fc*CC)/denom
   a=sqrt(a2)
   b=sqrt(b2)
   c=sqrt(c2)
   bfaceinfo(13,ibface)=sqrt(a2*b2*c2/ &
        ((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)))
   bfaceinfo(14:16,ibface)=AA
   bfaceinfo(17:19,ibface)=BB
   bfaceinfo(20:22,ibface)=CC
end do
end subroutine Tribfdata
