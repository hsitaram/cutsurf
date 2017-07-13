!==================================================================================
subroutine Distance2(d2,P1,P2)
implicit none
real*8, intent(in) :: P1(3),P2(3)
real*8, intent(out) :: d2
real*8 :: a,b,c
a=P1(1)-P2(1)
b=P1(2)-P2(2)
c=P1(3)-P2(3)
d2=a*a+b*b+c*c
return
end subroutine Distance2
!===================================================================================
subroutine CrossProd(v,v1,v2)
implicit none
real*8, intent(in) :: v1(3),v2(3)
real*8, intent(out) :: v(3)
v(1)=v1(2)*v2(3)-v1(3)*v2(2)
v(2)=v1(3)*v2(1)-v1(1)*v2(3)
v(3)=v1(1)*v2(2)-v1(2)*v2(1)
end subroutine CrossProd
!===================================================================================
subroutine DotProd(dp,v1,v2)
implicit none
real*8, intent(in) :: v1(3),v2(3)
real*8, intent(out) :: dp
dp=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
end subroutine DotProd
!===================================================================================
     




























































































































































