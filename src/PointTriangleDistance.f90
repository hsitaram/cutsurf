subroutine ComputePointTriangleDist(AA,BB,CC,Mdata,P,dist2,CP,s,t)
!
implicit none

real*8, intent(in) :: AA(3),BB(3),CC(3),Mdata(9)
real*8, intent(in) :: P(3)
real*8, intent(out) :: dist2
real*8, intent(out) :: CP(3),s,t
!
! Local variables
!
real*8 :: V(3),DD(3),E0(3),E1(3),n(3)
real*8 :: a,b,c,d,e,det,detinv,denom
real*8 :: numer,tmp0,tmp1,dp
integer :: BF(3),iregion
!
! Begin
!
n(1)=Mdata(7)
n(2)=Mdata(8)
n(3)=Mdata(9)
DD=BB-P
call DotProd(dp,n,DD)
if (dp>1e-10 .and. 1.eq.0) then
   dist2=1e15
   CP=0.0
   s=0.0
   t=0.0
   return
end if

a=Mdata(1)
b=Mdata(2)
c=Mdata(3)
det=Mdata(4)
detinv=Mdata(5)
denom=Mdata(6)
!
E0=CC-BB
E1=AA-BB
!
d=E0(1)*DD(1)+E0(2)*DD(2)+E0(3)*DD(3)
e=E1(1)*DD(1)+E1(2)*DD(2)+E1(3)*DD(3)
!
s=b*e-c*d;
t=b*d-a*e;

if ((s+t)<=det) then
   if (s<0) then
      if (t<0) then
         iregion=4;
      else
         iregion=3;
      end if
   else if (t<0) then
      iregion=5;
   else
      iregion=0;
   end if
else
   if (s<0) then
      iregion=2;
   else if (t<0) then
      iregion=6;
   else
      iregion=1;
   end if
end if

select case(iregion)
   
case(0)
   s=s*detinv;
   t=t*detinv;
   
case(1) 
   numer=c+e-b-d;
   if (numer<=0) then
      s=0;
   else
      if (numer>=denom) then
         s=1;
      else
         s=numer/denom;
      end if
   end if
   t=1-s;

case(2)
   tmp0=b+d;
   tmp1=c+e;
   if (tmp1>tmp0) then
      numer=tmp1-tmp0;
      if (numer>=denom) then
         s=1;
      else
         s=numer/denom;
      end if
      t=1-s;
   else
      s=0;
      if (tmp1<=0) then
         t=1;
      elseif (e>=0) then
         t=0;
      else
         t=-e/c;
      end if
   end if

case(3)
   s=0;
   if (e>=0) then
      t=0;
   else if (c+e<=0) then
      t=1;
   else
      t=-e/c;
   end if

case(4)
   if (d<0) then
      t=0;
      if (a+d<=0) then
         s=1;
      else
         s=-d/a;
      end if
   else
      s=0;
      if (c+e<=0) then
         t=1;
      else if (e>=0) then
         t=0;
      else
         t=-e/c;
      end if
   end if

case(5)
   t=0;
   if (d>=0) then
      s=0;
   else if (a+d<=0) then
      s=1;
   else
      s=-d/a;
   end if
       
case(6)
   tmp0=b+e;
   tmp1=a+d;
   if (tmp1>tmp0) then
      numer=tmp1-tmp0;
      if (numer>=denom) then
         t=1;
      else
         t=numer/denom;
      end if
      s=1-t;
   else
      t=0;
      if (tmp1<=0) then
         s=1;
      else if (d>=0) then
         s=0;
      else
         s=-d/a;
      end if
   end if

case DEFAULT
   write(*,*) 'ERROR in find triangle-pt distance'

end select

CP=BB+s*E0+t*E1
V=CP-P
dist2=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)
!
return
end subroutine ComputePointTriangleDist
