subroutine mlest (x,beta,mmin,mmax)
! Maximum likelihood estimator for the parameters of a Pareto law
implicit none
double precision, dimension(:), intent(in) :: x
double precision, dimension(size(x)) :: lgx
double precision, intent(out) :: beta,mmin,mmax
double precision :: a,sigma
double precision :: t,y,z,n

lgx=log(x)
n=float(size(lgx))
t=sum(lgx)
y=minval(x)
z=maxval(x)

 call mlebis(n,t,y,z,a)

beta=a+1.
sigma=(n/(n-2))**2 *a**2/(n-3)
sigma=sqrt(sigma)
mmin=z
mmax=z

contains

	subroutine mlebis (n,t,y,z,a)
	! Auxiliary funcion for the subroutine mlest
	! Calculate the minimum of the log max-likelihood functn
	! for data from a truncated Pareto distribution
	implicit none
	double precision, intent(in) :: n,t,y,z
	double precision, intent(out) :: a
	double precision :: f,fmid,xmid
	double precision :: rtbis,dx
	double precision :: left,right
	integer :: j
	
	left=-10.
	right=10.01
	
	f=mle(n,t,y,z,left)
	fmid=mle(n,t,y,z,right)
	do 
	if (f*fmid < 0) exit
	if (f*fmid > 0) then
!	print *,"Mod ML error in functn",f,fmid
	left=left*3
	right=right*3
	f=mle(n,t,y,z,left)
	fmid=mle(n,t,y,z,right)
	end if
	if ( abs(right) > 100000) exit
	end do
	if (abs(right) >= 100000) then
		print *,"some problem in ML estimator" 
		print *,"press any key to continue"
		read(*,*)
	end if
	if (f*fmid > 0) then
		print *,"still error in ML estimator"
		print *,"press any key to continue"
		read(*,*)
	end if
	
	if (f<0) then
	rtbis=left
	dx=right-left
	else
	rtbis=right
	dx=left-right
	end if
	j=1
	do 
	dx=dx*0.5
	xmid=rtbis+dx
	fmid=mle(n,t,y,z,xmid)
	if (fmid <= 0) rtbis=xmid
	if (abs(dx) < 0.0000001 .or. fmid==0) exit
	j=j+1
	!if (j>10000) pause "bisection"
	end do
	a=rtbis
	end subroutine mlebis

	function mle (n,t,y,z,beta)
	implicit none
	double precision :: mle,n,t,y,z,beta
	mle=n/beta + n*(z/y)**beta*log(z/y)/(1-(z/y)**beta) + n*log(z) - t
	end function mle
	
end subroutine mlest
