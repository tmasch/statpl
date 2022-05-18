program statpl
use statfunmod
use mffunmod
use mfestmod
use statpl_mod
use utilmod
implicit none

integer :: i,j
double precision :: tmp,tmpx,tmpy,step
character(len=20) :: chartmp

double precision :: infinity

double precision, dimension(:), allocatable :: masses

! input parameters
character (len=50) :: datafile
integer :: calcstddev
integer :: montecarlosize
integer :: nbin
integer :: makeplots
double precision :: significance
double precision, dimension(nestimators) :: est_mask



! Results
integer :: ndata
double precision :: beta,mmin,mmax,stddev
double precision :: dsp_inf,dsp_trunc


double precision, dimension(ntests) :: critval_inf,critval_trunc,power
double precision, dimension(ntests) :: statistics_data_inf,statistics_data_trunc
double precision, dimension(ntests,3,2) :: statval_data
double precision, dimension(:,:), allocatable :: stats_inf,stats_trunc
character (len=50), dimension(ntests) :: stat_names
character (len=50), dimension(nestimators) :: est_names
double precision, dimension(nestimators,3) :: estout
double precision, dimension(nestimators,9) :: eststat



!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
!!!!!!! Initialise
!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
infinity=1.e30
est_names=""
call name_estimators (est_names)
stat_names=""
call name_statistics (stat_names)

eststat=0.
statistics_data_inf=0.
statistics_data_trunc=0.

open(1000,file="statmfout.txt",action="write")



!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
!!!!!!! Read the parameters
!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
! Name of the data files
write(*,*) "Data file:"
read(*,*) datafile
write(*,*) "using ",datafile
write(1000,*) "Data file: "
write(1000,*) datafile

! Number of bins for binning (-1 = auto calculate)
write(*,*) "Input number of bins (-1 = auto):"
read(*,*) nbin
write(*,*) "using ",nbin
write(1000,*) "Number of bins (-1 = auto):"
write(1000,*) nbin


! Switches for the individual estimators
write(*,*) "Which estimators should be used"
write(1000,*) "Which estimators should be used"
do i=1,nestimators
write(*,*) est_names(i)
read(*,*) est_mask(i)
write(1000,*) est_names(i),est_mask(i)
write(*,*) "using ",est_names(i),est_mask(i)
end do

! Switch whether standard deviation should be calculated
write(*,*) "Calculate Standard deviation etc."
read(*,*) calcstddev
write(*,*) "using ",calcstddev
write(1000,*) "Calculate Standard deviation etc."
write(1000,*) calcstddev

! Size of the Monte-Carlo experiment to calculcate std. dev. and tests
write(*,*) "Monte-Carlo sampling size:"
read(*,*) montecarlosize
write(*,*) "using ",montecarlosize
write(1000,*) "Monte-Carlo sampling size:"
write(1000,*) montecarlosize

! Significance for the tests
write(*,*) "Significance Level of the tests"
read(*,*) significance
write(*,*) "using ",significance
write(1000,*) "Significance Level of the tests"
write(1000,*) significance

! Check whether plots should be made
write(*,*) "Make plots?"
read(*,*) makeplots
write(*,*) makeplots
write(1000,*) "Make plots?"
write(1000,*) makeplots


allocate(stats_inf(montecarlosize,ntests))
allocate(stats_trunc(montecarlosize,ntests))


!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
!!!!!!! Read the data
!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!

ndata=numberoflines(datafile)
write(1000,*) "Number of data"
write(1000,*) ndata
write(*,*) "Number of data"
write(*,*) ndata
allocate(masses(ndata))

open(1,file=datafile,action="read")
do i=1,ndata
read (1,*) masses(i)
end do
close(1)

call sort(masses)

write(1000,*) "Minimum of data"
write(1000,*) masses(1)
write(*,*) "Minimum of data"
write(*,*) masses(1)
write(1000,*) "Maximum of data"
write(1000,*) masses(ndata)
write(*,*) "Maximum of data"
write(*,*) masses(ndata)

if (nbin<0) nbin=floor(size(masses)**(2./5.))

write(*,*) "Used number of bins:"
write(*,*) nbin
write(1000,*) "Used number of bins:"
write(1000,*) nbin


!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
!!!!!!! Call the Estimators
!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!


estout=0

write(1000,*) "Results of the estimators (Binning: ",nbin," bins)"
write(*,*) "Results of the estimators (Binning: ",nbin," bins)"

call estimate_parameters(est_mask,masses,nbin,estout)

! Use the ModModML-Estimates in the following
beta=estout(imodml,1)
mmin=estout(imodml,2)
mmax=estout(imodml,3)

!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
!!!!!!! Calculate the statistics
!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!

if (calcstddev > 0) then
call calceststat (montecarlosize,size(masses),beta,mmin,mmax,&
nbin,est_mask,eststat,significance,critval_inf,critval_trunc,power,stats_inf,stats_trunc)
call calculate_statistics(masses,beta,mmin,infinity,statistics_data_inf)
call calculate_statistics(masses,beta,mmin,mmax,statistics_data_trunc)
end if

write(*,*) "ALPHA:              : Estimate   : Average    : Std. Dev.  : Bias       :"
write(1000,*) "ALPHA:              : Estimate   : Average    : Std. Dev.  : Bias       :"
do i=1,nestimators
write(1000,'(A20," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ")') &
est_names(i),estout(i,1),eststat(i,1),eststat(i,2),eststat(i,3)
write(*,'(A20," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ")') &
est_names(i),estout(i,1),eststat(i,1),eststat(i,2),eststat(i,3)
end do
write(*,*) "X MAX:              : Estimate   : Average    : Std. Dev.  : Bias       :"
write(1000,*) "X MAX:              : Estimate   : Average    : Std. Dev.  : Bias       :"
do i=1,nestimators
write(1000,'(A20," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ")') &
est_names(i),estout(i,3),eststat(i,7),eststat(i,8),eststat(i,9)
write(*,'(A20," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ")') &
est_names(i),estout(i,3),eststat(i,7),eststat(i,8),eststat(i,9)
end do

stddev=eststat(imodml,2)


!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
!!!!!!! Some Output
!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
! Estimated Values
write(1000,*) "Modified Maximum Likelihood"
write(1000,*) "Estimated exponent"
write(1000,*) beta
write(1000,*) "Standard deviation of the estimated exponent"
write(1000,*) stddev
write(1000,*) "Estimated lower limit"
write(1000,*) mmin
write(1000,*) "Estimated upper limit"
write(1000,*) mmax

write(*,*) "Modified Maximum Likelihood results, used furtheron"
write(*,*) "Estimated exponent"
write(*,*) beta
write(*,*) "Standard deviation of the estimated exponent"
write(*,*) stddev
write(*,*) "Estimated lower limit"
write(*,*) mmin
write(*,*) "Estimated upper limit"
write(*,*) mmax

! Results of the statistics
write(1000,*) "TESTS FOR NOT TRUNCATED POWER LAW, significance ",significance
write(*,*) "TESTS FOR NOT TRUNCATED POWER LAW, significance ",significance
write(1000,*) "Statistic           : data       : CritValInf : Perc. Sig. : 1-P.S.     : Power      : Conclusion"
write(*,*) "Statistic           : data       : CritValInf : Perc. Sig. : 1-P.S.     : Power      : Conclusion"
statval_data=0.
do i=1,ntests
do j=2,montecarlosize
if (statistics_data_inf(i) <= stats_inf(j,i) .AND. statval_data(i,1,1) < 1. ) then
	statval_data(i,1,1)=float(j-1)
	statval_data(i,1,2)=stats_inf(j-1,i)
	statval_data(i,2,1)=float(j)
	statval_data(i,2,2)=stats_inf(j,i)
	tmp= ( statistics_data_inf(i) - stats_inf(j-1,i) ) / ( stats_inf(j,i) - stats_inf(j-1,i))
	statval_data(i,3,1)=(tmp + float(j-1))/float(montecarlosize)*100
	statval_data(i,3,2)=statistics_data_inf(i)
end if
end do
if (statistics_data_inf(i) <= stats_inf(1,i)) then
statval_data(i,1,1)=1.
statval_data(i,1,2)=stats_inf(1,i)
statval_data(i,2,1)=float(montecarlosize)
statval_data(i,2,2)=stats_inf(montecarlosize,i)
statval_data(i,3,1)=0.
statval_data(i,3,2)=statistics_data_inf(i)
end if
if (statistics_data_inf(i) >= stats_inf(montecarlosize,i)) then
j=montecarlosize
statval_data(i,1,1)=1.
statval_data(i,1,2)=stats_inf(1,i)
statval_data(i,2,1)=float(montecarlosize)
statval_data(i,2,2)=stats_inf(montecarlosize,i)
statval_data(i,3,1)=100.
statval_data(i,3,2)=statistics_data_inf(i)
end if

chartmp="not truncated"
if (statistics_data_inf(i) > critval_inf(i)) chartmp="*** TRUNCATED ***"
write(*,'(A20," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ",A20)') &
stat_names(i),statistics_data_inf(i),critval_inf(i),statval_data(i,3,1),100.-statval_data(i,3,1),power(i),chartmp
write(1000,'(A20," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ",A20)')  &
stat_names(i),statistics_data_inf(i),critval_inf(i),statval_data(i,3,1),100.-statval_data(i,3,1),power(i),chartmp
end do




statval_data=0

write(*,*) "TESTS FOR TRUNCATED POWER LAW"
write(1000,*) "TESTS FOR TRUNCATED POWER LAW"
write(1000,*) "Statistic           : data       : CrValTrunc : Perc. Sig. : 1-P.S.     : Conclusion"
write(*,*) "Statistic           : data       : CrValTrunc : Perc. Sig. : 1-P.S.     : Conclusion"
	do i=1,ntests-4

! Determine first the percentiles for the data stemming from the hypothesis truncated
do j=2,montecarlosize
if (statistics_data_trunc(i) <= stats_trunc(j,i) .AND. statval_data(i,1,1) < 1. ) then
	statval_data(i,1,1)=float(j-1)
	statval_data(i,1,2)=stats_trunc(j-1,i)
	statval_data(i,2,1)=float(j)
	statval_data(i,2,2)=stats_trunc(j,i)
	tmp= ( statistics_data_trunc(i) - stats_trunc(j-1,i) ) / ( stats_trunc(j,i) - stats_trunc(j-1,i))
	statval_data(i,3,1)=(tmp + float(j-1))/float(montecarlosize)*100
	statval_data(i,3,2)=statistics_data_trunc(i)
end if
end do
if (statistics_data_trunc(i) <= stats_trunc(1,i)) then
statval_data(i,1,1)=1.
statval_data(i,1,2)=stats_trunc(1,i)
statval_data(i,2,1)=float(montecarlosize)
statval_data(i,2,2)=stats_trunc(montecarlosize,i)
statval_data(i,3,1)=0.
statval_data(i,3,2)=statistics_data_trunc(i)
end if
if (statistics_data_trunc(i) >= stats_trunc(montecarlosize,i)) then
j=montecarlosize
statval_data(i,1,1)=1.
statval_data(i,1,2)=stats_trunc(1,i)
statval_data(i,2,1)=float(montecarlosize)
statval_data(i,2,2)=stats_trunc(montecarlosize,i)
statval_data(i,3,1)=100.
statval_data(i,3,2)=statistics_data_trunc(i)
end if


! Check if data are consistent with a truncated power law and output
chartmp="*** CONSISTENT ***"
if (statistics_data_trunc(i) > critval_trunc(i)) chartmp="not consistent"
write(*,'(A20," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ",A20)') &
stat_names(i),statistics_data_trunc(i),critval_trunc(i),statval_data(i,3,1),100.-statval_data(i,3,1),chartmp
write(1000,'(A20," : ",F10.3," : ",F10.3," : ",F10.3," : ",F10.3," : ",A20)') &
stat_names(i),statistics_data_trunc(i),critval_trunc(i),statval_data(i,3,1),100.-statval_data(i,3,1),chartmp
end do


dsp_inf=0.001
dsp_trunc=0.001

! Changing and setting the variables for better readability below
dsp_inf=critval_inf(4)
dsp_trunc=critval_trunc(4)

if (makeplots == 1) then
	
	!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
	!!!!!!! Plot the stabilized probability (SPP) plot
	!!!!!!! Null hypothesis (diagonal) = not truncated/infinite PL
	!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
	
	open(1,file="sppinf_data.txt",action="write")
	open(2,file="sppinf_estimated.txt",action="write")
	open(3,file="sppinf_infinity.txt",action="write")
	do i=1,size(masses)
		tmpx=intpdfpl(masses(i),beta,mmin,infinity)
	! Data
		tmpy=(float(i)-0.5)/float(size(masses))
		write(1,*) stabquant(tmpx),stabquant(tmpy)
	end do
	tmp=0.9999
	do
		tmpy=tmp
		tmp=invintpdfpl(tmp,beta,mmin,infinity)
	! Theory, estimated 
		tmpx=intpdfpl(tmp,beta,mmin,mmax)
		if (tmp<mmax) write(2,*) stabquant(tmpy),stabquant(tmpx)
	! Theory, infinite
		tmpx=intpdfpl(tmp,beta,mmin,infinity)
		write(3,*) stabquant(tmpy),stabquant(tmpx)	
	! The following is for an approximately equal spacing of the plotting points	
		tmp=tmpy
		if(tmp>0.85) tmp=tmp-0.005
		if(tmp<=0.85 .AND. TMP >0.1) tmp=tmp-0.01
		if(tmp < 0.1) tmp=tmp-0.005
	if (tmp < 0.) exit
	end do
	close(1);close(2);close(3);
	
	! Plot the reference line (=diagonal from 0,0 to 1,1) and the 
	! acceptance region (d)
	open(1,file="sppinf_dplus.txt",action="write")
	write(1,*) 0,dsp_inf
	write(1,*) 1-dsp_inf,1
	close(1)
	open(1,file="sppinf_dminus.txt",action="write")
	write(1,*) dsp_inf,0
	write(1,*) 1,1-dsp_inf
	close(1)
	
	! Write some ticks and labels for the plot 
	! Calculate the masses corresponding to the theoretical quantiles
	! (used for plotting)
	tmp=log10(maxval(masses))-log10(minval(masses))
	tmp=tmp*10
	tmp=floor(tmp)
	step=tmp/100.
	step=0.1
	tmp=ceiling(10.*log10(minval(masses)))/10.
	open(1,file="sppinf_ticks.txt",action="write")
	do 
	tmpx=intpdfpl(10.**tmp,beta,mmin,infinity)
	tmpx=stabquant(tmpx)
	write(1,*) tmp,tmpx
	tmp=tmp+step
	if (tmp >= log10(mmax) ) exit
	end do
	close(1);
	
	!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
	!!!!!!! Plot the stabilized probability (SPP) plot
	!!!!!!! Null hypothesis (diagonal) = truncated PL
	!!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!! !!!!!!!
	open(1,file="spptrunc_data.txt",action="write")
	open(2,file="spptrunc_estimated.txt",action="write")
	open(3,file="spptrunc_infinity.txt",action="write")
	do i=1,size(masses)
		tmpx=intpdfpl(masses(i),beta,mmin,mmax)
	! Data
		tmpy=(float(i)-0.5)/float(size(masses))
		write(1,*) stabquant(tmpx),stabquant(tmpy)
	end do
	tmp=0.9999
	do
		tmpy=tmp
		tmp=invintpdfpl(tmp,beta,mmin,mmax)
	! Theory, estimated 
		tmpx=intpdfpl(tmp,beta,mmin,mmax)
		if (tmp<mmax) write(2,*) stabquant(tmpy),stabquant(tmpx)
	! Theory, infinite
		tmpx=intpdfpl(tmp,beta,mmin,infinity)
		write(3,*) stabquant(tmpy),stabquant(tmpx)
		
	! The following is for an approximately equal spacing of the plotting points	
		tmp=tmpy
		if(tmp>0.85) tmp=tmp-0.005
		if(tmp<=0.85 .AND. TMP >0.1) tmp=tmp-0.01
		if(tmp < 0.1) tmp=tmp-0.005
	if (tmp < 0.) exit
	end do
	close(1);close(2);close(3);
	! Plot the reference line (=diagonal from 0,0 to 1,1) and the 
	! acceptance region (d)
	open(1,file="spptrunc_dplus.txt",action="write")
	write(1,*) 0,dsp_trunc
	write(1,*) 1-dsp_trunc,1
	close(1)
	open(1,file="spptrunc_dminus.txt",action="write")
	write(1,*) dsp_trunc,0
	write(1,*) 1,1-dsp_trunc
	close(1)
	
	! Write some ticks and labels for the plot 
	! Calculate the masses corresponding to the theoretical quantiles
	! (used for plotting)
	tmp=log10(maxval(masses))-log10(minval(masses))
	tmp=tmp*10
	tmp=floor(tmp)
	step=tmp/100.
	step=0.1
	tmp=ceiling(10.*log10(minval(masses)))/10.
	
	open(1,file="spptrunc_ticks.txt",action="write")
	do 
	tmpx=intpdfpl(10.**tmp,beta,mmin,mmax)
	tmpx=stabquant(tmpx)
	write(1,*) tmp,tmpx
	tmp=tmp+step
	if (tmp >= log10(mmax) ) exit
	end do
	close(1);
	
	
end if


end program statpl
