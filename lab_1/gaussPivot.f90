! Solving a system of Linear Equations using Gauss Elimination Method
! This also includes Pivoting (Swapping the rows in order to make the diagonal elements far from 0)

program gaussPivot
	implicit none
! Initialisation of variables
	integer*8,parameter::len=3
	real*8,dimension(len,len)::Coeff
	real*8,dimension(len)::Constants,Xvalues
	integer*8::i,j,k
	real*8::summation,lamda,temp

! Reading values from file "input.txt"
	do i=1,len
		read (*,*)(Coeff(i,j),j=1,len)
	end do
	do i=1,len
		read (*,*)Constants(i)
	end do

! Printing the matrix and vector
	print *, "The matrix A and vector B :"
	print *, "---------------------------"
	do i=1,len
		write (*,*)(Coeff(i,j),j=1,len)
	end do
	print *, ""
	do i=1,len
		write (*,*)Constants(i)
	end do

! Gauss Elimination with partial pivoting
	do j=1,len-1
		if ( Coeff(j,j) < Coeff(j+1,j) ) then
			do i = 1,len
				temp = Coeff(j,i)
				Coeff(j,i) = Coeff(j+1,i)
				Coeff(j+1,i) = temp
			end do
		end if 

		do i=j+1,len
			lamda = Coeff(i,j)/Coeff(j,j)
			Coeff(i,j) = 0.0
			Constants(i) = Constants(i) - lamda*Constants(j)
			do k=j+1,len
				Coeff(i,k) = Coeff(i,k) - lamda*Coeff(j,k)
			end do
		end do
	end do

! Printing the matrix and vector
	print *, "The matrix A and vector B after elimination :"
	print *, "---------------------------------------------"
	do i=1,len
		write (*,*)(Coeff(i,j),j=1,len)
	end do
	print *, ""
	do i=1,len
		write (*,*)Constants(i)
	end do

! Getting the values of X
	Xvalues(len) = Constants(len)/Coeff(len,len)
	do i=len-1,1,-1
		summation = 0.0
		do j=1,len
			summation = summation + Coeff(i,j)*Xvalues(j)
		end do
		Xvalues(i) = (Constants(i) - summation)/Coeff(i,i)
	end do

! Printing values of X
	print *, "Solutions :"
	print *, "-----------"
	do i=1,len
		write (*,*)Xvalues(i)
	end do


end program gaussPivot
