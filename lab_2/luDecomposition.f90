! Solving a system of linear equations using LU Decomposition method
! A = LU
! LUX = B -> Taking UX = Y
! LY = B -> Solve using Forward Substitution -> Find Y
! Find X by solving UX = Y using Backward Substitution

program luDecomposition
    implicit none
! Initialisation of variables
    integer*8::len
    real*8,dimension(:,:),allocatable::Coeff,Lmat
    real*8,dimension(:),allocatable::Constants,Xvalues,Yvalues
    integer*8::i,j,k
    real*8::summation,lamda,temp,temp2

    print*,"Enter the no of variables (size of coeff matrix) followed by Coeff and Constants:"
    read*,len

    allocate(Coeff(len,len))
    allocate(Lmat(len,len))
    allocate(Constants(len))
    allocate(Xvalues(len))
    allocate(Yvalues(len))

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

! Making diagonal elements 1 in L matrix
    do i=1,len
        do j=1,len
            Lmat(i,j) = 0
            if(i.eq.j)Lmat(i,j) = 1
        end do
    end do

! Gauss Elimination with partial pivoting (The Coeff matrix becomes U after a series of row operations)
    do j=1,len-1
        if ( Coeff(j,j) < Coeff(j+1,j) .or. (Coeff(j,j).eq.0)) then
            do i = 1,len
                temp = Coeff(j,i)
                Coeff(j,i) = Coeff(j+1,i)
                Coeff(j+1,i) = temp
            end do
            temp2 = Constants(j)
            Constants(j) = Constants(j+1)
            Constants(j+1) = temp2
        end if 

        do i=j+1,len
            lamda = Coeff(i,j)/Coeff(j,j)
            Coeff(i,j) = 0.0
            Lmat(i,j) = lamda
            do k=j+1,len
                Coeff(i,k) = Coeff(i,k) - lamda*Coeff(j,k)
            end do
        end do
    end do

! Printing the matrix and vector
    print *, "The matrix L, U and vector B :"
    print *, "-----------------------------"
    do i=1,len
        write (*,*)(Lmat(i,j),j=1,len)
    end do
    print *, ""
    do i=1,len
        write (*,*)(Coeff(i,j),j=1,len)
    end do
    print *, ""
    do i=1,len
        write (*,*)Constants(i)
    end do

! Getting the values of Y matrix as LY=Constants
    Yvalues(1) = Constants(1)/Lmat(1,1)
    do i=2,len
        summation = 0.0
        do j=1,i-1
            summation = summation + Lmat(i,j)*Yvalues(j)
        end do
        Yvalues(i) = (Constants(i) -summation)/Lmat(i,i)
    end do

! Getting the values of X from Y matrix as UX=Y
    Xvalues(len) = Yvalues(len)/Coeff(len,len)
    do i=len-1,1,-1
        summation = 0.0
        do j=len,i+1,-1
            summation = summation + Coeff(i,j)*Xvalues(j)
        end do
        Xvalues(i) = (Yvalues(i) - summation)/Coeff(i,i)
    end do

! Printing values of X
    print *, "Solutions :"
    print *, "-----------"
    do i=1,len
        write (*,*)Xvalues(i)
    end do


end program luDecomposition
