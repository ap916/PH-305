! Finding inverse of a matrix using LU Decomposition 
! By solving for all columns of AX = I using LU Decomposition


program gaussPivot
    implicit none
! Initialisation of variables
    integer*8::len
    real*8,dimension(:,:),allocatable::Coeff,Lmat,Constants,Xvalues
    real*8,dimension(:),allocatable::Yvalues
    integer*8::i,j,k
    real*8::summation,lamda,temp


    print*,"Enter the no of variables (size of coeff matrix) followed by the matrix:"
    read*,len

    allocate(Coeff(len,len))
    allocate(Lmat(len,len))
    allocate(Constants(len,len))
    allocate(Xvalues(len,len))
    allocate(Yvalues(len))

! Reading values from file "input.txt"
    do i=1,len
        read (*,*)(Coeff(i,j),j=1,len)
    end do

! Printing the matrix and vector
    print *, "The matrix A  :"
    print *, "--------------"
    do i=1,len
        write (*,*)(Coeff(i,j),j=1,len)
    end do
    print *, ""


! Making L and Constants Identity Matrices
    do i=1,len
        do j=1,len
            Lmat(i,j) = 0
            Constants(i,j) = 0
            if(i.eq.j) then 
                Lmat(i,j) = 1
                Constants(i,j) = 1
            end if
        end do
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
            Lmat(i,j) = lamda
            do k=j+1,len
                Coeff(i,k) = Coeff(i,k) - lamda*Coeff(j,k)
            end do
        end do
    end do

! Printing the matrix and vector
    print *, "The matrix L, U and Constant matrix :"
    print *, "------------------------------------"
    do i=1,len
        write (*,*)(Lmat(i,j),j=1,len)
    end do
    print *, ""
    do i=1,len
        write (*,*)(Coeff(i,j),j=1,len)
    end do
    print *, ""
    do i=1,len
        write (*,*)(Constants(i,j),j=1,len)
    end do

    do k=1,len
! Getting the values of Y matrix as LY=Constants
        Yvalues(1) = Constants(1,k)/Lmat(1,1)
        do i=2,len
            summation = 0.0
            do j=1,i-1
                summation = summation + Lmat(i,j)*Yvalues(j)
            end do
            Yvalues(i) = (Constants(i,k) -summation)/Lmat(i,i)
        end do

! Getting the values of X from Y matrix as UX=Y
        Xvalues(len,k) = Yvalues(len)/Coeff(len,len)
        do i=len-1,1,-1
            summation = 0.0
            do j=len,i+1,-1
                summation = summation + Coeff(i,j)*Xvalues(j,k)
            end do
            Xvalues(i,k) = (Yvalues(i) - summation)/Coeff(i,i)
        end do
    end do


! Printing values of X
    print *, "Solutions :"
    print *, "-----------"
    do i=1,len
        write (*,*)(Xvalues(i,j),j=1,len)
    end do


end program gaussPivot
