! Cholesky Decomposition Method ( Similar to LU )

program Cholesky
    implicit none
! Initialisation of variables
    integer*8::len
    real*8,dimension(:,:),allocatable::Coeff,Lmat,Lmattrans,res
    integer*8::i,j,k
    real*8::summation,lamda,temp,temp2

    print*,"Enter the no of variables (size of coeff matrix) followed by Coeff and Constants:"
    read*,len

    allocate(Coeff(len,len))
    allocate(Lmat(len,len))
    allocate(Lmattrans(len,len))
    allocate(res(len,len))

! Reading values from file "input.txt"
    do i=1,len
        read (*,*)(Coeff(i,j),j=1,len)
    end do
! Making elements of Lmat 0
    do i=1,len
        do j=1,len
            Lmat(i,j) = 0
            Lmattrans(i,j) = 0
        end do
    end do

! Printing the matrix and vector
    print *, "The matrix A (Input) :"
    print *, "---------------------------"
    do i=1,len
        write (*,*)(Coeff(i,j),j=1,len)
    end do
    print *, ""

! Cholesky's decomposition
    do i=1,len
        do j=1,i

            if(i.eq.j) then
                temp = 0.0
                do k=1,i-1
                    temp = temp + Lmat(i,k)*Lmat(i,k)
                end do
                Lmat(i,i) = sqrt(Coeff(i,i) - temp)
                Lmattrans(i,i) = Lmat(i,i)
            else 
                temp = 0.0
                do k=1,i-1
                    temp = temp + Lmat(i,k)*Lmat(j,k)
                end do
                Lmat(i,j) = (Coeff(i,j) - temp)/Lmat(j,j)
                Lmattrans(j,i) = Lmat(i,j)
            end if

        end do
    end do

    res = matmul(Lmat,Lmattrans)

! Printing the matrix and vector
    print *, "The matrix L  :"
    print *, "---------------------------"
    do i=1,len
        write (*,*)(Lmat(i,j),j=1,len)
    end do
    print *, ""
    print *, "The transpose of matrix L  :"
    print *, "---------------------------"
    do i=1,len
        write (*,*)(Lmattrans(i,j),j=1,len)
    end do
    print *, ""
    print *, "The multiplication of L(transpose) and L  :"
    print *, "---------------------------"
    do i=1,len
        write (*,*)(res(i,j),j=1,len)
    end do
    print *, ""


end program Cholesky
