! Finding Eigenvalues and Eigenvectors of a matrix using Jacobi method

program main
    implicit none
    real*8,dimension(:,:),allocatable::A,R
    real*8::t,theta,c,s,maxerror,maxElement
    integer*8::len,i,j,p,q

    print*,"Enter the size of the matrix followed by the matrix:"
    read*,len

    allocate(A(len,len))
    allocate(R(len,len))

    maxerror = 1.0e-06

! Read values from user
    do i=1,len
        read(*,*)(A(i,j),j=1,len)
    end do
    
    print*, "The entered matrix is :"
    do i=1,len
        write(*,*)(A(i,j),j=1,len)
    end do
    print*,""

    call jacobi(A,R,maxerror,len)

! Printing eigenvalues
    print*, "The obtained eigenvalues :"
    write(*,*)(A(i,i),i=1,len)
    print*,""
! Printing eigenvectors
    print*, "Corresponding Eigenvectors : "
    do i=1,len
        write(*,*)(R(i,j),j=1,len)
    end do

end program main


! Subroutine for finding eigenvalues using Jacobi method
subroutine jacobi(A,R,abserror,len)
    implicit none
    integer*8::len
    real*8::abserror
    real*8,dimension(len,len)::A,R
    real*8::sqtotal,avgtotal,theta,t,c,s,cs,sc,argSign
    integer*8::i,j,k,l,iter

    ! Used in sign() function, -1 if theta<0, 1 otherwise
    argSign = 1.0 

    R=0.0
    do i=1,len
        R(i,i) = 1.0
    end do

    ! Finding sum of square of off diagonal elements
    sqtotal = 0.0
    do i=1,len
      do j=1,len
        if (i.ne.j) sqtotal = sqtotal + A(i,j)**2
      end do
    end do

    ! If Matrix is already diagonal, return
    if(sqtotal .le. abserror) return 

    avgtotal = 0.5*sqtotal/float(len*len)
    iter = 1
    do while(sqtotal .gt. abserror)
        do i=1,len-1
            do j=i+1,len
                if (A(j,i)**2 .le. avgtotal) cycle
                sqtotal = sqtotal - 2.0*A(j,i)**2
                avgtotal = 0.5*sqtotal/float(len*len)

            ! Calculating values for the givens matrix
                theta = (A(j,j)-A(i,i))/(2.0*A(i,j))
                ! t = sign(argSign,theta)/(abs(theta)+sqrt(theta**2+1.0))
                ! c = 1.0/sqrt(t**2+1.0)
                ! s = c*t
                t = 0.5*theta/sqrt(1.0+theta**2)
                s = sqrt(max(0.5+t,0.0))
                c = sqrt(max(0.5-t,0.0))

            ! Pre multiplication by K', Columns p,q are altered 
                do k=1,len
                    cs = c*A(i,k)+s*A(j,k)
                    sc = -s*A(i,k)+c*A(j,k)
                    A(i,k) = -1*sc
                    A(j,k) = cs
                end do
            ! Post multiplication by K, Only Columns i,j are altered
            ! R matrix contains the eigenvectors of the correcponding eigenvalues
                do k=1,len
                    cs =  c*A(k,i)+s*A(k,j)
                    sc = -s*A(k,i)+c*A(k,j)
                    A(k,i) = -1*sc
                    A(k,j) = cs
                    cs =  c*R(k,i)+s*R(k,j)
                    sc = -s*R(k,i)+c*R(k,j)
                    R(k,i) = -1*sc
                    R(k,j) = cs
                end do
            ! Printing Stuff     
                print *, "Iteration no : ", iter
                iter = iter + 1
                print*, "A matrix : "
                do k=1,len
                    write(*,*)(A(k,l),l=1,len)
                end do
                print*,""
                print*, "R matrix : "
                do k=1,len
                    write(*,*)(R(k,l),l=1,len)
                end do
                print*,""
                
            end do
        end do
    end do
    return
end subroutine jacobi
