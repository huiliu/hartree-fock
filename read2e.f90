!USAGE:
!       ./read2e basisNum x2e.int
!

subroutine Read2E(Nb, FileName)

    implicit none
    character(len=10) :: FileName
    integer,intent(in) :: NB
    integer :: NB2 
    integer :: NB4
    real*8,allocatable :: Ggf(:)
    integer :: i, j, k, l
    integer :: ii
    integer,external :: lab
    real*8 :: GG_us(NB, Nb, nb, nb)

    NB2 = (NB+1)*NB/2
    NB4 = (NB2+1)*NB2/2
    allocate(Ggf(NB4))
!   print *,"NB4",NB4

!   Write(0,*) 'Reading 2-e Integrals...'
!   Open(17,File='x2e.int',Form='Binary',Status='old')
    Open(17,File=FileName,Form='Binary',Status='old')
    Read(17) (Ggf(I),I=1,Nb4)
!   print *, Ggf
!   Write(0,*) 'Done'

    ii = 1
    do i = 1, Nb
        do j = 1, Nb
            do k = 1, Nb
                do l = 1, Nb
                    GG_us(i, j, k, l) = Ggf(lab(lab(i,j),lab(k,l)))
                     if (dabs(GG_us(i,j,k,l)) > 1.0E-8) then
                         write (*,"(I1,I1,I1,I1,F15.9)") i-1,j-1,k-1,l-1, GG_us(i,j,k,l)
                     end if
                    ii = ii + 1
                end do
            end do
        end do
    end do
      
End subroutine Read2E

Function Lab(I,J)
    Implicit Real*8 (A-H,O-Z)
    If (I.Gt.J) Then
    Lab=I*(I-1)/2+J
    Else
    Lab=J*(J-1)/2+I
    Endif
End Function

program main
    character(len=5) :: arg
    character(len=10) :: filename
    integer :: n

    call getarg(1, arg)
    call getarg(2, filename)
    read(arg, '(I5)') n

    call Read2E(n, filename)
end program main
