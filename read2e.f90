
subroutine Read2E(Nb)

    implicit none
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
    Open(17,File='x2e.int',Form='Binary',Status='old')
    Read(17) (Ggf(I),I=1,Nb4)
!   print *, Ggf
!   Write(0,*) 'Done'

    ii = 1
    do i = 1, Nb
        do j = 1, Nb
            do k = 1, Nb
                do l = 1, Nb
                    GG_us(i, j, k, l) = Ggf(lab(lab(i,j),lab(k,l)))
                    ii = ii + 1
                end do
            end do
        end do
    end do
      
    print *, GG_us
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
    integer :: n
    call Read2E(4)
    call getarg(1, arg)
    read(arg, '(I5)') n
    write(*, *) n
end program main
