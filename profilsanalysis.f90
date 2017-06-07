program pfa

implicit none

integer :: istart , iend ,di
character :: shearrate*100 , temperature*100, Z*100
logical :: there
integer :: n, i , k
real , allocatable :: shearrateD(:), temperatureD(:), ZD(:)
real :: tr

istart = 300
iend = 300
di = 100

do i = istart , iend , di

write(shearrate,'(I3.3,a)')I,'/profils/profilshearrate.txt'
write(temperature,'(I3.3,a)')I,'/profils/proftemperature.txt'
write(Z,'(I3.3,a)')I,'/profils/profZ.txt'

inquire(file=shearrate,exist=there)
if(.NOT.there) cycle

inquire(file=temperature,exist=there)
if(.NOT.there) cycle

inquire(file=Z,exist=there)
if(.NOT.there) cycle

print*,"Calcul"
!Ils ont tous la meme taille normalement

n = sizefile(shearrate)

allocate(shearrateD(n))
allocate(temperatureD(n))
allocate(ZD(n))

print*,shearrate
open(1,file=shearrate)
do k=1,n
read(1,*)tr,shearrateD(k)
enddo
close(1)

open(1,file=temperature)
do k=1,n
read(1,*)tr,temperatureD(k)
enddo
close(1)

open(1,file=Z)
do k=1,n
read(1,*)tr,ZD(k)
enddo
close(1)


!On fait la moyenne des deux branches:
!Symetrique par rapport a n/2

do k=1,n/2
print*,k,n-k-1
shearrateD(k)=0.5 * (shearrateD(k)+shearrateD(n-k-1))
enddo

do k=1,n/2
print*,k,n-k-1
ZD(k)=0.5 * (ZD(k)+ZD(n-k-1))
enddo

do k=1,n/2
print*,k,n-k-1
temperatureD(k)=0.5 * (temperatureD(k)+temperatureD(n-k-1))
enddo

open(2,file='testshearsym.txt')
do k=1,n/2
write(2,*)shearrateD(k),temperatureD(k),ZD(k)
enddo
close(2)


deallocate(shearrateD,temperatureD,ZD)


end do








contains


integer function sizefile(namefile) result(nl)
implicit none
character (len=*), intent(in) :: namefile
character :: command2*150
write(command2,'(a3,a50,a6)') 'wc ',namefile,' > bla'
call system (command2)
open(1,file='bla',status='OLD')
read(1,*)nl
close(1)
call system ("rm bla")
end function sizefile



end program
