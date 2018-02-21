program pf_sol_gen

!
!  written by Reed Maxwell (maxwell5@llnl.gov)
!  10-30-04
! Modified by Basile HECTOR (1st May 2017)
! This version reads a mask, identify perimeter cells (in order)
! then draws triangles for 3 patches. Last patch is the simplest:
! the perimeter: 2 triangles per 4 nodes
! top & bottom patches are similar only the definition of 
! triangles change as they should have a normal direction 
! pointing outside the domain.
! The big difference with previous version is that top&bottom
! triangles are defined accross each individual latitude band,
! and not for any single cell within the domain
! make sure you start at the correct longitude to not fall on an
! island at first...

!real, allocatable:: DEM(:,:), bottom_elev(:,:),points(:,:),&
!real, allocatable:: enter_ind(:,:), exit_ind(:,:),perim_reordered(:,:)
integer, allocatable:: points(:,:), triangles(:,:), patches(:,:), num_tri_patches(:),& 
mask(:,:), mask2(:,:), num_conti_zones(:),highdim(:), enter_ind(:,:), exit_ind(:,:),&
perim_reordered(:,:)

real datum, xmax, ymax, mask_read, dist, mdist
integer*8 num_solid, nt_y0, nt_ymax, nt_x0, nt_xmax, nt_surf, &
perimeter(60000,2),pfind(9,2), iperim,&
checkedge, punsort(60000,2), npunsort, imin, jmin, istart, jstart, &
istartsave,min_lowdim,max_lowdim,tmp1(17),tmp2(17),tmp3(17),&
single_block, enter_exist, exit_exist,diag_lat1,diag_lat2, &
num_points, num_triangles,num_patches,jj,ii,nx_dem,ny_dem,&
kk,i,j,k,num_tri_topface,num_per_pts,x0,y0,dx,dy,option_number,&
enterind, exitind, next_entereind, next_exitind, &
enterindtmp,exitindtmp,next_enterindtmp, next_exitindtmp

character*40 dem_filename, pfsol_filename, mask_filename, bottom_filename

! code creates a triangulated pf solid file
! it reads in a DEM file for the shape of the top surface, but assumes that
! the rest of the domain is rectangularly shaped (and is the same size as the
! extents of the DEM)
!
!NOTES:
! we assume the spatial discretization is the SAME for bottom and top of domain
! we assume rectangular sides
! we assume one (1) solid
! we assume six (6) patches
!
! input block, should be read in and made more general

!nx_dem = 1116   ! num dem pnts in x
!ny_dem = 660   ! num dem pnts in y
nx_dem = 5120   ! num dem pnts in x
ny_dem = 3840   ! num dem pnts in y

!latitude range for diagnostic printing
diag_lat1=10348
diag_lat2=10350
diag_lat1=1805
diag_lat2=1809

! now we need to read in the mask before we know how many triangles
! we will have.
! mask 1 = inactive
!      2 = active
allocate (mask(0:nx_dem+1,0:ny_dem+1))
allocate (mask2(0:nx_dem+1,0:ny_dem+1))

mask = 1
print*, "mask allocated"
print*, "read mask"
!mask_filename = 'mask.CONUS.txt'
mask_filename = 'mask.NA_1km.txt'
open(10, file = trim(mask_filename))
read(10,*)
do j = 1, ny_dem
do i =1, nx_dem
read(10,*) mask_read 
mask(i,j) = nint(mask_read)
!print*, mask(i,j)
!mask(i,j)=2
end do
end do
close(10)
print*, "read mask: done"

!mask2 = mask
print*, "clean up mask for thin land bands:"
do ii = 1,4
do j = 1, ny_dem
do i =1, nx_dem
checkedge = mask(i-1,j) *mask(i+1,j)
if((mask(i,j) == 2).and.(mask(i,j-1)==2).and.&
(mask(i,j+1)==2).and.(checkedge == 1)) mask(i,j) = 1
checkedge = mask(i,j-1) *mask(i,j+1)
if((mask(i,j) == 2).and.(mask(i-1,j)==2).and.&
(mask(i+1,j)==2).and.(checkedge == 1)) mask(i,j) = 1
end do
end do
end do
print*, "clean up mask for thin land bands: done"

! clean up mask for single cells
! BH: not clear: smoothing: seems to remov all cells that have 3
! empty nearby cells. Do it 5 times (why?)
print*, "clean up mask for single cells"
do ii = 1,5
do j = 1, ny_dem
do i =1, nx_dem
checkedge = mask(i-1,j) *mask(i,j-1)*mask(i+1,j)*mask(i,j+1)
if((mask(i,j) == 2).and.(checkedge <= 3)) mask(i,j) = 1

end do
end do
end do
print*, "clean up mask for single cells: done"


mask2 = mask
!BH TMP!! TO be remOVED!
!totarea=0
!do j = 1, ny_dem
!do i =1, nx_dem
!mask2(i,j)=mask2(i,j)-1
!totarea=totarea+mask2(i,j)
!end do
!end do
!print*, "total area :  ", totarea



!fill interior cells
! BH: ok same idea: set as domain cells inland empty cells
! but only for the temporary mask2 mask, doesn't alter real mask
print*, "fill interior cells"
do ii = 1, 10
do j = 1, ny_dem
do i =1, nx_dem
checkedge = mask2(i-1,j) *mask2(i,j-1)*mask2(i+1,j)*mask2(i,j+1)
if((mask2(i,j) == 1).and.(checkedge > 7)) mask2(i,j) = 2
end do
end do
end do
print*, "fill interior cells: done"



!x0 = 5000.    ! origin coord, y
!y0 = 5000.    ! origin coord, x
!dx = 10000.0   ! dx
!dy = 10000.0   ! dy
!x0 = 125.    ! origin coord, y
!y0 = 125.    ! origin coord, x
!dx = 250.0   ! dx
!dy = 250.0   ! dy
x0 = 500    ! origin coord, y
y0 = 500    ! origin coord, x
dx = 1000   ! dx
dy = 1000   ! dy
xmax = float(nx_dem)*dx + x0
ymax = float(ny_dem)*dy + y0

!! determine all perimeter cells, but not in order
open (999,file='perimeter_NA_1km.curve')
open (199,file='perimeter.ij_NA_1km.txt')


ii = 1

write(999,*) '#NA_1km'
do i = 1, nx_dem
do j = 1, ny_dem
! check if edge node
!BH: this allows no inside corner (i.e.) when checkedge=8 but mask2(i+1,j+1)=1
!for instance.
checkedge = mask2(i-1,j)+mask2(i,j-1)+mask2(i+1,j)+mask2(i,j+1)
!BH: new:
kk = mask2(i-1,j-1)+mask2(i-1,j+1)+mask2(i+1,j-1)+mask2(i+1,j+1)
if ((mask2(i,j) == 2) .and.( checkedge < 8).and.(checkedge > 4))  then
punsort(ii,1) = i
punsort(ii,2) = j
!print*, punsort(ii,1:2)
write(999,*) float(i-1)*dx+x0, float(j-1)*dy+y0
write(199,*) i,j
ii = ii + 1
end if
!BH: new:
if ((mask2(i,j) == 2) .and.( checkedge == 8).and.(kk < 8))  then
! this should help to remove these case (if i remember correctly)
!     *   * - *
!     |  /    |
! * - * /     * . . .
! |    /
! *   *
! |   |
! * - *
!
! But mostly makes corners:
! *               *
!  \              |
!   *  becomes:   * - * 
!   |                 |
!   *                 *

punsort(ii,1) = i
punsort(ii,2) = j
!print*, punsort(ii,1:2)
write(999,*) float(i-1)*dx+x0, float(j-1)*dy+y0
write(199,*) i,j
ii = ii + 1
end if

end do ! j
end do ! i
print*, "check if edge node: done"

close(999)
close(199)

npunsort = ii - 1
print*, "number of edge nodes", npunsort


!BH: Now find min max index of edges in smallest dimension (here y)
! before we alterate punsort. Could be used later, but on perimeter...!
min_lowdim=minval(punsort(1:npunsort,2),1)
max_lowdim=maxval(punsort(1:npunsort,2),1)
print*,'miny',min_lowdim
print*,'maxy',max_lowdim
allocate(num_conti_zones(ny_dem))
allocate(highdim(nx_dem))
num_conti_zones=0

!!!!!!!!!!!!!!!!!!!!! Perimeter section:
! finds perimeter cells in order (CCW)

print*, npunsort

iperim = 2

if (iperim == 2) then

 !determine perimeter nodes
! i=1, increase j until until mask = 2, then
! we sweep in a CC direction until mask goes from 1 to 2, then repeat
! ii = perimeter counter
ii = 1
!BH seems like i is the first (line?) (starting from north? South?)
! actually j = 1:ny, so j loops over lines.... i over col? yes
! i is the index for longitudes!
!i = 298 !corr BH
!i = 20000 !corr BH shoul dcheck why it was 298 ! with i=1, first j = end of col
!i = 10240 !(20480/2): to start in the middle of the X direction; not sure why
! BH: first point has to be on the mainland (no unconnected island)
i = 5000
!i = 1250
!i = 13000
do j = 1, ny_dem
!if (mask(i,j) == 2) exit
!BH:
if (mask2(i,j) == 2) exit
end do

perimeter(ii, 1) = i
perimeter(ii, 2) = j
istart = i
jstart = j

! remove the intial point

imin = perimeter(1,1)
jmin = perimeter(1,2)

do jj = 1, npunsort
if ((punsort(jj,1)==imin).and.(punsort(jj,2)==jmin)) then
punsort(jj,1) = -999999
punsort(jj,2) = -999999
istartsave = jj
exit
end if
end do ! npunsort

ii = ii + 1
! remove the second point
perimeter(ii, 1) = i+1
perimeter(ii, 2) = j

imin = perimeter(2,1)
jmin = perimeter(2,2)

do jj = 1, npunsort
if ((punsort(jj,1)==imin).and.(punsort(jj,2)==jmin)) then
punsort(jj,1) = -999999
punsort(jj,2) = -999999
exit
end if
end do ! npunsort

print*, 1, perimeter(1,1),perimeter(1,2)

!do ii = 3, 6000  !npunsort
!do ii = 3, 545651 !modif BH
!do ii = 3, 545398 ! BH: number ofedge node given above? Not sure if this is what! sweep around peremiter
!do ii = 3, 543087 ! BH: number ofedge node given above? Not sure if this is what! sweep around peremiter
!do ii = 3, 655495 ! BH: number ofedge node given above? Not sure if this is what! sweep around peremiter
do ii = 3, 57393 ! BH: number ofedge node given above? Not sure if this is what! sweep around peremiter
!do ii = 3, 40000 ! BH: number ofedge node given above? Not sure if this is what! sweep around peremiter




!go to 9090

mdist = 1000000.0
imin = 9999999
jmin = 9999999

option_number=0
!do jj = 1, 50000 
!do jj = 1, 545651 !modif BH
!do jj = 1, 545398 !modif BH
!do jj = 1, 543087 !modif BH
!do jj = 1, 655495 !modif BH
do jj = 1, 57393 !modif BH
!do jj = 1, 40000 !modif BH
! locate the closest cell

!if (jj==100) then
!punsort(istartsave,1) = istart
!punsort(istartsave,2) = jstart
!end if

!if ((punsort(jj,1).ne.perimeter(ii-1,1)).and.(punsort(jj,2).ne.perimeter(ii-1,2))) then

dist = sqrt( (float(punsort(jj,1))-float(perimeter(ii-1,1)))**2 + &
         (float(punsort(jj,2))-float(perimeter(ii-1,2)))**2 )
!print*, ii,jj, dist,punsort(jj,1), punsort(jj,2), perimeter(ii-1,1), perimeter(ii-1,2)
if (dist < mdist) then
imin = punsort(jj,1)
jmin = punsort(jj,2)
mdist = dist
end if
!end if

!if (dist <= sqrt(2.) ) then
!option_number=option_number+1
!end if

end do ! npunsort (jj?)
!if (option_number>0) then
!print*,'option number > 0: ',option_number
!end if
!option_number=0

if (mdist < 10) then ! if the closest cell has been located
!BH attempt:
!if (mdist <=9) then ! if the closest cell has been located
perimeter(ii,1) = imin
perimeter(ii,2) = jmin

do jj = 1, npunsort
if ((punsort(jj,1)==imin).and.(punsort(jj,2)==jmin)) then
punsort(jj,1) = -999999
punsort(jj,2) = -999999
end if
end do ! npunsort

!print*, ii, perimeter(ii,1),perimeter(ii,2)

!BH: commented:
9090 continue !similar to end do(?)

if ((perimeter(ii,1) == istart).and.(perimeter(ii,2) == jstart)) exit

else
perimeter(ii,1) = istart
perimeter(ii,2) = jstart
exit
end if

end do  ! ii


num_per_pts = ii

else

print*,'have to read perimieteer file'
!BH:
print*, "get into the else..."
!BH:
!num_per_pts = 5802
!num_per_pts = 545651
!num_per_pts = 545398
!num_per_pts = 543087
!num_per_pts = 655495
num_per_pts = 57393
!num_per_pts = 37996
!open (99,file='CONUS-perimeter.txt')
open (99,file='CONUS-perimeter.txt')

read(99,*)

do ii = 1, num_per_pts
!print*, ii, perimeter(ii,1), perimeter(ii,2)
read(99,*) perimeter(ii,1), perimeter(ii,2)
end do
end if

print*,'npunsort: ',npunsort,'num_per_points: ',num_per_pts




do jj = 1, npunsort
if ((punsort(jj,1)/=-999999).and.(punsort(jj,2)/=-999999 )) then
!print*, 'yes, found one'
!print*,punsort(jj,1),punsort(jj,2)
mask2(punsort(jj,1),punsort(jj,2))=1
end if
end do ! npunsort




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of perimeter section




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section starting here is made to detect, for each latitude !
! band, the number of contiguous blocks (mask =2), and store      !
! the index for entering or leaving each of these blocks.         !
! they will be further used to draw the triangles                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! EXAMPLE !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are examples for testing the fortran
! function to do this enter/exit locations:

!four different cases regarding the E-W edges
!tmp1=[1,1,2,2,2,1,1,2,2,2,2,1,2,2,2,1,1]
!tmp1=[2,2,2,2,2,1,1,2,2,2,2,1,2,2,2,1,1]
!tmp1=[2,2,2,2,2,1,1,2,2,2,2,1,2,2,2,2,2]
!tmp1=[1,1,2,2,2,1,1,2,2,2,2,1,2,2,2,2,2]
! then 3 cases for 1 single contiguous zone:
!tmp1=[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1]
!tmp1=[1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
!tmp1=[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
! case where enter ind = exit ind
tmp1=[1,1,2,2,2,1,1,1,1,2,1,1,2,2,2,1,1]

!!!! Identify number of contiguous zones and create vectors for 
!!!! index entering & exiting:
print*,'example           :',tmp1
tmp2=eoshift(tmp1,1)
print*,'eoshift(example,1):',tmp2
tmp3=tmp1-tmp2
print*,'(1) - (2)         :',tmp3
print*,'abs( (1) - (2) )  :',abs(tmp3)
print*,'sum :',sum(abs(tmp3))
!1 for the last element anyway, of eoshift
! +1 = 2 to enter first, +1 = 3 to exit first
! +1 =4 to re-enter, +1 =5 to exit again etc
! in this simple case (sum odd), number of contiguous area is (sum -1)/2 =3
! in the second case (sum even), you start by being in, total is sum/2=3
! in the 3rd case (sum even) you start by being in and exit being in: total is sum/2 = 3
! in the 4th case (sum odd) you start out and ends in: total is (sum-1)/2 =3

!!!! Show how to find 1st & 2nd enter/exit indexes (for case 1):
!in the example i is the number of contiguous areas:
i=floor(float(sum(abs(tmp3))/2))
print*, 'number of contiguous zones: ',i
!work on tmp2 instead
tmp2=tmp3
tmp2(17)=0
print*,'first entering index:'
print*,'minloc((1)-(2))   :',minloc(tmp2)
print*,'first exit index:'
print*,'maxloc((1)-(2))   :',maxloc(tmp2)
print*,'index removed:'
tmp2(minloc(tmp2))=0
tmp2(maxloc(tmp2))=0
print*,'second entering index:'
print*,'minloc((1)-(2))   :',minloc(tmp2)
print*,'second exit index:'
print*,'maxloc((1)-(2))   :',maxloc(tmp2)

!!!!!!now back on tmp3 for the real looping (still example)
!loop over each contiguous block in the other dimension:
tmp3(17)=0
!first last index:
ii=1
!populate entr & exit indices
! seems easier for the diff cases to separate entering & exit
! 'been tired, probably smarter ways to write that!
!!!!!entering first:
! if first col is domain:
if (tmp1(1)==2) then
k=0
print*,'entering 1 time: index ',k
tmp3(k+1)=0
ii=2
end if
! if only one zone and first col is not domain
if ((i==1).and.(ii==1)) then
k=minloc(tmp3,1)
print*,'entering 1 time: index ',k
! all other cases:
else
do j = ii,i
k=minloc(tmp3,1)
print*,'entering ',j,' time: index: ',k
tmp3(k)=0
end do
end if

ii=i
!!!!! exit indices:
! if only one zone and last col is domain
if ((i==1).and.(tmp1(17)==2)) then
k=17
print*,'exiting  1 time: index ',k
!if more than one zone and last col is domain
elseif ((i>=2).and.(tmp1(17)==2)) then
ii=ii-1
do j = 1,ii
k=maxloc(tmp3,1)
print*,'exiting  ',j,' time: index: ',k
tmp3(k)=0
end do
k = 17
print*,'exiting  ',i,' time: index: ',k
! other cases: one zone last col isnot domain
else
do j = 1,ii
k=maxloc(tmp3,1)
print*,'exiting  ',j,' time: index: ',k
tmp3(k)=0
end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! end EXAMPLE !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Now apply the algo to the real data:
do i=min_lowdim,max_lowdim
num_conti_zones(i)=1
highdim=mask2(1:nx_dem,i)
num_conti_zones(i)=floor(float(sum(abs(highdim-eoshift(highdim,1)))/2))
!highdim-eoshift(highdim,1)
end do
!print*,num_conti_zones
print*,'max number of contiguous zone per -lowest dimension i.e. lat or lon&
-band: ',maxval(num_conti_zones,1)

allocate(enter_ind(ny_dem,maxval(num_conti_zones,1)),&
exit_ind(ny_dem,maxval(num_conti_zones,1)))
enter_ind=0
exit_ind=0

!loop over lowest dimension within mask
do i = min_lowdim,max_lowdim
highdim=mask2(1:nx_dem,i) - eoshift(mask2(1:nx_dem,i),1)
! remove the last one, or it will be max value
highdim(nx_dem)=0

!first last index:
ii=1
!populate entr & exit indices
!!!!!entering first:
! if first col is domain:
if (mask2(1,i)==2) then
k=0
!enter_ind(i,1)=k*dx+dx/2.+x0
enter_ind(i,1)=k*dx+dx/2+x0
highdim(k+1)=0
ii=2
end if
! if only one zone and first col is not domain
if ((num_conti_zones(i)==1).and.(ii==1)) then
k=minloc(highdim,1)
!enter_ind(i,1)=k*dx+dx/2.+x0
enter_ind(i,1)=k*dx+dx/2+x0
! all other cases:
else
do j = ii,num_conti_zones(i)
k=minloc(highdim,1)
if (k == maxloc(highdim,1)-1) then
enter_ind(i,j)=0
else
!enter_ind(i,j)=k*dx+dx/2.+x0
enter_ind(i,j)=k*dx+dx/2+x0
end if
highdim(k)=0
end do
end if

!!!!! exit indices:
ii=num_conti_zones(i)
! if only one zone and last col is domain
if ((num_conti_zones(i)==1).and.(mask2(nx_dem,i)==2)) then
k=nx_dem
!exit_ind(i,1)=(k-1)*dx+dx/2.+x0
exit_ind(i,1)=(k-1)*dx+dx/2+x0
!if more than one zone and last col is domain
elseif ((num_conti_zones(i)>=2).and.(mask2(nx_dem,i)==2)) then
ii=ii-1
do j = 1,ii
if (enter_ind(i,j)==0) then
exit_ind(i,j)=0
highdim(k)=0
else
k=maxloc(highdim,1)
!exit_ind(i,j)=(k-1)*dx+dx/2.+x0
exit_ind(i,j)=(k-1)*dx+dx/2+x0
highdim(k)=0
end if
end do
exit_ind(i,num_conti_zones(i))=(nx_dem-1)*dx+dx/2+x0
k = nx_dem
! other cases: more than1 zone last col isnot domain
else
do j = 1,ii
if (enter_ind(i,j)==0) then
exit_ind(i,j)=0
else
k=maxloc(highdim,1)
!exit_ind(i,j)=(k-1)*dx+dx/2.+x0
exit_ind(i,j)=(k-1)*dx+dx/2+x0
highdim(k)=0
end if

end do
end if
end do

! check if enter_ind = exit_ind, if so, remove this case
do i = min_lowdim,max_lowdim
k=0
do j =1, num_conti_zones(i)
if ((enter_ind(i,j)==0).and.(j<=num_conti_zones(i))) then
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'enter ind(i,j): ',enter_ind(i,j),'enter ind(i,j+1): ',enter_ind(i,j+1)
print*,'exit ind(i,j) : ',exit_ind(i,j),'exit ind(i,j+1)  : ',exit_ind(i,j+1)
end if
do kk=j,num_conti_zones(i)-1
enter_ind(i,kk)=enter_ind(i,kk+1)
exit_ind(i,kk)=exit_ind(i,kk+1)
end do
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'enter ind(i,j): ',enter_ind(i,j),'enter ind(i,j+1): ',enter_ind(i,j+1)
print*,'exit ind(i,j) : ',exit_ind(i,j),'exit ind(i,j+1)  : ',exit_ind(i,j+1)
end if
num_conti_zones(i)=num_conti_zones(i)-1
print*,'latitude band: ',i,'remove block: ',j
end if
end do
end do 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 End of the first section                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!x0 = 5000.    ! origin coord, y
!y0 = 5000.    ! origin coord, x
!dx = 10000.0   ! dx
!dy = 10000.0   ! dy
!x0 = 125.    ! origin coord, y
!y0 = 125.    ! origin coord, x
!dx = 250.0   ! dx
!dy = 250.0   ! dy
x0 = 500    ! origin coord, y
y0 = 500    ! origin coord, x
dx = 1000   ! dx
dy = 1000   ! dy
xmax = float(nx_dem)*dx + x0
ymax = float(ny_dem)*dy + y0
datum = 0.0

!open (999,file='perimeter.loop_NA_250m_Vkeepsinglecells.curve')
!open (999,file='perimeter.loop_NA_250m_Vdelband2.curve')
!open (999,file='perimeter.loop_NA_250m_Vcurr.curve')
!open (999,file='perimeter.loop_NA_250m_puns.curve')
open (999,file='perimeter.loop_NA_1km.curve')
!open (999,file='perimeter.loop_NA_250m_Vtmp.curve')
write(999,*) '#NA_250m'
do ii = 1, num_per_pts
!write(999,*) float(perimeter(ii,1)-1)*dx+x0, float(perimeter(ii,2)-1)*dy+y0
!BH:
!write(999,*) float(perimeter(ii,1)-1)*dx+dx/2.+x0, float(perimeter(ii,2)-1)*dy+dy/2.+y0
write(999,*) (perimeter(ii,1)-1)*dx+dx/2+x0, (perimeter(ii,2)-1)*dy+dy/2+y0
end do
close(999)

!
num_solid = 1

! now we count up surface, L, R, F, B number of triangles and points
! to make the indexing easier, we dump out all points, so the 
! num points does not depend on the mask
! calculate total number of points
!num_points =  2*nx_dem*ny_dem 
!BH only perimeter:
num_points = 2 * (num_per_pts-1)
!Because perimeter(1,1)=perimeter(end,1)
! currently we are gridding 2 triangles per box of four nodes, top and bottom

!BH: mindim_range is the number of rectangles elongated along the largest dim
!BH: perimeter + topx2 (top + bottom)
!BH: the following is a way to calculate the number of triangles fot
!non-too-fancy geometries. just look at how many contiguous zones per latitude
!bands.
!ii=0
!do i = min_lowdim,max_lowdim - 1
!ii=ii+num_conti_zones(i)
!end do
!num_tri_topface=2*ii
!num_triangles = 2 * (num_per_pts-1) + 2*num_tri_topface
!BH: However, sometimes it's not that simple, and one contiguous zone in
!latitude band i matches more than one continguous zones in latitude band i+1.
! in those specific cases, you should move up to latitude band i+1, and look
! back south toward latitude band i, to define the triangles that link i & i+1
! So further we will quickly loop over the whole domain to allocate space for
! triangles.

! currently hard wired for 6 patches
! we also allocate the max size for patches- num_triangles
!num_patches = 6
num_patches = 3

!allocate(DEM(nx_dem,ny_dem),points(num_points,3),triangles(num_triangles,3),  &
!patches(num_patches,num_triangles),num_tri_patches(num_triangles),bottom_elev(nx_dem,ny_dem), &
!perim_reordered(num_per_pts-1,4))
!allocate(DEM(nx_dem,ny_dem),points(num_points,3), &
!bottom_elev(nx_dem,ny_dem), &
!perim_reordered(num_per_pts-1,4))
allocate(points(num_points,3),perim_reordered(num_per_pts-1,4))

!bottom_elev = 0.  ! elevation of the bottom of the domain, init to 0.
!dem = 1000.  ! elevation of the top of the domain, init to 0.

dem_filename = 'demcol.txt'
!pfsol_filename = 'CONUS_new_edge_test.pfsol'
!pfsol_filename = 'CONUS_new_edge_test_lighter_V2.pfsol'
!pfsol_filename = 'CONUS_new_edge_test_lighter_V3.pfsol'
!pfsol_filename = 'CONUS_new_edge_test_lighter_V4.pfsol'
!pfsol_filename = 'NA_250m_new_edge_test_lighter_Vkeepsinglecells.pfsol'
!pfsol_filename = 'NA_250m_lighter_Vdlband2.pfsol'
!pfsol_filename = 'NA_250m_lighter_Vcurr.pfsol'
!pfsol_filename = 'NA_250m_lighter_puns.pfsol'
pfsol_filename = 'NA_1km.pfsol'
!pfsol_filename = 'NA_250m_lighter_Vtmp.pfsol'

print*, 'flag 1'


!bottom_elev(i,j) = 0.0


!constant bottom, dem can be read from file
!do j = 1, ny_dem
!do i =1, nx_dem

 !bottom_elev(i,j) = 0.0
!end do
!end do


! loop over points and assign, jj = point counter
! first the DEM/upper surface
! we assume points are in the center of the digital elevation tile 

!jj = 1
!do i = 1, nx_dem
!do j = 1, ny_dem
!points(jj,1) = float(i-1)*dx + dx/2. + x0
!points(jj,2) = float(j-1)*dy + dy/2. + y0
!points(jj,3) = dem(i,j)
!jj = jj + 1
!end do
!end do

! then the lower surface; currently a flat surface - i.e. datum
!do i = 1, nx_dem
!do j = 1, ny_dem
!points(jj,1) = float(i-1)*dx + dx/2. + x0
!points(jj,2) = float(j-1)*dy + dy/2. + y0
!points(jj,3) = bottom_elev(i,j)
!jj = jj + 1
!end do
!end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! Fill up the vertice matrix (points)

points=0.0
!BH:
jj=1
!stop one before because perimeter (1,1) = perimeter(end,1)
do kk = 1, num_per_pts-1
!points(jj,1)=float(perimeter(kk,1)-1)*dx+x0
!points(jj,2)=float(perimeter(kk,2)-1)*dy+y0
!points(jj,1)=float(perimeter(kk,1)-1)*dx+dx/2.+x0
!points(jj,2)=float(perimeter(kk,2)-1)*dy+dy/2.+y0
!points(jj,3)=1000.0
points(jj,1)=(perimeter(kk,1)-1)*dx+dx/2+x0
points(jj,2)=(perimeter(kk,2)-1)*dy+dy/2+y0
points(jj,3)=1000
jj = jj + 1
end do
do kk = 1, num_per_pts-1
!points(jj,1)=float(perimeter(kk,1)-1)*dx+x0
!points(jj,2)=float(perimeter(kk,2)-1)*dy+y0
!points(jj,1)=float(perimeter(kk,1)-1)*dx+dx/2.+x0
!points(jj,2)=float(perimeter(kk,2)-1)*dy+dy/2.+y0
!points(jj,3)=0.0
points(jj,1)=(perimeter(kk,1)-1)*dx+dx/2+x0
points(jj,2)=(perimeter(kk,2)-1)*dy+dy/2+y0
points(jj,3)=0
jj = jj+ 1
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do kk = 1, num_per_pts-1
!perim_reordered(kk,1)=float(perimeter(kk,2)-1)*dy+y0
!perim_reordered(kk,2)=float(perimeter(kk,1)-1)*dx+x0
!perim_reordered(kk,1)=float(perimeter(kk,2)-1)*dy+dy/2.+y0
!perim_reordered(kk,2)=float(perimeter(kk,1)-1)*dx+dx/2.+x0
!perim_reordered(kk,3)= kk
!perim_reordered(kk,4)=perimeter(kk,2)
perim_reordered(kk,1)=(perimeter(kk,2)-1)*dy+dy/2+y0
perim_reordered(kk,2)=(perimeter(kk,1)-1)*dx+dx/2+x0
perim_reordered(kk,3)= kk
perim_reordered(kk,4)=perimeter(kk,2)
end do

!test: check in the case where exit ind is not found on perimeter and enter ind
!of next block neither: then remove that block;
do i = min_lowdim,max_lowdim
k=0
do j =1, num_conti_zones(i)

!!!! locate perimeter cell for enter & exit index
enter_exist=0
exit_exist=0
do k = 1,num_per_pts-1
if ((perim_reordered(k,4)==i).and.(perim_reordered(k,2)==enter_ind(i,j))) then
enterind=perim_reordered(k,3)
enter_exist=1
end if
if ((perim_reordered(k,4)==i).and.(perim_reordered(k,2)==exit_ind(i,j))) then
exitind=perim_reordered(k,3)
exit_exist=1
end if
end do

if ((enter_exist==1).and.(exit_exist==0).and.(j<num_conti_zones(i))) then

!identify the next matching exit:
do jj=j+1,num_conti_zones(i)

do k = 1,num_per_pts-1
if ((perim_reordered(k,4)==i).and.(perim_reordered(k,2)==exit_ind(i,jj))) then
exitind=perim_reordered(k,3)
exit_exist=1
end if
end do
if (exit_exist==1) exit
end do
!jj becomes the index location of the next matching exit
!say current j=5 , num_conti_zones(i) = 13 and jj=8
!=> from jj=8 to 13 all index (enter and exit) have to be shifted to
! 5 to 10 and new num conti zones(i)=10

exit_ind(i,j)=exit_ind(i,jj)
!exit_ind(i,j)=exit_ind(i,j+1)

if (j<num_conti_zones(i)-1) then
!do kk=j,num_conti_zones(i)-1
do kk=j,j+num_conti_zones(i)-jj
!enter_ind(i,kk+1)=enter_ind(i,kk+2)
!exit_ind(i,kk+1)=exit_ind(i,kk+2)
enter_ind(i,kk+1)=enter_ind(i,kk+(jj-j)+1)
exit_ind(i,kk+1)=exit_ind(i,kk+(jj-j)+1)
end do
end if
num_conti_zones(i)=num_conti_zones(i)-(jj-j)
!num_conti_zones(i)=num_conti_zones(i)-1

print*,'latitude band: ',i,'remove blocks (bcse not matching exit): ',j+1,'to',jj
end if
end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Here: loop over all the domain to count how many triangles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! to allocate space
num_tri_topface=0

!!!! Loop1:  biggest loop: look at each latitude band:
do i =min_lowdim,max_lowdim-1
!print*,i
exit_exist=1
!!!! Loop2: for each latitude, loop over contiguous zones:
if (i==perimeter(1,2)) then
print*,i
print*,num_conti_zones(i)
print*,num_conti_zones(i+1)
end if
if (i==perimeter(1,2)-1) then
print*,i
print*,num_conti_zones(i)
print*,num_conti_zones(i+1)
end if
!if (i>perimeter(1,2)) then
!print*,i
!end if

do j=1,num_conti_zones(i)

!!!! locate perimeter cell for enter & exit index
enter_exist=0
exit_exist=0
do k = 1,num_per_pts-1
if ((perim_reordered(k,4)==i).and.(perim_reordered(k,2)==enter_ind(i,j))) then
enterind=perim_reordered(k,3)
enter_exist=1
end if
if ((perim_reordered(k,4)==i).and.(perim_reordered(k,2)==exit_ind(i,j))) then
exitind=perim_reordered(k,3)
enter_exist=1
exit_exist=1
end if
end do

!if (i==perimeter(1,2)) then
!print*,'block number: ',j
!print*,'enter exist? :', enter_exist
!print*,' exit exist? :', exit_exist
!end if
!if (i==perimeter(1,2)-1) then
!print*,'block number: ',j
!print*,'enter exist? :', enter_exist
!print*,' exit exist? :', exit_exist
!end if

!!!! IF 1: check that enter & exit index are actually found on the perimeter,
!because perimeter is defined after enter & exit (above) using RM's code section
if ((enter_exist==1).and.(exit_exist==1)) then

! Identify the actual triangle vertice for enter & exit.
! for instance, if enter is like
!     *
! * * *
! ^   ^
! then dont take the first arrow (the identified enter index
! but instead take the last one before latitude change along the perimeter:

! this should be corrected close to the end of the perimeter list...
do k = 1,enterind
if (perim_reordered(enterind-k,4)>i) then
enterind=enterind-k+1
EXIT
end if
!on the opposite, if we get lower in latitude, its not good, it means we are in
! 
!* * *                     * * *
!    * *       * * * *     *   * 
!      * * * * *     * * * *   *
!              ^ or  ^         *
if (perim_reordered(enterind-k,4)<i) then
! in this case, exit
EXIT
end if
end do

do k = 1,exitind
if (perim_reordered(exitind+k,4)>i) then
exitind=exitind+k-1
EXIT
end if
if (perim_reordered(exitind+k,4)<i) then
EXIT
end if
end do

!if (i==perimeter(1,2)) then
!print*,'exitind: ',exitind
!print*,'enterind:',enterind
!end if
!if ((i>diag_lat1).and.(i<diag_lat2)) then
!if (i==perimeter(1,2)) then
!print*,'enterind: ',enterind,'lon: ',&
!(perimeter(enterind,1)-1)*dx + dx/2. +x0,'lat: ',&
!(perimeter(enterind,2)-1)*dy + dy/2. +y0
!print*,'exitind:  ',exitind,'lon: ',&
!(perimeter(exitind,1)-1)*dx + dx/2. +x0,'lat: ',&
!(perimeter(exitind,2)-1)*dy + dy/2. +y0
!end if

!if (i==perimeter(1,2)-1) then
!print*,'enterind: ',enterind,'lon: ',&
!(perimeter(enterind,1)-1)*dx + dx/2. +x0,'lat: ',&
!(perimeter(enterind,2)-1)*dy + dy/2. +y0
!print*,'exitind:  ',exitind,'lon: ',&
!(perimeter(exitind,1)-1)*dx + dx/2. +x0,'lat: ',&
!(perimeter(exitind,2)-1)*dy + dy/2. +y0
!end if


! here it's about the same idea when looking for next 
! enter /exit index. In the case mentioned above, set ii to 1
! and this will mean that no triangle should be drawn there
ii=0
do k = 1,enterind
if (perim_reordered(enterind-k,4)>i) then
next_enterind=enterind-k
EXIT
end if
if (perim_reordered(enterind-k,4)<i) then
ii=1
EXIT
end if
end do
do k = 1,exitind
if (perim_reordered(exitind+k,4)>i) then
next_exitind=exitind+k
EXIT
end if
if (perim_reordered(exitind+k,4)<i) then
ii=1
EXIT
end if
end do


!if (i==perimeter(1,2)) then
!print*,'next exitind: ',next_exitind
!print*,'next enterind:',next_enterind
!print*,'nextenter ',next_enterind,'lon: ',&
!(perimeter(next_enterind,1)-1)*dx + dx/2. +x0,'lat: ',&
!(perimeter(next_enterind,2)-1)*dy + dy/2. +y0
!print*,'nextexit: ',next_exitind,'lon: ',&
!(perimeter(next_exitind,1)-1)*dx + dx/2. +x0,'lat: ',&
!(perimeter(next_exitind,2)-1)*dy + dy/2. +y0
!end if

!if (i==perimeter(1,2)-1) then
!print*,'next exitind: ',next_exitind
!print*,'next enterind:',next_enterind
!print*,'nextenter ',next_enterind,'lon: ',&
!(perimeter(next_enterind,1)-1)*dx + dx/2. +x0,'lat: ',&
!(perimeter(next_enterind,2)-1)*dy + dy/2. +y0
!print*,'nextexit: ',next_exitind,'lon: ',&
!(perimeter(next_exitind,1)-1)*dx + dx/2. +x0,'lat: ',&
!(perimeter(next_exitind,2)-1)*dy + dy/2. +y0
!end if

!particular cases for first and ast perieter index
if (exitind==num_per_pts-1) then
next_exitind = 1
end if
if (enterind==1) then
next_enterind=num_per_pts-1
end if

if (i==perimeter(1,2)) then
print*,'next exitind: ',next_exitind
print*,'next enterind:',next_enterind
end if
single_block=0
!Check if these index (or any previous enter index/next exit index)
!correspond to one single block in the next zone (latitude i+1)
if (i==perimeter(1,2)) then
print*,'enterindtmp: ',enterindtmp
print*,'exitindtmp:  ',exitindtmp
print*,'enterindtmp: ',enterindtmp,'lon: ',&
(perimeter(enterindtmp,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(enterindtmp,2)-1)*dy + dy/2. +y0
print*,'exitindtmp:  ',exitindtmp,'lon: ',&
(perimeter(exitindtmp,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(exitindtmp,2)-1)*dy + dy/2. +y0
end if

do kk = 1,num_conti_zones(i+1)
do k = 1,num_per_pts-1
if ((perim_reordered(k,4)==i+1).and.(perim_reordered(k,2)==enter_ind(i+1,kk))) then
enterindtmp=perim_reordered(k,3)
end if
if ((perim_reordered(k,4)==i+1).and.(perim_reordered(k,2)==exit_ind(i+1,kk))) then
exitindtmp=perim_reordered(k,3)
end if
end do

if (i==perimeter(1,2)-1) then
print*,'enterindtmp: ',enterindtmp
print*,'exitindtmp:  ',exitindtmp
print*,'enterindtmp: ',enterindtmp,'lon: ',&
(perimeter(enterindtmp,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(enterindtmp,2)-1)*dy + dy/2. +y0
print*,'exitindtmp:  ',exitindtmp,'lon: ',&
(perimeter(exitindtmp,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(exitindtmp,2)-1)*dy + dy/2. +y0
end if


if ((perim_reordered(enterindtmp,2)<=perim_reordered(next_enterind,2)).and.&
(perim_reordered(exitindtmp,2)>=perim_reordered(next_exitind,2)).and.&
(ii==0)) then
single_block=1
end if
end do

if (i==perimeter(1,2)) then
print*,'single block? ',single_block
end if
!!!! IF 2: 
! if current block correspond to one single block in the next zone
if (single_block==1) then
num_tri_topface=num_tri_topface+2

else !if single_block==0
do kk = 1, num_conti_zones(i+1)
if ((exit_ind(i+1,kk)>=enter_ind(i,j)).and.&
(enter_ind(i+1,kk)<=exit_ind(i,j))) then

num_tri_topface=num_tri_topface+2

end if! end if the zone of next band is inked to current band
end do
end if!  END IF2 (single block)
end if ! END IF1
end do ! END Loop2
end do ! END Loop1

num_triangles = 2 * (num_per_pts-1) + 2*num_tri_topface

!BH: compare actual number of triangles vs the simple case:
ii=0
do k = min_lowdim,max_lowdim - 1
ii=ii+num_conti_zones(k)
end do
print*, 'number tri top old: ',2*ii
print*, 'number tri top new: ',num_tri_topface



!allocate(triangles(num_triangles,3),  &
!patches(num_patches,num_triangles),num_tri_patches(num_triangles))
allocate(triangles(num_triangles,3),  &
patches(num_patches,num_triangles),num_tri_patches(num_patches))
!num_tri_topface triangles for top patch:
num_tri_patches(1)=num_tri_topface
!num_tri_topface triangles for bottom patch:
num_tri_patches(2)=num_tri_topface
!2*(num_perimeter _pts -1) triangles for perimeter patch:
num_tri_patches(3)=2*(num_per_pts-1)

! end triangle allocation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!print*, 'points: ',kk - 1, num_points


! loop over triangles and assign, jj = triangle counter; we also do the patches at the same time
! note: point looping starts at zero but we will start at 1, then sub later
! first the DEM/upper surface
! notes:
! patch(1,:) = upper surface (z = zmax)
! patch(2,:) = lower surface (z = z0)
! patch(3,:) = x = x0
! patch(4,:) = x = xmax
! patch(5,:) = y = y0
! patch(6,:) = y = ymax
! jj = triangle counter
! ii = patch counter; ii restarts w/ every face
!
! NOTE: the order in which the triangle points verticies are specified is critical -
! they must be specified so the normal vector points out of the domain always (ie RHR
!  counterclockwise for the top face, etc)
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Top & Bottom surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BH: top surface first
! recall:
! min_lowdim, max_lowdim: min latitude zone, max latitude zone
! num_conti_zones(ny_dem): number of contiguous zones per latitude band
! enter_ind(ny_dem,max_numberofcontizone): enter index for each conti zone
! exit_ind, same for exit index

jj=1 !(triangle counter)

!!!! Loop1:  biggest loop: look at each latitude band:
do i =min_lowdim,max_lowdim-1

! print for checking:
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'process latitude band: ',i
print*,'number of contiguous zones: ',num_conti_zones(i)
end if

!!!! Loop2: for each latitude, loop over contiguous zones:
do j=1,num_conti_zones(i)

!print for checking:
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'zone: ',j
print*,'enterind: ',enter_ind(i,j),'exitind: ',exit_ind(i,j)
end if

enter_exist=0
exit_exist=0

!!!! locate perimeter cell for enter & exit index
do k = 1,num_per_pts-1
if ((perim_reordered(k,4)==i).and.(perim_reordered(k,2)==enter_ind(i,j))) then
enterind=perim_reordered(k,3)
!print for checking:
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'yes,found one'
end if
enter_exist=1
end if

if ((perim_reordered(k,4)==i).and.(perim_reordered(k,2)==exit_ind(i,j))) then
exitind=perim_reordered(k,3)
!print for checking:
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'yes,found one'
end if
enter_exist=1
exit_exist=1
end if
end do

!!!! IF 1: check that enter & exit index are actually found on the perimeter,
!because perimeter is defined after enter & exit (above) using RM's code section
! could be inconsistencies
if ((enter_exist==1).and.(exit_exist==1)) then

!print for chcking
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'block enter & exit match perimeter cells: ok'
end if

! Identify the actual triangle vertice for enter & exit.
! for instance, if enter is like
!     *
! * * *
! ^   ^
! then dont take the first arrow (the identified enter index
! but instead take the last one before latitude change along the perimeter:

! this should be corrected close to the end of the perimeter list...
do k = 1,enterind
if (perim_reordered(enterind-k,4)>i) then
enterind=enterind-k+1
EXIT
end if
!on the opposite, if we get lower in latitude, its not good, it means we are in
! 
!* * *                     * * *
!    * *       * * * *     *   * 
!      * * * * *     * * * *   *
!              ^ or  ^         *
if (perim_reordered(enterind-k,4)<i) then
! in this case, exit
EXIT
end if
end do
do k = 1,exitind
if (perim_reordered(exitind+k,4)>i) then
exitind=exitind+k-1
EXIT
end if
if (perim_reordered(exitind+k,4)<i) then
EXIT
end if
end do

!particular cases for first and ast perieter index
if (exitind==num_per_pts-1) then
next_exitind = 1
end if
if (enterind==1) then
next_enterind=num_per_pts-1
end if

! Identify the triangle vertice of the i+1 latitude band.
! could simply be:
!next_enterind = enterind - 1
!next_exitind = exitind + 1
! ... that is, the previous perimeter node with respect to the enter index and
! the next perimeter node with respect to the 'exit' index.
! However, it's better to take instead the first of the previous perimeter node
! that actually belong to latitude band i+1 (otherwise risk to stay on same lat)

! this should be corrected close to the end of the perimeter list...

! here it's about the same idea when looking for next 
! enter /exit index. In the case mentioned above (when <i), set ii to 1
! and this will mean that no triangle should be drawn there
ii=0
do k = 1,enterind
if (perim_reordered(enterind-k,4)>i) then
next_enterind=enterind-k
EXIT
end if
if (perim_reordered(enterind-k,4)<i) then
ii=1
EXIT
end if
end do
do k = 1,exitind
if (perim_reordered(exitind+k,4)>i) then
next_exitind=exitind+k
EXIT
end if
if (perim_reordered(exitind+k,4)<i) then
ii=1
EXIT
end if
end do

!particular cases for first and ast perieter index
if (exitind==num_per_pts-1) then
next_exitind = 1
end if
if (enterind==1) then
next_enterind=num_per_pts-1
end if


single_block=0
!Check if these index (or any previous enter index/next exit index)
!correspond to one single block in the next zone

do kk = 1,num_conti_zones(i+1)
do k = 1,num_per_pts-1
if ((perim_reordered(k,4)==i+1).and.(perim_reordered(k,2)==enter_ind(i+1,kk))) then
enterindtmp=perim_reordered(k,3)
end if
if ((perim_reordered(k,4)==i+1).and.(perim_reordered(k,2)==exit_ind(i+1,kk))) then
exitindtmp=perim_reordered(k,3)
end if
end do

!if ((enterindtmp<=next_enterind).and.(exitindiiiiiiiiiiiiiitmp>=next_exitind)) then
if ((perim_reordered(enterindtmp,2)<=perim_reordered(next_enterind,2)).and.&
(perim_reordered(exitindtmp,2)>=perim_reordered(next_exitind,2)).and.&
(ii==0)) then
! if any of the block of latitude band i+1 has an enter index (x -longitude-
! value) lower or equal to expected enter index from the current block in lat
! band i AND an exit index higher or equal to the expected exit index from
! currnt block in latitude band i, then the two blocks at i and i+1 match
! directly. Otherwise, should move up to lat band i+1 and lood down south...!

if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'single block? Yes...'
end if
single_block=1
end if
end do
if (single_block==0) then
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'single block? No...'
end if
end if

! print for checking:
!if (i<4+min_lowdim) then
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'enterind: ',enterind,'lon: ',&
(perimeter(enterind,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(enterind,2)-1)*dy + dy/2. +y0
print*,'exitind:  ',exitind,'lon: ',&
(perimeter(exitind,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(exitind,2)-1)*dy + dy/2. +y0
print*,'nextenter ',next_enterind,'lon: ',&
(perimeter(next_enterind,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(next_enterind,2)-1)*dy + dy/2. +y0
print*,'nextexit: ',next_exitind,'lon: ',&
(perimeter(next_exitind,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(next_exitind,2)-1)*dy + dy/2. +y0
end if


!!!! IF 2: 
! if current block correspond to one single block in the next zone
if (single_block==1) then

! draw 2 triangles per block(W/E), for top & bot:
! upper West (top) (order 1-2-3, CCW for normal going out):
! 1- 3  -> latitude band i + 1, pts 1 & 3 are previous/next pts of 2/4
! | /
! 2/    -> current latitude band i, point '2' is the enter index for this
! current block
triangles(jj,1)=next_enterind
triangles(jj,2)=enterind
triangles(jj,3)=next_exitind
patches(1,jj)=jj
! upper West (bottom) (order 1-2-3, CW for normal going out below):
! 2- 3
! | /
! 1/
! 
triangles(jj+num_tri_topface,1)=enterind+num_per_pts-1
triangles(jj+num_tri_topface,2)=next_enterind+num_per_pts-1
triangles(jj+num_tri_topface,3)=next_exitind+num_per_pts-1
patches(2,jj)=jj+num_tri_topface
jj=jj+1

! lower East (top):
!   /3
!  / |
! 1- 2
triangles(jj,1)=enterind
triangles(jj,2)=exitind
triangles(jj,3)=next_exitind
patches(1,jj)=jj
! lower East (bottom):
!   /2
!  / |
! 1- 3
triangles(jj+num_tri_topface,1)=enterind+num_per_pts-1
triangles(jj+num_tri_topface,2)=next_exitind+num_per_pts-1
triangles(jj+num_tri_topface,3)=exitind+num_per_pts-1
patches(2,jj)=jj+num_tri_topface
jj=jj+1

else !if single_block==0
! in this case should look from the point of view of next latitude band:
! but just for this current zone...

! Identify which zones of the next band to take into account
! all thse having an exit > to the enter of current zone and
! and enter < to the exit of current zone
do kk = 1, num_conti_zones(i+1)
if ((exit_ind(i+1,kk)>=enter_ind(i,j)).and.&
(enter_ind(i+1,kk)<=exit_ind(i,j))) then
! then if next_enterind previously calculated from the current zone
! is > enter ind, should take this one (next_enterind) instead
! same for exit if next_exitind< exit ind should take this one instead.
! and should always identify the new 'next' perimeter node looking toward 
! the south now...

!locate perimeter  cell: (loop over num_per_points instead, as RM?)
do k = 1,num_per_pts-1
if ((perim_reordered(k,4)==i+1).and.(perim_reordered(k,2)==enter_ind(i+1,kk))) then
enterindtmp=perim_reordered(k,3)
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'yes,found one'
end if
end if

if ((perim_reordered(k,4)==i+1).and.(perim_reordered(k,2)==exit_ind(i+1,kk))) then
exitindtmp=perim_reordered(k,3)
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'yes,found one'
end if
end if
end do

do k = 1,enterindtmp
if (perim_reordered(enterindtmp+k,4)==i) then
enterindtmp=enterindtmp+k-1
EXIT
end if
end do
do k = 1,exitindtmp
if (perim_reordered(exitindtmp-k,4)==i) then
exitindtmp=exitindtmp-k+1
EXIT
end if
end do

!looking toward the south, so looking the other way round now (check consistency
!with upper part of the domain!!!)
!next_enterindtmp = enterindtmp + 1
!next_exitindtmp = exitindtmp - 1

do k = 1,enterindtmp
if (perim_reordered(enterindtmp+k,4)==i) then
next_enterindtmp=enterindtmp+k
EXIT
end if
end do
do k = 1,exitindtmp
if (perim_reordered(exitindtmp-k,4)==i) then
next_exitindtmp=exitindtmp-k
EXIT
end if
end do

!particular cases for first and ast perieter index
if (enterindtmp==num_per_pts-1) then
next_enterindtmp = 1
end if
if (exitindtmp==1) then
next_exitindtmp=num_per_pts-1
end if

! what is this for?
!if (enterindtmp<=next_enterind) then
!enterindtmp=next_enterind
!end if
!if (exitindtmp>=next_exitind) then
!exitindtmp=next_exitind
!end if

! print for checking:
!if (i<4+min_lowdim) then
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'enterind: ',enterindtmp,'lon: ',&
(perimeter(enterindtmp,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(enterindtmp,2)-1)*dy + dy/2. +y0
print*,'exitind:  ',exitindtmp,'lon: ',&
(perimeter(exitindtmp,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(exitindtmp,2)-1)*dy + dy/2. +y0
print*,'nextenter ',next_enterindtmp,'lon: ',&
(perimeter(next_enterindtmp,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(next_enterindtmp,2)-1)*dy + dy/2. +y0
print*,'nextexit: ',next_exitindtmp,'lon: ',&
(perimeter(next_exitindtmp,1)-1)*dx + dx/2. +x0,'lat: ',&
(perimeter(next_exitindtmp,2)-1)*dy + dy/2. +y0
end if


! draw 2 triangles per block(W/E), for top & bot:
! upper West (top):
! 1- 3
! | /
! 2/
triangles(jj,1)=enterindtmp
triangles(jj,2)=next_enterindtmp
triangles(jj,3)=exitindtmp
patches(1,jj)=jj
! upper West (bottom):
! 2- 3
! | /
! 1/
triangles(jj+num_tri_topface,1)=next_enterindtmp+num_per_pts-1
triangles(jj+num_tri_topface,2)=enterindtmp+num_per_pts-1
triangles(jj+num_tri_topface,3)=exitindtmp+num_per_pts-1
patches(2,jj)=jj+num_tri_topface
jj=jj+1

! lower East (top):
!   /3
!  / |
! 1- 2
triangles(jj,1)=next_enterindtmp
triangles(jj,2)=next_exitindtmp
triangles(jj,3)=exitindtmp
patches(1,jj)=jj
! lower East (bottom):
!   /2
!  / |
! 1- 3
triangles(jj+num_tri_topface,1)=next_enterindtmp+num_per_pts-1
triangles(jj+num_tri_topface,2)=exitindtmp+num_per_pts-1
triangles(jj+num_tri_topface,3)=next_exitindtmp+num_per_pts-1
patches(2,jj)=jj+num_tri_topface
jj=jj+1

!
end if! end if the zone of next band is inked to current band

end do
end if! END IF 2 ( single block)

!print for checking:
!if (i<4+min_lowdim) then
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'triangles Upper West (top):'
print*,triangles(jj-2,1),triangles(jj-2,2),triangles(jj-2,3)
print*,'triangles Upper West (bottom):'
print*,triangles(jj-2+num_tri_topface,1),triangles(jj-2+num_tri_topface,2),&
triangles(jj-2+num_tri_topface,3)
print*,'triangles Lower East (top):'
print*,triangles(jj-1,1),triangles(jj-1,2),triangles(jj-1,3)
print*,'triangles Lower East (bottom):'
print*,triangles(jj+num_tri_topface-1,1),triangles(jj+num_tri_topface-1,2),&
triangles(jj+num_tri_topface-1,3)
end if

else
if ((i>diag_lat1).and.(i<diag_lat2)) then
print*,'block enter & exit dont match perimeter cells: aie'
end if ! END IF 2 single block

end if ! END IF 1 
end do ! END LOOP 2
end do ! END LOOP 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Perimeter surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!loop around the perimeter
jj=2*num_tri_topface+1
ii=1
do kk = 1, num_per_pts-2
!upper triangle
! first points stored are higher elevation ones:
! 1- 3
!  \ |
!   \2
triangles(jj,1)=kk
triangles(jj,2)=kk+num_per_pts-1+1
triangles(jj,3)=kk+1
patches(3,ii)=jj
jj=jj+1
ii=ii+1
!lower triangle
! 1
! | \ 
! 2- 3
triangles(jj,1)=kk
triangles(jj,2)=kk+num_per_pts-1
triangles(jj,3)=kk+num_per_pts-1+1
patches(3,ii)=jj
jj=jj+1
ii=ii+1
end do
!last triangles:
!upper triangle
triangles(jj,1)=num_per_pts-1
triangles(jj,2)=num_per_pts
triangles(jj,3)=1
patches(3,ii)=jj
jj=jj+1
ii=ii+1
!lower triangle
triangles(jj,1)=num_per_pts-1
triangles(jj,2)=num_per_pts-1+num_per_pts-1
triangles(jj,3)=num_per_pts
patches(3,ii)=jj

!num_tri_patches(3) = ii - 1
ii = 1


!print*, 'triangles: ',jj - 1, num_triangles
!num_triangles = jj -1
!BH:
!print*, 'triangles: ',jj - 1, num_triangles
!num_triangles = jj -1

open(21,file= trim(pfsol_filename))
! write version
write(21,'(i1)') 1
! write num vertices/points
!write(20,'(i8)') num_points
!BH:
!write(21,'(i8)') (num_per_pts-1)*2
write(21,'(i9)') (num_per_pts-1)*2
! write points
!do i = 1, num_points
!BH
do i = 1, (num_per_pts-1)*2
!write(21,'(3(f15.4,2x))') points(i,1), points(i,2), points(i,3)
write(21,'(3(f15.4,2x))') float(points(i,1)), float(points(i,2)), float(points(i,3))
end do
! currently we are assuming 1 solid 
write(21,'(i4)') num_solid
do k = 1, num_solid

! write num trianglesw
!write(21,'(i8)') num_triangles
write(21,'(i9)') num_triangles
! write triangles
do i = 1, num_triangles
!BH:
! IMPORTANT: in the .pfsol, first point is number 0, and called as such within
! the triangles, that's why we subtract 1 here (same in patches)
!write(21,'(3(i8,2x))') triangles(i,1) -1, triangles(i,2) -1, triangles(i,3) -1
write(21,'(3(i9,2x))') triangles(i,1) -1, triangles(i,2) -1, triangles(i,3) -1
end do ! num_triangles

!write number of patches
write(21,'(i3)') num_patches
do i = 1, num_patches
! write num triangles in each patch
!write(21,'(1x,i8)') num_tri_patches(i)
write(21,'(1x,i9)') num_tri_patches(i)
do j = 1, num_tri_patches(i)
!write(21,'(i8)') patches(i,j) - 1
!BH
!write(21,'(i8)') patches(i,j) -1
write(21,'(i9)') patches(i,j) -1
!write(20,'(<num_tri_patches(i)>(i8,2x))') patches(i,1:num_tri_patches(i)) - 1
!write(20,*) patches(i,1:num_tri_patches(i)) - 1
end do  ! j, num_tri_patch
end do  ! i, num_patch

end do ! k, num solids
close (21)

! write dots to check w/ chunk
open(30,file= 'chunk_check.pt.txt')
! write points
!write(30,'(5(f15.4,2x))') points(1,1), points(1,2), points(1,3),1.,0.
write(30,'(5(f15.4,2x))') float(points(1,1)), float(points(1,2)), float(points(1,3)),1.,0.
do i = 2, num_points
!write(30,'(5(f15.4,2x))') points(i,1), points(i,2), points(i,3),1.,1.
write(30,'(5(f15.4,2x))') float(points(i,1)), float(points(i,2)), float(points(i,3)),1.,1.
end do
close(30)

! write 3 line triangles segments to check w/ chunk
open(30,file= 'chunk_check.tri.txt')
! write triangles
do i = 1, num_patches 
do j =  1  , num_tri_patches(i) 
!if (i>nx_dem*ny_dem) ii = 2
!if (i>2*nx_dem*ny_dem) ii = 3
!if (i>2*nx_dem*ny_dem+ 2*(nx_dem-1)) ii = 4
!if (i>2*nx_dem*ny_dem+ 4*(nx_dem-1)) ii = 5
!if (i>2*nx_dem*ny_dem+ 4*(nx_dem-1)+2*(ny_dem-1)) ii = 6
!print*, i,j,patches(i,j), triangles(patches(i,j),1), triangles(patches(i,j),2), triangles(patches(i,j),3)
!write(30,'(5(f15.4,2x))') points(triangles(patches(i,j),1),1), &
!  points(triangles(patches(i,j),1),2), points(triangles(patches(i,j),1),3),float(i),0.
!write(30,'(5(f15.4,2x))') points(triangles(patches(i,j),2),1), &
!  points(triangles(patches(i,j),2),2), points(triangles(patches(i,j),2),3),float(i),float(i)
!write(30,'(5(f15.4,2x))') points(triangles(patches(i,j),3),1), &
!  points(triangles(patches(i,j),3),2), points(triangles(patches(i,j),3),3),float(i),float(i)
!write(30,'(5(f15.4,2x))') points(triangles(patches(i,j),1),1), &
!  points(triangles(patches(i,j),1),2), points(triangles(patches(i,j),1),3),float(i),float(i)
write(30,'(5(f15.4,2x))') float(points(triangles(patches(i,j),1),1)), &
  float(points(triangles(patches(i,j),1),2)), float(points(triangles(patches(i,j),1),3)),float(i),0.
write(30,'(5(f15.4,2x))') float(points(triangles(patches(i,j),2),1)), &
  float(points(triangles(patches(i,j),2),2)), float(points(triangles(patches(i,j),2),3)),float(i),float(i)
write(30,'(5(f15.4,2x))') float(points(triangles(patches(i,j),3),1)), &
  float(points(triangles(patches(i,j),3),2)), float(points(triangles(patches(i,j),3),3)),float(i),float(i)
write(30,'(5(f15.4,2x))') float(points(triangles(patches(i,j),1),1)), &
  float(points(triangles(patches(i,j),1),2)), float(points(triangles(patches(i,j),1),3)),float(i),float(i)
end do
end do
close(30)

!open(40, file = 'chunk.mask.check.asc')
!write(40,*) nx_dem, ny_dem, 1
!do j = 1, ny_dem
!do i =1, nx_dem
!write(40,*) float(mask2(i,j))
!end do
!end do
!close(40)


end program pf_sol_gen


