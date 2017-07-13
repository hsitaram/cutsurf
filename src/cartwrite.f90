module cartwrite

implicit none
integer, parameter :: MAXSTRSIZE=500
integer, parameter :: DOUBLESIZE=8
integer, parameter :: INTSIZE=4
contains
!===============================================================================================
subroutine Writevtrfile(solnvars,nscalars,scalarnames,lx,ly,lz,dx,dy,dz,offx,offy,offz)

integer, intent(in) :: lx,ly,lz
integer, intent(in) :: nscalars
character(LEN=*), intent(in) :: scalarnames(nscalars)
real*8, intent(in) :: solnvars(1:lx,1:ly,1:lz,nscalars)
real*8, intent(in) :: dx,dy,dz
real*8, intent(in) :: offx,offy,offz

character(LEN=MAXSTRSIZE) :: fname
character(LEN=MAXSTRSIZE) :: str
character(LEN=8) :: offset
character(LEN=8) :: xmin_str
character(LEN=8) :: ymin_str
character(LEN=8) :: zmin_str
character(LEN=8) :: xmax_str
character(LEN=8) :: ymax_str
character(LEN=8) :: zmax_str
character(LEN=8) :: proc_str
character(LEN=MAXSTRSIZE) :: global_extent
character(LEN=MAXSTRSIZE) :: local_extent
character :: slashn

integer :: ierr
integer :: filenum
integer :: i,j,k
integer :: byte_offset
integer :: scalars_size

real*8,allocatable :: local_xc(:),local_yc(:),local_zc(:)

filenum = 9
slashn  = char(10)

fname = 'file.vtr'

allocate(local_xc(lx+1))
allocate(local_yc(ly+1))
allocate(local_zc(lz+1))
	
do i=1,lx+1
	local_xc(i) = offx + (i-1)*dx
enddo

do i=1,ly+1
	local_yc(i) = offy + (i-1)*dy
enddo

do i=1,lz+1
	local_zc(i) = offz + (i-1)*dz
enddo

scalars_size   = lx*ly*lz * DOUBLESIZE
	
write(xmin_str(1:8),'(i8)') 0
write(ymin_str(1:8),'(i8)') 0
write(zmin_str(1:8),'(i8)') 0
write(xmax_str(1:8),'(i8)') lx
write(ymax_str(1:8),'(i8)') ly
write(zmax_str(1:8),'(i8)') lz

local_extent = trim(xmin_str)//' '//trim(xmax_str)//' '//trim(ymin_str)//' '&
		//trim(ymax_str)//trim(zmin_str)//' '//trim(zmax_str)

byte_offset = 0                             
open(unit=filenum,file=fname,form='unformatted',access='stream')

str = '<?xml version="1.0"?>'//slashn                                                                                                   
write(filenum) trim(str)

str = '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'//slashn		                               
write(filenum) trim(str)

str = '  <RectilinearGrid WholeExtent="'//trim(local_extent)//'">'//slashn                                       
write(filenum) trim(str)

str = '    <Piece Extent="'//trim(local_extent)//'">'//slashn
write(filenum) trim(str)

str = '      <PointData>  </PointData>'//slashn
write(filenum) trim(str)

str = '      <CellData> '//slashn                                                       
write(filenum) trim(str)

do i=1,nscalars

	write(offset(1:8),'(i8)') byte_offset
	str = '         <DataArray type="Float64" Name="'//trim(scalarnames(i))//&
		'" format="appended" offset="'//offset//'"           />'//slashn
	write(filenum) trim(str)
	byte_offset = byte_offset + intsize + scalars_size

enddo

str = '      </CellData>'//slashn            
write(filenum) trim(str)

str = '      <Coordinates>'//slashn                                                                                   
write(filenum) trim(str)

write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="X"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + intsize + (lx+1)*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="Y"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

byte_offset = byte_offset + intsize + (ly+1)*DOUBLESIZE
write(offset(1:8),'(i8)') byte_offset
str = '<DataArray type="Float64" Name="Z"  format="appended" offset="'//offset//'" />'//slashn
write(filenum) trim(str)

str = '      </Coordinates>'//slashn       
write(filenum) trim(str)

str = '    </Piece>'//slashn                                                                                                  
write(filenum) trim(str)

str = '  </RectilinearGrid>'//slashn                                                                                                   
write(filenum) trim(str)

str = '  <AppendedData encoding="raw">'//slashn                                                                                         
write(filenum) trim(str)

str = '_'                                                                                                                           
write(filenum) trim(str)

do i=1,nscalars
write(filenum) scalars_size  ,  solnvars(1:lx,1:ly,1:lz,i)
enddo

write(filenum) (lx+1)*DOUBLESIZE  , local_xc(:)
write(filenum) (ly+1)*DOUBLESIZE  , local_yc(:)
write(filenum) (lz+1)*DOUBLESIZE  , local_zc(:)

str = slashn//'  </AppendedData>'//slashn                                                                             
write(filenum) trim(str)

str = '</VTKFile>'//slashn                                                                                                             
write(filenum) trim(str)

close(filenum)

end subroutine Writevtrfile
!===============================================================================================

end module cartwrite

