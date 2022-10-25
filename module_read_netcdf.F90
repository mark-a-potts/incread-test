#define ESMF_ERR_RETURN(rc) \
    if (ESMF_LogFoundError(rc, msg="Breaking out of subroutine", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

#define NC_ERR_STOP(status) \
    if (status /= nf90_noerr) write(0,*) "file: ", __FILE__, " line: ", __LINE__, trim(nf90_strerror(status)); \
    if (status /= nf90_noerr) call ESMF_Finalize(endflag=ESMF_END_ABORT)

module module_read_netcdf

  use esmf
  use netcdf
  use atmosphere_mod,     only: Atm
  use tracer_manager_mod, only: get_tracer_index
  use field_manager_mod,  only: MODEL_ATMOS
  use mpi

  implicit none
  public read_netcdf

  logical :: par, testing


  contains

!----------------------------------------------------------------------------------------
  subroutine read_netcdf(filename, &
                          use_parallel_netcdf, mpi_comm, mype, &
                          grid_id,rc)
!
!   type(ESMF_FieldBundle), intent(in) :: readfb
    character(*), intent(in)           :: filename
    logical, intent(in)                :: use_parallel_netcdf
    integer, intent(in)                :: mpi_comm
    integer, intent(in)                :: mype
    integer, intent(in)                :: grid_id
    integer, optional,intent(out)      :: rc
!
!** local vars
    integer :: i,j,t, istart,iend,jstart,jend
    integer :: im, jm, lm

    integer, dimension(:), allocatable              :: fldlev

    real(ESMF_KIND_R4), dimension(:,:), pointer     :: array_r4
    real(ESMF_KIND_R4), dimension(:,:,:), pointer   :: array_r4_cube
    real(ESMF_KIND_R4), dimension(:,:,:), pointer   :: array_r4_3d
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer :: array_r4_3d_cube

    real(ESMF_KIND_R8), dimension(:), pointer     :: tileval
    real(ESMF_KIND_R8), dimension(:), pointer     :: timeval
    real(ESMF_KIND_R8), dimension(:), pointer     :: phalf
    real(ESMF_KIND_R8), dimension(:), pointer     :: pfull
    real(ESMF_KIND_R8), dimension(:,:), pointer     :: array_r8
    real(ESMF_KIND_R8), dimension(:,:,:), pointer   :: array_r8_cube
    real(ESMF_KIND_R8), dimension(:,:,:), pointer   :: array_r8_3d
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: array_r8_3d_cube
    real(ESMF_KIND_R8), dimension(:,:,:,:,:), pointer :: array_r8_3d_tiled
    real(ESMF_KIND_R8), dimension(:,:,:,:,:), pointer :: ice_wat
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: lon

    real(8), dimension(:), allocatable :: x,y
    integer, dimension(:), allocatable :: xt,yt
    integer :: fieldCount, fieldDimCount, gridDimCount
    integer, dimension(:), allocatable   :: ungriddedLBound, ungriddedUBound
    integer, dimension(:), allocatable   :: start_idx
    integer, dimension(5)                :: start,endpos

    type(ESMF_Field), allocatable        :: fcstField(:)
    type(ESMF_TypeKind_Flag)             :: typekind
    type(ESMF_TypeKind_Flag)             :: attTypeKind
    type(ESMF_Grid)                      :: readgrid
    type(ESMF_Array)                     :: array
    type(ESMF_DistGrid)                  :: distgrid

    integer :: attCount
    character(len=ESMF_MAXSTR) :: attName, fldName

    integer :: varival
    real(4) :: varr4val, dataMin, dataMax
    real(4), allocatable, dimension(:) :: compress_err
    real(8) :: varr8val
    character(len=ESMF_MAXSTR) :: varcval

    integer :: ncerr,ierr
    integer :: ncid
    integer :: oldMode
    integer :: im_dimid, jm_dimid, tile_dimid, pfull_dimid, phalf_dimid, time_dimid, ch_dimid
    integer :: tm,tl,pf,ph
    integer :: im_varid, jm_varid, tile_varid, lon_varid, lat_varid, timeiso_varid
    integer :: time_varid, pfull_varid, phalf_varid, icewat_varid 
    integer :: liqwat_varid, sphum_varid, o3mr_varid
    integer :: ugrd_varid, vgrd_varid, dpres_varid, tmp_varid 
    integer, dimension(:), allocatable :: dimids_2d, dimids_3d
    logical shuffle

    logical :: is_cubed_sphere
    integer :: rank, deCount, localDeCount, dimCount, tileCount
    integer :: my_tile, start_i, start_j
    integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe
    integer, dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
    integer, dimension(:), allocatable :: deToTileMap, localDeToDeMap
    logical :: do_io
    integer :: par_access, nvar, xtype, ndims, nAtts, numdims, rhvarid, numatts
    integer :: sphum_idx, liq_wat_idx, ice_wat_idx, o3mr_idx
    integer, dimension(:), allocatable :: varids
    character(len=NF90_MAX_NAME) :: varname
    integer, dimension(nf90_max_var_dims) :: rhDimIds, dims
!
    testing = .true.
    is_cubed_sphere = .true.
    tileCount = 0
    my_tile = 0
    start_i = -10000000
    start_j = -10000000

    par = use_parallel_netcdf

    ! create netcdf file and enter define mode

    if (par) then
       ! not implemented yet
       ncerr = nf90_open(trim(filename),&
               mode=nf90_nowrite, &
               comm=mpi_comm, info = MPI_INFO_NULL, ncid=ncid); NC_ERR_STOP(ncerr)
    else
       ncerr = nf90_open(trim(filename),&
               mode=nf90_nowrite, &
               ncid=ncid); NC_ERR_STOP(ncerr)
    end if
    ncerr = nf90_inquire(ncid, nvariables = nvar)
    allocate(varids(nvar))
    ncerr = nf90_inq_varids(ncid, nvar, varids)
    do i=1,nvar
      ncerr = nf90_inquire_variable(ncid, varids(i), name=varname, xtype=xtype, ndims=ndims, nAtts=nAtts)
      write(6,*) 'Name is ',trim(varname), varids(i),ndims
    enddo
    varname = "grid_xt"
    ncerr = nf90_inq_dimid(ncid, "grid_xt", im_dimid) 
    ncerr = nf90_inquire_dimension(ncid,im_dimid,len=im)
    ncerr = nf90_inq_varid(ncid,varname,im_varid)
    varname = "grid_yt"
    ncerr = nf90_inq_dimid(ncid, varname, jm_dimid) 
    ncerr = nf90_inq_varid(ncid,varname,jm_varid)
    ncerr = nf90_inquire_dimension(ncid,jm_dimid,len=jm)
    varname = "time"
    ncerr = nf90_inq_dimid(ncid, varname, time_dimid) 
    ncerr = nf90_inquire_dimension(ncid,time_dimid,len=tm)
    ncerr = nf90_inq_varid(ncid,varname,time_varid)
    varname = "tile"
    ncerr = nf90_inq_dimid(ncid, varname, tile_dimid) 
    ncerr = nf90_inquire_dimension(ncid,tile_dimid,len=tileCount)
    write(6,*) 'HEY, tilecount is, ',tileCount 
    ncerr = nf90_inq_varid(ncid,varname,tile_varid)
    varname = "pfull"
    ncerr = nf90_inq_dimid(ncid, varname, pfull_dimid) 
    ncerr = nf90_inquire_dimension(ncid,pfull_dimid,len=pf)
    ncerr = nf90_inq_varid(ncid,varname,pfull_varid)
    varname = "phalf"
    ncerr = nf90_inq_dimid(ncid, varname, phalf_dimid) 
    ncerr = nf90_inquire_dimension(ncid,phalf_dimid,len=ph)
    ncerr = nf90_inq_varid(ncid,varname,phalf_varid)
    varname = "ice_wat"
    ncerr = nf90_inq_varid(ncid,varname,icewat_varid)
    varname = "liq_wat"
    ncerr = nf90_inq_varid(ncid,varname,liqwat_varid)
    varname = "sphum"
    ncerr = nf90_inq_varid(ncid,varname,sphum_varid)
    varname = "o3m4"
    ncerr = nf90_inq_varid(ncid,varname,o3mr_varid)
    varname = "ugrd"
    ncerr = nf90_inq_varid(ncid,varname,ugrd_varid)
    varname = "vgrd"
    ncerr = nf90_inq_varid(ncid,varname,vgrd_varid)
    varname = "dpres"
    ncerr = nf90_inq_varid(ncid,varname,dpres_varid)
    varname = "tmp"
    ncerr = nf90_inq_varid(ncid,varname,tmp_varid)
    allocate(x(im))
    allocate(y(jm))
    allocate(timeval(tm))
    allocate(tileval(tileCount))
    allocate(phalf(ph))
    allocate(pfull(pf))
    allocate(lon(im,jm,tileCount))
    allocate(ice_wat(im,jm,pf,tileCount,tm))
    if(testing) then
      ! allocate 6 tiles for Atm
      allocate(Atm(6))
      ! assign dummy indices
      sphum_idx   = 1
      liq_wat_idx = 2
      ice_wat_idx = 3
      o3mr_idx    = 4
      ! Allocate space in Atm for testing 
      do i=1,tileCount
        allocate(Atm(i)%u(im,jm,pf))
        allocate(Atm(i)%v(im,jm,pf))
        allocate(Atm(i)%pt(im,jm,pf))
        allocate(Atm(i)%delp(im,jm,pf))
        allocate(Atm(i)%q(im,jm,pf,4)) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
      enddo
    else
      sphum_idx   = get_tracer_index(MODEL_ATMOS, 'sphum')
      liq_wat_idx = get_tracer_index(MODEL_ATMOS, 'liq_wat')
      ice_wat_idx = get_tracer_index(MODEL_ATMOS, 'ice_wat')
      o3mr_idx    = get_tracer_index(MODEL_ATMOS, 'o3mr')
    endif
    allocate(array_r8_3d_tiled(im,jm,pf,tileCount,tm))
    ncerr = nf90_get_var(ncid, jm_varid, y) 
    ncerr = nf90_get_var(ncid, im_varid, x) 
    ncerr = nf90_get_var(ncid, time_varid, timeval) 
    ncerr = nf90_get_var(ncid, tile_varid, tileval) 
    ncerr = nf90_get_var(ncid, pfull_varid, pfull) 
    ncerr = nf90_get_var(ncid, phalf_varid, phalf) 
    ncerr = nf90_get_var(ncid, lon_varid, lon) 
    ncerr = nf90_get_var(ncid, icewat_varid, ice_wat) 

!   start = [2,2,2,1,1]
!   endpos = [4,4,4,6,1]
!   array_r8_3d_tiled(:,:,:,:,:) = 0.0
!   ncerr = nf90_get_var(ncid, ugrd_varid, array_r8_3d_tiled, start=start, count=endpos) 
!   write(6,*) 'ugrd is ',array_r8_3d_tiled(1:6,1,1,1,1)
    ncerr = nf90_get_var(ncid, ugrd_varid, array_r8_3d_tiled) 
    do i=1,tileCount
       Atm(i)%u(:,:,:) = array_r8_3d_tiled(:,:,:,i,1)
    enddo 
    ncerr = nf90_get_var(ncid, vgrd_varid, array_r8_3d_tiled) 
    do i=1,tileCount
       Atm(i)%v(:,:,:) = array_r8_3d_tiled(:,:,:,i,1)
    enddo 
    ncerr = nf90_get_var(ncid, tmp_varid, array_r8_3d_tiled) 
    do i=1,tileCount
       Atm(i)%pt(:,:,:) = array_r8_3d_tiled(:,:,:,i,1)
    enddo 
    ncerr = nf90_get_var(ncid, dpres_varid, array_r8_3d_tiled) 
    do i=1,tileCount
       Atm(i)%delp(:,:,:) = array_r8_3d_tiled(:,:,:,i,1)
    enddo 
    ncerr = nf90_get_var(ncid, sphum_varid, array_r8_3d_tiled) 
    do i=1,tileCount
       Atm(i)%q(:,:,:,sphum_idx) = array_r8_3d_tiled(:,:,:,i,1)
    enddo 
    ncerr = nf90_get_var(ncid, icewat_varid, array_r8_3d_tiled) 
    do i=1,tileCount
       Atm(i)%q(:,:,:,ice_wat_idx) = array_r8_3d_tiled(:,:,:,i,1)
    enddo 
    ncerr = nf90_get_var(ncid, liqwat_varid, array_r8_3d_tiled) 
    do i=1,tileCount
       Atm(i)%q(:,:,:,liq_wat_idx) = array_r8_3d_tiled(:,:,:,i,1)
    enddo 
    ncerr = nf90_get_var(ncid, o3mr_varid, array_r8_3d_tiled) 
    do i=1,tileCount
       Atm(i)%q(:,:,:,o3mr_idx) = array_r8_3d_tiled(:,:,:,i,1)
    enddo 
    if(testing) then
      do i=1,tileCount
        deallocate(Atm(i)%u)
        deallocate(Atm(i)%v)
        deallocate(Atm(i)%pt)
        deallocate(Atm(i)%delp)
        deallocate(Atm(i)%q) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
      enddo
      deallocate(Atm) ! currently, just need specific humdity, liq_wat, ice_wat, and o3mr
    endif 
    deallocate(array_r8_3d_tiled)

  end subroutine read_netcdf

!----------------------------------------------------------------------------------------
end module module_read_netcdf
