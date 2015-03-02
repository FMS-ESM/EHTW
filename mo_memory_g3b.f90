#undef G3BTEST
MODULE mo_memory_g3b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add => add_stream_element, &
                            default_stream_setting,                   &
                            ABOVESUR2, ABOVESUR10, BELOWSUR, HYBRID_H,&
                            HYBRID, SNOWICE, SNOWICE2, WATER, WATERF, &
                            OCEAN, OCEANF, OGAUSSIAN
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g3b ! routine to construct the g3b table
  PUBLIC :: destruct_g3b  ! routine to destruct  the g3b table
  PUBLIC :: g3b           ! the g3b table

  ! declaration of predefined fields within this module

  REAL(dp), POINTER, PUBLIC :: geosp(:,:)
  REAL(dp), POINTER, PUBLIC :: tsl(:,:)
  REAL(dp), POINTER, PUBLIC :: ws(:,:)
  REAL(dp), POINTER, PUBLIC :: wl(:,:)
  REAL(dp), POINTER, PUBLIC :: sn(:,:)
  REAL(dp), POINTER, PUBLIC :: slm(:,:)
  REAL(dp), POINTER, PUBLIC :: az0(:,:)
  REAL(dp), POINTER, PUBLIC :: alb(:,:)
  REAL(dp), POINTER, PUBLIC :: forest(:,:)
  REAL(dp), POINTER, PUBLIC :: vgrat(:,:)
  REAL(dp), POINTER, PUBLIC :: vlt(:,:)
  REAL(dp), POINTER, PUBLIC :: wsmx(:,:)
  REAL(dp), POINTER, PUBLIC :: fao(:,:)
  REAL(dp), POINTER, PUBLIC :: aps(:,:)
  REAL(dp), POINTER, PUBLIC :: aprl(:,:)
  REAL(dp), POINTER, PUBLIC :: aprc(:,:)
  REAL(dp), POINTER, PUBLIC :: aprs(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrgw(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrgw(:,:)
  REAL(dp), POINTER, PUBLIC :: vdisgw(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcov(:,:)
  REAL(dp), POINTER, PUBLIC :: temp2(:,:)
  REAL(dp), POINTER, PUBLIC :: dew2(:,:)
  REAL(dp), POINTER, PUBLIC :: wind10(:,:)
  REAL(dp), POINTER, PUBLIC :: swnir(:,:)
  REAL(dp), POINTER, PUBLIC :: swdifnir(:,:)
  REAL(dp), POINTER, PUBLIC :: swvis(:,:)
  REAL(dp), POINTER, PUBLIC :: swdifvis(:,:)
  REAL(dp), POINTER, PUBLIC :: swnirac(:,:)
  REAL(dp), POINTER, PUBLIC :: swdifnirac(:,:)
  REAL(dp), POINTER, PUBLIC :: swvisac(:,:)
  REAL(dp), POINTER, PUBLIC :: swdifvisac(:,:)
  REAL(dp), POINTER, PUBLIC :: u10(:,:)
  REAL(dp), POINTER, PUBLIC :: v10(:,:)
  REAL(dp), POINTER, PUBLIC :: srads(:,:)
  REAL(dp), POINTER, PUBLIC :: trads(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0(:,:)
  REAL(dp), POINTER, PUBLIC :: trad0(:,:)
  REAL(dp), POINTER, PUBLIC :: vdis(:,:)
  REAL(dp), POINTER, PUBLIC :: ustr(:,:)
  REAL(dp), POINTER, PUBLIC :: vstr(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfs(:,:)
  REAL(dp), POINTER, PUBLIC :: evap(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfl(:,:)
  REAL(dp), POINTER, PUBLIC :: tslm(:,:)
  REAL(dp), POINTER, PUBLIC :: tslm1(:,:)
  REAL(dp), POINTER, PUBLIC :: emter(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsol(:,:,:)
  REAL(dp), POINTER, PUBLIC :: runoff(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0u(:,:)
  REAL(dp), POINTER, PUBLIC :: sradsu(:,:)
  REAL(dp), POINTER, PUBLIC :: tradsu(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo(:,:)
  REAL(dp), POINTER, PUBLIC :: tsurf(:,:)
  REAL(dp), POINTER, PUBLIC :: seaice(:,:)
  REAL(dp), POINTER, PUBLIC :: obsseaice(:,:)
  REAL(dp), POINTER, PUBLIC :: siced(:,:)
  REAL(dp), POINTER, PUBLIC :: relhum(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wind10w(:,:)
  REAL(dp), POINTER, PUBLIC :: glac(:,:)
  REAL(dp), POINTER, PUBLIC :: gld(:,:)
  REAL(dp), POINTER, PUBLIC :: aclc(:,:,:)
  REAL(dp), POINTER, PUBLIC :: aclcac(:,:,:)
  REAL(dp), POINTER, PUBLIC :: snmel(:,:)
  REAL(dp), POINTER, PUBLIC :: runtoc(:,:)
  REAL(dp), POINTER, PUBLIC :: apmegl(:,:)
  REAL(dp), POINTER, PUBLIC :: t2max(:,:)
  REAL(dp), POINTER, PUBLIC :: t2min(:,:)
  REAL(dp), POINTER, PUBLIC :: wimax(:,:)
  REAL(dp), POINTER, PUBLIC :: topmax(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcv(:,:)
  REAL(dp), POINTER, PUBLIC :: qvi(:,:)
  REAL(dp), POINTER, PUBLIC :: xlvi(:,:)
  REAL(dp), POINTER, PUBLIC :: xivi(:,:)
  REAL(dp), POINTER, PUBLIC :: runlnd(:,:)
  REAL(dp), POINTER, PUBLIC :: rgcgn(:,:)
  REAL(dp), POINTER, PUBLIC :: sodif(:,:)
  REAL(dp), POINTER, PUBLIC :: srafs(:,:)
  REAL(dp), POINTER, PUBLIC :: trafs(:,:)
  REAL(dp), POINTER, PUBLIC :: sraf0(:,:)
  REAL(dp), POINTER, PUBLIC :: traf0(:,:)
  REAL(dp), POINTER, PUBLIC :: emtef(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsof(:,:,:)
  REAL(dp), POINTER, PUBLIC :: emtef0(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsof0(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tke(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tkem(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tkem1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xvar(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xskew(:,:,:)
  REAL(dp), POINTER, PUBLIC :: drain(:,:)
  REAL(dp), POINTER, PUBLIC :: grndcapc(:,:)
  REAL(dp), POINTER, PUBLIC :: grndhflx(:,:)
  REAL(dp), POINTER, PUBLIC :: grndflux(:,:)
  REAL(dp), POINTER, PUBLIC :: tsoil(:,:,:)
  REAL(dp), POINTER, PUBLIC :: grndc(:,:,:)
  REAL(dp), POINTER, PUBLIC :: grndd(:,:,:)
  REAL(dp), POINTER, PUBLIC :: srad0d(:,:)
  REAL(dp), POINTER, PUBLIC :: acdnc(:,:,:)
  REAL(dp), POINTER, PUBLIC :: snacl(:,:)
  REAL(dp), POINTER, PUBLIC :: rogl(:,:)
  REAL(dp), POINTER, PUBLIC :: alake(:,:)
  REAL(dp), POINTER, PUBLIC :: aprflux(:,:)
  REAL(dp), POINTER, PUBLIC :: acvtype(:,:)
  REAL(dp), POINTER, PUBLIC :: xtec(:,:,:)
  REAL(dp), POINTER, PUBLIC :: slf(:,:)
  REAL(dp), POINTER, PUBLIC :: snc(:,:)
  REAL(dp), POINTER, PUBLIC :: rtype(:,:)
  REAL(dp), POINTER, PUBLIC :: rintop(:,:)
  REAL(dp), POINTER, PUBLIC :: apmeb(:,:)
  REAL(dp), POINTER, PUBLIC :: apmebco(:,:)
  REAL(dp), POINTER, PUBLIC :: rain(:,:)
  REAL(dp), POINTER, PUBLIC :: qtnew(:,:)
  REAL(dp), POINTER, PUBLIC :: abso4(:,:)
  REAL(dp), POINTER, PUBLIC :: so4nat(:,:,:)
  REAL(dp), POINTER, PUBLIC :: so4all(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ao3(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tropo(:,:)
  !
  !  variables for fractional surface coverage
  !
  REAL(dp), POINTER, PUBLIC :: tsi(:,:)
  REAL(dp), POINTER, PUBLIC :: tsw(:,:)
  REAL(dp), POINTER, PUBLIC :: sni(:,:)
  REAL(dp), POINTER, PUBLIC :: ustri(:,:)
  REAL(dp), POINTER, PUBLIC :: vstri(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrw(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrw(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrl(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrl(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfli(:,:)
  REAL(dp), POINTER, PUBLIC :: ahflw(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfll(:,:)
  REAL(dp), POINTER, PUBLIC :: az0i(:,:)
  REAL(dp), POINTER, PUBLIC :: az0w(:,:)
  REAL(dp), POINTER, PUBLIC :: az0l(:,:)
  REAL(dp), POINTER, PUBLIC :: alsoi(:,:)
  REAL(dp), POINTER, PUBLIC :: alsow(:,:)
  REAL(dp), POINTER, PUBLIC :: alsol(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfice(:,:)
  REAL(dp), POINTER, PUBLIC :: qres(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfcon(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfres(:,:)
  REAL(dp), POINTER, PUBLIC :: fluxres(:,:)
  !
  REAL(dp), POINTER, PUBLIC :: ahfsiac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfswac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfslac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfliac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahflwac(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfllac(:,:)
  REAL(dp), POINTER, PUBLIC :: evapiac(:,:)
  REAL(dp), POINTER, PUBLIC :: evapwac(:,:)
  REAL(dp), POINTER, PUBLIC :: evaplac(:,:)
  REAL(dp), POINTER, PUBLIC :: trfllac(:,:)
  REAL(dp), POINTER, PUBLIC :: trflwac(:,:)
  REAL(dp), POINTER, PUBLIC :: trfliac(:,:)
  REAL(dp), POINTER, PUBLIC :: sofllac(:,:)
  REAL(dp), POINTER, PUBLIC :: soflwac(:,:)
  REAL(dp), POINTER, PUBLIC :: sofliac(:,:)
  REAL(dp), POINTER, PUBLIC :: friac(:,:)
  !
  !  variables for ocean coupling only
  !
  REAL(dp), POINTER, PUBLIC :: awhea(:,:)
  REAL(dp), POINTER, PUBLIC :: awsol(:,:)
  REAL(dp), POINTER, PUBLIC :: awfre(:,:)
  REAL(dp), POINTER, PUBLIC :: awust(:,:)
  REAL(dp), POINTER, PUBLIC :: awvst(:,:)
  REAL(dp), POINTER, PUBLIC :: awsta(:,:)
  REAL(dp), POINTER, PUBLIC :: aicon(:,:)
  REAL(dp), POINTER, PUBLIC :: aiqre(:,:)
  REAL(dp), POINTER, PUBLIC :: aifre(:,:)
  REAL(dp), POINTER, PUBLIC :: aiust(:,:)
  REAL(dp), POINTER, PUBLIC :: aivst(:,:)
  REAL(dp), POINTER, PUBLIC :: ocu(:,:)
  REAL(dp), POINTER, PUBLIC :: ocv(:,:)
  !
  !  new surface variables proposed by bjt
  !
  REAL(dp), POINTER, PUBLIC :: lclass(:,:)
  REAL(dp), POINTER, PUBLIC :: tinertia(:,:)
  REAL(dp), POINTER, PUBLIC :: bathy(:,:)
  !
  !  variables for sit coupling only
  !
  REAL(dp), POINTER, PUBLIC :: awufl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: awvfl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: awtfl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: awsfl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: awtkefl(:,:,:)

  REAL(dp), POINTER, PUBLIC :: wtfn(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wsfn(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wtfns(:,:)
  REAL(dp), POINTER, PUBLIC :: wsfns(:,:)
  REAL(dp), POINTER, PUBLIC :: sitmask(:,:)
  REAL(dp), POINTER, PUBLIC :: oceanid(:,:)
  REAL(dp), POINTER, PUBLIC :: ctfreez2(:,:)

  REAL(dp), POINTER, PUBLIC :: obstsw(:,:)
  REAL(dp), POINTER, PUBLIC :: obswsb(:,:)
  REAL(dp), POINTER, PUBLIC :: obswt(:,:,:)
  REAL(dp), POINTER, PUBLIC :: obsws(:,:,:)
  REAL(dp), POINTER, PUBLIC :: obswu(:,:,:)
  REAL(dp), POINTER, PUBLIC :: obswv(:,:,:)

  REAL(dp), POINTER, PUBLIC :: sitzsi(:,:,:)    
  REAL(dp), POINTER, PUBLIC :: sitsilw(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sittsi(:,:,:)  
  REAL(dp), POINTER, PUBLIC :: sitwt(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitwu(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitwv(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitww(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitwp(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitws(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitwtke(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitwlvl(:,:)

  REAL(dp), POINTER, PUBLIC :: sitwtb(:,:)
  REAL(dp), POINTER, PUBLIC :: sitwub(:,:)
  REAL(dp), POINTER, PUBLIC :: sitwvb(:,:)
  REAL(dp), POINTER, PUBLIC :: sitwsb(:,:)

  REAL(dp), POINTER, PUBLIC :: sitwlmx(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitwldisp(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitwkm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sitwkh(:,:,:)
  REAL(dp), POINTER, PUBLIC :: wrho(:,:,:)
  !
  ! Instantaneous variables required for calculation of skin temperature (SIT - model )
  !
!!!  REAL(dp), POINTER, PUBLIC :: fluxs(:,:)
  REAL(dp), POINTER, PUBLIC :: dfluxs(:,:)
  REAL(dp), POINTER, PUBLIC :: fluxiw(:,:)

  REAL(dp), POINTER, PUBLIC :: cc(:,:)
  REAL(dp), POINTER, PUBLIC :: hc(:,:)
  REAL(dp), POINTER, PUBLIC :: engwac(:,:)
!!!  REAL(dp), POINTER, PUBLIC :: engw(:,:)
!!!  REAL(dp), POINTER, PUBLIC :: engw2(:,:)

  REAL(dp), POINTER, PUBLIC :: pme(:,:)
  REAL(dp), POINTER, PUBLIC :: pme2(:,:)
  REAL(dp), POINTER, PUBLIC :: sc(:,:)
  REAL(dp), POINTER, PUBLIC :: saltwac(:,:)

  REAL(dp), POINTER, PUBLIC :: subfluxw(:,:)
  REAL(dp), POINTER, PUBLIC :: wsubsal(:,:)

  !
  !  variables for embedded ocean (ocean) coupling only
  !
  REAL(dp), POINTER, PUBLIC :: ocn_oromea(:,:),ocn_wf(:,:),ocn_bsl_oromea(:,:),ocn_divzmea(:,:),    &
                               ocn_divzmin(:,:),ocn_divzmax(:,:),ocn_bsl_oropic(:,:),               &
                               ocn_bsl_oroval(:,:),ocn_por(:,:,:),ocn_porx(:,:,:),ocn_pory(:,:,:)
 
  REAL(dp), POINTER, PUBLIC :: obox_mask(:,:)
  REAL(dp), POINTER, PUBLIC :: ocnmask(:,:)
  REAL(dp), POINTER, PUBLIC :: ocnp(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ocnt(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ocnu(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ocnv(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ocnw(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ocns(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ocnkvm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ocnkvh(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ocnkhm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: ocnkhh(:,:,:)

  REAL(dp), POINTER, PUBLIC :: afluxs(:,:)
  REAL(dp), POINTER, PUBLIC :: apme(:,:)  
  REAL(dp), POINTER, PUBLIC :: awust2(:,:)
  REAL(dp), POINTER, PUBLIC :: awvst2(:,:)  
  REAL(dp), POINTER, PUBLIC :: aicon2(:,:)
  REAL(dp), POINTER, PUBLIC :: aiqre2(:,:)  
  REAL(dp), POINTER, PUBLIC :: aiust2(:,:)  
  REAL(dp), POINTER, PUBLIC :: aivst2(:,:)
  REAL(dp), POINTER, PUBLIC :: awsol2(:,:)  
  REAL(dp), POINTER, PUBLIC :: awsta2(:,:)

  REAL(dp), POINTER, PUBLIC :: afluxiw(:,:)  
  REAL(dp), POINTER, PUBLIC :: apme2(:,:)  
  REAL(dp), POINTER, PUBLIC :: asubfluxw(:,:)  
  REAL(dp), POINTER, PUBLIC :: awsubsal(:,:)  
  !
  !  variables for coupling with HD-model and calving model only
  !
  REAL(dp), POINTER, PUBLIC :: aros(:,:)
  REAL(dp), POINTER, PUBLIC :: adrain(:,:)
  REAL(dp), POINTER, PUBLIC :: disch(:,:)
  REAL(dp), POINTER, PUBLIC :: apmecal(:,:)
  !
  !  variables for sso parametrization
  !
  REAL(dp), POINTER, PUBLIC :: oromea(:,:)
  REAL(dp), POINTER, PUBLIC :: orostd(:,:)
  REAL(dp), POINTER, PUBLIC :: orosig(:,:)
  REAL(dp), POINTER, PUBLIC :: orogam(:,:)
  REAL(dp), POINTER, PUBLIC :: orothe(:,:)
  REAL(dp), POINTER, PUBLIC :: oropic(:,:)
  REAL(dp), POINTER, PUBLIC :: oroval(:,:)
!
!  REAL(dp), POINTER, PUBLIC :: trflw(:,:)
!  REAL(dp), POINTER, PUBLIC :: soflw(:,:)
!  REAL(dp), POINTER, PUBLIC :: ahfsw(:,:)
  !
  !
  !  variables for mixed layer ocean only
  !
  REAL(dp), POINTER, PUBLIC :: amlcorr(:,:)
  REAL(dp), POINTER, PUBLIC :: amlcorac(:,:)
  REAL(dp), POINTER, PUBLIC :: amlheatac(:,:)
  !
  !  variables for 200mb radiation
  !
  REAL(dp), POINTER, PUBLIC :: tradl(:,:)
  REAL(dp), POINTER, PUBLIC :: sradl(:,:)
  REAL(dp), POINTER, PUBLIC :: trafl(:,:)
  REAL(dp), POINTER, PUBLIC :: srafl(:,:)

  ! declaration of table with 2d- and 3d-field entries

  TYPE (t_stream), POINTER     :: g3b

CONTAINS

!!! Reserved code number by after burner
!!! code= 34  	low_cld  	 low cloud  	 single  	 223
!!! code= 35 	        mid_cld 	mid cloud 	single 	223
!!! code= 36 	        hih_cld 	high cloud 	single 	223
!!! code= 131 	u 	u-velocity 	atm (ml+pl) 	138, 155
!!! code= 132 	v 	v-velocity 	atm (ml+pl) 	138, 155
!!! code= 135 	omega 	vertical velocity 	atm (ml+pl) 	138, 152, 155
!!! code= 148 	stream 	streamfunction 	atm (ml+pl) 	131, 132
!!! code= 149 	velopot 	velocity potential 	atm (ml+pl) 	131, 132
!!! code= 151 	slp 	mean sea level pressure 	surface 	129, 130, 152
!!! code= 156 	geopoth 	geopotential height 	atm (ml+pl) 	129, 130, 133, 152
!!! code= 157 	rhumidity 	relative humidity 	atm (ml+pl) 	130, 133, 152
!!! code= 189 	sclfs 	surface solar cloud forcing 	surface 	176-185
!!! code= 190 	tclfs 	surface thermal cloud forcing 	surface 	177-186
!!! code= 191 	sclf0 	top solar cloud forcing 	surface 	178-187
!!! code= 192 	tclf0 	top thermal cloud forcing 	surface 	179-188
!!! code= 259 	windspeed 	windspeed 	atm (ml+pl) 	sqrt(u2+v2)
!!! code= 260 	precip 	total precipitation 	surface 	142+143
!!! code= 261 	net_top 	total top radiation 	surface 	178+179 	Not available for ECHAM5
!!! code= 262 	net_bot 	total surface radiation 	surface 	176+177 	Not available for ECHAM5
!!! code= 263 	net_heat 	net surface heat flux 	surface 	146+147+176+177-220-C*218    C=Lf*RhoH2O
!!!                Lf: Latent heat of fusion    RhoH2O: Density of water 	Not available for ECHAM5
!!! code= 264 	net_water 	total surface water 	surface 	142+143+182-160-221 	Not available for ECHAM5 

  SUBROUTINE construct_g3b

    USE mo_control,    ONLY: lmidatm,lcouple,lhd,lsit,locaf,lasia,lgodas,lasia,locn, &
                             lwarning_msg,locn_msg
    IMPLICIT NONE                             
    LOGICAL :: lwarnGE2,lwarnGE2orlocn_msg
    
    ! set default attributes for the g3b stream
    lwarnGE2=(lwarning_msg.GE.2)
    lwarnGE2orlocn_msg=((lwarning_msg.GE.2).OR.locn_msg)
    
    CALL default_stream_setting (g3b               &
                                ,lrerun=.TRUE.     &
                                ,lpost=.TRUE.      &
                                ,table=128 ,bits=16)

    ! Add fields to the g3b stream.
    ! despite some 3-d fields (lpost=.FALSE.) these fields are written out by default

    CALL add (g3b,'qtnew',    qtnew    ,lpost=.FALSE.,contnorest=.true.)
    CALL add (g3b,'swnir',    swnir    ,lpost=.FALSE., longname='net surface NIR'          &
         ,units='W/m**2' ,contnorest=.true.)
    CALL add (g3b,'swdifnir', swdifnir ,lpost=.FALSE., longname='fraction of diffuse NIR'  &
         ,units='' ,contnorest=.true.)
    CALL add (g3b,'swvis',    swvis    ,lpost=.FALSE., longname='net surface visible'      &
         ,units='W/m**2' ,contnorest=.true.)
    CALL add (g3b,'swdifvis', swdifvis ,lpost=.FALSE., longname='fraction of diffuse visible' &
      ,units='' ,contnorest=.true.)
    CALL add (g3b,'swnirac',  swnirac,code= 79,laccu=.TRUE. ,longname='net surface NIR flux acc. ' &
            ,units='W/m**2'    ,contnorest=.true.)
    CALL add (g3b,'swdifnirac',swdifnirac,code= 80,laccu=.TRUE. ,longname='fraction of diffuse NIR acc. '    &
  ,units='W/m**2'    ,contnorest=.true.)
    CALL add (g3b,'swvisac',  swvisac,code= 81,laccu=.TRUE. ,longname='net surface visible flux acc. '      &
  ,units='W/m**2'    ,contnorest=.true.)
    CALL add (g3b,'swdifvisac',swdifvisac,code= 82,laccu=.TRUE. ,longname='fraction of diffuse visibleacc. '&
  ,units='W/m**2'    ,contnorest=.true.)

    CALL add (g3b,'trfliac',  trfliac  ,code= 91,laccu=.TRUE. ,lpost=lwarnGE2,longname='LW flux over ice (+ downward)'          ,units='W/m**2'   )
    CALL add (g3b,'trflwac',  trflwac  ,code= 92,laccu=.TRUE. ,lpost=lwarnGE2,longname='LW flux over water  (+ downward)'       ,units='W/m**2'   )
    CALL add (g3b,'trfllac',  trfllac  ,code= 93,laccu=.TRUE. ,lpost=lwarnGE2,longname='LW flux over land (+ downward)'         ,units='W/m**2'   )
    CALL add (g3b,'sofliac',  sofliac  ,code= 94,laccu=.TRUE. ,lpost=lwarnGE2,longname='SW flux over ice (+ downward)'          ,units='W/m**2'   )
    CALL add (g3b,'soflwac',  soflwac  ,code= 95,laccu=.TRUE. ,lpost=lwarnGE2,longname='SW flux over water (+ downward)'        ,units='W/m**2'   )
    CALL add (g3b,'sofllac',  sofllac  ,code= 96,laccu=.TRUE. ,lpost=lwarnGE2,longname='SW flux over land (+ downward)'         ,units='W/m**2'   )
    CALL add (g3b,'friac',    friac    ,code= 97,laccu=.TRUE. ,longname='ice cover (fraction of grid box)'                         )
    CALL add (g3b,'tsi',      tsi      ,code=102              ,longname='surface temperature of ice'             ,units='K'        )
    CALL add (g3b,'tsw',      tsw      ,code=103,lmiss=.TRUE. ,longname='surface temperature of water'           ,units='K'        )
    CALL add (g3b,'ustri',    ustri    ,code=104              ,lpost=lwarnGE2,longname='zonal      wind stress over ice'        ,units='Pa'       )
    CALL add (g3b,'vstri',    vstri    ,code=105              ,lpost=lwarnGE2,longname='meridional wind stress over ice'        ,units='Pa'       )
    CALL add (g3b,'ustrw',    ustrw    ,code=106              ,longname='zonal      wind stress over water'      ,units='Pa'       )
    CALL add (g3b,'vstrw',    vstrw    ,code=107              ,longname='meridional wind stress over water'      ,units='Pa'       )
    CALL add (g3b,'ustrl',    ustrl    ,code=108              ,lpost=lwarnGE2,longname='zonal      wind stress over land'       ,units='Pa'       )
    CALL add (g3b,'vstrl',    vstrl    ,code=109              ,lpost=lwarnGE2,longname='meridional wind stress over land'       ,units='Pa'       )
    CALL add (g3b,'ahfliac',  ahfliac  ,code=110,laccu=.TRUE. ,lpost=lwarnGE2,longname='latent heat flux over ice (+ downward)' ,units='W/m**2'   )
    CALL add (g3b,'ahflwac',  ahflwac  ,code=111,laccu=.TRUE. ,lpost=lwarnGE2,longname='latent heat flux over water (+ downward)'  ,units='W/m**2'   )
    CALL add (g3b,'ahfllac',  ahfllac  ,code=112,laccu=.TRUE. ,lpost=lwarnGE2,longname='latent heat flux over land (+ downward)',units='W/m**2'   )
    CALL add (g3b,'evapiac',  evapiac  ,code=113,laccu=.TRUE. ,lpost=lwarnGE2,longname='evaporation over ice'                   ,units='kg/m**2s' )
    CALL add (g3b,'evapwac',  evapwac  ,code=114,laccu=.TRUE. ,lpost=lwarnGE2,longname='evaporation over water'                 ,units='kg/m**2s' )
    CALL add (g3b,'evaplac',  evaplac  ,code=115,laccu=.TRUE. ,lpost=lwarnGE2,longname='evaporation over land'                  ,units='kg/m**2s' )
    CALL add (g3b,'az0i',     az0i     ,code=116              ,lpost=lwarnGE2,longname='roughness length over ice'              ,units='m'        )
    CALL add (g3b,'az0w',     az0w     ,code=117              ,lpost=lwarnGE2,longname='roughness length over water'            ,units='m'        )
    CALL add (g3b,'az0l',     az0l     ,code=118              ,lpost=lwarnGE2,longname='roughness length over land'             ,units='m'        )
    CALL add (g3b,'ahfsiac',  ahfsiac  ,code=119,laccu=.TRUE. ,lpost=lwarnGE2,longname='sensible heat flux over ice (+ downward)'  ,units='W/m**2'   )
    CALL add (g3b,'ahfswac',  ahfswac  ,code=120,laccu=.TRUE. ,lpost=lwarnGE2,longname='sensible heat flux over water (+ downward)',units='W/m**2'   )
    CALL add (g3b,'ahfslac',  ahfslac  ,code=121,laccu=.TRUE. ,lpost=lwarnGE2,longname='sensible heat flux over land (+ downward)' ,units='W/m**2'   )
    CALL add (g3b,'alsoi',    alsoi    ,code=122              ,lpost=lwarnGE2,longname='albedo of ice'                                            )
    CALL add (g3b,'alsow',    alsow    ,code=123              ,lpost=lwarnGE2,longname='albedo of water'                                          )
    CALL add (g3b,'alsol',    alsol    ,code=124              ,lpost=lwarnGE2,longname='albedo of land'                                           )
    CALL add (g3b,'ahfice',   ahfice   ,code=125              ,lpost=lwarnGE2,longname='conductive heat flux'                   ,units='W/m**2'   )
    CALL add (g3b,'qres',     qres     ,code=126              ,lpost=lwarnGE2,longname='residual heat flux for melting sea ice' ,units='W/m**2'   )
    CALL add (g3b,'alake',    alake    ,code=127,lpost=.FALSE.,longname='lake fraction of grid box'                                )
    CALL add (g3b,'rintop',   rintop   ,code=128,lpost=.FALSE.,longname='low level inversion      '                                )
    CALL add (g3b,'geosp',    geosp    ,code=129              ,longname='surface geopotential (orography)'       ,units='m**2/s**2')
    !               stp                      130                         temperature                                     K
    !                                        131                         u-velocity                                      m/s
    !                                        132                         v-velocity                                      m/s
    !                                        133                         specific humidity                               kg/kg
    CALL add (g3b,'aps',      aps      ,code=134              ,longname='surface pressure'                       ,units='Pa'       )
    !                                        135                         vertical velocity                               Pa/s
    CALL add (g3b,'acdnc',    acdnc    ,code=136,lpost=.FALSE.,longname='cloud droplet number concentration'     ,units='1/m**3'   )
    CALL add (g3b,'apmeb',    apmeb    ,code=137,laccu=.TRUE. ,longname='vert.integr.tendencies of water',bits=24,units='kg/m**2s' )
    !                         svo            138                         vorticity                                       1/s
    CALL add (g3b,'tslm1',    tslm1    ,code=139,lpost=lwarnGE2,longname='surface temperature of land'            ,units='K'        )
    CALL add (g3b,'ws',       ws       ,code=140              ,longname='soil wetness'                           ,units='m'        )
    CALL add (g3b,'sn',       sn       ,code=141              ,longname='snow depth'                             ,units='m'        )
    CALL add (g3b,'aprl',     aprl     ,code=142,laccu=.TRUE. ,longname='large scale precipitation'      ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'aprc',     aprc     ,code=143,laccu=.TRUE. ,longname='convective  precipitation'      ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'aprs',     aprs     ,code=144,laccu=.TRUE. ,longname='snow fall'                      ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'vdis',     vdis     ,code=145,laccu=.TRUE. ,longname='boundary layer dissipation'             ,units='W/m**2'   )
    CALL add (g3b,'ahfs',     ahfs     ,code=146,laccu=.TRUE. ,longname='sensible heat flux (+ downward)'        ,units='W/m**2'   )
    CALL add (g3b,'ahfl',     ahfl     ,code=147,laccu=.TRUE. ,longname='latent heat flux (+ downward)'          ,units='W/m**2'   )
    !                                        148                         streamfunction                                  m**2/s
    !                                        149                         velocity potential                              m**2/s
    CALL add (g3b,'xivi',     xivi     ,code=150,laccu=.TRUE. ,longname='vertically integrated cloud ice'        ,units='kg/m**2'  )
    !                                        151                         mean sea level pressure                         Pa
    !                         stp(20)        152                         log surface pressure
    !                         xl             153                         cloud water                                     kg/kg
    !                         xi             154                         cloud ice                                       kg/kg
    !                         sd             155                         divergence                                      1/s
    !                                        156                         geopotential height                             gpm
    CALL add (g3b,'relhum',   relhum   ,code=157              ,longname='relative humidity'                   )
    !                                        158                         tendency of surface pressure                    Pa/s
    CALL add (g3b,'wind10w', wind10w   ,code=159,lpost=.FALSE.,longname='10m windspeed over water'               ,units='m/s'      )
    CALL add (g3b,'runoff',  runoff    ,code=160,laccu=.TRUE. ,longname='surface runoff and drainage'    ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'drain',   drain     ,code=161,laccu=.TRUE. ,longname='drainage'                       ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'aclc',    aclc      ,code=162,lpost=.FALSE.,longname='cloud cover'                                              )
    CALL add (g3b,'aclcv',   aclcv     ,code=163,lpost=.FALSE.,longname='total cloud cover'                                        )
    CALL add (g3b,'aclcov',  aclcov    ,code=164,laccu=.TRUE. ,longname='total cloud cover'                                        )
    CALL add (g3b,'u10',     u10       ,code=165              ,longname='10m u-velocity'    ,units='m/s'      )
    CALL add (g3b,'v10',     v10       ,code=166              ,longname='10m v-velocity'    ,units='m/s'      )
    CALL add (g3b,'temp2',   temp2     ,code=167              ,longname='2m temperature'    ,units='K'        )
    CALL add (g3b,'dew2',    dew2      ,code=168       ,longname='2m dew point temperature' ,units='K'        )
    CALL add (g3b,'tsurf',   tsurf     ,code=169,laccu=.TRUE. ,longname='surface temperature'                    ,units='K'        )
    CALL add (g3b,'xvar',    xvar      ,code=170,lpost=.FALSE.,longname='variance of total water amount qv+qi+ql',units='kg/kg'    )
    CALL add (g3b,'wind10',  wind10    ,code=171,laccu=.TRUE. ,longname='10m windspeed'     ,units='m/s'      )
    CALL add (g3b,'slm',     slm       ,code=172              ,longname='land sea mask (1. = land, 0. = sea/lakes)'                )
    CALL add (g3b,'az0',     az0       ,code=173,lpost=.FALSE.,longname='roughness length'                       ,units='m'        )
    CALL add (g3b,'alb',     alb       ,code=174,lpost=.FALSE.,longname='surface background albedo'                                )
    CALL add (g3b,'albedo',  albedo    ,code=175              ,longname='surface albedo'                                           )
    CALL add (g3b,'srads',   srads     ,code=176,laccu=.TRUE. ,longname='net surface solar radiation (+ downward)'   ,units='W/m**2'   )
    CALL add (g3b,'trads',   trads     ,code=177,laccu=.TRUE. ,longname='net surface thermal radiation (+ downward)' ,units='W/m**2'   )
    CALL add (g3b,'srad0',   srad0     ,code=178,laccu=.TRUE. ,longname='net top solar radiation (+ downward)'       ,units='W/m**2'   )
    CALL add (g3b,'trad0',   trad0     ,code=179,laccu=.TRUE. ,longname='top thermal radiation (OLR) (+ downward)'   ,units='W/m**2'   )
    CALL add (g3b,'ustr',    ustr      ,code=180,laccu=.TRUE. ,longname='u-stress'                               ,units='Pa'       )
    CALL add (g3b,'vstr',    vstr      ,code=181,laccu=.TRUE. ,longname='v-stress'                               ,units='Pa'       )
    CALL add (g3b,'evap',    evap      ,code=182,laccu=.TRUE. ,longname='evaporation'                    ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'xskew',   xskew     ,code=183,lpost=.FALSE.,longname='skewness of total water amount qv+qi+ql'                  )
    CALL add (g3b,'srad0d',  srad0d    ,code=184,laccu=.TRUE. ,longname='top incoming solar radiation'           ,units='W/m**2'   )
    CALL add (g3b,'srafs',   srafs     ,code=185,laccu=.TRUE. ,longname='net surf. solar radiation   (clear sky)',units='W/m**2'   )
    CALL add (g3b,'trafs',   trafs     ,code=186,laccu=.TRUE. ,longname='net surf. thermal radiation (clear sky)',units='W/m**2'   )
    CALL add (g3b,'sraf0',   sraf0     ,code=187,laccu=.TRUE. ,longname='net top solar radiation     (clear sky)',units='W/m**2'   )
    CALL add (g3b,'traf0',   traf0     ,code=188,laccu=.TRUE. ,longname='net top thermal radiation   (clear sky)',units='W/m**2'   )
    !                                        189                         surface solar cloud forcing                     W/m**2
    !                                        190                         surface thermal cloud forcing                   W/m**2
    !                                        191                         SW top cloud forcing (178-187)                  W/m**2
    !                                        192                         LW top cloud forcing (179-188)                  W/m**2
    CALL add (g3b,'wl',      wl        ,code=193              ,longname='skin reservoir content'                 ,units='m'        )
    CALL add (g3b,'slf',     slf       ,code=194,lpost=.FALSE.,longname='sea land fraction (1. = land, 0. = sea/lakes)'            )
    CALL add (g3b,'ustrgw',  ustrgw    ,code=195,lpost=.FALSE.,longname='u-gravity wave stress'                  ,units='Pa'       )
    CALL add (g3b,'vstrgw',  vstrgw    ,code=196,lpost=.FALSE.,longname='v-gravity wave stress'                  ,units='Pa'       )
    CALL add (g3b,'vdisgw',  vdisgw    ,code=197              ,longname='gravity wave dissipation'               ,units='W/m**2'   )
    CALL add (g3b,'vgrat',   vgrat     ,code=198,lpost=.FALSE.,longname='vegetation ratio'                                         )
    CALL add (g3b,'orostd',  orostd    ,code=199,lpost=.FALSE.,longname='orographic standard deviation'          ,units='m'        )
    CALL add (g3b,'vlt',     vlt       ,code=200,lpost=.FALSE.,longname='leaf area index'                                          )
    CALL add (g3b,'t2max',   t2max     ,code=201,reset=-99.0_dp   ,longname='maximum 2m temperature',units='K'     )
    CALL add (g3b,'t2min',   t2min     ,code=202,reset=999.0_dp   ,longname='minimum 2m temperature',units='K'     )
    CALL add (g3b,'srad0u',  srad0u    ,code=203,laccu=.TRUE. ,longname='top solar radiation upward'             ,units='W/m**2'   )
    CALL add (g3b,'sradsu',  sradsu    ,code=204,laccu=.TRUE. ,longname='surface solar radiation upward'         ,units='W/m**2'   )
    CALL add (g3b,'tradsu',  tradsu    ,code=205,laccu=.TRUE. ,longname='surface thermal radiation upward'       ,units='W/m**2'   )
    CALL add (g3b,'grndflux',grndflux  ,code=206,laccu=.TRUE. ,longname='surface ground heat flux (+ upward)'    ,units='W/m**2'   )
    CALL add (g3b,'tsoil',   tsoil     ,code=207              ,longname='deep soil temperatures',leveltype=BELOWSUR,units='K'      )
    CALL add (g3b,'ahfcon',  ahfcon    ,code=208,laccu=.TRUE. ,longname='conductive heat flux through ice'       ,units='W/m**2'   )
    CALL add (g3b,'ahfres',  ahfres    ,code=209,laccu=.TRUE. ,longname='melting of ice'                         ,units='W/m**2'   )
    CALL add (g3b,'seaice',  seaice    ,code=210              ,longname='ice cover (fraction of 1-SLM)'                            )
    CALL add (g3b,'siced',   siced     ,code=211              ,longname='ice depth'                              ,units='m'        )
    CALL add (g3b,'forest',  forest    ,code=212,lpost=.FALSE.,longname='forest fraction'                                          )
    CALL add (g3b,'gld',     gld       ,code=213              ,longname='glacier depth'                          ,units='m'        )
    CALL add (g3b,'sni',     sni       ,code=214              ,longname='water equivalent of snow on ice'        ,units='m'        )
    CALL add (g3b,'rogl',    rogl      ,code=215,laccu=.TRUE.,lpost=.FALSE. ,longname='glacier runoff'           ,units='kg/m**2s' )
    CALL add (g3b,'wimax',   wimax     ,code=216,reset=-99.0_dp   ,longname='maximum 10m-wind speed',units='m/s'  )
    CALL add (g3b,'topmax',  topmax    ,code=217,reset=99999.0_dp ,longname='maximum height of convective cloud tops',units='Pa'   )
    CALL add (g3b,'snmel',   snmel     ,code=218,laccu=.TRUE. ,longname='snow melt'                              ,units='kg/m**2s' )
!!!    CALL add (g3b,'runtoc',  runtoc    ,code=219,lpost=lwarnGE2,laccu=.TRUE.,               &
!!!                 longname='surface runoff into ocean',bits=24,units='kg/m**2s' )
!!!    CALL add (g3b,'runlnd',  runlnd    ,code=220,lpost=lwarnGE2,                            &
!!!                 longname='surface runnof not running into ocean'  ,units='kg/m**2s' )
    CALL add (g3b,'runtoc',  runtoc    ,code=219,lpost=.TRUE.,laccu=.TRUE.,               &
                 longname='surface runoff into ocean',bits=24,units='kg/m**2s' )
!!!    CALL add (g3b,'runlnd',  runlnd    ,code=220,lpost=.TRUE.,                            &
!!!                 longname='surface runnof not running into ocean'  ,units='kg/m**2s' )
    CALL add (g3b,'runlnd',  runlnd    ,lpost=.FALSE.,                            &
                 longname='surface runnof not running into ocean'  ,units='kg/m**2s' )
    CALL add (g3b,'apmegl',  apmegl    ,code=221,laccu=.TRUE. ,longname='P-E over land ice'              ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'snacl',   snacl     ,code=222,laccu=.TRUE. ,longname='snow accumulation over land'            ,units='kg/m**2s' )
    CALL add (g3b,'aclcac',  aclcac    ,code=223,laccu=.TRUE. ,longname='cloud cover'                                              )
    CALL add (g3b,'tke',     tke       ,code=224,lpost=.false.,longname='turbulent kinetic energy'               ,units='m**2/s**2')
    CALL add (g3b,'tkem1',   tkem1     ,code=225,lpost=.false.,longname='turbulent kinetic energy (t-1)'         ,units='m**2/s**2')

    CALL add (g3b,'fao',     fao       ,code=226,lpost=lwarnGE2,longname='FAO data set (soil data flags 0...5.)'                    )
    CALL add (g3b,'rgcgn',   rgcgn     ,code=227,lpost=lwarnGE2,longname='volumetric heat capacity of soil and land ice',            &
                   units='J/(m**3*K)')
    CALL add (g3b,'sodif',   sodif     ,code=228,lpost=lwarnGE2,longname='diffusivity  of soil and land ice'      ,units='m**2/s'   )
!!!    CALL add (g3b,'wsmx',    wsmx      ,code=229              ,longname='field capacity of soil'               ,units='m'        )
    CALL add (g3b,'wsmx',    wsmx      ,code=229              ,longname='rooting depth of soil'                   ,units='m'        )
    CALL add (g3b,'qvi',     qvi       ,code=230,laccu=.TRUE. ,longname='vertically integrated water vapor'      ,units='kg/m**2'  )
    CALL add (g3b,'xlvi',    xlvi      ,code=231,laccu=.TRUE. ,longname='vertically integrated cloud water'      ,units='kg/m**2'  )
    CALL add (g3b,'glac',    glac      ,code=232              ,longname='fraction of land covered by glaciers'                     )
    CALL add (g3b,'snc',     snc       ,code=233              ,longname='snow depth at the canopy'               ,units='m'        )
    CALL add (g3b,'rtype',   rtype     ,code=234,lpost=.FALSE.,longname='type of convection 0...3.'                                )
    !                                        259                         windspeed (sqrt(u**2+v**2))
    !                                        260                         total precipitation (142+143)
    !                                        261                         total top radiation (178+179)
    !                                        262                         total surface radiation (176+177)
    !                                        263                         net surface heat flux 
    !                                                                    (146+147+176+177-C*218-208*fice-209); C=3.345E5*fland
    !                                        264                         total surface water (142+143+182-160-221)
    !
    CALL add (g3b,'abso4',    abso4    ,code=235,laccu=.TRUE. ,                  &
         longname='anthropogenic sulfur burden'  ,units='kg/m**2' ,contnorest=.true.)
    CALL add (g3b,'so4nat',   so4nat   ,lpost=.FALSE.,leveltype=HYBRID,          &
         longname='natural sulfate'  ,units='kg/kg' ,contnorest=.true.)
    CALL add (g3b,'so4all',   so4all   ,lpost=.FALSE.,leveltype=HYBRID,          &
         longname='total sulfur'  ,units='kg/m**2' ,contnorest=.true.)
    CALL add (g3b,'ao3',      ao3      ,code=236,leveltype=HYBRID,               &
         longname='ipcc ozone'  ,units='kg/kg' ,contnorest=.true.)
    CALL add (g3b,'tropo',    tropo    ,code=237 ,contnorest=.true.,             &
         longname='WMO defined tropopause height'          ,units='Pa'       )

    CALL add (g3b,'obstsw',   obstsw,       code=238 ,lmiss=.TRUE.,            &
         longname='bulk surface temperature of water',units='K')
         
    CALL add (g3b,'obswsb',   obswsb,     code=239 ,lpost=lwarnGE2,lmiss=.TRUE.,              &
         longname='bulk salinity of water',units='PSU')
    CALL add (g3b,'lclass', lclass, code=254,lpost=lwarnGE2,lmiss=.TRUE.,          &
         longname='land cover class, 1: land, 2: ocean, 3: lake, 4: glacier')
    IF (lasia) THEN
      CALL add (g3b,'tinertia', tinertia, code=15,lpost=lwarnGE2,lmiss=.TRUE.,      &
           longname='thermal inertia',units='J m**-2 K**-1 s**-0.5 (tiu)')
    END IF           
    !
    !  variables for fractional surface coverage
    !
    CALL add (g3b,'ahfli',    ahfli    ,lpost=.FALSE.,longname='LE over ice (+ downward)',units='W/m**2')
    CALL add (g3b,'ahflw',    ahflw    ,lpost=.FALSE.,longname='LE over water (+ downward)',units='W/m**2')
    CALL add (g3b,'ahfll',    ahfll    ,lpost=.FALSE.,longname='LE over land  (+ downward)',units='W/m**2')
    CALL add (g3b,'obox_mask', obox_mask,   code=16,lpost=lwarnGE2,longname='ocean boxes nudging mask',lmiss=.TRUE.)
    CALL add (g3b,'ocnmask', ocnmask,   code=17,lpost=lwarnGE2,longname='3-D ocean mask',lmiss=.TRUE.)
    CALL add (g3b,'sitmask', sitmask,   code=18,lpost=lwarnGE2,longname='sit mask')
    CALL add (g3b,'fluxres',  fluxres  ,code=19,lpost=lwarnGE2,longname='residual flux for melting (+ downward)',units='W/m**2')
    CALL add (g3b,'bathy', bathy, lpost=.FALSE.,                             &
         longname='orography (mountain ridge + ocean bathymetry)',           &
         units='m')                                                                          
!!!    CALL add (g3b,'bathy', bathy, lpost=.FALSE.,lmiss=.TRUE.,                &
!!!         longname='orography (mountain ridge + ocean bathymetry)',           &
!!!         units='m')                                                                          
    CALL add (g3b,'sitwlvl', sitwlvl, lpost=.FALSE.,lmiss=.TRUE.,            &
         longname='current water level (ice/water interface)' ,              &
         units='m')
    CALL add (g3b,'oceanid', oceanid,lpost=.FALSE.,lmiss=.TRUE.,             &
         longname='ocean basin ID')
!!!
!!!  new ocn porosity model
!!!
    CALL add (g3b,'ocn_oromea', ocn_oromea,lpost=.FALSE.,lmiss=.TRUE.,longname='ocean Mean orography',units='m')
    CALL add (g3b,'ocn_wf', ocn_wf,lpost=.FALSE.,lmiss=.TRUE.,longname='water fraction',units='fractional, [0,1]')
    CALL add (g3b,'ocn_bsl_oromea', ocn_bsl_oromea,lpost=.FALSE.,lmiss=.TRUE.,longname='mean depth over water fraction',units='m, + upward')
    CALL add (g3b,'ocn_divzmea', ocn_divzmea,lpost=.FALSE.,lmiss=.TRUE.,longname='mean divergence over water fraction',units='m/gd2')
    CALL add (g3b,'ocn_divzmin', ocn_divzmin,lpost=.FALSE.,lmiss=.TRUE.,longname='min divergence over water fraction',units='m/gd2')
    CALL add (g3b,'ocn_divzmax', ocn_divzmax,lpost=.FALSE.,lmiss=.TRUE.,longname='max divergence over water fraction',units='m/gd2')
    CALL add (g3b,'ocn_bsl_oropic', ocn_bsl_oropic,lpost=.FALSE.,lmiss=.TRUE.,longname='Orographic peak elevation over water fraction',units='m')
    CALL add (g3b,'ocn_bsl_oroval', ocn_bsl_oroval,lpost=.FALSE.,lmiss=.TRUE.,longname='Orographic valley elevation over water fraction',units='m')
    CALL add (g3b,'ocn_por', ocn_por,lpost=.FALSE.,lmiss=.TRUE.,longname='porosity of each ocn cell for water',leveltype=OCEAN,units='fractional, [0,1]')
    CALL add (g3b,'ocn_porx', ocn_porx,lpost=.FALSE.,lmiss=.TRUE.,longname='porosity of each ocn cell for water going in x dir',leveltype=OCEAN,units='fractional, [0,1]')
    CALL add (g3b,'ocn_pory', ocn_pory,lpost=.FALSE.,lmiss=.TRUE.,longname='porosity of each ocn cell for water going in y dir',leveltype=OCEAN,units='fractional, [0,1]')

    IF (lsit.OR.locn) THEN
  !
  !  prognostic variables 
      CALL add (g3b,'sitzsi', sitzsi,   code=42 ,lmiss=.TRUE.,               &
           longname='water equivalent of snow and ice over water' ,          &
           leveltype=SNOWICE,units='m')
      CALL add (g3b,'sitsilw', sitsilw, code=43 ,lmiss=.TRUE.,               &
           longname='liquid water in snow and ice'    ,                      &
           leveltype=SNOWICE,units='m')  
      CALL add (g3b,'sittsi', sittsi,   code=44 ,lmiss=.TRUE.,               &
           longname='temperature of snow and ice'     ,                      &
           leveltype=SNOWICE2,units='K')
      CALL add (g3b,'sitwt', sitwt,     code=45 ,lmiss=.TRUE.,               &
           longname='potential water temperature' ,leveltype=WATER,units='K')     
      CALL add (g3b,'sitws', sitws,     code=46 ,lmiss=.TRUE.,               &
           longname='salinity'         ,leveltype=WATER,units='0/00') 
      CALL add (g3b,'sitwu', sitwu,     code=47 ,lmiss=.TRUE.,               &
           longname='water current U'    ,leveltype=WATER,units='m/s') 
      CALL add (g3b,'sitwv', sitwv,     code=48 ,lmiss=.TRUE.,               &
           longname='water current V'    ,leveltype=WATER,units='m/s')   
      CALL add (g3b,'sitww', sitww,     code=49 ,lpost=locn,lmiss=.TRUE.,               &
           longname='water subsidence W (down plus)'    ,leveltype=WATER,units='m/s')   
      CALL add (g3b,'sitwp', sitwp,     code=50 ,lpost=locn,lmiss=.TRUE.,               &
           longname='ocean water pressure'    ,leveltype=WATER,units='Pa')
      CALL add (g3b,'sitwtke', sitwtke, code=39 ,lpost=.TRUE.,lmiss=.TRUE.,               &
           longname='water TKE'          ,                                   &
           leveltype=WATERF,units='m**2/s**2') 
  !  diagnostic variables
      CALL add (g3b,'sitwlmx', sitwlmx, code=40,lpost=.TRUE., lmiss=.TRUE.,         &
           longname='mixing length in water'       ,                         &
           leveltype=WATERF,units='m')  
      CALL add (g3b,'sitwldisp', sitwldisp, code=41,lpost=.TRUE., lmiss=.TRUE.,      &
           longname='dissipation length in water' ,                          &
           leveltype=WATERF,units='m') 
      CALL add (g3b,'sitwkm', sitwkm, lpost=.FALSE.,lmiss=.TRUE.,            &
           longname='momentum diffusivity in water'     ,                    &
           leveltype=WATERF,units='m**2/s') 
      CALL add (g3b,'sitwkh', sitwkh, code=51,lmiss=.TRUE.,                  &
           longname='heat diffusivity in water'         ,                    &
           leveltype=WATERF,units='m**2/s') 
!      CALL add (g3b,'wrho', wrho, code=52,lpost=lwarnGE2,lmiss=.TRUE.,                      &
      CALL add (g3b,'wrho', wrho, code=52,lmiss=.TRUE.,                      &
           longname='1000-m potential water density',                              &
           leveltype=WATER,units='kg/m**3')   
!
      CALL add (g3b,'wtfn', wtfn, code=53,laccu=.TRUE.,lmiss=.TRUE.,             &
           longname='mean temperature flux nudging at each level (positive into ocean)' ,  &
           leveltype=WATER, units='W/m2')
      CALL add (g3b,'wsfn', wsfn, code=54,laccu=.TRUE.,lmiss=.TRUE.,             &
           longname='mean salinity flux nudging at each level (positive into ocean)' ,     &
           leveltype=WATER, units='PSU*m/s')
      CALL add (g3b,'wtfns', wtfns, code=55,laccu=.TRUE.,lmiss=.TRUE.,           &
           longname='mean temperature flux nudging for an entire water column (positive into ocean)' ,  &
           units='W/m2')
      CALL add (g3b,'wsfns', wsfns, code=56,laccu=.TRUE.,lmiss=.TRUE.,           &
           longname='mean salinity flux nudging for an entire water column (positive into ocean)' ,  &
           units='PSU*m/s')
      CALL add (g3b,'awtfl', awtfl, code=57, laccu=.FALSE.,lmiss=.TRUE.,  &
           longname='advected temperature flux at each level (positive into ocean)' ,  &
           leveltype=WATER, units='W/m2')
      CALL add (g3b,'awsfl', awsfl, code=58, laccu=.FALSE.,lmiss=.TRUE.,  &
           longname='advected salinity flux at each level (positive into ocean)' , &
           leveltype=WATER, units='PSU*m/s')
      CALL add (g3b,'awufl', awufl,   laccu=.FALSE.,lpost=.FALSE.,lmiss=.TRUE.,  &
           longname='advected u flux at each level (positive into ocean)' ,  &
           leveltype=WATER, units='m2/s2')
      CALL add (g3b,'awvfl', awvfl,   laccu=.FALSE.,lpost=.FALSE.,lmiss=.TRUE.,  &
           longname='advected v flux at each level (positive into ocean)' ,  &
           leveltype=WATER, units='m2/s2')
      CALL add (g3b,'awtkefl', awtkefl, laccu=.FALSE.,lpost=.FALSE.,lmiss=.TRUE.,  &
           longname='advected tke flux at each level (positive into ocean)' ,    &
           leveltype=WATERF, units='m3/s3')                                       
      CALL add (g3b,'ctfreez2', ctfreez2,lpost=.FALSE.,lmiss=.TRUE.,                         &
           longname='ref water freezing temperature',units='K')
!
      CALL add (g3b,'fluxiw',fluxiw, code=59, lpost=.TRUE.,lmiss=.TRUE.,                       &
           longname='below ice, over water surface heat flux (+ upward)',units='W/m2')         
!
      CALL add (g3b,'cc', cc, code=60, lmiss=.TRUE.,                         &
           longname='cold content (energy need to melt snow and ice) of snow/ice column',  &
           bits=24,units='J/m**2')
      CALL add (g3b,'hc', hc, code=61, lmiss=.TRUE.,                                       &
           bits=24,longname='heat content of water column above tmelt' ,units='J/m**2')
!!!      CALL add (g3b,'engw', engw, lmiss=.TRUE.,laccu=.TRUE.,lpost=.FALSE.,                        &
!!!           longname='mean net heat flux over water fraction (+ downward)',units='W/m**2')
!!!      CALL add (g3b,'engw2', engw2, lmiss=.TRUE.,laccu=.TRUE.,lpost=.FALSE.,                      &
!!!           longname='mean snowfall-corrected net heat flux over water fraction (+ downward)', &
!!!           units='W/m**2')
      CALL add (g3b,'engwac', engwac, code=62, lmiss=.TRUE.,                               &
           longname='accumulated energy over water fraction (+ downward)',units='J/m**2')

      CALL add (g3b,'pme',pme, lpost=.FALSE.,lmiss=.TRUE.,                             &
           longname='precipitation minus evaporation (+ downward)',units='m/s')                
      CALL add (g3b,'pme2',pme2, code=63, lpost=.TRUE.,lmiss=.TRUE.,                       &
           longname='ice-corrected net fresh water into ocean (+ downward)',units='m/s')       
      CALL add (g3b,'saltwac',saltwac, code=64, lmiss=.TRUE.,                              &
           longname='accumulated salt into water (+ downward)',units='PSU*m')       
      CALL add (g3b,'sc', sc, code=65, lmiss=.TRUE.,                                       &
           longname='salinity content of water column' ,bits=24,units='PSU*m')
  !  Note tha code number might not exceed 255 bjt
  !
  !  variables required for skin SST and sit calculation
  !
      IF (lgodas) THEN     
        CALL add (g3b,'obswt', obswt, code=66, lmiss=.TRUE.,                       &
             longname='obs potential water temperature' ,leveltype=WATER,units='K')
        CALL add (g3b,'obsws', obsws, code=67, lmiss=.TRUE.,                       &
             longname='obs salinity' ,leveltype=WATER,units='0/00')
        CALL add (g3b,'obswu', obswu, code=68, lmiss=.TRUE.,                       &
             longname='obs u current' ,leveltype=WATER,units='m/s')
        CALL add (g3b,'obswv', obswv, code=69, lmiss=.TRUE.,                       &
             longname='obs v current' ,leveltype=WATER,units='m/s')
      END IF
    END IF
#if defined (G3BTEST)            
    IF (locn) THEN
      CALL add (g3b,'ocnp', ocnp, code=70, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,                         &
           longname='ocean water pressure' ,leveltype=OCEAN,repr=OGAUSSIAN,units='Pa')                         
      CALL add (g3b,'ocnt', ocnt, code=71, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='potential water temperature' ,leveltype=OCEAN,repr=OGAUSSIAN,units='K')                            
      CALL add (g3b,'ocns', ocns, code=72, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='salinity'         ,leveltype=OCEAN,repr=OGAUSSIAN,units='0/00')                          
      CALL add (g3b,'ocnu', ocnu, code=73, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='water current U'    ,leveltype=OCEAN,repr=OGAUSSIAN,units='m/s')                         
      CALL add (g3b,'ocnv', ocnv, code=74, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='water current V'    ,leveltype=OCEAN,repr=OGAUSSIAN,units='m/s')                         
      CALL add (g3b,'ocnw', ocnw, code=75, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='water subsidence W (down plus)',                                          &
           leveltype=OCEANF,repr=OGAUSSIAN,units='m/s')
      CALL add (g3b,'ocnkvm', ocnkvm, code=76, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,           &
           longname='vertical momentum diffusivity in water',                                  &
           leveltype=OCEANF,repr=OGAUSSIAN,units='m**2/s')                                                    
      CALL add (g3b,'ocnkvh', ocnkvh, lpost=.FALSE., lmiss=.TRUE.,                               &
           longname='vertical heat diffusivity in water',                                      &
           leveltype=OCEANF,repr=OGAUSSIAN,units='m**2/s')                                                    
      !!! CALL add (g3b,'ocnkhm', ocnkhm, code=77, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,           &
      !!!      longname='horizontal momentum diffusivity in water',                                &
      !!!      leveltype=OCEAN,repr=OGAUSSIAN,units='m**2/s')                                                    
      CALL add (g3b,'ocnkhm', ocnkhm, code=77, lpost=.TRUE., lmiss=.TRUE.,           &
           longname='horizontal momentum diffusivity in water',                                &
           leveltype=OCEAN,repr=OGAUSSIAN,units='m**2/s')                                                    
      CALL add (g3b,'ocnkhh', ocnkhh, lpost=.FALSE., lmiss=.TRUE.,                               &
           longname='horizontal heat diffusivity in water',                                    &
           leveltype=OCEAN,repr=OGAUSSIAN,units='m**2/s')                                                    
    END IF
#else
    IF (locn) THEN
      CALL add (g3b,'ocnp', ocnp, code=70, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,                         &
           longname='ocean water pressure' ,leveltype=OCEAN,units='Pa')                         
      CALL add (g3b,'ocnt', ocnt, code=71, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='potential water temperature' ,leveltype=OCEAN,units='K')                            
      CALL add (g3b,'ocns', ocns, code=72, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='salinity'         ,leveltype=OCEAN,units='0/00')                          
      CALL add (g3b,'ocnu', ocnu, code=73, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='water current U'    ,leveltype=OCEAN,units='m/s')                         
      CALL add (g3b,'ocnv', ocnv, code=74, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='water current V'    ,leveltype=OCEAN,units='m/s')                         
      CALL add (g3b,'ocnw', ocnw, code=75, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,             &
           longname='water subsidence W (down plus)',                                          &
           leveltype=OCEANF,units='m/s')
      CALL add (g3b,'ocnkvm', ocnkvm, code=76, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,           &
           longname='vertical momentum diffusivity in water',                                  &
           leveltype=OCEANF,units='m**2/s')                                                    
      CALL add (g3b,'ocnkvh', ocnkvh, lpost=.FALSE., lmiss=.TRUE.,                               &
           longname='vertical heat diffusivity in water',                                      &
           leveltype=OCEANF,units='m**2/s')                                                    
      !!! CALL add (g3b,'ocnkhm', ocnkhm, code=77, lpost=lwarnGE2orlocn_msg, lmiss=.TRUE.,           &
      !!!      longname='horizontal momentum diffusivity in water',                                &
      !!!      leveltype=OCEAN,units='m**2/s')                                                    
      CALL add (g3b,'ocnkhm', ocnkhm, code=77, lpost=.TRUE., lmiss=.TRUE.,           &
           longname='horizontal momentum diffusivity in water',                                &
           leveltype=OCEAN,units='m**2/s')                                                    
      CALL add (g3b,'ocnkhh', ocnkhh, lpost=.FALSE., lmiss=.TRUE.,                               &
           longname='horizontal heat diffusivity in water',                                    &
           leveltype=OCEAN,units='m**2/s')                                                    
    END IF
#endif    
                                                                                         
    IF (lsit.AND.locn) THEN                                                                    
      CALL add (g3b,'subfluxw',subfluxw, code=78, lpost=.TRUE.,lmiss=.TRUE.,                   &
           longname='water subsurface heat flux (+ upward)',units='W/m2')                      
      CALL add (g3b,'wsubsal',wsubsal, lpost=.FALSE.,lmiss=.TRUE.,                             &
           longname='water subsurface salinity flux (+ upward)',units='m*PSU/s')               
      CALL add (g3b,'sitwtb', sitwtb, lpost=.FALSE., lmiss=.TRUE.,                             &
           longname='acc. bulk surface water temperature' ,units='K*s')                        
      CALL add (g3b,'sitwub', sitwub, lpost=.FALSE., lmiss=.TRUE.,                             &
           longname='acc. bulk surface water current U'    ,units='m/s*s')                     
      CALL add (g3b,'sitwvb', sitwvb, lpost=.FALSE., lmiss=.TRUE.,                             &
           longname='bulk surface water current V'    ,units='m/s*s')                          
      CALL add (g3b,'sitwsb', sitwsb, lpost=.FALSE., lmiss=.TRUE.,                             &
           longname='bulk surface salinity'         ,units='PSU*s')
!            
      CALL add (g3b,'afluxiw',afluxiw,lpost=.FALSE.,lmiss=.TRUE.,       &
           longname='acc. below ice, surface ocean water heat flux (+ upward)',units='W/m2*s')
      CALL add (g3b,'apme2',apme2,lpost=.FALSE.,lmiss=.TRUE.,       &
           longname='acc. ice-corrected net fresh water into ocean (+ downward)',units='m/s*s')
      CALL add (g3b,'asubfluxw',asubfluxw,lpost=.FALSE.,lmiss=.TRUE.,       &
           longname='acc. subsurface ocean heat flux (+ upward)',units='W/m2*s')
      CALL add (g3b,'awsubsal',awsubsal,lpost=.FALSE.,lmiss=.TRUE.,       &
           longname='acc. subsurface ocean salinity flux (+ upward)',units='m*PSU/s*s')
    END IF
    
    ! Add fields not written to the output stream
    !
    CALL add (g3b,'dfluxs', dfluxs , lpost=.FALSE.,longname='d(flux)/dT over water (+ upward)',units='W/m**2/K')
!    CALL add (g3b,'trflw',  trflw  , lpost=.FALSE.,longname='net LW flux over water instantaneous',units='W/m**2')    
!    CALL add (g3b,'soflw',  soflw  , lpost=.FALSE.,longname='net SW flux over water instantaneous',units='W/m**2')
!    CALL add (g3b,'ahfsw',  ahfsw  , lpost=.FALSE.,longname='sensible heat flux over water instantaneous',units='W/m**2')


    CALL add (g3b,'tsl',     tsl,code=20, lpost=lwarnGE2,lmiss=.TRUE.,                         &
           longname='surface temperature of land'            ,units='K')
    CALL add (g3b,'tslm',    tslm,code=21, lpost=lwarnGE2,lmiss=.TRUE.,                        &
           longname='surface temperature of land'            ,units='K')

    CALL add (g3b,'emter',   emter   ,lpost=.FALSE. ,leveltype=HYBRID_H)
    CALL add (g3b,'trsol',   trsol   ,lpost=.FALSE. ,leveltype=HYBRID_H)
    CALL add (g3b,'emtef0',  emtef0  ,lpost=.FALSE. ,leveltype=HYBRID_H,contnorest=.true.)
    CALL add (g3b,'trsof0',  trsof0  ,lpost=.FALSE. ,leveltype=HYBRID_H,contnorest=.true.)
    CALL add (g3b,'emtef',   emtef   ,lpost=.FALSE. ,klev=2)
    CALL add (g3b,'trsof',   trsof   ,lpost=.FALSE. ,klev=2)
    CALL add (g3b,'tkem',    tkem    ,lpost=.FALSE.)
!   
    CALL add (g3b,'grndcapc',grndcapc,code=22, lpost=lwarnGE2,lmiss=.TRUE.,                    &
         longname='areal heat capacity of the uppermost ground layer (snow/ice/water/soil)',units='j/m**2/K')
    CALL add (g3b,'grndhflx',grndhflx,code=23, lpost=lwarnGE2,lmiss=.TRUE.,                    &
         longname='subsurface ground heat flux (+ upward)',units='W/m2')

    CALL add (g3b,'obsseaice',obsseaice,code=24, lpost=.TRUE.,lmiss=.TRUE.,longname='obs ice cover (fraction of 1-SLM)')    
    CALL add (g3b,'grndc',   grndc   ,lpost=.FALSE. ,leveltype=BELOWSUR)
    CALL add (g3b,'grndd',   grndd   ,lpost=.FALSE. ,leveltype=BELOWSUR)
    CALL add (g3b,'acvtype', acvtype ,lpost=.FALSE.)
    CALL add (g3b,'xtec',    xtec    ,lpost=.FALSE.)
    !
    ! variables for middle atmosphere only
    !
    IF (lmidatm) &
      CALL add (g3b,'aprflux',  aprflux ,lpost=.FALSE.)
    !
    !  variables for mixed layer ocean only
    !
    CALL add (g3b,'amlcorr',  amlcorr  ,lpost=.FALSE.)
    CALL add (g3b,'amlcorac', amlcorac ,code= 89,laccu=.TRUE.)
    CALL add (g3b,'amlheatac',amlheatac,code= 90,laccu=.TRUE.,lpost=.FALSE.)
    !
    !  variables for ocean coupling only
    !
    CALL add (g3b,'apmebco',  apmebco  ,lpost=.FALSE.,contnorest=.true.)
    CALL add (g3b,'rain',     rain     ,lpost=.FALSE.,contnorest=.true.)
    
    !
    !
    IF (lcouple) THEN
      CALL add (g3b,'awhea',    awhea   ,lpost=.FALSE.,longname='net surface heat flux into water (+ downward)',units='W/m**2*s?')
      CALL add (g3b,'awsol',    awsol   ,lpost=.FALSE.,longname='net solar radiation into water (+ downward)',units='W/m**2*s?')
      CALL add (g3b,'awfre',    awfre   ,lpost=.FALSE.,longname='net fresh water (P-E) into water (+ downward)',units='m/s*s?')
      CALL add (g3b,'awust',    awust   ,lpost=.FALSE.,longname='ustrw',units='Pa*s')
      CALL add (g3b,'awvst',    awvst   ,lpost=.FALSE.,longname='vstrw',units='Pa*s')
      CALL add (g3b,'awsta',    awsta   ,lpost=.FALSE.,longname='wind10',units='Pa*s')
      CALL add (g3b,'aicon',    aicon   ,lpost=.FALSE.,longname='conductive heat flux' ,units='W/m**2*s?'   )
      CALL add (g3b,'aiqre',    aiqre   ,lpost=.FALSE.,longname='residual heat flux for melting sea ice' ,units='W/m**2*s?'   )
      CALL add (g3b,'aifre',    aifre   ,lpost=.FALSE.,longname='net snowfall (P-E) over ice (+ downward)',units='m/s*s?')
      CALL add (g3b,'aiust',    aiust   ,lpost=.FALSE.,longname='ustri',units='Pa*s')
      CALL add (g3b,'aivst',    aivst   ,lpost=.FALSE.,longname='vstri',units='Pa*s')
    END IF
    
    IF (locn) THEN
      CALL add (g3b,'afluxs',afluxs,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. net surface heat flux into water (+ downward)',units='W/m**2*s')
      CALL add (g3b,'awsol2',awsol2,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. net solar radiation into water (+ downward)',units='W/m**2*s')
      CALL add (g3b,'apme',apme,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. net fresh water (P-E) (+ downward)',units='m/s*s')
      CALL add (g3b,'awust2',awust2,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. ustrw',units='Pa*s')
      CALL add (g3b,'awvst2',awvst2,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. vstrw',units='Pa*s')
      CALL add (g3b,'awsta2',awsta2,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. wind10',units='m/s*s')
      CALL add (g3b,'aicon2',aicon2,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. conductive heat flux' ,units='W/m**2*s'   )
      CALL add (g3b,'aiqre2',aiqre2,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. residual heat flux for melting sea ice' ,units='W/m**2*s'   )
      CALL add (g3b,'aiust2',aiust2,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. ustri',units='Pa*s')
      CALL add (g3b,'aivst2',aivst2,lpost=.FALSE.,lmiss=.TRUE.,longname='acc. vstri',units='Pa*s')
    END IF 

    CALL add (g3b,'ocu',  ocu    ,code= 83,                                             &
         longname='ocean eastw. velocity',bits=24,units='m/s',contnorest=.true. )
    CALL add (g3b,'ocv',  ocv    ,code= 84,                                             &
         longname='ocean northw. velocity',bits=24,units='m/s',contnorest=.true. )
    !
    !  variables for 200mb radiation
    !
    CALL add (g3b,'tradl',   tradl     ,code= 85,laccu=.TRUE.,units='W/m**2',contnorest=.true.   )
    CALL add (g3b,'sradl',   sradl     ,code= 86,laccu=.TRUE.,units='W/m**2',contnorest=.true.   )
    CALL add (g3b,'trafl',   trafl     ,code= 87,laccu=.TRUE.,units='W/m**2',contnorest=.true.   )
    CALL add (g3b,'srafl',   srafl     ,code= 88,laccu=.TRUE.,units='W/m**2',contnorest=.true.   )
    !
    !  variables for coupling with HD-model only
    !
    IF (lhd) THEN
      CALL add (g3b,'aros',     aros    ,lpost=.FALSE.)
      CALL add (g3b,'adrain',   adrain  ,lpost=.FALSE.)
!!!      CALL add (g3b,'disch',    disch   ,lpost=.FALSE.)
      CALL add (g3b,'apmecal',  apmecal ,lpost=.FALSE.)
    END IF
    CALL add (g3b,'disch',  disch    ,code=220,lpost=.TRUE.,                            &
             longname='surface runoff running into ocean'  ,units='m/s' )
    !
    !  variables for sso parametrization
    !
    CALL add (g3b,'oromea', oromea,lpost=.FALSE.)
    CALL add (g3b,'orosig', orosig,lpost=.FALSE.)
    CALL add (g3b,'orogam', orogam,lpost=.FALSE.)
    CALL add (g3b,'orothe', orothe,lpost=.FALSE.)
    CALL add (g3b,'oropic', oropic,lpost=.FALSE.)
    CALL add (g3b,'oroval', oroval,lpost=.FALSE.)

  END SUBROUTINE construct_g3b

  SUBROUTINE destruct_g3b

    CALL delete_stream (g3b)

  END SUBROUTINE destruct_g3b

END MODULE mo_memory_g3b
