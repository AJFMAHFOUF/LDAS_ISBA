SUBROUTINE SOLAR_ANGLE(IDAT,NSSSSS,PLAT,PLON,PMU0,PMU0M)
!-----------------------------------------------------------------------
!
!     Computation of solar zenith angle 
!     ---------------------------------
!
!     INPUT PARAMETERS :
!
!     IDAT    : DATE in the following form => YYYYMMDD
!     NSSSSS  : TIME of the day in seconds
!     PLAT    : LATITUDE  in Degrees
!     PLON    : LONGITUDE in Degrees
!
!     OUTPUT PARAMETERS :
!
!     PMU0    : Cosine of solar zenith angle
!     PMU0M   : Cosine of solar zenith angle (mean value)
!
!
!     J.F. Mahfouf (4/12/97) from IFS/ARPEGE routines
!
!-----------------------------------------------------------------------
!
! - Astronomical functions
!   you will find the description in the annex 1 of the documentation
!   RRS is the distance Sun-Earth
!   RDS is the declination of the Earth
!   RET is the equation of time
!
!   Orbit of the earth
!     
      RTETA(PTIME)=PTIME/(RDAY*365.25)
      REL(PTETA)=1.7535+6.283076*PTETA
      REM(PTETA)=6.240075+6.283020*PTETA
      RRS(PTETA)=REA*(1.0001-0.0163*SIN(REL(PTETA)) &
     &           +0.0037*COS(REL(PTETA)))
!   Relative movement Sun/Earth
      RLLS(PTETA)=4.8951+6.283076*PTETA
      RLLLS(PTETA)=4.8952+6.283320*PTETA-0.0075*SIN(REL(PTETA)) & 
     &          -0.0326*COS(REL(PTETA))-0.0003*SIN(2.*REL(PTETA)) &
     &          +0.0002*COS(2.*REL(PTETA))
      RDS(PTETA)=ASIN(SIN(REPSM)*SIN(RLLLS(PTETA)))
      RET(PTETA)=591.8*SIN(2.*RLLS(PTETA))-459.4*SIN(REM(PTETA)) &
     &   +39.5*SIN(REM(PTETA))*COS(2.*RLLS(PTETA)) &
     &   -12.7*SIN(4.*RLLS(PTETA))-4.8*SIN(2.*REM(PTETA))
!     
!    -------------------------------------------------------------
!
! - Time functions
!   the descriptions are in the annex 1 of the documentation
!
!   TIME
!
!   NDD   : extraxt dd from ccaammdd
!   NMM   : extract mm from ccaammdd
!   NAA   : extract aa from ccaammdd
!   NCCAA : extract ccaa from ccaammdd
!   NAMD  : extract aammdd from ccaammdd
!   NCTH  : turn seconds into hours
!   RTIME : returns the time of the model (in seconds of course!)
!
      NDD(KGRDAT)  =MOD(KGRDAT,100)
      NMM(KGRDAT)  =MOD((KGRDAT-NDD(KGRDAT))/100,100)
      NCCAA(KGRDAT)=KGRDAT/10000
      NAA(KGRDAT)=MOD(NCCAA(KGRDAT),100)
      NAMD(KGRDAT)=MOD(KGRDAT,1000000)
      NCTH(KSEC)=KSEC/3600
!
      NZZAA(KAAAA,KMM)=KAAAA-( (1-ISIGN(1,KMM-3))/2 )
      NZZMM(KMM)=KMM+6*(1-ISIGN(1,KMM-3))
      RJUDAT(KAAAA,KMM,KDD)=1720994.5 + FLOAT( &
     &   2-NZZAA(KAAAA,KMM)/100 + (NZZAA(KAAAA,KMM)/100)/4 &
     & + INT(365.25*FLOAT(NZZAA(KAAAA,KMM))) &
     & + INT(30.601*FLOAT(NZZMM(KMM)+1)) &
     & + KDD)
      RTIME(KAAAA,KMM,KDD,KSS)=(RJUDAT(KAAAA,KMM,KDD)-2451545.) &
     &     *RDAY+FLOAT(KSS)
!
!     Set-up constants
!     
      RPI=2.*ASIN(1.)
      RDAY=86400.
      REPSM=0.409093                ! obliquity
!
!     Angle conversions 
!      
      PGEMU=SIN(PLAT*RPI/180.)      ! sinus of latitude
      PGELAM=PLON*RPI/180.          ! longitude
!      
      ID=NDD(IDAT)
      IM=NMM(IDAT)
      IA=NCCAA(IDAT)
      RTIMTR=RTIME(IA,IM,ID,NSSSSS)
      ZTETA=RTETA(RTIMTR)
      RDECLI=RDS(ZTETA)             ! declinaison
      REQTIM=RET(ZTETA)
      RHGMT=REAL(MOD(NSSSSS,NINT(RDAY)))
      RSOVR =REQTIM+RHGMT    
      RWSOVR=RSOVR*2.*RPI/RDAY      ! hour angle
!      
      RCODEC=COS(RDECLI)
      RSIDEC=SIN(RDECLI)
!
      RCOVSR=COS(RWSOVR)
      RSIVSR=SIN(RWSOVR)
!      
      PMU0=MAX( RSIDEC*PGEMU  &
     &         -RCODEC*RCOVSR*SQRT(1.-PGEMU**2)*COS(PGELAM) &
     &         +RCODEC*RSIVSR*SQRT(1.-PGEMU**2)*SIN(PGELAM) , 0.)
      IF (PMU0.GT.0.) THEN
        PMU0=SQRT(1224.*PMU0*PMU0 +1.)/35. ! Magnification factor
      ENDIF  
!
!     Mean angle over the previous 6 hours
!     ------------------------------------
!      
      ITRAD=10800
      RTIMTRM=RTIME(IA,IM,ID,NSSSSS-ITRAD)
      ZTETAM=RTETA(RTIMTRM)
      RDECLIM=RDS(ZTETAM)             ! declinaison
      REQTIMM=RET(ZTETAM)
      IF ((NSSSSS-ITRAD).LT.0) THEN
        NSSSSS=NSSSSS+86400
      ENDIF  
      RHGMTM=REAL(MOD(NSSSSS-ITRAD,NINT(RDAY)))
      RSOVRM =REQTIMM+RHGMTM    
      RWSOVRM=RSOVRM*2.*RPI/RDAY      ! hour angle
!      
      RCODECM=COS(RDECLIM)
      RSIDECM=SIN(RDECLIM)
!
      RCOVSRM=COS(RWSOVRM)
      RSIVSRM=SIN(RWSOVRM)
!     
      PMU0M=MAX( RSIDECM*PGEMU &
     &         -RCODECM*RCOVSRM*SQRT(1.-PGEMU**2)*COS(PGELAM) &
     &         +RCODECM*RSIVSRM*SQRT(1.-PGEMU**2)*SIN(PGELAM) , 0.)
      IF (PMU0M.GT.0.) THEN
        PMU0M=SQRT(1224.*PMU0M*PMU0M +1.)/35. ! Magnification factor
      ENDIF  
!
END SUBROUTINE SOLAR_ANGLE
