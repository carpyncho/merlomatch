      parameter (ncmax=80,nrmax=2000000)
      real*4 rec(nrmax,ncmax)
      real*4 magzpt,magzrr,maglim,magobj,magerr
      character*80 infile,argv,errmsg,dateobs,utstart
      character*72 comment
      character*68 value
      character*16 ttype(ncmax),tunit(ncmax),tform(ncmax)
      character*11 filter
      character*1 csign,czero
      integer status,blocksize
      common /rd/ tpa,tpd,a,b,c,d,e,f,projp1,projp3,projp5,tand,secd,
     1            tpi,iwcstype
c *** FITSIO_CAT_LIST - simple program to read in cats fix fluxes etc....
c *** ncmax - max no. of columns in table
c *** nrmax - max no. of rows in table
      li=3
      lr=5
      lt=6
      pi=4.0*atan(1.0)
      tpi=8.0*atan(1.0)
      pio2=pi/2.0
      czero='0'

c *** try for presence of command line arguments
      m=iargc()
c *** if correct number of arguments try and interpret them
      if(m.eq.1)then
        call getarg(1,argv)
        infile=argv
      else
        print *,' '
        print *,'Arguments: filename'
        print *,' '
        stop
      endif
c *** open fits file and read in catalogue info
      iunit=1
      status=0
      call ftopen(iunit,infile,1,blocksize,status)
      if(status.ne.0)then
        call ftgerr(status,errmsg)
        print *,'FITSIO Error status =',status,': ',errmsg
        ier=1
        status=0
        call ftclos(iunit,status)
      endif
      call ftthdu(iunit,nhdus,status)
      nmefs=max(1,nhdus-1)

c *** clobber output file if exists
      call clobber(3,'catalogue.asc')    
      open(unit=3,file='catalogue.asc',status='new')

c *** see if VIRCAM and get some info from PHU
      ivircam=0
      call ftgkey(iunit,'INSTRUME',value,comment,status)
      if(status.eq.202)then
        status=0
      else
        if(index(value,'VIRCAM').gt.0)then
          ivircam=1
          call ftgkys(iunit,'ESO INS FILT1 NAME',filter,comment,
     1                status)
          call ftgkye(iunit,'ESO TEL AIRM END',airmass,comment,
     1                status)
          call ftgcrd(iunit,'DATE-OBS',dateobs,status)
          utstart='UTST    ='
        endif
      endif

c *** now go through each extension
      do mef=1,nmefs
      call fitsio_rtable(rec,iunit,mef,nrows,ncols,apcor,delapcor,
     1  flim,filter,exptime,utstart,dateobs,percorr,magzpt,magzrr,
     2  airmass,extinct,skyvar,gain,saturate,pltscl,theta,ivircam)

c *** some magnitude info
      magzpt=magzpt-(airmass-1.0)*extinct+2.5*log10(exptime)
      maglim=magzpt-2.5*log10(flim)                       ! 5-sigma mag limit

c *** print out some stuff as an example
      do k=1,nrows
      if(ncols.eq.32)then
        xcord=rec(k,5)
        ycord=rec(k,6)
        flux=max(0.1,rec(k,4))                            ! apcor flux
        fluxerr=sqrt(flux/gain)+skyvar                    ! approx error
        peak=rec(k,10)
        ellipt=rec(k,8)
        pa=rec(k,9)-theta                                 ! correct to sky PA
        pa=90.0-pa                                        ! don't ask
c        if(pa.gt.180.0)pa=pa-180.0
        icls=rec(k,25)
c *** flag possibly saturated objects and update core fluxes
        if(peak.gt.saturate)then
          icls=-9                                         ! saturated ?
          fixed=(rec(k,21)-rec(k,4))/(delapcor-1.0) 
          flux=max(fixed,flux)                            ! cunning sat corr
        endif
c *** flag objects containing bad pixels
        if(rec(k,30).gt.0.0)icls=-7
      else
        xcord=rec(k,3)
        ycord=rec(k,5)
        flux=max(0.1,rec(k,24))
        fluxerr=max(1.0,rec(k,25))
        peak=rec(k,18)
        ellipt=rec(k,8)
        pa=rec(k,9)-theta                                 ! correct to sky PA
        pa=90.0-pa                                        ! don't ask
c        if(pa.gt.180.0)pa=pa-180.0
        icls=rec(k,61)
c *** flag possibly saturated objects and update core fluxes
        if(peak.gt.saturate)then
          icls=-9                                         ! saturated ?
          fixed=(rec(k,28)-rec(k,24))/(delapcor-1.0) 
          flux=max(fixed,flux)                            ! cunning sat corr
        endif
c *** flag objects containing bad pixels
        if(rec(k,55).gt.0.0)icls=-7
      endif

c *** a bit of astrometry
      xx=xcord
      yy=ycord
      call radeczp(xx,yy)
      call rahour(xx,ih,im,ss)
      call radeg(yy,id,in,rr)
      if(yy.lt.0)then
        csign='-'
        id=abs(id)
        in=abs(in)
        rr=abs(rr)
      else
        csign='+'
      endif

c *** a bit of flux fixing
      call distort(xcord,ycord,distortcorr)            ! for field distortion
      flux=flux*apcor*percorr*distortcorr              ! + aperture correction
      fluxerr=fluxerr*apcor*percorr*distortcorr
      magobj=magzpt-2.5*log10(flux)                    ! object magnitude
      magerr=2.5*log10(1.0+fluxerr/flux)
      magerr=max(0.01,magerr)                          ! allow for systematics

c *** and some dumping
      if(id.lt.10)then
        if(flux.gt.0.5)then                           ! only detected objects
        write(3,'(2i3,f7.3,1x,2a,i1,i3,f6.2,2f9.2,2f8.3,2i3,2f7.2)') 
     1  ih,im,ss,csign,czero,id,in,rr,xcord,ycord,magobj,magerr,mef,
     2  icls,ellipt,pa
        endif
      else
        if(flux.ge.0.5)then
        write(3,'(2i3,f7.3,1x,a,i2,i3,f6.2,2f9.2,2f8.3,2i3,2f7.2)') 
     1  ih,im,ss,csign,id,in,rr,xcord,ycord,magobj,magerr,mef,icls,
     2  ellipt,pa
        endif
      endif

      enddo                                            ! nrows loop
      enddo                                            ! mef loop

c *** and finish off
      call ftclos(iunit,status)
      if(status.ne.0)then
        call ftgerr(status,errmsg)
        print *,'FITSIO Error status =',status,': ',errmsg
      endif
      end

c *** -----------------------------------------------------------------------

      subroutine fitsio_rtable(rec,iunit,mef,nrows,ncols,apcor,
     1 delapcor,flim,filter,exptime,utstart,dateobs,percorr,magzpt,
     2 magzrr,airmass,extinct,skyvar,gain,saturate,pltscl,theta,ivircam)
      parameter (ncmax=80,nrmax=2000000)
      real*4 rec(nrmax,ncmax)
      integer status,hdutype,naxes(2)
      real*4 enullval,magzpt,magzrr
      character*80 errmsg,utstart,dateobs
      character*72 comment
      character*68 value
      character*16 ttype(ncmax),tunit(ncmax),tform(ncmax)
      character*11 filter
      logical anynull
      common /rd/ tpa,tpd,a,b,c,d,e,f,projp1,projp3,projp5,tand,secd,
     1            tpi,iwcstype
c *** FITSIO_RTABLE - reads fits binary tables
      pi=4.0*atan(1.0)
      anynull=.false.
      enullval=0.0
c *** find out what sort of table extension it is
      status=0
      call ftmahd(iunit,mef+1,hdutype,status)
      if(hdutype.eq.1)then
        print *,'Extension is an ASCII table'
        stop ' - not implemented yet'
      else if(hdutype.eq.2)then
        print *,'Extension is a BINARY table'
      endif
c *** read assorted keywords to find table dimensions etc.
      call ftgknj(iunit,'NAXIS',1,2,naxes,nfound,status)
      if(nfound.ne.2)stop 'Not enough NAXIS keywords'
      ncols=naxes(1)/4                                    ! bytes to r*4
      nrows=naxes(2)
      print *,' '
      print *,'No. of columns =',ncols
      print *,'No. of rows    =',nrows
      if(nrows.gt.nrmax)then
        print *, ' '
        print *, 'nrmax limit exceeded'
        print *, ' '
        stop 
      endif
      if(ivircam.eq.0)then
        call ftgkys(iunit,'WFFBAND',filter,comment,status)
        if(status.eq.202)then
          status=0
          call ftgkys(iunit,'FILTER',filter,comment,status)
          if(status.eq.202)then
            status=0
            call ftgkys(iunit,'HIERARCH ESO INS FILT1 NAME',
     1                  filter,comment,status)
          endif
        endif
        if(status.eq.202)then
          status=0
          filter='         '
        endif
      endif

      call ftgkye(iunit,'EXPTIME',exptime,comment,status)
      if(status.eq.202)then
        status=0
        call ftgkye(iunit,'EXPOSED',exptime,comment,status)
        if(status.eq.202)then
          status=0
          call ftgkye(iunit,'EXP_TIME',exptime,comment,status)
        endif
      endif
      if(status.eq.202)then
        status=0
        exptime=0.0
      else
        exptime=abs(exptime)
        exptime=max(1.0,exptime)
      endif

      call ftgkye(iunit,'SKYLEVEL',skylevel,comment,status)
      if(status.eq.202)then
        status=0.0
        call ftgkye(iunit,'ESO DRS SKYLEVEL',skylevel,comment,status)
        if(status.eq.202)then
          status=0
          skylevel=0.0
        endif
      endif

      call ftgkye(iunit,'SKYNOISE',skynoise,comment,status)
      if(status.eq.202)then
        status=0
        call ftgkye(iunit,'ESO DRS SKYNOISE',skynoise,comment,status)
        if(status.eq.202)then
          status=0
          skynoise=1.0
        endif
      endif

      call ftgkye(iunit,'RCORE',rcore,comment,status)
      if(status.eq.202)then
        status=0
        call ftgkye(iunit,'ESO DRS RCORE',rcore,comment,status)
        if(status.eq.202)then
          status=0
          rcore=1.0
        endif
      endif

      flim=5.0*sqrt(pi*rcore**2)*skynoise           ! 5-sigma flux limit
      skyvar=(pi*rcore**2)*skynoise**2              ! for error analysis

      call ftgkye(iunit,'GAIN',gain,comment,status)
      if(status.eq.202)then
        status=0
        gain=2.0                                    ! default
      endif
      if(ivircam.eq.1)gain=4.0

      call ftgkye(iunit,'SATURATE',saturate,comment,status)
      if(status.eq.202)then
        status=0
        saturate=30000.0                            ! safe guess
      endif
      saturate=0.9*(saturate-skylevel)              ! relative to sky

      if(status.ne.0)then
        call ftgerr(status,errmsg)
        print *,'FITSIO Error status =',status,': ',errmsg
      endif

      if(ncols.eq.32)then
        call ftgkye(iunit,'APCOR',apcor,comment,status)
        if(status.eq.202)then
          status=0
          apcor=0.25
        endif
      else
        call ftgkye(iunit,'APCOR3',apcor,comment,status)
        if(status.eq.202)then
          status=0
          apcor=0.25                                ! sensible default
        endif
      endif
      if(ncols.eq.32)then
        call ftgkye(iunit,'APCOR3',apcor3,comment,status)
        if(status.eq.202)then
          status=0
          apcor3=0.0
        endif
      else
        call ftgkye(iunit,'APCOR5',apcor3,comment,status)
        if(status.eq.202)then
          status=0
          apcor3=0.0
        endif
      endif
      delapcor=10.0**(0.4*(apcor-apcor3))
      apcor=10.0**(0.4*apcor)
 
      call ftgkye(iunit,'PERCORR',percorr,comment,status)
      if(status.eq.202)then
        status=0
        percorr=0.0
      endif
      percorr=10.0**(0.4*percorr)

      flim=flim*apcor*percorr                    ! including ap and per corr

      if(ivircam.eq.0)then
        call ftgkye(iunit,'AIRMASS',airmass,comment,status)
        if(status.eq.202)then
          status=0
          call ftgkye(iunit,'AMSTART',airmass1,comment,status)
          call ftgkye(iunit,'AMEND',airmass2,comment,status)
          if(status.eq.202)then
            status=0
            airmass=1.0
          else
            airmass=0.5*(airmass1+airmass2)
          endif
        endif
      endif

      call ftgkye(iunit,'MAGZPT',magzpt,comment,status)
      if(status.eq.202)then
        status=0
        magzpt=-1.0
      endif

      call ftgkye(iunit,'MAGZRR',magzrr,comment,status)
      if(status.eq.202)then
        status=0
        magzrr=-1.0
      endif

c *** get extinction value - if not present use default
      call ftgkye(iunit,'EXTINCT',extinct,comment,status)
      if(status.eq.202)then
        status=0
        extinct=0.05
      endif

c *** check if WCS information present
      call ftgkey(iunit,'CTYPE1',value,comment,status)
      if(status.eq.202)then
        status=0
        call ftgkey(iunit,'TCTYP3',value,comment,status)
      endif
      if(status.ne.202)then
        if(index(value,'RA---').gt.0)then           ! need a better test
          iwcs=1                  
        else
          iwcs=0
        endif
      else
        iwcs=0
        status=0
      endif
      if(iwcs.eq.1)then
        print *,'Valid WCS present'
      else
        print *,'No WCS present'
      endif
c *** now get celestial coordinate information the hard way 
      if(iwcs.eq.1)then
        iwcstype=1                                             ! ZPN default
        if(ivircam.eq.0)then
          call ftgkys(iunit,'CTYPE1',value,comment,status)
          if(index(value,'RA---TAN').gt.0)iwcstype=0           ! TAN only other
          call ftgkye(iunit,'CRPIX1',c,comment,status)
          call ftgkye(iunit,'CRPIX2',f,comment,status)
          call ftgkye(iunit,'CRVAL1',tpa,comment,status)
          call ftgkye(iunit,'CRVAL2',tpd,comment,status)
        else
          call ftgkys(iunit,'TCTYP3',value,comment,status)
          if(index(value,'RA---TAN').gt.0)iwcstype=0           ! TAN only other
          call ftgkye(iunit,'TCRPX3',c,comment,status)
          call ftgkye(iunit,'TCRPX5',f,comment,status)
          call ftgkye(iunit,'TCRVL3',tpa,comment,status)
          call ftgkye(iunit,'TCRVL5',tpd,comment,status)
        endif
        tpa=tpa*pi/180.0
        tpd=tpd*pi/180.0
        tand=tan(tpd)
        secd=1.0/cos(tpd)
        call ftgkye(iunit,'CD1_1',a,comment,status)
        if(status.eq.202)then
          status=0
          call ftgkye(iunit,'TC3_3',a,comment,status)
        endif  
        a=a*pi/180.0
        call ftgkye(iunit,'CD2_2',e,comment,status)
        if(status.eq.202)then
          status=0
          call ftgkye(iunit,'TC5_5',e,comment,status)
        endif
        e=e*pi/180.0
        call ftgkye(iunit,'CD1_2',b,comment,status)
        if(status.eq.202)then
          status=0
          call ftgkye(iunit,'TC3_5',b,comment,status)
          if(status.eq.202)then
            status=0
            b=0.0
          endif
        endif
        b=b*pi/180.0
        call ftgkye(iunit,'CD2_1',d,comment,status)
        if(status.eq.202)then
          status=0
          call ftgkye(iunit,'TC5_3',d,comment,status)
          if(status.eq.202)then
            status=0
            d=0.0
          endif
        endif
        d=d*pi/180.0
        if(iverbose.eq.1)then
          print *,' '
          print *,'Frame transform constants:'
          print *,a,b,c
          print *,d,e,f
          print *,'Tangent points:',tpa,tpd
        endif
c *** pixel scale and orientation
        pltscl=sqrt(0.5*(a**2+b**2+d**2+e**2))*180.0/pi               ! deg/pix
        theta=0.5*(atan(abs(b)/abs(a))+atan(abs(d)/abs(e)))*180.0/pi  ! deg
c *** read zp constants from header
        projp1=1.0
        call ftgkye(iunit,'PV2_3',projp3,comment,status)
        if(status.eq.202)then
          status=0
          call ftgkye(iunit,'PROJP3',projp3,comment,status)
          if(status.eq.202)then
            status=0
            call ftgkye(iunit,'TV5_3',projp3,comment,status)
            if(status.eq.202)then
              status=0
              projp3=0.0
            endif
          endif
        endif
        call ftgkye(iunit,'PV2_5',projp5,comment,status)
        if(status.eq.202)then
          status=0
          call ftgkye(iunit,'TV5_5',projp5,comment,status)
          if(status.eq.202)then
            status=0
            projp5=0.0
          endif 
        endif
      endif
c *** observation time and date
      if(ivircam.eq.0)then
        call ftgcrd(iunit,'UTSTART',utstart,status)
        if(status.eq.202)then
          status=0
          call ftgcrd(iunit,'UTC-OBS',utstart,status)
        endif
        if(status.eq.202)then
          utstart='UTST    ='
          status=0
        endif
        call ftgcrd(iunit,'DATE-OBS',dateobs,status)
        if(status.eq.202)then
          dateobs='DATE    ='
          status=0
        endif
      endif
      if(status.ne.0)then
        call ftgerr(status,errmsg)
        print *,'FITSIO Error status =',status,': ',errmsg
      endif
c *** now get rest of table info
      call ftgkns(iunit,'TTYPE',1,ncmax,ttype,nfound,status)
      call ftgkns(iunit,'TUNIT',1,ncmax,tunit,nfound,status)
      call ftgkns(iunit,'TFORM',1,ncmax,tform,nfound,status)
      if(ncols.ne.nfound)stop 'Something screwy with header'
      if(iverbose.eq.1)then
        print *,' '
        print *,'Parameters in table:'
        print *,(ttype(i),i=1,ncols)
        print *,' '
      endif
c *** now read data column by column (yuk!)
      do j=1,ncols
      call ftgcve(iunit,j,1,1,nrows,enullval,rec(1,j),anynull,status)
      enddo
c *** check for errors
      if(status.ne.0)then
        call ftgerr(status,errmsg)
        print *,'FITSIO Error status =',status,': ',errmsg
      endif
      return
      end

c *** -----------------------------------------------------------------------

      subroutine radeczp(x,y)
      common /rd/ tpa,tpd,a,b,c,d,e,f,projp1,projp3,projp5,tand,secd,
     1            tpi,iwcstype
c *** RADECZP converts x,y coordinates to ra and dec using ZP wrt optical axis
c *** NB. This is an extension of the ARC projection using Zenithal polynomials
      x=x-c
      y=y-f
      xi=a*x+b*y
      xn=d*x+e*y
      r=sqrt(xi**2+xn**2)
      if(iwcstype.eq.1)then
        rfac=projp1+projp3*r**2+projp5*r**4    ! NB this is a 1st order approx
        r=r/rfac
        rfac=projp1+projp3*r**2+projp5*r**4    ! now 2nd order correction
      else
        if(r.eq.0.0)then
          rfac=1.0
        else
          rfac=r/tan(r)
        endif
      endif
      xi=xi/rfac 
      xn=xn/rfac
      aa=atan(xi*secd/(1.0-xn*tand))
      alpha=aa+tpa
      delta=atan((xn+tand)*sin(aa)/(xi*secd))
      x=alpha
      y=delta
      if(x.gt.tpi)x=x-tpi
      if(x.lt.0.0)x=x+tpi
      return
      end

c *** -----------------------------------------------------------------------

      subroutine distort(x,y,distortcorr)
      common /rd/ tpa,tpd,a,b,c,d,e,f,projp1,projp3,projp5,tand,secd,
     1            tpi,iwcstype
c *** DISTORT converts x,y coordinates to standard coordinates 
c ***         and works out flux distortion factor
c *** NB. This is an extension of the ARC projection using Zenithal polynomials
      xi=a*(x-c)+b*(y-f)
      xn=d*(x-c)+e*(y-f)
      r=sqrt(xi**2+xn**2)
      if(iwcstype.eq.1)then
        distortcorr=1.0+3.0*projp3*r**2/projp1+
     1                  5.0*projp5*r**4/projp1 ! this is only a 1st order app
        distortcorr=distortcorr*(1.0+(projp3*r**2+projp5*r**4)/projp1)
      else
        if(r.eq.0.0)then
          distortcorr=1.0
        else
          distortcorr=1.0/(1.0+tan(r)**2)
          distortcorr=distortcorr*atan(r)/r
        endif
      endif
      distortcorr=1.0/distortcorr              ! is this the right sign
      return                                   ! r.dr.dtheta --> r'dr'dtheta
      end

c *** -----------------------------------------------------------------------

      subroutine clobber(iunit,name)
c *** CLOBBER checks to see if file name exists and then deletes it
      integer ier,stat,statb(13),status,system
      character*(*) name
      character*80 filename
      character*3 del /'rm '/
      character*1 space /' '/
      filename=name
      i=1
      do while (filename(i:i).ne.' '.and.i.lt.80)
      i=i+1
      enddo
      open(unit=iunit,file=filename)
      ier=stat(filename,statb)
      if(statb(8).eq.0.or.ier.ne.0)then
        status=system(del//filename)
      else
        print *,'File: '//filename(1:i)//' already exists - clobbering'
        status=system(del//filename)
      endif
      close(unit=iunit)
      return
      end

c *** ------------------------------------------------------------------------

      subroutine rahour(radian,ihour,min,sec)
      character*6 dum
c *** RAHOUR  converts radians into hours,mins,secs
      tpi=8.0*atan(1.0)
      if(radian.lt.0.0)radian=radian+tpi
c *** convert to secs
      temp=10800.0*abs(radian)/atan(1.0)
      ir=int(temp)
      ihour=ir/3600
      min=(ir-3600*ihour)/60
      sec=(temp-ihour*3600-min*60)

c *** check for illegal sexagesimals
      write(dum,'(f6.2)') sec
      if(index(dum,'60.').ne.0)then
        sec=0.0
        min=min+1
      endif
      if(min.eq.60)then
        min=0
        ihour=ihour+1
        if(ihour.eq.24)ihour=0
      endif

      return
      end

c *** ------------------------------------------------------------------------

      subroutine radeg(radian,ideg,min,sec)
      character*5 dum
c *** RADEG  converts radians into degs,mins,secs
      iflag=1
      if(radian.lt.0.0)iflag=-1
c *** convert to secs
      temp=162000.0*abs(radian)/atan(1.0)
      ir=int(temp)
      ideg=ir/3600
      min=(ir-3600*ideg)/60
      sec=(temp-ideg*3600-min*60)*iflag

c *** check for illegal sexagesimals
      write(dum,'(f5.1)') sec
      if(index(dum,'60.').ne.0)then
        sec=0.0
        min=min+1
      endif
      if(min.eq.60)then
        min=0
        ideg=ideg+1
      endif

      min=min*iflag
      ideg=ideg*iflag
      return
      end
