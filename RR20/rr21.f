c rr21.f combine  mm and time loop
c negative ions removed
c rr20a.f has abs electric field printed in record
c rr20a.f is for ele calculations for eelectron density after current zero
c rr20.f attempt to include space charge at base in pois in initial mm iterationss
c r19.f has modificaions in z velocities to allow for reversed velocities after current zero, see sub. elec
c r19.f omits electron diffusion to avoid complexity.
c r18.f hs large component of jetair8d.f for solution of electron density
c rr17 is a working copy of rr16.f
c rr16.f has positive ion and negative ion densities added
c rr15.f has attachment just for edge of arc from diffusion of electrons
c rr14.f has electrin convection added to rr11.f
c rr11.f power supply impedance line 493
c rr11.f adds time calculation  line 612 nrrvtest 0
c rr10.f has account of power supply impedance included; line 488 - not executed
c rr10.f has thermal ionization term sc made zero for volts<10 line 3995
c rr10.f is rr9ff.f with correct p(z) scaling factors
c rr9f.f has p(z)*p(z) in attach1
c rr9e.f has p(z) increase of attachment
c rr9d.f has ambipolar diffusion coefficient of de voto replaced by D = (1/3)lambda*velocity
c rr9c has i less pont for zz for matlab
c rr9.f has axial vwlocities at outer radius calculated from (1/2) rho v**2 = pressure drop
c rr9.f has 10cm long arc and pressures at outer radius from 2000A BBC arc
c rrrv8 has solution for 10 cm long arc
c rrrv8 is initial try for sonic flow; from rrrv7.f
c rrrv7 has lte option; if lte = 1, lte results with lte ne.
c rrrv6 has fluc functions removed
c rrrv5 has first read of matf as attachment rate coefficient and temp above which attachment is zero
c rrrv5 has attachment; argon recombination; attach coef const below given temperature 
c Tungsten cathode also after current zero
c rrrv4 has time dependent electron continuity equation, didt before current zero and rrrv after current zero 
c rrrv3 has attachment
c rrrv2 has code for didt of circuit breakeri; input controlled by actest
c recombination coefficients are from Hinnov and herschberg, for argon
c input is kappa gas rther than kappa/cp
c rrrv1.f starts with flux5.f
c flux5 changes: (1) common cath statement put in readmatf to transfer dtdh to other subroutines
c (2) in readmatf   dtdh statement changed to  dtdh = 1700.d0/(ah(3) - ah(1)) 
c (3) in reaadmatf changes to cazet,anzet : cazet(j) = cazet(j)/acp(j);
c     anzet(j) = anzet(j)/acp(j)
c (4) introduce kappa arrays in read statement to preserve zet for kappa/cp
c (5) Ref pressure in liquid aref =p(1,nza)-p(1,nza-1) in mom so p arc = p liquid at r=0 and nza
c (6) Energy balance printed
c Term dsdt(j,i)*dtdh*h(j,i) is included on both sides of the energy balance equation by Haidar
c    to assist numerical stability for points on the surface of the electrodes. dsdt = ds/dT
c    dtdh = dT/dh. S includes electron heating or cooling + ion heating + Black Body radiation cooling.
c    dsdt(j,i) is evaluated in subroutine cash.
c insulating layer for j.gt.npr
c arcfn shortened
c fux material functions in anode spot;  jjl
c weld pool included;  MATLAB plot files included; jjl
c initial velocities and pressure in the anode can be set to zero line 298
c radial velocity at top of weldpool zero, line 2879
c jXB can be removed, line 3285, mom
c includes Marangoni surface tension, line 3287, mom
c gravity  can be removed, line 3341, fillv
c arcfn does not have temperatures, only enthalpy
c-----------------------------------------------------------------------
c	30 / 11 / 1999
c      version of laurent sansonnens modified by jjl for Macintosh FORTRAN 77
c      this version of the code solve the electron continuity equation
c      for the whole plasma, contains the electron diffusion current in
c      the current continuity equation and uses a non-lte treatment for the
c      electrical conductivity.
c
c      this version of the code can handle:
c      1- straight polarity 
c      3- various electrode configuration (conical and rectangular)
c      4- time dependent calculation (fixed electrode configuration)
c         varying arc current. not tested
c      5- rho's cp's solid are used for the electrodes.
c     	the upper electrode is designated by array ica
c       ica=1 in the upper electrode and ica=0 elsewhere. 
c     	the lower electrode is designated by array ian
c       ian=1 in the electrode and ian=0 elsewhere.

c jjl trial printout
c  568  format(' Trial Function')
c  569  format(' i =',i4, ' z=',f8.5, ' nr =',i4, ' nz =',i4)
c 572  format(1p8e10.2)
c      do 780 i=79,83
c         write(9,569) i,z(i),nr,nz
c         write(9,572) (bjz(j,i),j=1,20)
c 780  continue
c      do 781 i=79,83
c         write(9,569) i,z(i),nr,nz
c         write(9,572) (bjr(j,i),j=1,20)
c 781  continue
c       pause
c
c jjl trial printout
c 568  format(' Trial Function')
c 569  format(' i =',i4, ' z=',f8.5, ' nr =',i4, ' nz =',i4)
c 572  format(1p8e10.2)
c          write(9,568)
c 780  continue
c      do 781 i=78,83
c         write(9,569) i,z(i),nr,nz
c         write(9,572) (dc(j,i),j=1,20)
c 781  continue
c      do 782 i=78,83
c         write(9,569) i,z(i),nr,nz
c         write(9,572) (ann(j,i),j=1,20)
c 782  continue
c      do 783 i=78,83
c         write(9,569) i,z(i),nr,nz
c         write(9,572) (ans(j,i),j=1,20)
c 783  continue
c      do 784 i=78,83
c         write(9,569) i,z(i),nr,nz
c         write(9,572) (anw(j,i),j=1,20)
c 784  continue
c      do 785 i=78,83
c         write(9,569) i,z(i),nr,nz
c         write(9,572) (ane(j,i),j=1,20)
c 785  continue
c       pause
c subroutine cont produces plot output for pv-wave.
c  this default setting is done in line after data instruction
c energy eq.  solves for enthalpy
c      0 brouwer   main
c      1 initialise initialise the common variables to zero
c      2 initmp    generates initial t when stated from scratch
c      3 bound     tranlates  temp. at boundaries into h
c      4 intpol    converts t to h
c      5 amf       material functions lte sigma, ne and mobility
c      6 ener      call fille, solve for enthalpie and tranlates to t
c      7 fille     fill coefficient for energy eq.
c      8 htemp     translates enthalphy into temperature
c      9 cpsolid
c     10 rvpol     reverse polarity
c     11 stpol     straight polarity
c     12 ansh      energy at anode surface 
c     13 cash      energy at cathode surface
c     15 mom       call fillv, solve for vz, vr and p
c     17 fillv     fill coefficient for vz, vr
c     18 mbal
c     19 zerocoef  
c     20 poiss     current continuity eq. with diffusion current
c     21 fill      fill coefficient for poiss
c     22 elec      electron cont. equation
c     23 fillde    fill coefficient for elec if using Patankar method.         
c     25 solvne    solver
c     26 solve     solver
c     27 solvpp    solver
c     28 solvie    solver
c     29 residu    
c     30 tdma      
c     31 meshls    generates mesh from data in gridls
c     32 multcoef  used in meshls
c     33 volume    used in meshls
c     34 cellpos   used in meshls 
c     35 readmatf  read input mat. funct. from matf
c     36 readarcin reads input parameters from arcin
c     37 readfn reads intial set and writes final solution
c     38 cont      writes to plot files if plot is set true
c     39 ppp       
c     40 meshrec   
c
c-----------------------------------------------------------------------
      program  brouwer
      implicit double precision (a-h,o-z)
      logical rvpol,actest,plot
c-----------------------------------------------------------------------
      character*40 stfile,fname
      parameter(mr=200, mz=100, nsh=21)
      dimension cd0(0:mr),hetr(0:mr)
      common/itoor/ i itoor
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/visc/
     r  eta(0:mr,0:mz),etm(0:mr,0:mz),etr(0:mr,0:mz)  
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
      common/liqu/  asigliqu(41),tmla,rholiqu,
     1  rhocoeff,cpliqu,cpcoeff,npr 
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr  
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/pol/
     r  rvpol
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz)     
      common/cath/
     r  dtdh, dsdt(0:mr, 0:mz),rc
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi
      common/vcoef/
     r   aen(mr,mz),aes(mr,mz),aee(mr,mz),aew(mr,mz),
     1  aec(mr,mz),bsr(mr,mz),ann(mr,mz),ans(mr,mz),ane(mr,mz),
     2  anw(mr,mz),anc(mr,mz),bsz(mr,mz)
      common/mag/
     r   bth(0:mr,0:mz),bjz(0:mr,0:mz),bjr(0:mr,0:mz)
      common/tcatht/
     r  setjz(0:mz),richjz(0:mz),setjr(0:mr),richjr(0:mr),caions
      common/currpoi/
     r  zi(0:mz),cathi
      common/hbound/
     r  hamb,hc,he,hanod(0:mr),hwall
      common/emiss/
     r  emissc,emissa
      common/surface/
     r  rsf(500),zsf(500),rsfo(0:mr,0:mz),
     1  zsfo(0:mr,0:mz),pstr(0:mr,0:mz),pstz(0:mr,0:mz)
      common/plot/ iplot(50),jplot(50),iplotn,jplotn,ib,jend
      dimension tt(0:mr,0:mz),vvz(0:mr,0:mz),phirr(0:mr,0:mz),
     r    vvr(0:mr,0:mz),err(0:mr,0:mz),ezz(0:mr,0:mz),rrr(0:mr)
c---------------------------------------------------------------------
      data pid,pi/6.283185307d+00,3.141592654d+00/
c---------------------------------------------------------------------
        time=0.0
        e = 1.60206d-19
        epsilon = 8.85d-14
        icurr=0
        irs=0
        irun=0
        iarcin=0
        stfile='arcst'
        mctoor=0
        itoor = 0
c---------------------------------------------------------------------
        resmax = 1.0
        resphi = 1.0
        resvr = 1.0
        resvz = 1.0
        resp = 1.0
        rese = 1.0
        resne = 1.0
c-----------------------------------------------------------------
c initialise the variables to zero
        call initialise(cd0,hetr)
c------------------------------------------------------------------
c reads the material function in files 'matf' and 'matliquid'.
      call readmatf
c reads the input parameters in file 'arcin'
        call readarcin(currm,flowin,rc,twall,changf,ifread,
     1  te,tc,rlxphi,rlxp,rlxvr,rlxvz,rlxt,rlxne,deltel,didt,delt,
     2  p0,jscreen,rlamt,tamb,rvpol,rrrv,rsdl,rldl,ntimesteps,mm,
     4  ielec,nra,nzci,istore,lte,nnn,npr,nnel,
     5  actest)
c        call readarcin(currm,flowin,rc,twall,changf,actest,delt,
c     1  te,tc,rlxphi,rlxp,rlxvr,rlxvz,rlxt,rlxne,deltel,didt,
c     2  p0,jscreen,rlamt,tamb,rvpol,rrrv,rsdl,rldl,ntimesteps,mm,
c     4  ielec,nra,nzci,istore,lte,nnn,npr,
c     5  ifread)
	
c------------------------------------------------------------------
c generates the mesh from the data read in file 'gridls'
        call meshls
        write(*,*)'mesh generation completed'
c-----------------------------------------------------------------------
      curr = currm
      delt0=delt
      rvpol=.false.
c      if(curr.gt.0)rvpol=.true.
      do 5 i=1,41
        au(i) = 2.*pid*au(i)
   5  continue  
c--------------------
      if(jscreen.eq.1)then
       write(*,10) curr,nz,nr,nzc,nza,nra
       write(*,20) twall,te,tc,p0,nrspc,flowin
       write(*,25) rlxphi,rlxp,rlxvr,rlxvz,rlxt,rlxne,mm
       if (rvpol) write(*,255)
       if (actest) write(*,256) currm,delt,didt,rrrv 
       if (ifread.eq.0) write(*,26)
       if (ifread.eq.1) write(*,27)
      endif
c-----------------------------------------------------------------
c     initial conditions:
c if ifread = 1 strat from file arcfn
c if ifread = 0 strat from scratch.
c----------------------------------------------------------------
      if (ifread .eq. 1) then
        open(2,file='arcfn')
        call readfn(1,2,time,curr,0,nrp,0,nzp)
        close(2)
        call htemp(itoor,itor,itoz,htoor)
        curr=currm
c make phi values positive, not negative, as applies after current zero
         do 480 i=0,nzp
         do 481 j=0,nrp
            phi(j,i) = abs(phi(j,i))
  481 continue
  480 continue
      else
        call initmp(t,twall)
         do 100 i=0,nzp
         do 98 j=0,nrp
             s = t(j,i)/1000.
c---linear interpolaton t->h
	    call intpol(s,val)
	    h(j,i)=val
            rho(j,i) = arho(19)
            vz(j,i) = 0.d0
            vr(j,i) = 0.d0
            p(j,i) = 0.d0
  98     continue
 100     continue

c---linear interoplaton t->h
            s=te/1000.
	    call intpol(s,val)
            he=val 

         do 101 i=0,nzc
         do  97 j=0,nrca(i)
            if(ica(j,i).eq.1)then
              s=t(j,i)/1000
	      call intpol(s,val)
	      h(j,i) =val
	    endif
  97     continue
 101     continue
        endif
        open(305,file='pzfile')
        rewind 305
c        do 776 i=1,3
	    read(305,*) (ppp(i),i=1,nz)
c  776   continue
        close(305)
c------- end initial condition ------------------------------------

c------------------------------------------
c     set up backplane boundary
c
        do 152 j=0,nr
          vz(j,0)=vz(j,1)
          vr(j,0)=vr(j,1)
  152   continue
c      if (flowin .ne. 0.) then
c        router = r(nrspc) + 0.5*dr(nrspc)
c        rinner = r(nrca(1))+ 0.5*dr(nrca(1))
c        fac = rinner/router
c        facsq = fac*fac
c        alofac = dlog(router/rinner)
c        term2 = (1.-facsq)/alofac
c        annul = pi*(router*router-rinner*rinner)
c        flrm = flowin/annul
c        denom = (1+facsq-term2)
c        sum = 0.d0
c        do 151 j=nrca(1)+1,nrspc
c           rr = r(j)/router
c           vz(j,0) = 2.*flrm*(1.-rr*rr+term2*dlog(rr))/denom
c           sum = sum + vz(j,0)*an(j)
c 151    continue
c      else
c        do 152 j=0,nrp
c          vz(j,0)=0.d0
c  152   continue
c      end if
c        flowout = pid*sum
c        if (jscreen.eq.1) write(*,1050) flowin, flowout
c---------------------------------------------------------------
c     setting up boundary conditions.
c
c--- get values for boundaries in enthalpy 
         call bound(twall,tamb,te,tc)

      do 110 i=0,nzp
         t(nrp,i) = twall
	 h(nrp,i)= hwall
         vz(nrp,i) = 0.d0
         t(0,i) = t(1,i)
	 h(0,i)=h(1,i)
c         vz(0,i) = vz(1,i)
c         vr(0,i) = -vr(1,i)
 110  continue
      do 130 j=0,nrca(1)
        if(irca(1).gt.1.and.j.lt.irca(1))then
         t(j,0) = twall
         h(j,0) = hwall
        else
         t(j,0) = te
	 h(j,0) =he
         vz(j,0) = 0.
        endif
 130  continue
      do 140 j=nrca(1)+1,nrp
         t(j,0) = twall
         h(j,0) = hwall
 140  continue
      do 150 j=0,nrp
         t(j,nzp) = twall
         h(j,nzp) = hwall
 150  continue

c----------------------------------------------------------------------
      time = 0.d0
      nrrvtest = 0
      curr = currm
      open(305,file='timevalues')
        rewind 305
c  initial densities of positive and negative ions
      if (time.eq.0.0) then
      do 151 j=0,nrp
      do 153 i=0,nzp
           rnp(j,i) = rne(j,i)
           rni(j,i) = 0.0
  153 continue
  151 continue
      endif
c----------------------------------------------------------------
c     time iterations start here. 
c         Over all calculation "do 6000" with n total time steps
c           within this do loop are three separate sub do loops
c           (1) initial calculation of dc current; "do 200", m to mm
c           (2) time dependent caculation n to 6000; ntimesteps
c           (3) calculation of electron densiy at end of each step of (2)
c               by iterations within subroutine "elec".
c----------------------------------------------------------------------
c      do 6000 n = 0,ntimesteps
c  nnn is number of time steps to current zero
          nnn = abs(currm/(didt*delt))
          timezero = nnn*delt
c Boundary condition for phi for z = 0 cm
c    For toimes after current zero voltages positive; set voltcz = 900.0
         voltcz = -900.0
            nrrvtest = 1
c          if (time.gt.timezero) then 
c            nrrvtest = 1
c	    voltcz = (time - timezero)*rrrv
c          else
c             nrrvtest = 0
c             do 301 j=0,nrp
c                 phi(j,1) = voltcz*(time - timezero)/timezero
c 301         continue
c          endif
c if nrrvtest = 1 applies rrrv after current zero
c     remember the old values
        do 745 i=0,nzp
        do 746 j=0,nrp
            s = t(j,i)/1000.
            l = int(s)
            l1 = l+1
            l2 = l1+1
            frac = s-l
            rfr = 1.-frac
           rho1(j,i) = arho(l2)*frac + arho(l1)*rfr
 	 if(ica(j,i).eq.1)rho1(j,i)=carho
 	 if(ian(j,i).eq.1)rho1(j,i)=anrho
 746     continue
 745     continue
c
         do 765 i=0,nz
         do 766 j=0,nr
            rhmr1(j,i) = f(j)*rho1(j+1,i) + (1.-f(j))*rho1(j,i)
            rhmz1(j,i) = g(i)*rho1(j,i+1) + (1.-g(i))*rho1(j,i)
  766    continue
  765    continue
         do 770 i=0,nz
            rhmr1(nrp,i) = rho1(nrp,i)
            rhmz1(nrp,i) = rho1(nrp,i)
 770     continue
         do 721 i=0,nz
         do 722 j=0,nr
            t1(j,i) = t(j,i)
            h1(j,i) = h(j,i)
            vz1(j,i) = vz(j,i)
            vr1(j,i) = vr(j,i)
 722  continue
 721  continue
c
c     setting up cathode current density.
c
        rcath=r(nrca(1))+0.5*dr(nrca(1))
        rcathi=0.
        if(irca(1).gt.1) rcathi=r(irca(1))-0.5*dr(irca(1))
          cd=1.*curr/(pi*(rcath**2-rcathi**2))
          itirc=0
          if(irca(1).gt.1)itirc=irca(1)
          do 600 j=itirc,nrp
              cd0(j) = cd
  600     continue
      cop = 0.d0
      jmax=0
      do 46 j=1,nrp
         cop = cop + pid*an(j)*cd0(j)
	if (dabs(cop) .gt. dabs(curr)) go to 47
  46  continue
      go to 48
  47  jmax = j
      cd0(jmax) = cd0(jmax) - (cop-curr)/(pid*an(jmax))
      if (dabs(cd0(jmax)) .lt. 1.0d-6) cd0(jmax) = 0.d0
      do 49 j=jmax+1,nrp
         cd0(j) = 0.d0
  49  continue
  48  if (jscreen.eq.1) write(*,801) cop
      if (jscreen.eq.1) write(*,802) curr
      if (jscreen.eq.1) write(*,1060) jmax
      if (jscreen.eq.1) write(*,1070)
      do 120 j=0,nrp
            cdz(j,0) = cd0(j)
 120  continue
      timeelec = 0.0
      nrrvtest =  1
c Start initial calculaions with poiss
      volt = 0.0
      do 200 m=1,mm
         timeelec = timeelec + deltel
c        call amf(lte)
c        if(itoor.eq.1)then
c             write(*,*)' 1             temp>30,000, stop'
c             go to 461
c        endif
          if (jscreen.eq.1) write(*,990)
 990    format(1x,' <solving for phi>')
         volt = timeelec*rrrv
         call poiss(curr,volt,nrrvtest)
         call mom(flowin,irst)
	 call elec

      do 812 i=1,nz
      do 813 j=1,nr+1
        ezzz = (ez(j,i-1) + ez(j,i))/2
	errr = (er(j-1,i) + er(j,i))/2
        eabs(j,i) = sqrt(ezzz**2 + errr**2)
  813  continue
  812  continue

c          if (jscreen.eq.1) write(*,991)
 991    format(1x,' <solving for t>')
c         call ener(itoor,itor,itoz,htoor,irst)
c        if(itoor.eq.1) then
c        go to 461
c        endif
c         volt = phi(1,1)
c--- get current  difference
        diffcm=0.0
        idiffcm=0
      do 1100 i=1,nz
        diffc=dabs(zi(i))-dabs(curr)
        if (dabs(diffc) .gt. dabs(diffcm) ) then
                diffcm=diffc
                idiffcm=i
        endif
 1100 continue
      resmax = dmax1(resp, resvr, resvz, rese,resphi,resne)
      relmax = dmax1(relp, relvr, relvz, rele,relphi,relne)
      if (m.eq.mm) then
        write(*,850) m, resmax,relmax
        write(*,859)
        write(*,861) resphi,jph,iph,phi(jph,iph),relphi,jjph,
     1        iiph,phi(jjph,iiph)
        write(*,862) rese,je,ie,t(je,ie),rele,jje,iie,t(jje,iie)
c        write(*,863) resvz,jvz,ivz,vz(jvz,ivz),relvz,jjvz,iivz,
c     1        vz(jjvz,iivz)
c        write(*,864) resvr,jvr,ivr,vr(jvr,ivr),relvr,jjvr,iivr,
c     1        vr(jjvr,iivr)
        write(*,865) resp,jp,ip,p(jp,ip),relp,jjp,iip,p(jjp,iip)
	write(*,866) resne,jne,ine,rne(jne,ine),relne,jjne,iine,
     1        rne(jjne,iine)
        endif

        nzcp=nzc+nrca(nzc)-1
        nzc1=nzcp-nzc
        nt1=nzc+nrca(nzc)-irca(nzc)
        nt2=nrca(nzc)-irca(nzc)
c      if (nrca(nzc).gt.1) write(*,502)(richjz(nt1-i),i=0,nt2)
c      if (jscreen.eq.1) write(*, 502) (richjz( nzc + 1 - i), i = 1, nzc)
c      if (nrca(nzc).gt.1) write(*,502) (t(j,nzc),j=irca(nzc),nrca(nzc))
	write(463,*)'not conv. at time  ',time
c     converged?
	if (resmax.le.rsdl)  then
           if (m.gt.1) go to 210
        end if
 200  continue
c
c     end internal iteration with m; (1).
c

 210	continue
 
       if(actest) write(*,*)irs,ncrs
      open(304,file='arcfnout')
        rewind 304
        call readfn(2,304,time,curr,0,nrp,0,nzp)
      close(304)
c       if (.not. actest) go to 6001
         call amf(lte)
        if(itoor.eq.1)then
        write(*,*)' 2             temp>30,000, stop'
        go to 461
        endif
      do 220 i=0,nzp
      do 221 j=0,nrp
         t1(j,i) = t(j,i)
         h1(j,i) = h(j,i)
         vz1(j,i) = vz(j,i)
         vr1(j,i) = vr(j,i)
         rho1(j,i) = rho(j,i)
         rhmr1(j,i) = rhmr(j,i)
         rhmz1(j,i) = rhmz(j,i)
	 rne1(j,i) = rne(j,i)
 221  continue
 220  continue
c values of current and voltage with time
      if (time.le.timezero) then
          voltp = phi(1,1)
      else
          voltp = voltcz
      endif	  
      write(*,226) n, time, curr, volt
      write(*,181) curr,time,n,volt
      write(305,181) curr,time,n,volt
       if(time.le.timezero)   curr = curr + didt*delt
          time = time + delt
c 6000 continue
      curr = curr + didt*delt
      time = time + delt
c End time iterations with n; (3)
         close(305)

c Electron density calculation (3)
        timeel = 0.0
        do 222 iii = 1, ntimesteps
             call elec
             curr = currm
             nrrvtest = 1
            call poiss(curr,volt,nrrvtest)
             timeel = timeel + deltel
 222    continue
        irs=irs+1
	if(incrs.gt.0)then
	ncrs=ncrs+incrs
	incrs=0
	endif
      write(*,226) mm, timeel, curr, volt
c  end (3)

c	call cellpos
c
         call amf(lte)
        if(itoor.eq.1)then
        write(*,*)' 3             temp>30,000, stop'
        go to 461
	endif
c---calculating total radiated power.
          radpow=0.d0
        do 719 i=1,nz
        do 718 j=1,nr
          radpow= radpow+rad(j,i)*vol(j,i)
  718   continue
  719   continue
          radpow=radpow*pid
c--- fraction of radiation and total power  in percent
c    radiatedt power fraction
         radpof= radpow/curr/phi(1,1)*100.
         if (jscreen.eq.1) write(*,517)radpow,radpof
         if (jscreen.eq.1) write(*,514)
c      if (nrca(nzc).gt.1) write(*,500) (t(j,nzc),j=irca(nzc),nrca(nzc))
c-----------
c
c     calculating heat input into anode
c
      do 720 j=1,nr
         hetr(j) = zet(j,nz)*(t(j,nz)-t(j,nzp))/(0.5*dz(nz))
 720  continue
c
c     calculating heat input into anode and cathode
c
      nzcp =  nzc + 1
      ahanode = 0.d0
      bhanode = 0.d0
      zhanode = 0.d0
      hcathode = 0.d0
      nza1 = nza+1
      nza0 = nza-1
       do 730 j=1,nr
          akapp = cpz(j,nza0)*zemz(j,nza0)
          dhanode = 0.d0-akapp*(t(j,nza)-t(j,nza0))*an(j)/
     1         (z(nza)-z(nza0))
          ahanode = ahanode + dhanode
          bkapp = cpz(j,nza)*zemz(j,nza)
          dhanode = 0.d0-bkapp*(t(j,nza1)-t(j,nza))*an(j)/
     1         (z(nza1)-z(nza))
          bhanode = bhanode + dhanode
          zkapp = cpz(j,nz)*zemz(j,nz)
          dhanode = 0.d0-zkapp*(t(j,nz)-t(j,nz-1))*an(j)/
     1         (z(nz)-z(nz-1))
          zhanode = zhanode + dhanode
          zkapp = cpz(j,1)*zemz(j,1)
          dcathode = zkapp*(t(j,2)-t(j,1))*an(j)/
     1         (z(2)-z(1))
          hcathode = hcathode + dcathode
 730  continue
      ahanode = 2.*pi*ahanode
      bhanode = 2.*pi*bhanode
      zhanode = 2.*pi*zhanode
      hcathode = 2.*pi*hcathode
c
c   calculating heat input into electrode sides
c
      hside = 0.d0
      do 740 i=1,nzc
         nrcn = nrca(i)
         nrcp = nrcn + 1
         hside = hside + 2.*zemr(nrcn,i)*(t(nrcp,i)-te)*ae(nrcn,i)
     1          /(dr(nrcn)+dr(nrcp))
         if (nrcn .lt. nrc) then
            hside = hside + 2.*zemz(nrcp,i-1)*(t(nrcp,i)-te)*an(nrcp)
     1              /(dz(i)+dz(i-1))
         end if

        if(irca(i).gt.1)then
         ircn = irca(i)
         ircp = ircn - 1
         hside = hside + 2.*zemr(ircp,i)*(t(ircp,i)-te)*ae(ircp,i)
     1          /(dr(ircn)+dr(ircp))
         if (irca(i) .gt. irca(i-1)) then
            hside = hside + 2.*zemz(ircp,i-1)*(t(ircp,i)-te)*an(ircp)
     1              /(dz(i)+dz(i-1))
         end if
        endif
 740  continue
      hside = 2.*pi*hside
      htot =  hanode +  hside
         volt = phi(irca(1),1)
      pow = volt*curr
 978    format(i3.3)
 909    format(1x,i3,6(2pe12.4))
 910    format(7x,'xne',10x,'sne',9x,'pos.',8x,'temp.',8x,
     &  'diff',8x,'s-funct.')
 463    format(2x,5(3x,f6.3),3x,i5)
 496     format(2e15.5)
 995     format(1p6e13.5)
 519     format (3e24.8)
  1       format(1x,' enthalpie out of range  node j= '
     &    ,i2,' i= ',i2,' h= ',e12.4 )
  10  format(1x,'curr =',f7.1,'a, nz=',i3,', nr=',i3,', nzc=',i3,
     1    ', nza=',i3,', nra=',i3)
   8  format(' Uses lte sigma not from non-lte ne/')
   9  format(' Laurent code with ne solution instead of sigma for arc'/)
  20  format(1x,'temps: wall, electrde, cathde, (k) =',f10.2,
     1     f10.2,f10.2,/1x,'ambient pres. =',f5.1,
     2     ', nrspc =',i3,', flowin =',f10.1,' c.c/s')
  25   format(1x,'rlxn. factors =',6(f6.2,','),' no. iterations =',i5)
  26  format(1x,'initial t, vz, vr, phi fields generated by program')
  27  format(1x,'initial t, vz, vr and phi fields from <arcst>')
 918    format(2(1p,i4),2(2p,f12.4))
 920    format((1p,i4),2(2p,f12.4))
  30  format(1x,'nodal positions - axial direction'/(9f8.5))
  34  format(2(e11.4))
  37  format(/,' cathode suface temperatures')
  38  format(/,' cathode z current densities')
  39  format(/,' cathode rich. current densities')
  40  format(1x,'nodal positions - radial direction'/(9f8.5))
  41  format(/,' cathode r current densities')
  42  format(/,' ion heating of cathode',f10.2,' W')
  43   format(/,' Energy Balance',f10.2,' %')
 495     format(2e15.5)
  50     format(1p6e13.5)
  53    format(5(e11.4))
 181    format('current = ',f10.4,5x,' time= ',e10.4, ' Time No =',i4,
     &      ' voltcz =',f10.4)
 182    format('current = ',f10.4,5x,' time= ',e10.4,2x,'nb'2x,i4)
 226  format('Time Itn =',i4,1x,'time =',0p,e10.4' curr =',0p,f10.2,
     & ' voltage =',f10.3)
 227  format(1x,'curr =',0p,f10.2,' voltage   =',f10.3)
 255  format(1x,'reverse polarity arc')
 256  format(1x,'ac current =',1x,e12.4,' delt =',1x,e12.4,' didt =',
     & 1x,e12.4,' rrrv =',1x,e12.4)
 500     format(8f9.0)
 501     format(/)
 502     format(8f9.0)
 503     format('  z  velocity')
 504     format('  r  velocity')
 505     format('  temperature')
 506  format('  pressure')
 507  format('  potential')
 508  format('  z  current  density')
 509  format('  r  current  density')
 510  format('  i =',i4,'  z(i)  =',e14.6)
 511  format(14f7.2)
 512  format(14f6.0)
 513  format(11f8.0)
 514  format(' cathode surface temperatures ')
 516  format(9f7.0)
 517  format ( ' total radiated power  ',f6.0,' [w]',
     &             ' percentage of total power',f6.3,' % ')
 518  format (3e24.8)
 520  format (f9.3,2x, f10.5)
 521  format (1pe10.2,1pe10.2,1pe10.2,1pe10.2,1pe10.2,1pe10.2,
     1       1pe10.2,1pe10.2,1pe10.2,1pe10.2)
 556  format(1x,'heat transfer to anode below surface (w)',1p7e11.3)
 559  format(1x,'heat transfer to anode above surface (w)',1p7e11.3)
 560  format(1x,'heat transfer to anode bottom surface (w)',1p7e11.3)
 561  format(1x,'heat transfer from cathode top surface (w)',1p7e11.3)
 562  format(1x,'heat from anode-cathode work functions (w)',1p7e11.3)
 572  format(1p6e12.4)
 801    format(1x,'calculated current = ', f10.4, ' amps')
 802    format(1x,'current = ', f10.4, ' amps')
 850  format(1x/,'iteration no.',i4,
     * '   maximum residual =', 1pe9.2,'   max relerror =',e9.2/)
 851  format(' max. res. radial momentum  =',e12.3,' at node',2i4,e12.3)
 852  format(' max. res. axial momentum   =',e12.3,' at node',2i4,e12.3)
 853  format(' max. res. enthalpie        =',e12.3,' at node',2i4,f8.0)
 854  format(' max. res. phi   (0.1 %)    =',e12.3,' at node',2i4,e12.3)
 855  format (' max. curr diff. ',e12.3,' [a] ','at row i= ',i3)
 856  format (' current emerging from cathode ',e12.3,'  [a] ')
 859  format('          res    je   ie      phi        rel      jje  iie
     1   phi')
 860  format(' max. residual mass balance =',e12.3,' at node',2i4)
 861  format(' phi',1pe12.2,2i4,1pe12.2,1pe12.2,2i4,e12.2)
 862  format(' h  ',1pe12.2,2i4,1pe12.2,1pe12.2,2i4,e12.2)
 863  format(' vz ',1pe12.2,2i4,1pe12.2,1pe12.2,2i4,e12.2)
 864  format(' vr ',1pe12.2,2i4,1pe12.2,1pe12.2,2i4,e12.2)
 865  format(' p  ',1pe12.2,2i4,1pe12.2,1pe12.2,2i4,e12.2)
 866  format(' ne ',1pe12.2,2i4,1pe12.2,1pe12.2,2i4,e12.2)
 899     format(1x,'iteration no. ',i4)
 1050   format(1x,'flowin =',f10.2,', flowout =',f10.2,' c.c/sec')
 1060  format( ' jmax =',i4)
 1070  format(/'  insufficient j(r) for input current'/)
 1072  format('  deltel =',2e12.3)
c
      write(*,1072) deltel
      write(*,226) n, time, curr, volt
      if (resmax.le.rsdl)  write(*,226)n, time, curr, volt
        fname='record'
      open(5,file=fname)
      rewind 5
      if (lte.eq.1) then
         write(5,8)
      else
         write(5,9)
      endif
      if(rvpol)then
	write(5,*)' ***********reverse polarity***********'
      else
	write(5,*)' ***********straight polarity***********'
      endif
      write(5,10) curr,nz,nr,nzc,nza,nra
      write(5,20) twall,te,tc,p0,nrspc,flowin
      write(5,25) rlxphi,rlxp,rlxvr,rlxvz,rlxt,rlxne,mm
      if (rvpol)  write(5,255)
      if (actest) write(5,256) currm,delt,didt,rrrv
      if (ifread .eq. 0) write(5,26)
      if (ifread .eq. 1) write(5,27)
        write(5,30) (z(i),i=0,nzp)
        write(5,40) (r(j),j=0,nrp)
        write(5,1050) flowin, flowout
c        write(5,801) cop
        write(5,1060) jmax
        write(5,850) mc, resmax,relmax
        write(5,859)
        write(5,861) resphi,jph,iph,phi(jph,iph),relphi,jjph,
     1        iiph,phi(jjph,iiph)
        write(5,862) rese,je,ie,t(je,ie),rele,jje,iie,t(jje,iie)
c        write(5,863) resvz,jvz,ivz,vz(jvz,ivz),relvz,jjvz,iivz,
c     1        vz(jjvz,iivz)
c        write(5,864) resvr,jvr,ivr,vr(jvr,ivr),relvr,jjvr,iivr,
c     1        vr(jjvr,iivr)
        write(5,865) resp,jp,ip,p(jp,ip),relp,jjp,iip,p(jjp,iip)
	write(5,866) resne,jne,ine,rne(jne,ine),relne,jjne,iine,
     1        rne(jjne,iine)
        write(5,855) diffcm,idiffcm
        write(5,856) cathi
        write(5,227) curr,phi(1,1)
      write(5,555) volt,pow
      write(5,517) radpow,radpof
      write(5,556) bhanode
      write(5,559) ahanode
      write(5,560) zhanode
      write(5,561) hcathode
      workf =0d0- curr*(aworkf-cworke)
      write(5,562) workf
      write(5,42) caions
      enbalance = 100*(bhanode+radpow+hcathode-pow-workf-caions)/pow
      write(5,43) enbalance
         write(5,501)
         write(5,514)
      nzcp=nzc+nrca(nzc)-1
      nzc1=nzcp-nzc
      if (nrca(nzc).gt.1) write(5,500) (t(j,nzc),j=irca(nzc),nrca(nzc))
      write(5, 500) (t(nrca(nzc + 1 -i), nzc + 1 - i), i = 1, nzc)
      write(5,38)
      write(5,502) (setjz(nntip+i),i=1,nuni)
      write(5, 502) (setjz( nntip+ 1 - i), i = 1, nzc)
      write(5,39)
      write(5,502) (richjz(nntip+i),i=1,nuni)
      write(5, 502) (richjz( nntip+ 1 - i), i = 1, nzc)
      write(5,41)
      write(5, 502) (setjr( nzc + 1 - i), i = 1, nzc)
      write(5,39)
      write(5, 502) (richjr( nzc + 1 - i), i = 1, nzc)
         write(5,501)
         write(5,503)
c      do 775 i=0,nza-1
c         write(5,511) (vz(j,i)*0.01,j=0,nrp)
c 775  continue
c      do 789 i=nza,nzp
c      do 789 i=nza,nzp
      do 789 i=0,nzp
         write(5,510) i,z(i)
         write(5,521) (vz(j,i),j=0,nrp)
 789  continue
         write(5,501)
         write(5,504)
c      do 776 i=1,nza-1
c         write(5,510) i,z(i)
c         write(5,511) (vr(j,i)*0.01,j=0,nrp)
c 776  continue
c      do 788 i=nza,nzp
      do 788 i=0,nzp
         write(5,510) i,z(i)
         write(5,521) (vr(j,i),j=0,nrp)
 788  continue
         write(5,501)
         write(5,505)
      do 777 i=0,nzp
         write(5,510) i,z(i)
         write(5,512) (t(j,i),j=0,nrp)
 777  continue
         write(5,501)
         write(5,506)
      do 782 i=0,nzp
         write(5,510) i,z(i)
         write(5,521) (p(j,i),j=0,nrp)
 782  continue
         write(5,501)
         write(5,507)
      do 778 i=0,nzp
         write(5,510) i,z(i)
c         write(5,511) (phi(j,i),j=0,nrp)
         write(5,521) (phi(j,i),j=0,nrp)
 778  continue
         write(5,501)
         write(5,508)
      do 779 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (cdz(j,i),j=0,nrp)
 779  continue
         write(5,501)
         write(5,509)
      do 781 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (cdr(j,i),j=0,nrp)
 781  continue
c         write(5,501)
c         write(5,*) 'z diffusion current'
c      do 784 i=0,nzp
c         write(5,510) i,z(i)
c         write(5,572) (cdez(j,i),j=0,nrp)
c  784  continue
c         write(5,501)
c         write(5,*) 'r diffusion current'
c      do 783 i=0,nzp
c         write(5,510) i,z(i)
c         write(5,572) (cder(j,i),j=0,nrp)
c 783  continue 
         write(5,501)
         write(5,*) 'Electron density'
      do 785 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (rne(j,i),j=0,nrp)
 785  continue
         write(5,501)
         write(5,*) 'LTE electron density'
      do 797 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (rneth(j,i),j=0,nrp)
 797  continue  
         write(5,501)
         write(5,*) 'Positive Ion Density'
      do 811 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (rnp(j,i),j=0,nrp)
 811  continue  
         write(5,501)
         write(5,*) 'Negative Ion  Density'
      do 790 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (rni(j,i),j=0,nrp)
 790  continue  
         write(5,501)
         write(5,*) 'Net Charge Density'
      do 814 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (rnet(j,i),j=0,nrp)
 814  continue  
         write(5,501)
         write(5,*) 'Electron increments'
      do 791 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (dne(j,i),j=0,nrp)
 791  continue  
         write(5,501)
         write(5,*) 'Negative Ion increments'
      do 795 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (dni(j,i),j=0,nrp)
 795  continue  
         write(5,501)
         write(5,*) 'Positive Ion Increments'
      do 794 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (dni(j,i),j=0,nrp)
 794  continue  
         write(5,501)
         write(5,*) 'Absolute Electric Field'
      do 796 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (eabs(j,i),j=0,nrp)
 796  continue  
         write(5,501)
         write(5,*) 'current profile'
      do 786 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (rtest(j,i),j=0,nrp)
 786  continue
       if (jscreen.eq.1) write(*,555) volt,pow
      write(5,555) volt,pow
      close(5)
 555  format(1x,'arc voltage =',f10.2,'v arc power =',f10.1,'w')
      if (jscreen.eq.1) write(*,556) (hetr(j),j=1,nr)
      if (jscreen.eq.1) write(*,557) hanode, hside,  htot
 557  format(1x,'heat entering anode tip, side and total anode heat'/
     1        1p3e12.3)
     
c plot files for MATLAB or IDL
      open(20,file='r.dat')
      write(20,*) (r(j),j=1,jend)
      close(20)
       open(21,file='z.dat')
      write(21,*) (z(i),i=ib,nzm)
c      write(21,*) (z(i),i=ib,nz)
      close(21)
       open(22,file='t.dat')
       do 793 i=ib,nz
          write(22,*) (t(j,i),j=1,jend)
 793  continue
       close(22)
       open(25,file='vz.dat')
c      do 794 i=ib,nz
c          write(25,*) (vz(j,i),j=1,jend)
c 794  continue
      close(25)
       open(26,file='vr.dat')
c      do 795 i=ib,nz
c          write(26,*) (vr(j,i),j=1,jend)
c 795  continue
      close(26)
      open(27,file='rv.dat')
          write(27,*) (r(jplot(j)),j=1,jplotn)
      close(27)
c tt and rr files for full plot from r = -R to r = +R; jjl
      jjend = 2*jend
      do 142 j=1,jjend
          jjj = jend - j + 1
          if (jjj.le.0)  then 
              jjj = j - jend
              rrr(j) = r(jjj)
          else
              rrr(j) = -r(jjj)
          endif
 142   continue
      do 798 i = ib,nz
      do 799 j = 1,jjend
          jjj = jend - j + 1
          if (jjj.le.0) jjj = j - jend
          tt(j,i) = t(jjj,i)
          phirr(j,i) = phi(jjj,i)
          ezz(j,i) = ez(jjj,i)
          err(j,i) = er(jjj,i)
 799  continue
 798  continue      
       open(23,file='rr.dat')
      write(23,*) (rrr(j),j=1,jjend)
      close(23)
       open(24,file='tt.dat')
       do 792 i=ib,nz
          write(24,*) (tt(j,i),j=1,jjend)
 792  continue
      close(24)
       open(25,file='phirr.dat')
       do 815 i=ib,nz
          write(25,*) (phirr(j,i),j=1,jjend)
 815  continue
      close(25)
       open(26,file='ezz.dat')
       do 816 i=ib,nz
          write(26,*) (ezz(j,i),j=1,jjend)
 816  continue
      close(26)
       open(27,file='err.dat')
       do 817 i=ib,nz
          write(27,*) (err(j,i),j=1,jjend)
 817  continue
      close(27)
      do 803 i=ib,nz
      do 804 j=1,jjend
          jjj = jend-j+1
          if (jjj.le.0) jjj = j - jend
          rrne(j,i) = rne(jjj,i)
          if (rrne(j,i).le.1.0D-7) rrne(j,i)=1.0D-7
          rlogne(j,i) = log10(rrne(j,i))
 804  continue
 803  continue 
       open(28,file='logne.dat')
       do 800 i=ib,nz
          write(28,*) (rlogne(j,i),j=1,jjend)
 800  continue
      close(28)
      do 805 i=ib,nz
      do 806 j=1,jjend
          jjj = jend-j+1
          if (jjj.le.0) jjj = j - jend
          rrni(j,i) = rni(jjj,i)
          if (rrni(j,i).le.1.0D-7) rrni(j,i)=1.0D-7
          rlogni(j,i) = log10(rrni(j,i))
 806  continue
 805  continue 
       open(29,file='logni.dat')
       do 807 i=ib,nz
          write(29,*) (rlogni(j,i),j=1,jjend)
 807  continue
      close(29)
      do 808 i=ib,nz
      do 809 j=1,jjend
          jjj = jend-j+1
          if (jjj.le.0) jjj = j - jend
          rrnp(j,i) = rnp(jjj,i)
          if (rrnp(j,i).le.1.0D-7) rrnp(j,i)=1.0D-7
          rlognp(j,i) = log10(rrnp(j,i))
 809  continue
 808  continue 
       open(30,file='lognp.dat')
       do 810 i=ib,nz
          write(30,*) (rlognp(j,i),j=1,jjend)
 810  continue
      close(30)


  461   continue
      end 
c-------------------------------------------------
c  end main
c--------------------------------------------------------------------------
c----------------------------------------------------------------
c     subroutine initialise
c     v.1 ls 16.9.98
c----------------------------------------------------------------

        subroutine initialise(cd0,hetr)
      implicit double precision (a-h,o-z)
      logical rvpol,actest
      parameter(mr=200, mz=100, nsh=21)
      dimension cd0(0:mr),hetr(0:mr)
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/visc/
     r  eta(0:mr,0:mz),etm(0:mr,0:mz),etr(0:mr,0:mz)  
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr    
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)   
      common/pol/
     r  rvpol
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz) 
      common/cath/
     r  dtdh, dsdt(0:mr, 0:mz),rc
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi
      common/vcoef/
     r   aen(mr,mz),aes(mr,mz),aee(mr,mz),aew(mr,mz),
     1  aec(mr,mz),bsr(mr,mz),ann(mr,mz),ans(mr,mz),ane(mr,mz),
     2  anw(mr,mz),anc(mr,mz),bsz(mr,mz)
      common/mag/
     r   bth(0:mr,0:mz),bjz(0:mr,0:mz),bjr(0:mr,0:mz)
      common/tcatht/
     r  setjz(0:mz),richjz(0:mz),setjr(0:mr),richjr(0:mr),caions
      common/currpoi/
     r  zi(0:mz),cathi
      common/hbound/
     r  hamb,hc,he,hanod(0:mr),hwall 
      common/emiss/
     r  emissc,emissa
      common/surface/
     r  rsf(500),zsf(500),rsfo(0:mr,0:mz),
     1  zsfo(0:mr,0:mz),pstr(0:mr,0:mz),pstz(0:mr,0:mz)

        do 605 i=0,mz
        do 604 j=0,mr
        phi(j,i) = 1.d-06
        eohm(j,i)=0.0
 604    continue
 605    continue
        do 606 j=0,mr
        cd0(j)=0.d0
        hetr(j)=0.d0
 606    continue
        do 200 i=1,500
           rsf(i)=0.
           zsf(i)=0.
 200    continue
      do 490 i=0,mz
      do 491 j=0,mr
          vol(j,i)=0.
          pstr(j,i)=0.
          pstz(j,i)=0.
          rsfo(j,i)=0.
          zsfo(j,i)=0.
          ica(j,i)=0.0
          ian(j,i)=0.0
          ae(j,i)=0.0
          bth(j,i) = 0.
          bjz(j,i) = 0.
          bjr(j,i)= 0.
	  cdz(j,i)=0.0
          cdr(j,i)=0.0
          cp(j,i)=0.
          cdez(j,i)=0.0
          cder(j,i)=0.0
          dni(j,i)=0.0
          dnp(j,i)=0.0
          dne(j,i)=0.0
          dsdt(j,i) = 0.0
          diffa(j,i)=0.0
          diffar(j,i)=0.0
	  diffaz(j,i)=0.0
          eta(j,i)=0.
          etm(j,i)=0.
          etr(j,i)=0.
          eohm(j,i)=0.
          h(j,i)=0.0
	  h1(j,i)=0.0
          p(j,i)=0.0
          rho(j,i)=0.
          rhmr(j,i)=0.
          rhmz(j,i)=0.
          rjnpr(j,i) = 0.0
          rjnpz(j,i) = 0.0
          rjniz(j,i) = 0.0
          rjnir(j,i) = 0.0
          rjnez(j,i) = 0.0
          rjner(j,i) = 0.0
          rlogne(j,i) = 0.0
          rlogni(j,i) = 0.0
          rlognp(j,i) = 0.0
          rrne(j,i) = 0.0
          rrni(j,i) = 0.0
          rrnp(j,i) = 0.0
          rne(j,i)=0.0
          rnp(j,i)=0.0
          rni(j,i)=0.0
          rneth(j,i)=0.0
          rnet(j,i)=0.0
          rtest(j,i)=0.0
          ran(j,i)=0.0
          rad(j,i)=0.
          sig(j,i)=0.
          SS(j,i) = 0.0
          t(j,i)=0.0
	  t1(j,i)=0.0
          vr1(j,i)=0.0
          vz1(j,i)=0.0
          vr1(j,i)=0.0
          vz1(j,i)=0.0
          zet(j,i)=0.
          zemr(j,i)=0.
          zemz(j,i)=0.
  491 continue
  490 continue
      do 492 i=1,mz
      do 493 j=1,mr
          aen(j,i)=0.0
          aes(j,i)=0.0
          aee(j,i)=0.0
          aew(j,i)=0.0
          aec(j,i)=0.0
          bsr(j,i)=0.0
          ann(j,i)=0.0
          ans(j,i)=0.0
          ane(j,i)=0.0
          anw(j,i)=0.0
          anc(j,i)=0.0
          bsz(j,i)=0.0
  493 continue
  492 continue
      do 494 i=1,mz
          a(i)=0.
          b(i)=0.
          c(i)=0.
          d(i)=0.
  494 continue
      do 2 j=0,mr
          dr(j)=0.
          f(j)=0.
          hanod(j)=0.
          r(j)=0.
  2    continue
      do 3 i=0,mz
          dz(i)=0.
          g(i)=0.
          nrca(i)=0.
          iran(i)=0
          nran(i)=0
          irca(i)=0
          richjz(i)=0.
          setjz(i)=0.
          z(i)=0.
          zi(i)=0.
  3   continue
        end

c-------------------------------------------------
c end of subroutine initialise
c-------------------------------------------------

c-------------------------------------------------
c     subroutine initmp
c     v.1 ls 21.9.98
c
c     generates an initial temperature distribution
c---------------------------------------------------

        subroutine initmp(t,twall)
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      dimension t(0:mr,0:mz)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)

        tcath1=twall
        tcath2=10000.
        tcstep=(tcath2-tcath1)/nzc
cls        tarc1=20000.
cls        tarc2=20000
        tarc1=10000.
        tarc2=10000.
        tarstep=(tarc1-tarc2)/(nza-1-(nzc+1))
        tano1=10000.
        tano2=twall
        tanstep=(tano1-tano2)/(nzp-nza)

        do 700 j=0,nrp
        t(j,0)=twall
        t(j,nzp)=twall
  700   continue

        do 701 i=0,nzp
        t(nrp,i)=twall
  701   continue

        t(0,0)=tcath1
        do 699 i=1,nzc
           t(0,i)=tcath1+i*tcstep
           tcast1=(t(0,i)-t(nrp,i))/nrp
        do 702 j=1,nrp
           t(j,i)=t(0,i)-j*tcast1
  702   continue
  699   continue


        t(0,nzc+1)=tarc1
        do 703 i=nzc+2,nza-1
           t(0,i)=t(0,i-1)-tarstep
  703   continue

        t(0,nza)=tano1
        do 704 i=nza+1,nzp
        t(0,i)=t(0,nza)-(i-nza)*tanstep
  704   continue


        do 709 i=nza,nzp
        tarst1=(t(0,i)-t(nrp,i))/nrp
        do 705 j=1,nrp
        t(j,i)=t(0,i)-j*tarst1
  705   continue
  709   continue
        do 808 i=nza,nz
           tstp=(t(nran(nz)+1,i)-t(0,i))/nran(nz)
        do 809 j=1,nran(nz)
           t(j,i)=t(0,i)+tstp
 809    continue
 808    continue


        do 708 i=nzc+1,nza-1
        tarst1=(t(0,i)-t(nrp,i))/nrp
        do 710 j=1,nrp
c            t(j,i) = 1.0e+04*(1.-(r(j)/r(nrp))**0.7) + twall
        t(j,i)=t(0,i)-j*tarst1
  710   continue
  708   continue


        do 707 i=1,nzc
        tcast1=(t(nrca(i)+1,i)-t(nrp,i))/(nrp-nrca(i)-1)
        do 711 j=nrca(i)+2,nrp
        t(j,i)=t(j-1,i)-tcast1
  711   continue
  707   continue

        do 806 j=nran(nz)+1,nrp
           tstp=(t(j,nza-1)-t(j,nz))/(nz-nza+1)
        do 807 i=nza,nzp
           t(j,i)=t(j,i-1)-tstp
 807    continue
 806    continue

	do 878 i=0,nz
	do 879 j=0,nr
	if(ica(j,i).gt.ica(j+1,i).and.i.gt.nzc-5)then
	t(j+1,i)=tarc1
	endif
	if(ian(j,i).lt.ian(j,i+1).and.i.le.nza)then
	t(j,i)=tarc2
	endif
	if(ica(j,i).eq.0.and.ian(j,i).eq.0.and.
     &	j.lt.nrca(1).and.i.le.nza)then
	t(j,i)=tarc1
	endif
 879	continue
 878	continue
        end

c----------------------------------------------------
c end of subroutine initmp
c----------------------------------------------------
c-------------------------------------------------
c subroutine bound
c     v.1 ls 21.9.98
c    
c     translates t into h at boundaries
c-------------------------------------------------

       subroutine bound(twall,tamb,te,tc)

       implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      common/hbound/
     r  hamb,hc,he,hanod(0:mr),hwall
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/tcatht/
     r  setjz(0:mz),richjz(0:mz),setjr(0:mr),richjr(0:mr),caions
	
	  s=twall/1000.
	  call intpol(s,hwall)
	  s=tamb/1000.
	  call intpol(s,hamb)

	  
	  s=tc/1000.
	  call intpol(s,hc)
	  s=te/1000.
	  call intpol(s,he)
          end

c-------------------------------------------------------------
c end of subroutine bound
c------------------------------------------------------------

c----------------------------------------------------------
c     subroutine intpol
c     v.1 ls 21.9.98
c---------------------------------------------------------

	subroutine intpol(s,val)
        implicit double precision (a-h,o-z)
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
	  l=int(s)
	  l1 =l+1
	  l2 =l+2
	  frac = s-l
	  rfr =1.-frac 
	  val= ah(l2)*frac +ah(l1)*rfr
        end

c-----------------------------------------------------------
c end of subroutine intpol
c----------------------------------------------------------
c---------------------------------------------------------------------
c     subroutine amf(lte)
c     v.2 ls 29.9.98
c
c     this module is used to calculate the material properties at
c     each node using the 4-points lagrangian interpolation technique
c     for non-uniform distribution of the x-variable.
c     modified in order to calculate the diffusion current in the plasma
c---------------------------------------------------------------------


      subroutine amf(lte)
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
c      dimension zmob(0:mr,0:mz)
      common/itoor/ i itoor
c      common/ac/
c     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
c     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),
c     2   rhmz1(0:mr,0:mz),actest
      common/liqu/  asigliqu(41),tmla,rholiqu,
     1  rhocoeff,cpliqu,cpcoeff,npr  
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr           
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/visc/
     r  eta(0:mr,0:mz),etm(0:mr,0:mz),etr(0:mr,0:mz)    
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi
      common/cath/
     r  dtdh, dsdt(0:mr, 0:mz),rc
      common/tcatht/
     r  setjz(0:mz),richjz(0:mz),setjr(0:mr),richjr(0:mr),caions
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/emiss/
     r  emissc,emissa
      data pi4/1.256637061e+01/
      electron  = 1.60206d-19
      boltz     = 1.38044d-23
      elecmas   = 9.1094d-31
      pi        =  3.141592654d+00
c--- linar interpolation of   coefficients vs. temperature
      do 45 i=0,nzp
      do 46 j=0,nrp
         s = t(j,i)/1000.
         l = int(s)
         l1 = l+1
         l2 = l1+1
         frac = s-l
         if (s.lt.1.d0) frac = (s-0.3)/0.7
         rfr = 1.-frac
         if (s.lt.0.3) then
c jjl modification from
c            frac = 1.0
c            rfr = 0.0
             frac = 0.0
             rfr = 1.0
         end if
         rho(j,i) = arho(l2)*frac + arho(l1)*rfr
c--zet thermal conductivity / specif heat 
         zet(j,i) = azet(l2)*frac + azet(l1)*rfr
         cp(j,i) = acp(l2)*frac + acp(l1)*rfr
         rad(j,i) = au(l2)*frac + au(l1)*rfr
         eta(j,i) = aet(l2)*frac + aet(l1)*rfr
cls v.5  schmit function at every mesh point
         schm(j,i) =aschm(l2)*frac + aschm(l1)*rfr
cls v.2  lte electron density from nasaplus data
         rneth(j,i) = aneth(l2)*frac + aneth(l1)*rfr
c SS(j,i) is gamma*rneth**2;thermal ionization in electron density equation
         SS(j,i) = gamma*rneth(j,i)**2
         ran(j,i) = aan(l2)*frac + aan(l1)*rfr
         zran = ran(j,i)+2*(rneth(j,i)-rne(j,i))
         if(zran.lt.1.d-10) zran=0.d0
cls v.3  ambipolar diffusion coefficient
c        formula from devoto
c         diffa(j,i) =8.4d-6*t(j,i)**1.62355
c         formula from (1/3)lambda*velocity); Q=1.0E-16
         diffa(j,i) =0.22*t(j,i)**1.5
cls v.4  electron mobility
c        me = e/(me*veth*(ne*qei+na*qea))
         if (lte.ne.1) then 
c           if(t(j,i).lt.2.0d4) then
c           zqea = (0.39-0.551d-4*t(j,i)+0.595d-8*t(j,i)**2)*1.d-15*1.2
c jjl modification
           if (rne(j,i).lt.1.d0) rne(j,i) = 1.d0
c           zqei=(1.95d-6/t(j,i)**2)*dlog(1.53d8*t(j,i)**3/rne(j,i))*0.71
c           veth = dsqrt(8*boltz*t(j,i)/(pi*elecmas))*1.d2
c           znue = veth*(rne(j,i)*zqei+zran*zqea)
c jjl modifi jjl modification
c           if (rne(j,i).lt.10.d0) rne(j,i) = 1.d0
c           zqei=(1.95d-6/t(j,i)**2)*dlog(1.53d8*t(j,i)**3/rne(j,i))*0.71
c           veth = dsqrt(8*boltz*t(j,i)/(pi*elecmas))*1.d2
c           znue = veth*(rne(j,i)*zqei+zran*zqea)
c jjl modification
c           if (znue.le.0.001) then
c           zmob(j,i) = (amue(l2)*frac + amue(l1)*rfr)
c           else
c           zmob(j,i)=electron/(elecmas*znue)*1.d4
c           endif
c         else   
c           zmob(j,i) = (amue(l2)*frac + amue(l1)*rfr)
c         endif
cation
c           if (znue.le.0.001) then
c           zmob(j,i) = (amue(l2)*frac + amue(l1)*rfr)
c           else
c           zmob(j,i)=electron/(elecmas*znue)*1.d4
c           endif
c         else   
           zmob(j,i) = (amue(l2)*frac + amue(l1)*rfr)
c         endif
c non-lte electrical conductivity sigma = ne*e*mu
             sig(j,i) = rne(j,i)*electron*zmob(j,i)
	 else
	     sig(j,i) = asig(l2)*frac + asig(l1)*rfr
         endif
  46  continue
  45  continue
c---handle cathode coefficients
c      nzam = nza - 1
      nza1 = nza + 1
      do 303 i=0,nzp
      do 304 j=0,nrp
          if(ica(j,i).eq.1)then	
         s = t(j,i)/1000.
         l = int(s)
         l1 = l+1
         l2 = l1+1
         frac = s-l
         rfr = 1.-frac
         sig(j,i) = casig(l2)*frac + casig(l1)*rfr
         sig(j,i) = sig(j,i) 
         zet(j,i) = cazet(l2)*frac + cazet(l1)*rfr
         eta(j,i) = 1.0e+10
         rad(j,i) = 0.0
         rho(j,i)=carho
          endif
c --- handle anode coefficients
      if(ian(j,i).eq.1)then	
         s = t(j,i)/1000.
         l = int(s)
         l1 = l+1
         l2 = l1+1
         frac = s-l
         rfr = 1.-frac
         sig(j,i) = ansig(l2)*frac + ansig(l1)*rfr
c  insert low resistance for insulating layer of flux;   jjl
c      if (i.eq.nza.and.j.ge.npr)  sig(j,i) = 0.0001
c      if (i.eq.nza1.and.j.ge.npr)  sig(j,i) = 0.0001  
         zet(j,i) = anzet(l2)*frac + anzet(l1)*rfr
         eta(j,i) = 1.0e+10
         rad(j,i) = 0.0
         rho(j,i)=anrho
c  coefficients for liquid metal
      if(t(j,i).gt.tmla.and.i.ge.nza) then
          sig(j,i) = asigliqu(l2)*frac + asigliqu(l1)*rfr
c  insert low resistance for insulating layer of flux;   jjl
      if (i.eq.nza.and.j.ge.npr)  sig(j,i) = 0.0001
      if (i.eq.nza1.and.j.ge.npr)  sig(j,i) = 0.0001  
          rho(j,i) = rholiqu + rhocoeff*(t(j,i) - tmla)
      endif
      endif
  304  continue
  303  continue
c
c---calculate reciprocal coefficents patankar
c-  num.heat transfer and fluid flow  page 45.

      do 65 i=0,nz
      do 66 j=0,nr

c---cpr and cpz only used in conjuntion with power law in fille
         cpr(j,i) = f(j)/cp(j,i) + (1.-f(j))/cp(j+1,i)
        if (cpr(j,i) .lt. 1.e-09) then
            cpr(j,i) = 1.0e+10
	else 
	   cpr(j,i)=1./cpr(j,i)

         end if

         cpz(j,i) = g(i)/cp(j,i) + (1.-g(i))/cp(j,i+1)
         if (cpz(j,i) .lt. 1.e-09) then
            cpz(j,i) = 1.0e+10
         else
            cpz(j,i) = 1./cpz(j,i)
         end if

         sigr(j,i) = f(j)/sig(j,i) + (1.-f(j))/sig(j+1,i)
         if (sigr(j,i) .lt. 1.e-09) then
            sigr(j,i) = 1.0e+10
         else
            sigr(j,i) = 1./sigr(j,i)
         end if

         sigz(j,i) = g(i)/sig(j,i) + (1.-g(i))/sig(j,i+1)
         if (sigz(j,i) .lt. 1.e-09) then
            sigz(j,i) = 1.0e+10
         else
            sigz(j,i) = 1./sigz(j,i)
         end if

cls calculation of the thermal conductivity
         if (h(j,i+1).ne.h(j,i)) then
           zemz(j,i) = (schm(j,i+1)-schm(j,i))/(h(j,i+1)-h(j,i))
         else
           zemz(j,i) = zet(j,i)
         endif   
         if (h(j+1,i).ne.h(j,i)) then
           zemr(j,i) = (schm(j+1,i)-schm(j,i))/(h(j+1,i)-h(j,i))
         else
           zemr(j,i)= zet(j,i)
         endif   
         
         if(ica(j,i).eq.1.or.ian(j,i).eq.1)then
             zemz(j,i) = g(i)/zet(j,i) + (1.-g(i))/zet(j,i+1)
          if (zemz(j,i) .lt. 1.e-09) then
             zemz(j,i) = 1.0e+10
          else
             zemz(j,i) = 1./zemz(j,i)
          endif
c
          zemr(j,i) = f(j)/zet(j,i) + (1.-f(j))/zet(j+1,i)
          if (zemr(j,i) .lt. 1.e-09) then
             zemr(j,i) = 1.0e+10
          else
             zemr(j,i) = 1./zemr(j,i)
          endif
         endif
cls end of thermal conductivity calculation

         rhmr(j,i) = f(j)*rho(j+1,i) + (1.-f(j))*rho(j,i)
c need to avoid averaging rho at interface of weldpool; jjl
c         if (i.ne.nzam)  then
c             rhmz(j,i) = g(i)*rho(j,i+1) + (1.-g(i))*rho(j,i)
c         else
c             rhmz(j,i) = rho(j,i)
c         endif

c--etm viscosity half a control volume shifted
         etm(j,i) = f(j)*g(i)/eta(j,i) + (1.-f(j))*g(i)/eta(j+1,i)
     1   + f(j)*(1.-g(i))/eta(j,i+1) + (1.-f(j))*(1.-g(i))/eta(j+1,i+1)
         if (etm(j,i) .lt. 1.e-09) then
            etm(j,i) = 1.0e+10
         else
            etm(j,i) = 1./etm(j,i)
         end if
c
         etr(j,i) = f(j)/eta(j,i) + (1.-f(j))/eta(j+1,i)
         if (etr(j,i) .lt. 1.e-09) then
            etr(j,i) = 1.0e+10
         else
            etr(j,i) = 1./etr(j,i)
         end if
   66    continue
   65    continue

      do 70 i=0,nz
         rhmr(nrp,i) = rho(nrp,i)
         rhmz(nrp,i) = rho(nrp,i)
  70  continue

cls v.2 calcul of the electron diffusion current on the face

      do 165 i=0,nz
      do 166 j=0,nr

        zder1= zmob(j,i)*t(j,i)
        zder2= zmob(j+1,i)*t(j+1,i)

        zder = f(j)/zder1 + (1.-f(j))/zder2
        if (zder .lt. 1.e-15) then
            zder = 1.0e+16*boltz
        else
            zder = 1./zder*boltz
        end if
        cder(j,i)=zder*(rne(j+1,i)-rne(j,i))*2/(dr(j)+dr(j+1))

        zdez1= zmob(j,i)*t(j,i)
        zdez2= zmob(j,i+1)*t(j,i+1)

        zdez = g(i)/zdez1 + (1.-g(i))/zdez2
        if (zdez .lt. 1.e-15) then
           zdez = 1.0e+16*boltz
        else
           zdez = 1./zdez*boltz
        end if
        cdez(j,i)=zdez*(rne(j,i+1)-rne(j,i))*2/(dz(i)+dz(i+1))

c      no diffusion current in the electrodes
       if (ica(j,i).eq.1.or.ian(j,i).eq.1) then
          if (ica(j+1,i).eq.1.or.ian(j+1,i).eq.1) then
            cder(j,i)=0.0
          endif
          if (ica(j,i+1).eq.1.or.ian(j,i+1).eq.1) then
           cdez(j,i)=0.0
          endif
       endif

c boundary conditions added by jjl
        if (j.le.1) cder(j,i) = 0.0
        if (j.le.1) cdez(j,i) = 0.0

cls v.3 ambipolar diffusion
       diffar(j,i) = f(j)/diffa(j,i) + (1.-f(j))/diffa(j+1,i)
         if (diffar(j,i) .lt. 1.e-09) then
            diffar(j,i) = 1.0e+10
         else
            diffar(j,i) = 1./diffar(j,i)
         end if

       diffaz(j,i) = g(i)/diffa(j,i) + (1.-g(i))/diffa(j,i+1)
         if (diffaz(j,i) .lt. 1.e-09) then
            diffaz(j,i) = 1.0e+10
         else
            diffaz(j,i) = 1./diffaz(j,i)
         end if
 166  continue
 165  continue

c corection of the kappa value and the electron diffusion current at the 
c interface between the plasma and the electrode

      do 777 i = 1, nz
      do 778 j=1,nr
          if(ica(j,i) .gt. ica(j-1,i).and.j .gt. 1) then
          t1 = t(j,i)
          h1 = h(j,i)
          schm1=schm(j,i)
          t2=t(j-1,i)
          h2=h(j-1,i)
          schm2=schm(j-1,i)
          if (h1 .ne. h2) then
             zkcp = (schm2-schm1)/(h2-h1)
          else
             zkcp = zet(j-1,i)
          endif   
             zemr(j-1,i) = f(j-1)/zkcp + (1.-f(j-1))/zet(j,i)
          if (zemr(j-1,i) .lt. 1.e-09) then
             zemr(j-1,i) = 1.0e+10
          else
             zemr(j-1,i) = 1./zemr(j-1,i)
          endif          
          cder(j-1,i)=cder(j-1,i)*(dr(j)+dr(j-1))/dr(j-1)
          diffar(j-1,i)=diffar(j-1,i)*(dr(j)+dr(j-1))/dr(j-1)
        endif

        if (ica(j,i).gt.ica(j+1,i)) then
          t1=t(j,i)
          h1=h(j,i)
          schm1=schm(j,i)
          t2=t(j+1,i)
          h2=h(j+1,i)
          schm2=schm(j+1,i)
          if (h1.ne.h2) then
             zkcp = (schm2-schm1)/(h2-h1)
          else
             zkcp = zet(j+1,i)
          endif   
             zemr(j,i) = f(j)/zet(j,i) + (1.-f(j))/zkcp
          if (zemr(j,i) .lt. 1.e-09) then
             zemr(j,i) = 1.0e+10
          else
             zemr(j,i) = 1./zemr(j,i)
          endif          
          cder(j,i)=cder(j,i)*(dr(j)+dr(j+1))/dr(j+1)
          diffar(j,i)=diffar(j,i)*(dr(j)+dr(j+1))/dr(j+1)
        endif

        if(ica(j,i).gt.ica(j,i-1))then
          t1=t(j,i)
          h1=h(j,i)
          schm1=schm(j,i)
          t2=t(j,i-1)
          h2=h(j,i-1)
          schm2=schm(j,i-1)
          if (h1.ne.h2) then
             zkcp = (schm2-schm1)/(h2-h1)
          else
             zkcp = zet(j,i-1)
	endif   
             zemz(j,i-1) = g(i-1)/zkcp + (1.-g(i-1))/zet(j,i)
          if (zemz(j,i-1) .lt. 1.e-09) then
             zemz(j,i-1) = 1.0e+10
          else
             zemz(j,i-1) = 1./zemz(j,i-1)
          endif          
          cdez(j,i-1)= cdez(j,i-1)*(dz(i-1)+dz(i))/dz(i-1)
          diffaz(j,i-1)= diffaz(j,i-1)*(dz(i-1)+dz(i))/dz(i-1)
        endif

        if(ica(j,i).gt.ica(j,i+1))then
          t1=t(j,i)
          h1=h(j,i)
          schm1=schm(j,i)
          t2=t(j,i+1)
          h2=h(j,i+1)
          schm2=schm(j,i+1)
          if (h1.ne.h2) then
             zkcp = (schm2-schm1)/(h2-h1)
          else
             zkcp = zet(j,i+1)
          endif   
             zemz(j,i) = g(i)/zet(j,i) + (1.-g(i))/zkcp
          if (zemz(j,i) .lt. 1.e-09) then
             zemz(j,i) = 1.0e+10
          else
             zemz(j,i) = 1./zemz(j,i)
          endif          
          cdez(j,i)= cdez(j,i)*(dz(i)+dz(i+1))/dz(i+1)
          diffaz(j,i)= diffaz(j,i)*(dz(i)+dz(i+1))/dz(i+1)
        endif


c points for anode interface

        if(ian(j,i).gt.ian(j,i-1))then
          t1=t(j,i)
          h1=h(j,i)
          schm1=schm(j,i)
          t2=t(j,i-1)
          h2=h(j,i-1)
          schm2=schm(j,i-1)
          if (h1.ne.h2) then
             zkcp = (schm2-schm1)/(h2-h1)
          else
             zkcp = zet(j,i-1)
          endif   
             zemz(j,i-1) = g(i-1)/zkcp + (1.-g(i-1))/zet(j,i)
          if (zemz(j,i-1) .lt. 1.e-09) then
             zemz(j,i-1) = 1.0e+10
          else
             zemz(j,i-1) = 1./zemz(j,i-1)
          endif          
          cdez(j,i-1)= cdez(j,i-1)*(dz(i-1)+dz(i))/dz(i-1)
          diffaz(j,i-1)= diffaz(j,i-1)*(dz(i-1)+dz(i))/dz(i-1)
        endif

        if(ian(j,i).gt.ian(j,i+1))then
          t1=t(j,i)
          h1=h(j,i)
          schm1=schm(j,i)
          t2=t(j,i+1)
          h2=h(j,i+1)
          schm2=schm(j,i+1)
          if (h1.ne.h2) then
             zkcp = (schm2-schm1)/(h2-h1)
          else
             zkcp = zet(j,i+1)
          endif   
             zemz(j,i) = g(i)/zet(j,i) + (1.-g(i))/zkcp
          if (zemz(j,i) .lt. 1.e-09) then
             zemz(j,i) = 1.0e+10
          else
             zemz(j,i) = 1./zemz(j,i)
          endif          
          cdez(j,i)= cdez(j,i)*(dz(i)+dz(i+1))/dz(i+1)
          diffaz(j,i)= diffaz(j,i)*(dz(i)+dz(i+1))/dz(i+1)
        endif

        if(ian(j,i).gt.ian(j+1,i))then
          t1=t(j,i)
          h1=h(j,i)
          schm1=schm(j,i)
          t2=t(j+1,i)
          h2=h(j+1,i)
          schm2=schm(j+1,i)
          if (h1.ne.h2) then
             zkcp = (schm2-schm1)/(h2-h1)
          else
             zkcp = zet(j+1,i)
          endif   
             zemr(j,i) = f(j)/zet(j,i) + (1.-f(j))/zkcp
          if (zemr(j,i) .lt. 1.e-09) then
             zemr(j,i) = 1.0e+10
          else
             zemr(j,i) = 1./zemr(j,i)
          endif          
          cder(j,i)=cder(j,i)*(dr(j)+dr(j+1))/dr(j+1)
          diffar(j,i)=diffar(j,i)*(dr(j)+dr(j+1))/dr(j+1)
        endif

        if(ian(j,i).gt.ian(j-1,i))then
          t1=t(j,i)
          h1=h(j,i)
          schm1=schm(j,i)
          t2=t(j-1,i)
          h2=h(j-1,i)
          schm2=schm(j-1,i)
          if (h1.ne.h2) then
             zkcp = (schm2-schm1)/(h2-h1)
          else
             zkcp = zet(j-1,i)
          endif   
             zemr(j-1,i) = f(j-1)/zkcp + (1.-f(j-1))/zet(j,i)
          if (zemr(j-1,i) .lt. 1.e-09) then
             zemr(j-1,i) = 1.0e+10
          else
             zemr(j-1,i) = 1./zemr(j-1,i)
          endif          
          cder(j-1,i)=cder(j-1,i)*(dr(j)+dr(j-1))/dr(j-1)
          diffar(j-1,i)=diffar(j-1,i)*(dr(j)+dr(j-1))/dr(j-1)
        endif

 778   continue
 777   continue
 789   continue

       do 791 j = 1,nr
          cdez(j,0)=0.0
          cder(j,0) = 0.0
 791   continue
       do 792 i = 1,nz
          cder(nr,i) = 0.0
 792   continue
      return
      end
c end subroutine amf

c------------------------------------------------------------------------
c     subroutine ener
c     v.1 ls 18.9.98
c     
c     ener calculates the enthalpy of the domain of integration
c-------------------------------------------------------------------------

      subroutine ener(itoor,itor,itoz,htoor,irst)
      implicit double precision (a-h,o-z)
      logical actest
      parameter(mr=200, mz=100, nsh=21)
      dimension cc(mr,mz),dfc(mr,mz)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)      
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/cath/
     r  dtdh, dsdt(0:mr, 0:mz),rc
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr    
      common/hbound/
     r  hamb,hc,he,hanod(0:mr),hwall
      common/tcatht/
     r  setjz(0:mz),richjz(0:mz),setjr(0:mr),richjr(0:mr),caions
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
c---------------------------------------------------------------------------
c
c     step 1: calculate coefficients in energy equation.
      call fille
c
c     introduction of relaxation factor.
c
      rc1=0.1*rc
      do 10 i=1,nz
      do 11 j=1,nr
        if(ica(j,i).eq.0.and.ian(j,i).eq.0)then
         cc(j,i) = dc(j,i)/rlxt
         dfc(j,i) = df(j,i) + (1.-rlxt)*cc(j,i)*h(j,i)
         cc(j,i) = cc(j,i) + abs(dfc(j,i)/(rc*h(j,i)))
         dfc(j,i) = dfc(j,i) + abs(dfc(j,i)/rc) 
        else
         cc(j,i) = dc(j,i)
         dfc(j,i) = df(j,i)
         cc(j,i) = cc(j,i) + abs(dfc(j,i)/(rc1*h(j,i)))
         dfc(j,i) = dfc(j,i) + abs(dfc(j,i)/rc1) 
        endif
 11   continue
 10   continue
c
c     step 2: solve the energy equation.
c
      call solve(1,nr,1,nz,de,dw,dn,ds,cc,dfc,h,1,1)
c
      do 30 i=1,nz
         h(0,i) = h(1,i)
         h(nrp,i) = hamb
  30  continue

      do 35 j=nrca(1)+1,nr
         h(j,0) = hamb
  35	continue

	if(irca(1).gt.1)then
      do 935 j=0,irca(1)-1
         h(j,0) = hamb
 935  continue
	endif

      rese = 0.d0
      rele = 0.d0
      do 40 i=1,nz
      do 41 j=1,nr
	  east= de(j,i)*h(j+1,i)/dc(j,i)
          west= dw(j,i)*h(j-1,i)/dc(j,i)
          rnorth= dn(j,i)*h(j,i+1)/dc(j,i)
          south= ds(j,i)*h(j,i-1)/dc(j,i)
          centre= df(j,i)/dc(j,i)
          res = east+west+south+rnorth+centre - h(j,i)
          east = dabs(east)
          west = dabs(west)
          rnorth = dabs(rnorth)
          south = dabs(south)
          centre = dabs(centre)
          ph = dabs(h(j,i))
          res = dabs(res)
          rel = res/dmax1(east,west,south,rnorth,centre,ph)

c--- res in 0.1 %
        if( abs(h(j,i)) .ge. 1.d-6 ) then 
		res=res/h(j,i)
        else
		res = res/1.d-6
        endif

        res = abs(res)*1000.d0

         if (res .gt. rese) then
            rese = res
            ie = i
            je = j
         end if
        if (rel .gt.rele) then
            rele = rel
            iie = i
            jje = j
        endif
  41  continue
  40  continue

c--- translate enthalpy into temperature
      call htemp(itoor,itor,itoz,htoor)

      return
      end

c--------------------------------------------------------------
c end of subroutine ener
c--------------------------------------------------------------

c--------------------------------------------------------------
c     subroutine fille
c     v.1 ls 18.9.98
c
c     fille calculates the coefficients needed for the solution
c     of the enthalpy equation
c--------------------------------------------------------------

      subroutine fille
      implicit double precision (a-h,o-z)
      logical actest,rvpol
      parameter(mr=200, mz=100, nsh=21)
      dimension cdrp(0:mr,0:mz),cdzp(0:mr,0:mz)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz) 
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/pol/
     r  rvpol
      common/cath/
     r  dtdh, dsdt(0:mr, 0:mz),rc
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi
      common/tcatht/
     r  setjz(0:mz),richjz(0:mz),setjr(0:mr),richjr(0:mr),caions
      common/hbound/
     r  hamb,hc,he,hanod(0:mr),hwall
      common/emiss/
     r  emissc,emissa
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai

      data bke25/2.154337d-04/
      stefbolt2 = 5.669d-12
c------------------------------------------------------------------------

c      emiss = 0.3
      nzcp = nzc + 1
      do 2 i=0,mz
      do 3 j=0,mr
 	dsdt(j,i)=0.
  3   continue
  2   continue
      do 5 i=1,nz
      do 6 j=1,nr
	 dc(j,i) = 0.d0
	 de(j,i) = 0.d0
         dw(j,i) = 0.d0
         dn(j,i) = 0.d0
         ds(j,i) = 0.d0
	 df(j,i)=  0.d0
         vol(j,i) = an(j)*dz(i)
         cdzp(j,i) = 0.5*(cdz(j,i)+cdz(j,i-1))
         cdrp(j,i) = 0.5*(cdr(j,i)+cdr(j-1,i))
  6   continue
  5   continue

c---thermal haet flux r-direction
      do 10 i=1,nz
      do 11 j=0,nr
         temp = 2.*ae(j,i)*zemr(j,i)/(dr(j)+dr(j+1))
         if (j .ne. 0) de(j,i) = temp
         dw(j+1,i) = temp
  11  continue
  10  continue

c---thermal haet flux z-direction
      do 20 i=0,nz
      do 21 j=1,nr
         temp = 2.*an(j)*zemz(j,i)/(dz(i)+dz(i+1))
         if (i .ne. 0) dn(j,i) = temp
         ds(j,i+1) = temp
  21  continue 
  20  continue 

c---convection r-diretion
         do 40 i=0,nz
         do 41 j=0,nr
          if(ica(j,i).eq.1.or.ian(j,i).eq.1)then
              call cpsolid(j,i,t(j,i),ian(j,i),ica(j,i),cpf,nza)
          else
              cpf=1.
          endif
          if(i.eq.0)go to 30
          temp= ae(j,i)*rhmr(j,i)*vr(j,i)
          if (j.ne.0) de(j,i) = de(j,i) + cpf*dmax1(0.d0,-temp)
          dw(j+1,i) = dw(j+1,i) + cpf*dmax1(0.d0,temp)

c---convection z-direction
  30     continue
	if(j.eq.0)go to 40
	  temp= an(j)*rhmz(j,i)*vz(j,i)
	  if (i.ne.0) dn(j,i)= dn(j,i)+ cpf*dmax1(0.d0,-temp)
	  ds(j,i+1)= ds(j,i+1) +cpf*dmax1(0.d0,temp)
  41	continue	
  40	continue	

c--- ohmic heating
        do 50 i=1,nz
        do 51 j=1,nr
        if (ica(j,i).eq.1.or.ian(j,i).eq.1) then   
         temp = cdzp(j,i)*cdzp(j,i)+cdrp(j,i)*cdrp(j,i)
         df(j,i) = temp*vol(j,i)/sig(j,i)
        else
         temp1=cdz(j,i)*2*(phi(j,i+1)-phi(j,i))/(dz(i+1)+dz(i))
         if(ica(j,i+1).eq.1.or.ian(j,i+1).eq.1)then
         temp1=cdz(j,i)*2*(phi(j,i+1)-phi(j,i))/(dz(i))
         endif

         temp2=cdz(j,i-1)*2*(phi(j,i)-phi(j,i-1))/(dz(i)+dz(i-1))
         if(ica(j,i-1).eq.1.or.ian(j,i-1).eq.1)then
         temp2=cdz(j,i-1)*2*(phi(j,i)-phi(j,i-1))/(dz(i))
         endif

         temp3=cdr(j,i)*2*(phi(j+1,i)-phi(j,i))/(dr(j+1)+dr(j))
         if(ica(j+1,i).eq.1.or.ian(j+1,i).eq.1)then
         temp3=cdr(j,i)*2*(phi(j+1,i)-phi(j,i))/(dr(j))
         endif

         temp4=cdr(j-1,i)*2*(phi(j,i)-phi(j-1,i))/(dr(j)+dr(j-1))
         if(ica(j-1,i).eq.1.or.ian(j-1,i).eq.1)then
         temp4=cdr(j-1,i)*2*(phi(j,i)-phi(j-1,i))/(dr(j))
         endif

         df(j,i)=-0.5*(temp1+temp2+temp3+temp4)*vol(j,i)
        endif 
        eohm(j,i)=df(j,i)/vol(j,i)
        if(actest)then
	call cpsolid(j,i,t(j,i),ian(j,i),ica(j,i),cpf,nza)
        df(j,i)=df(j,i)+rho1(j,i)*cpf*h1(j,i)*vol(j,i)/delt
	endif
  51    continue
  50    continue
c
c--- pfender term
c         do 60 i=1,nz
c        do 60 j=1,nr
c          temp = -bke25*cdr(j,i)*ae(j,i)/cpr(j,i)
c          de(j,i) = de(j,i) +dmax1(0.d0,-temp)
c           dw(j+1,i)=dw(j+1,i) +dmax1(0.d0,temp)
c   60   continue
c        do 70 i=1,nz
c        do 70 j=1,nr
c          temp = -bke25*cdz(j,i)*an(j)/cpz(j,i)
c          dn(j,i)=dn(j,i)+dmax1(0.d0,-temp)
c          ds(j,i+1)=ds(j,i+1)+dmax1(0.d0,temp)
c   70    continue
c	
crm source terms include ion heating and thermionic cooling at cathode. 
crm add blackbody radiation.
        if(rvpol)then
            call rvpo
        else
            call stpo
        endif
c-----------------------
      boltz     = 1.38044d-23

c--- radiation as  source term; cast to the left for stability
      do 80 i=1,nz
      do 81 j=1,nr
         dc(j,i) = dc(j,i)+de(j,i) + dw(j,i) + dn(j,i) + ds(j,i) +
     1        (rad(j,i)/h(j,i) - dsdt(j, i)*dtdh)*vol(j,i)
        if(actest)then
	call cpsolid(j,i,t(j,i),ian(j,i),ica(j,i),cpf,nza)
        dc(j,i)=dc(j,i)+rho1(j,i)*cpf*vol(j,i)/delt
        endif
   81 continue
   80 continue
c----------

c---boundary condition
      do 90 i=1,nz

c--dh/dr =0. at r= 0.
         dc(1,i) = dc(1,i) - dw(1,i)
         dw(1,i) = 0.d0
c--h=ho at boundary.
         df(nr,i) = df(nr,i) + de(nr,i)*h(nrp,i)
         de(nr,i) = 0.d0
   90 continue
      do 100 j=1,nr
         df(j,nz) = df(j,nz) + dn(j,nz)*h(j,nzp)
         dn(j,nz) = 0.d0
         df(j,1) = df(j,1) + ds(j,1)*h(j,0)
         ds(j,1) = 0.d0
  100 continue
 1234    format(i3,x,f10.3)
      return
      end
c end of subroutine fil
c------------------------------------------------------------

c     subroutine htemp
c     v.1 ls 18.9.98;    htemp converts the henthalpy in temperature

      subroutine htemp(itoor,itor,itoz,htoor)
      implicit double precision (a-h,o-z)
      logical actest
      parameter(mr=200, mz=100, nsh=21)
      dimension hah(0:50),idiv(0:15)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
c-------------------------------------------------------------

        hah(0)=ah(1)
     	do 5 i=1,41
	  	hah(i)=ah(i)
  5  	continue
     	do 10 i=42,50
	  hah(i)=ah(41)+ah(40)
  10	continue

          idiv(0)=1	
        do 15 i=1,10
	  idiv(i)=2*idiv(i-1) 
  15 	continue
c
c--- main loop convert all h(j,i) to t(j,i)
      	do 20 i=0,nz+1
      	do 21 j=0,nr+1

    	  if( h(j,i).lt. hah(1)  ) then
c	h(j,i)=hah(1)
	 if (jscreen.eq.1) write(*,1) j,i,h(j,i)
		write(*,1) j,i,h(j,i)
	itor=j
	itoz=i
	htoor=h(j,i)
	itoor=1
                return
	  else
	  itoor=0
  
          endif

  1       format(1x,' enthalpie out of range  node j= '
     &    ,i2,' i= ',i2,' h= ',e12.4 )

    	  if( h(j,i).gt. hah(41)  ) then
	  if (jscreen.eq.1) write(*,1) j,i,h(j,i)
		write(461,1) j,i,h(j,i)
		h(j,i)=ah(41)
	   endif

c---index loop: find index 1-41 for cp(lh) in matf  and base temp.
c*********
          lho=16
          lh=16
	 do 30 ih=1,6

	    diffh=h(j,i)-hah(lh)
c	    write(*,*)' ih = ',ih,' diffh= ',diffh

c--- check if plus or minus half interval
	    if (diffh .ge. 0.d0 ) then

c---upper interval			
c------------- already found ?
	      if(hah(lh+1)-hah(lh) .gt. diffh)  goto 99

c------------- half interval
	      lh=lh+lho/idiv(ih)	

	    else

c--- lower interval
c-------------already found ?
              if(hah(lh-1) -hah(lh) .lt. diffh) then
			lh=lh-1
			diffh=h(j,i)-hah(lh)
			goto 99
	      endif

c--------------half interval
	      lh=lh-lho/idiv(ih)	

            endif

  30  	  continue

c---end of index loop


  99      alpha=(acp(lh+1)-acp(lh))
  2     format(1x,' enthalpie  j= ',i2,' i= ',i2,' h= ',e12.4 )

	  dh=hah(lh+1)-hah(lh)
            t(j,i)=diffh/dh*1000. +(lh-1)*1000.
  21 	continue
  20 	continue


c---end of main loop


  100  format('  i = ',i4)
  110  format(14f6.0)	       
        end

c-----------------------------------------------------------------
c end of subroutine htemp
c-----------------------------------------------------------------

c-----------------------------------------------------------------
c     subroutine cpsolid
c     v.1 ls 18.9.98
c
c     cpsolid replaces plasma cp by electrode cp  for flow in electrodes
c-----------------------------------------------------------------
        subroutine cpsolid(j,i,t,ian,ica,cpf,nza)
c replaces plasma cp by electrode cp
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz) 
      common/liqu/  asigliqu(41),tmla,rholiqu,
     1  rhocoeff,cpliqu,cpcoeff,npr  

        s=t/1000.
        l=int(s)
        l1=l+1
        l2=l+2
        frac=1.-s
        rfr=1.-frac
        cpg = cp(j,i)
        cps = cp(j,i)
        if(ian.eq.1) cps = ancp(l2)*frac + ancp(l1)*rfr
        if(ica.eq.1) cps = cacp(l2)*frac + cacp(l1)*rfr
        if(t.gt.tmla.and.i.ge.nza) then  
            cps = cpliqu + cpcoeff*(t - tmla)      
        end if
        cpf=cps/cpg
        end

c-------------------------------------------------------------
c end of subroutine cpsolid
c--------------------------------------------------------------

c--------------------------------------------------------------------
c     subroutine rvpo
c     v.1 ls 18.9.98
c--------------------------------------------------------------------

      subroutine rvpo
      implicit double precision (a-h,o-z)
      logical actest,rvpol  
      parameter(mr=200, mz=100, nsh=21)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz) 
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/pol/
     r  rvpol
      common/cath/
     r  dtdh, dsdt(0:mr, 0:mz),rc
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi
      common/tcatht/
     r  setjz(0:mz),richjz(0:mz),setjr(0:mr),richjr(0:mr),caions
      common/hbound/
     r  hamb,hc,he,hanod(0:mr),hwall
      common/emiss/
     r  emissc,emissa
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai

      data bke25/2.154337d-04/
c----------------------------------------------------------------------
      stefbolt2 = 5.669d-12
c      emiss = 0.3

crm source terms include ion heating and thermionic cooling at cathode.
crm add blackbody radiation.

        iz=0

          do 51 i=1,nz
          do 52 j=1,nr
	if(ica(j,i).eq.1)then
        if(ica(j,i).gt.ica(j,i-1)) call ansh(t(j,i),df(j,i),
     &  emissc,cdz(j,i-1),vol(j,i),dz(i),cworkf)
        if(ica(j,i).gt.ica(j,i+1)) call ansh(t(j,i),df(j,i),
     &  emissc,cdz(j,i),vol(j,i),dz(i),cworkf)
        if(ica(j,i).gt.ica(j+1,i)) call ansh(t(j,i),df(j,i),
     &  emissc,cdr(j,i),vol(j,i),dr(j),cworkf)
        if(ica(j,i).gt.ica(j-1,i)) call ansh(t(j,i),df(j,i),
     &  emissc,cdr(j-1,i),vol(j,i),dr(j),cworkf)
	go to 51
	endif
crm
c anode  electron heating and radiation cooling

        sdsdt=0.d0
	if(ian(j,i).eq.1)then
        if(ian(j,i).gt.ian(j+1,i))then
        call cash(potion,emissa,anrich,aworkf,t(j,i),cdr(j,i),
     &  sdsdt,df(j,i),dr(j),vol(j,i),h(j,i),
     &  setj,richj,dtdh,dcaions)
        dsdt(j,i)=dsdt(j,i)+sdsdt
        endif

        if(ian(j,i).gt.ian(j,i+1))then
        call cash(potion,emissa,anrich,aworkf,t(j,i),cdz(j,i),
     &  sdsdt,df(j,i),dz(i),vol(j,i),h(j,i),
     &  setj,richj,dtdh,dcaions)
        dsdt(j,i)=dsdt(j,i)+sdsdt
        iz=iz+1
        setjz(iz)=setj
        richjz(iz)=richj
        endif

        if(ian(j,i).gt.ian(j-1,i).and.j.gt.2)then
        call cash(potion,emissa,anrich,aworkf,t(j,i),cdr(j-1,i),
     &  sdsdt,df(j,i),dr(j),vol(j,i),h(j,i),
     &  setj,richj,dtdh,dcaions)
        dsdt(j,i)=dsdt(j,i)+sdsdt
        endif

        if(ian(j,i).gt.ian(j,i-1).and.i.lt.nz-1)then
        call cash(potion,emissa,anrich,aworkf,t(j,i),cdz(j,i-1),
     &  sdsdt,df(j,i),dz(i),vol(j,i),h(j,i),
     &  setj,richj,dtdh,dcaions)
        dsdt(j,i)=dsdt(j,i)+sdsdt
        iz=iz+1
        setjz(iz)=setj
        richjz(iz)=richj
        endif

        dsdt(j,i)=-dmax1(0.d0,-dsdt(j,i))
        df(j,i)=df(j,i)-dsdt(j,i)*dtdh*h(j,i)*vol(j,i)
	endif
 52      continue
 51      continue
      return
      end

c--------------------------------------------------------------------------
c end of subroutine rvpo
c-------------------------------------------------------------------------

c---------------------------------------------------------------------------
c     subroutine stpo
c     v.1 ls 18.9.98
c---------------------------------------------------------------------------

      subroutine stpo
      implicit double precision (a-h,o-z)
      logical actest,rvpol
      parameter(mr=200, mz=100, nsh=21)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz) 
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/pol/
     r  rvpol
      common/cath/
     r  dtdh, dsdt(0:mr, 0:mz),rc
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi
      common/tcatht/
     r  setjz(0:mz),richjz(0:mz),setjr(0:mr),richjr(0:mr),caions
      common/hbound/
     r  hamb,hc,he,hanod(0:mr),hwall
      common/emiss/
     r  emissc,emissa
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai

      data bke25/2.154337d-04/
c----------------------------------------------------------------

      stefbolt2 = 5.669d-12
c      emiss = 0.3
      nzcp = nzc + 1

crm source terms include ion heating and thermionic cooling at cathode.
crm add blackbody radiation.
      electron  = 1.60206d-19
      boltz     = 1.38044d-23
      b0        = cworkf*electron/boltz
      caions = 0.
	iz=0
        iiz=0
      do 51 i=1,nz
      do 52 j=1,nr
        sdsdt=0.d0
	

	if(ica(j,i).eq.1)then

        if(ica(j,i).gt.ica(j+1,i))then
        call cash(potion,emissc,carich,cworkf,t(j,i),cdr(j,i),
     &  sdsdt,df(j,i),dr(j),vol(j,i),h(j,i),
     &  setj,richj,dtdh,dcaions)
        caions = caions + dcaions
        dsdt(j,i)=dsdt(j,i)+sdsdt
        iiz=iiz+1
        setjr(iiz)=setj
        richjr(iiz)=richj
        endif

        if(ica(j,i).gt.ica(j,i+1))then
        call cash(potion,emissc,carich,cworkf,t(j,i),cdz(j,i),
     &  sdsdt,df(j,i),dz(i),vol(j,i),h(j,i),
     &  setj,richj,dtdh,dcaions)
        caions = caions + dcaions
        dsdt(j,i)=dsdt(j,i)+sdsdt
        iz=iz+1
        setjz(iz)=setj
        richjz(iz)=richj
        endif

        if(ica(j,i).gt.ica(j-1,i).and.j.gt.2)then
        call cash(potion,emissc,carich,cworkf,t(j,i),cdr(j-1,i),
     &  sdsdt,df(j,i),dr(j),vol(j,i),h(j,i),
     &  setj,richj,dtdh,dcaions)
        dsdt(j,i)=dsdt(j,i)+sdsdt
        endif

        if(ica(j,i).gt.ica(j,i-1).and.i.gt.2)then
        call cash(potion,emissc,carich,cworkf,t(j,i),cdz(j,i-1),
     &  sdsdt,df(j,i),dz(i),vol(j,i),h(j,i),
     &  setj,richj,dtdh,dcaions)
        dsdt(j,i)=dsdt(j,i)+sdsdt
        iz=iz+1
        setjz(iz)=setj
        richjz(iz)=richj
        endif

        dsdt(j,i)=-dmax1(0.d0,-dsdt(j,i))
        df(j,i)=df(j,i)-dsdt(j,i)*dtdh*h(j,i)*vol(j,i)

	go to 51
	endif
crm
c anode  electron heating and radiation cooling
	if(ian(j,i).eq.1)then
        if(ian(j,i).gt.ian(j,i-1)) call ansh(t(j,i),df(j,i),
     &  emissa,cdz(j,i-1),vol(j,i),dz(i),aworkf)
        if(ian(j,i).gt.ian(j,i+1)) call ansh(t(j,i),df(j,i),
     &  emissa,cdz(j,i),vol(j,i),dz(i),aworkf)
        if(ian(j,i).gt.ian(j+1,i)) call ansh(t(j,i),df(j,i),
     &  emissa,cdr(j,i),vol(j,i),dr(j),aworkf)
        if(ian(j,i).gt.ian(j-1,i)) call ansh(t(j,i),df(j,i),
     &  emissa,cdr(j-1,i),vol(j,i),dr(j),aworkf)
	endif
 52      continue
 51      continue
      return
      end

c----------------------------------------------------------------------
c end of subroutine stpo
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c     subroutine ansh
c     v.1 18.9.98
c
c     energy balance for the anode sheath
c----------------------------------------------------------------------

        subroutine ansh(t,df,emissa,cdz,vol,dz,aworkf)

          implicit double precision (a-h,o-z)

          stefbolt2 = 5.669d-12
          coolstef  = emissa*stefbolt2*t**4
          setj      = dabs(cdz)
          df = df + (dabs(setj)*aworkf - coolstef)*vol/dz
        end
c----------------------------------------------------------------------
c end of subroutine ansh
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c     subroutine cash
c     v.1 18.9.98
c
c     energy balance for the cathode sheath
c----------------------------------------------------------------------

        subroutine cash(potion1,emissc,carich1,cworkf1,t2,cd,
     &  sdsdt,df,dl,vol,h,setj,richj,dtdh,dcaions)

      implicit double precision (a-h,o-z)
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi

      stefbolt2 = 5.669d-12
      electron  = 1.60206d-19
      boltz     = 1.38044d-23
      b0        = cworkf*electron/boltz
      coolstef  = emissc*stefbolt2*t2**4
      richj     = carich*t2*t2*dexp(-b0/t2)
      setj      = dabs(cd)
      currion   = setj*1.d-2*0.0
        sdsdt=0.0
      dcaions=0.
      if(richj.lt.setj) then
      xcurrion   = setj - richj
      if(xcurrion.gt.currion) currion = xcurrion
      sdsdt =  - ((2.*t2 + b0)*carich*
     *(potion + cworkf)*exp(-b0/t2) - 4*emissc*stefbolt2*t2**3)/dl
      endif
c      currion=0.0
      curelec   = setj - currion
      dcaions = currion*potion*vol/dl
        df = df
     * + (currion*potion - dabs(curelec)*cworke - coolstef)*vol/dl
c
        end
c end of subroutine cash


c------------------------------------------------------------------------------
c     subroutine mom
c     v.1 ls 18.9.97
c
c     mom calculates the axial and radial velocities and pressure
c     according to the simplec method of patankar
c     (s.v. patankar, num. heat transf., vol. 4, p409-425, 1981)
c     updated by van doormal and raithby (numerical heat transfer)
c------------------------------------------------------------------------------

      subroutine mom(flowin,irst)
      implicit double precision (a-h,o-z)
      logical actest
      parameter(mr=200, mz=100, nsh=21)
       dimension tn(mr,mz),tee(mr,mz),bal(mr,mz),pp(0:mr,0:mz)
      common/surface/
     r  rsf(500),zsf(500),rsfo(0:mr,0:mz),
     1  zsfo(0:mr,0:mz),pstr(0:mr,0:mz),pstz(0:mr,0:mz)
      common/liqu/  asigliqu(41),tmla,rholiqu,
     1  rhocoeff,cpliqu,cpcoeff,npr 
      common/vcoef/
     r   aen(mr,mz),aes(mr,mz),aee(mr,mz),aew(mr,mz),
     1  aec(mr,mz),bsr(mr,mz),ann(mr,mz),ans(mr,mz),ane(mr,mz),
     2  anw(mr,mz),anc(mr,mz),bsz(mr,mz)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm, 
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip,
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr  
	common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/emiss/
     r  emissc,emissa
c-------------------------------------------------------------------------------

      do 2 i=0,mz
      do 1 j=0,mr
         pp(j,i) = 0.
  1   continue
  2   continue
      do 4 i=1,mz
      do 3 j=1,mr
         bal(j,i) = 0.
         df(j,i) = 0.
         tn(j,i) = 0.
         tee(j,i) = 0.
  3   continue
  4   continue
c
c   step 1: calculate coefficients in momentum and pressure equations

c      call fillv
c
c     put in boundary conditions.
c     constant pressure at outer boundary and vr = 0 at centre
        do 5 i=1,nz
        ane(nr,i) = 0.
        aew(1,i) = 0.
        aec(nrm,i) = aec(nrm,i) - aee(nrm,i)*rhmr(nrm,i)*ae(nrm,i)/
     1              (rhmr(nr,i)*ae(nr,i))
        if (actest) bsr(nrm,i) = bsr(nrm,i) +
     1              aee(nrm,i)*(rho1(nr,i)-rho(nr,i))*vol(nr,i)
     2              /(rhmr(nr,i)*ae(nr,i)*delt)
        aee(nrm,i) = 0.
 5      continue
        do 6 j=1,nr
            aes(j,1) = 0.
            aen(j,nz) = 0.
            ann(j,nzm)=0.
            bsz(j,1) = bsz(j,1) + ans(j,1)*vz(j,0)
            ans(j,1) = 0.
 6      continue
c        do 9 i=1,nz
c        do 10 j=1,nrm
c            dc(j,i) = aec(j,i)/rlxvr
c            df(j,i) = bsr(j,i) + (1.-rlxvr)*dc(j,i)*vr(j,i)
c            df(j,i)=df(j,i)   + ae(j,i)*(p(j,i)- p(j+1,i))
c 10     continue
c  9     continue
c      do 40 i=1,nz 
c      do 44 j=1,nr
c        j1 = j+1
c        if  (t(j1,i).lt.tmla) vr(j,i) = 0.0
c        if  (t(j1,i).lt.tmla) then
c            call zerocoef(dc,df,aee,aew,aen,aes,0.d0,j,i)
c        end if
c  44  continue
c  40  continue
c  vr velocities at upper surface of cathode can be set to zero; jjl
c      do 44 j=1,nrc
c          call zerocoef(dc,df,aee,aew,aen,aes,0.d0,j,1)
c  44  continue
c test for viscous drag; jjl
c vr velocities above anode can be set to zero;  jjl
c      do 49 j=1,nr
c          call zerocoef(dc,df,aee,aew,aen,aes,0.d0,j,nza-1)
c 49   continue
c vr velocities at outer surfaces of cathode set to zero;  jjl
c      do 46 i=1,nzc
c          j = nrca(i)
c          call zerocoef(dc,df,aee,aew,aen,aes,0.d0,j,i)
c           vr(j,i) = 0.0
c 46   continue
c
c     step 3: solve the 2 momentum equations

c        if (jscreen.eq.1) write(*,991)
 991    format(1x,'<solving for vr>')
c        call solve(1,nrm,1,nz,aee,aew,aen,aes,dc,df,vr,0,1)
cccccccccccccccccccccccccccccccccccc
c        do 27 i=1,nzm
c        do 28 j=1,nr
c            dc(j,i) = anc(j,i)/rlxvz
c            df(j,i) = bsz(j,i) + (1.-rlxvz)*dc(j,i)*vz(j,i)
c            df(j,i)=df(j,i)  + an(j)*(p(j,i)-p(j,i+1))
c  28    continue
c  27    continue
c         nzam = nza-1
c        if (jscreen.eq.1) write(*,992)
 992    format(1x,'<solving for vz>')
c        call solve(1,nr,1,nzm,ane,anw,ann,ans,dc,df,vz,0,1)
c        do 70 i=1,nz
c        do 71 j=1,nr
c        if (i .ne. nz) then
c            tn(j,i) = an(j)/(dc(j,i)-ann(j,i)-ans(j,i)-anw(j,i)-
c     1              ane(j,i))
c        else
c            tn(j,i) = 0.
c        end if
c 71     continue
c 70     continue

c     step 4: calculate mass balance for pressure correction.
c
         call mbal(vr,vz,bal,1)
c
c     step 6: calculate pressure correction.
c
        if (actest) then
        do 110 i=1,nz
        do 111 j=1,nr
        df(j,i) = bal(j,i) + (rho1(j,i)-rho(j,i))*vol(j,i)/delt
 111    continue
 110    continue
        else
        do 120 i=1,nz
        do 121 j=1,nr
            df(j,i) = bal(j,i)
 121    continue
 120    continue
        end if
c
c     step 5: calculate coefficients for pressure correction equation.
c
c        do 130 i=1,nz
c        do 131 j=1,nr
c            temp1 = rhmr(j,i)*tee(j,i)*ae(j,i)
c            de(j,i) = temp1
c            dw(j+1,i) = temp1
c            temp2 = rhmz(j,i)*tn(j,i)*an(j)
c            dn(j,i) = temp2
c            ds(j,i+1) = temp2
c 131    continue
c 130    continue

        do 290 i=1,nz
c            vz(0,i) = vz(1,i)
c            vr(0,i) = -vr(1,i)
c            vz(nrp,i) = vz(nr,i)
c        vr(nr,i) = vr(nrm,i)*rhmr(nrm,i)*ae(nrm,i)/
c     1          (rhmr(nr,i)*ae(nr,i))
            p(nrp,i) = ppp(i)*1.0E6
 290    continue
        p(nrp,0) = p(nrp,1) + 0.1E7
        do 260 i=0,nz
        do 261 j=0,nr
            p(j,i) = p(nrp,i) 
 261    continue
 260    continue
 510  format('  i =',i4,'   z =',e13.5)
 511  format(10E9.2)
        open(7,file='rhoz1')
      do 775 i=0,nz
         write(7,510) i,z(i)
c         write(7,511) (rho(j,i),j=0,nrp)
         write(7,511) (p(j,i),j=0,nrp)
 775  continue
        close(7)
c vz velocities set to maximum value from (1/2)rho vz**2 = pressdure drop
c vz and vr velocities set to zero wihin cathode and anode
      do 64 i=0,nz
      do 65 j=1,nr
         vz(j,i) = sqrt(p(nr,0)-p(j,i))*2.0/rho(j,i)
         if (i.ge.nza.and.j.le.nra) vz(j,i) = 0.0
         if (i.le.nzc.and.j.le.nrca(i)) vz(j,i) = 0.0
         if (i.ge.nza.and.j.le.nra) vr(j,i) = 0.0
         if (i.le.nzc.and.j.le.nrca(i)) vr(j,i) = 0.0
  65  continue
  64  continue  
c radial velocities in sold cathode set to zero  
c      do 9 i=1,nzc 
c      do 10 j=1,nrca(j)
c          vr(j,i) = 0.0
c  10  continue
c   9  continue
c radial velocities in sold anode set to zero  
c      do 40 i=nza,nz 
c      do 44 j=1,nra
c          vr(j,i) = 0.0
c  44  continue
c  40  continue

c Boundary condition: velocities of input flow vz(j,0)
        do 101 j=1,nr
            vz(j,0) = vz(j,1)
            vr(j,0) = vr(j,1) 
            vz(j,nz+1) = vz(j,nz)
            vr(j,nz+1) = vr(j,nz)
 101    continue
        do 100 i=1,nz
            vz(0,i) = vz(1,i)
            vr(0,i) = 0.0
            vz(nrp,i) = vz(nr,i)
            vr(nrp,i) = vz(nr,i)
 100    continue
c vz velocities within cathode set to zero;  jjl
c      do 47 i=1,nzc
c      do 48 j=1,nrca(i)
c           vz(j,i) = 0.0
c 48   continue
c 47   continue
c axial velocities in solid anode set to zero
c      do 45 i=nza,nz 
c          i1 = i+1
c      do 48 j=1,nr
c       if (t(j,i1).lt.tmla) then
c           call zerocoef(dc,df,ane,anw,ann,ans,0.d0,j,i)
c           vz(j,i) = 0.0
c       end if
c  48  continue
c  45  continue
c
c        do 280 j=1,nr
c        vz(j,nz) = 0.
c        vr(j,nzp) = 0.
c        vr(j,0) = 0.
c        p(j,0) = p(j,1)
c 280    continue

c        if (actest) then
c        do 300 i=1,nz
c        vr(nr,i) = vr(nr,i) +(rho1(nr,i)-rho(nr,i))*an(nr)*dz(i)
c     1            /(rhmr(nr,i)*ae(nr,i)*delt)
c 300    continue
c        end if
        
c     calculate residuals for mass balance.
c
        resp = 0.d0
        relp=0.d0
        do 220 i=1,nz
        do 221 j=1,nr
          east= de(j,i)*pp(j+1,i)/dc(j,i)
          west= dw(j,i)*pp(j-1,i)/dc(j,i)
          rnorth= dn(j,i)*pp(j,i+1)/dc(j,i)
          south= ds(j,i)*pp(j,i-1)/dc(j,i)
          centre= df(j,i)/dc(j,i)
          res = east+west+south+rnorth+centre - pp(j,i)
          east = dabs(east)
          west = dabs(west)
          rnorth = dabs(rnorth)
          south = dabs(south)
	centre = dabs(centre)
	ph = dabs(pp(j,i))
          res = dabs(res)
          aaa=0.001
          if (ph.gt.aaa) rel=res/dmax1(east,west,south,rnorth,centre,ph)
         if (res .gt. resp.and.ica(j,i).eq.0) then
            resp = res
            ip = i
            jp = j
         end if
        if (rel .gt.relp) then
            relp = rel
            iip = i
            jjp = j
        endif
 221    continue
 220    continue
c
c     calculate residuals for the momentum equation.
c
      resvr = 0.
      relvr = 0.
      resvz = 0.
      relvz = 0.
      do 321 i=1,nzm
          if (i.eq.nza-1) go to 321
          i1 = i+1
      do 320 j=1,nr
          if(ica(j,i).eq.1)go to 320
          if (t(j,i1).lt.tmla) then 
              go to 320
          end if
          east= ane(j,i)*vz(j+1,i)/anc(j,i)
          west= anw(j,i)*vz(j-1,i)/anc(j,i)
          rnorth= ann(j,i)*vz(j,i+1)/anc(j,i)
          south= ans(j,i)*vz(j,i-1)/anc(j,i)
          centre= bsz(j,i)/anc(j,i)
          bss = an(j)*(p(j,i)-p(j,i+1))/anc(j,i)
          res = east+west+south+rnorth+centre+bss - vz(j,i)
         if (res .gt. resvz) then
            resvz = res
            ivz = i
            jvz = j
         end if
          east = dabs(east)
          west = dabs(west)
          rnorth = dabs(rnorth)
          south = dabs(south)
          centre = dabs(centre)
          bss = dabs(bss)
          ph = dabs(vz(j,i))
          res = dabs(res)
          rel = res/dmax1(east,west,south,rnorth,centre,bss,ph)
        if (rel .gt.relvz) then
            relvz = rel
            iivz = i
            jjvz = j
        endif
 320  continue
 321  continue
        do 330 i=1,nz
        do 331 j=1,nrm
          if(ica(j,i).eq.1)go to 330
          j1 = j+1
          if  (t(j1,i).lt.tmla) then
            go to 330
          end if
          east= aee(j,i)*vr(j+1,i)/aec(j,i)
          west= aew(j,i)*vr(j-1,i)/aec(j,i)
          rnorth= aen(j,i)*vr(j,i+1)/aec(j,i)
          south= aes(j,i)*vr(j,i-1)/aec(j,i)
          centre= bsr(j,i)/aec(j,i)
          bss = ae(j,i)*(p(j,i)-p(j+1,i))/aec(j,i)
          res = east+west+south+rnorth+centre+bss - vr(j,i)
          east = dabs(east)
          west = dabs(west)
          rnorth = dabs(rnorth)
          south = dabs(south)
          centre = dabs(centre)
          bss = dabs(bss)
          ph = dabs(vr(j,i))
          res = dabs(res)
          rel = res/dmax1(east,west,south,rnorth,centre,bss,ph)
         if (res .gt. resvr) then
         resvr = res
         ivr = i
         jvr = j
         end if
        if (rel .gt.relvr) then
         relvr = rel
         iivr = i
         jjvr = j
        endif
 331    continue
 330    continue
        return
      end

c-----------------------------------------------------------------------------
c end subroutine mom
c-----------------------------------------------------------------------------

c-----------------------------------------------------------------------------
c     subroutine fillv
c     v.1 ls 18.9.98
c
c     this subroutine calculates the coefficients needed for the solution
c     of axial and radial velocities and pressure
c------------------------------------------------------------------------------

      subroutine fillv
      implicit double precision (a-h,o-z)
      logical actest, rvpol
      parameter(mr=200, mz=100, nsh=21)
      dimension fe(0:mr,0:mz),fw(0:mr,0:mz),fn(0:mr,0:mz),fs(0:mr,0:mz)
      common/vcoef/
     r   aen(mr,mz),aes(mr,mz),aee(mr,mz),aew(mr,mz),
     1  aec(mr,mz),bsr(mr,mz),ann(mr,mz),ans(mr,mz),ane(mr,mz),
     2  anw(mr,mz),anc(mr,mz),bsz(mr,mz)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/visc/
     r  eta(0:mr,0:mz),etm(0:mr,0:mz),etr(0:mr,0:mz)  
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/mag/
     r   bth(0:mr,0:mz),bjz(0:mr,0:mz),bjr(0:mr,0:mz)
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/pol/
     r  rvpol
      common/liqu/  asigliqu(41),tmla,rholiqu,
     1  rhocoeff,cpliqu,cpcoeff,npr 
c--------------------------------------------------------------------------------
      do 15 i=0,mz
      do 16 j=0,mr
         fe(j,i) = 0.
         fs(j,i) = 0.
         fn(j,i) = 0.
         fw(j,i) = 0.
  16  continue
  15  continue
      do 20 i=1,nz
      do 22 j=0,nrm
         temp1 = 2.*r(j+1)*dz(i)*eta(j+1,i)/dr(j+1)
         temp2 = 0.5*(rhmr(j,i)*vr(j,i)*ae(j,i) +
     1               rhmr(j+1,i)*vr(j+1,i)*ae(j+1,i))
         if (j .ne. 0) then
            de(j,i) = temp1
            fe(j,i) = temp2
         end if
         dw(j+1,i) = temp1
         fw(j+1,i) = temp2
  22  continue
  20  continue
      do 21 i=1,nzc
         dw(1,i) = 2.*ae(1,i)*eta(1,i)/dr(1)
  21  continue
      do 30 i=0,nz
      do 31 j=1,nrm
         temp1 = (r(j+1)**2-r(j)**2)*etm(j,i)/(dz(i)+dz(i+1))
         temp2 = 0.5*(rhmz(j,i)*vz(j,i)*an(j) +
     1               rhmz(j+1,i)*vz(j+1,i)*an(j+1))
         if (j .eq. 1) temp2 = temp2 + 0.5*rhmz(1,i)*vz(1,i)*an(1)
         if (i .ne. 0) then
            dn(j,i) = temp1
            fn(j,i) = temp2
         end if
         ds(j,i+1) = temp1
         fs(j,i+1) = temp2

  31  continue
  30  continue
      do 45 i=1,nz
      do 46 j=1,nrm
	rint1=rhmr1(j,i)
         volr = dz(i)*0.5*(r(j+1)**2-r(j)**2)
         bsr(j,i) = -bjz(j,i)*volr + (etm(j,i)*(vz(j+1,i)-vz(j,i)) -
     1         etm(j,i-1)*(vz(j+1,i-1)-vz(j,i-1)))*(r(j)+r(j+1))*0.5
c removal of jXB
c         if (i.ge.nza) bsr(j,i) = (etm(j,i)*(vz(j+1,i)-vz(j,i)) -
c     1         etm(j,i-1)*(vz(j+1,i-1)-vz(j,i-1)))*(r(j)+r(j+1))*0.5
      if (actest) bsr(j,i)= bsr(j,i)+rint1*volr*vr1(j,i)/delt
         aee(j,i) = de(j,i) + dmax1(0.d0,-fe(j,i))
         aew(j,i) = dw(j,i) + dmax1(0.d0,fw(j,i))
         aen(j,i) = dn(j,i) + dmax1(0.d0,-fn(j,i))
         aes(j,i) = ds(j,i) + dmax1(0.d0,fs(j,i))

        aec(j,i)=aes(j,i)+aen(j,i)+aee(j,i)+aew(j,i)
     1   +  2.*etr(j,i)*volr/(r(j)+0.5*dr(j))**2
         if (actest) aec(j,i) = aec(j,i) + rint1*volr/delt 
  46  continue
  45  continue
c
      do 50 i=0,nzm
      do 51 j=1,nr
         temp1 = 2.*an(j)*eta(j,i+1)/dz(i+1)
         temp2 = 0.5*(rhmz(j,i)*vz(j,i)+rhmz(j,i+1)*vz(j,i+1))*an(j)
         if (i .ne. 0) then
            dn(j,i) = temp1
            fn(j,i) = temp2
         end if
         ds(j,i+1) = temp1
         fs(j,i+1) = temp2
  51  continue
  50  continue
      do 60 i=1,nzm
      do 61 j=0,nr
         temp1 = (ae(j,i+1)+ae(j,i))*etm(j,i)/(dr(j)+dr(j+1))
         temp2 = 0.5*(rhmr(j,i)*vr(j,i)*ae(j,i) +
     1               rhmr(j,i+1)*vr(j,i+1)*ae(j,i+1))
         if (j .ne. 0) then
            de(j,i) = temp1
            fe(j,i) = temp2
         end if
         dw(j+1,i) = temp1
         fw(j+1,i) = temp2

  61  continue
  60  continue
      do 66 i=1,nz
         fw(1,i) = 0.
  66  continue
      do 75 i=1,nzm
      do 76 j=1,nr
	rint=rhmz(j,i)
	rint1=rhmz1(j,i)
         volz = 0.5*(vol(j,i)+vol(j,i+1))
         temp1 = bjr(j,i)*volz
         if (j .eq. 1) then
            temp2 = 0.
         else
            temp2 = etm(j,i)*(vr(j,i+1)-vr(j,i))*(r(j)+0.5*dr(j))
     1      - etm(j-1,i)*(vr(j-1,i+1)-vr(j-1,i))*(r(j)-0.5*dr(j))
         end if
c gravity can be removed
         bsz(j,i) = temp1 + temp2 + rint*volz*980.65
c         bsz(j,i) = temp1 + temp2 
c         if (i.ge.nza) bsz(j,i) = temp2 + rint*volz*980.65
c         if (i.ge.nza) bsz(j,i) = temp1 + temp2
      if (actest) bsz(j,i)=bsz(j,i)+rint1*volz*vz1(j,i)/delt
         ann(j,i) = dn(j,i) + dmax1(0.d0,-fn(j,i))
         ans(j,i) = ds(j,i) + dmax1(0.d0,fs(j,i))
         ane(j,i) = de(j,i) + dmax1(0.d0,-fe(j,i))
         anw(j,i) = dw(j,i) + dmax1(0.d0,fw(j,i))
         anc(j,i)=ans(j,i)+ann(j,i)+ane(j,i)+anw(j,i)
         if (actest) anc(j,i) = anc(j,i) + rint1*volz/delt      
  76  continue
  75  continue
      return
      end
c end of subroutine fillv

c--------------------------------------------------------------------------
c     subroutine mbal
c     v.1 ls 18.9.98;   calculates net mass flow out of each cell

      subroutine mbal(vr,vz,bal,its)
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      dimension vr(0:mr,0:mz),vz(0:mr,0:mz),bal(mr,mz)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)  

      do 20 i=1,nz
      do 21 j=1,nr
         flw = rhmr(j-1,i)*vr(j-1,i)*ae(j-1,i)
         fle = rhmr(j,i)*vr(j,i)*ae(j,i)
         fls = rhmz(j,i-1)*vz(j,i-1)*an(j)
         fln = rhmz(j,i)*vz(j,i)*an(j)
         bal(j,i) = flw - fle + fls - fln
 21   continue
 20   continue
      return
      end
c end of subroutine mbal

c------------------------------------------------------
c     subroutine zerocoef
c     v.1 ls 18.9.98;   sets variable to dfval at j,i

        subroutine zerocoef(dc,df,de,dw,dn,ds,dfval,j,i) 
        implicit double precision (a-h,o-z)
        parameter(mr=200, mz=100, nsh=21)
        dimension  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1      df(mr,mz)
                dc(j,i)=1.
                df(j,i)=dfval
                de(j,i)=0.d0
                dw(j,i)=0.d0
                dn(j,i)=0.d0
                ds(j,i)=0.d0
        end
c end of subroutine zerocoef

c-------------------------------------------------------------------
c     subroutine poiss
c     v.2 ls 4.11.98
c
c     this subroutine solves the current continuity equation and calculates
c     the electrical potential (phi) for the arc
c     and also calculates the axial (cdz) and radial (cdr) current 
c     densities and the magnetic field (bth).
c     poiss  uses the subroutine fill to calculate the coeff for the phi 
c     equation and uses the subroutine solve to resolve the phi equation.
c     in this version the current in the plasma has two components, the 
c     conduction current and the elctron diffusion current.
c-----------------------------------------------------------------------
      subroutine poiss(curr,volt,nrrvtest)
      implicit double precision (a-h,o-z)
      logical rvpol
      parameter(mr=200, mz=100, nsh=21)
      dimension cd0(0:mr)
      dimension bjh(0:mr,0:mz),zdc(mr,mz),zdf(mr,mz)
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/mag/
     r   bth(0:mr,0:mz),bjz(0:mr,0:mz),bjr(0:mr,0:mz) 
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/pol/
     r  rvpol
      common/currpoi/
     r  zi(0:mz),cathi
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr  

      data dpi /6.283185307/
c--------------------------------------------------------------------
      do 1 j=0,mr
         cd0(j) = 0.
   1  continue
c
      do 5 i=0,nzp
      do 6 j=0,nrp
         bth(j,i) = 0.
         bjz(j,i) = 0.
         bjr(j,i) = 0.
         bjh(j,i) = 0.
   6  continue
   5  continue
c      do 10 j=1,nr
c            cd0(j) = cdz(j,0)
c  10  continue
      do 775 i=0,nz
         cdez(0,i) = cdez(1,i)
c         cdez(0,i) = 0.0
         cder(0,i) = 0.0
         cder(nr,i) = 0.0
  775  continue
c---calculate coefficients  
      call fill(cd0,volt,nrrvtest)

c---relaxation
      do 15 i=1,nz
      do 16 j=1,nr
         zdc(j,i) =dc(j,i)
         zdf(j,i)=df(j,i)
         dc(j,i) = dc(j,i)/rlxphi
         df(j,i) = df(j,i) + (1.d0-rlxphi)*dc(j,i)*phi(j,i)
  16  continue
  15  continue

c include space charge sheath field at lower electrode - for current continuity, not Poisson's equation
c      do 84 i=1,nz
c      do 85 j=0,nr
c          temp = 2.*ae(j,i)/(dr(j) + dr(j+1))
c          if (j.ne.0) de(j,i) = temp
c          dw(j+1,i) = temp
c  85  continue
c  84  continue
c       do 20 i=0,nz
c       do 21 j=1,nr
c         if (i.ge.60.and.i.le.71) then
c          temp = 2.*an(j)/(dz(i)+dz(i+1))
c          if (i .ne. 0) dn(j,i) = temp
c          ds(j,i+1) = temp
c electric field from space charge effects
c        if (i.ne.0) df(j,i) = (e/epsilon)*(rnp(j,i)-rne(j,i)-rni(j,i))*vol(j,i)
c          if (i.ne.0)  df(j,i) = 0.0
c        endif
c  21   continue
c  20   continue

c---solve matrix
  17      call solve(1,nr,1,nz,de,dw,dn,ds,dc,df,phi,1,1)
c--reset boundries for current evaluation
      do 377 j=0,nr
         phi(j,0) = phi(j,1)
         phi(j,nzp) = 0.
  377 continue
      do 30 i=0,nz
         phi(nrp,i) = phi(nr,i)
         phi(0,i)=phi(1,i)
  30  continue
      voltcz = phi(1,1)
c--calculate residual
      resphi = 0.d0
      relphi = 0.d0
      do 40 i=1,nz
      do 41 j=1,nr
          dc(j,i)=zdc(j,i)
          df(j,i)=zdf(j,i)
          east= de(j,i)*phi(j+1,i)/dc(j,i)
          west= dw(j,i)*phi(j-1,i)/dc(j,i)
          rnorth= dn(j,i)*phi(j,i+1)/dc(j,i)
          south= ds(j,i)*phi(j,i-1)/dc(j,i)
          centre= df(j,i)/dc(j,i)
          res = east+west+south+rnorth+centre- phi(j,i)
          east = dabs(east)
          west = dabs(west)
          rnorth = dabs(rnorth)
          south = dabs(south)
          centre = dabs(centre)
          ph = dabs(phi(j,i))
          res = dabs(res)
          rel = res/dmax1(east,west,south,rnorth,centre,ph)

c---res in 0.1 %
        if( abs(phi(j,i)) .ge. 1e-6 ) then 
		res=res/phi(j,i)
        else
		res = res/1e-6
        endif

        ares = abs(res)*1000.d0
        if (ares .gt. resphi) then
           resphi = ares
           iph = i
           jph = j
        end if
        if (rel .gt.relphi) then
            relphi = rel
            iiph = i
            jjph = j
        endif
c----end of residu calculation
 41   continue
 40   continue

c calculate electric fields
      do 165 i=0,nz
      do 166 j=0,nrp
        er(j,i)=(phi(j,i)-phi(j+1,i))*2/(dr(j)+dr(j+1))
        ez(j,i)=(phi(j,i)-phi(j,i+1))*2/(dz(i)+dz(i+1))
 166  continue
 165  continue

c--- calculate the axial current density
      do 70 i=1,nzm
      do 72 j=1,nr
         cdz(j,i) = dn(j,i)*(phi(j,i)-phi(j,i+1))/an(j)
     +             +cdez(j,i)

  72  continue
  70  continue
      do 71 j=1,nr
            cdz(j,nz) = 2.*phi(j,nz)*sig(j,nzp)/dz(nz)
  71  continue
c----calculate the radial current density
      do 80 i=1,nz
      do 81 j=1,nr
         cdr(j,i) = de(j,i)*(phi(j,i)-phi(j+1,i))/ae(j,i)
     +             +cder(j,i)

  81  continue
  80  continue
      do 75 i=0,nz
         cdz(0,i) = cdz(1,i)
         cdr(0,i) = 0.0
  75  continue
      
c---- calculate current emerging from cathode surface
      cathi=0.d0
      anodi=0.d0
        do 90 i=1,nz
        do 91 j=1,nr
c---- anode
        if(ian(j,i).gt.ian(j,i+1))anodi=anodi
     +           +(cdz(j,i))*an(j)
        if(ian(j,i).gt.ian(j,i-1))anodi=anodi
     +           +(cdz(j,i-1))*an(j)
        if(ian(j,i).gt.ian(j+1,i))anodi=anodi
     +           +(cdr(j,i))*ae(j,i)
        if(ian(j,i).gt.ian(j-1,i))anodi=anodi
     +           +(cdr(j-1,i))*ae(j-1,i)
c---- cathode
        if(ica(j,i).gt.ica(j,i+1))cathi=cathi
     +           +(cdz(j,i))*an(j)
        if(ica(j,i).gt.ica(j,i-1))cathi=cathi
     +           +(cdz(j,i-1))*an(j)
        if(ica(j,i).gt.ica(j+1,i))cathi=cathi
     +           +(cdr(j,i))*ae(j,i)
        if(ica(j,i).gt.ica(j-1,i))cathi=cathi
     +           +(cdr(j-1,i))*ae(j-1,i)

 91     continue
 90     continue
        cathi=cathi*dpi
        anodi=anodi*dpi
c        write(*,102) cathi,anodi
c        write(*,225) curr,phi(1,1)
 102  format(1x,' icath  =',1pe10.3,' ianode =',e10.3)
 225  format(1x,'curr =',0p,f10.2,' voltage   =',f10.3)
        zjmax=0.0
      do 100 i=1,nz
         cop = 0.
         do 105 j=1,nr
            cop = cop + (cdz(j,i))*an(j)
c---calculate current at any cross section 
            zi(i)=cop*dpi
            rtest(j,i)=cop*dpi/curr
c---- calculate the magnetic field
            bth(j,i) = 0.12566*cop/(r(j)+0.5*dr(j))
 105     continue
         zi(i)=cop*dpi
         zjj= dabs(curr-cop*dpi)
         if (zjj .gt. zjmax) then
            zjmax=zjj
            ii = i
         end if   
 100  continue
c         write(*,103)  zjj, ii
 103  format(1x,'max j = ',e10.2,i3)
c---- calculate the magnetic force	
      do 110 i=1,nz
      do 111 j=1,nr

          cdzm=0.25*(cdz(j,i)+cdz(j+1,i)+cdz(j,i-1)+cdz(j+1,i-1))
         bjz(j,i) = 0.5*(bth(j,i)+bth(j,i-1))*cdzm   
 111  continue
 110  continue
      do 120 i=1,nz-1
      do 121 j=1,nr
          cdrm=0.25*(cdr(j,i)+cdr(j,i+1)+cdr(j-1,i)+cdr(j-1,i+1))
         bjr(j,i) = 0.5*(bth(j,i)+bth(j-1,i))*cdrm         
 121  continue
 120  continue
      if (nrrvtest.eq.1) curr = zi(5)
cc	
      return
      end
c------------------------------------------------------------------
c end of subroutine poiss
c-------------------------------------------------------------------

c-------------------------------------------------------------------
c     subroutine fill
c     v.1 ls 18.9.98
c
c     fill calculates coefficients needed in the calculation 
c     of the electric potential 
c     fill is used by the subroutine poiss.
c----------------------------------------------------------------------
      subroutine fill(cd0,volt,nrrvtest)
      implicit double precision (a-h,o-z)
      logical rvpol
      parameter(mr=200, mz=100, nsh=21)
      dimension cd0(0:mr)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/pol/
     r  rvpol 
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
c----------------------------------------------------------------------
      e = 1.60206d-19
      epsilon = 8.85d-14
      do 10 i=1,nz
      do 11 j=1,nr
	 df(j,i)=(e/epsilon)*(rnp(j,i)-rne(j,i)-rni(j,i))
c	 df(j,i)=0.0
         dn(j,i) = 2/(dz(i)*(dz(i)+dz(i)))
         ds(j,i+1) = dn(j,i)
	 de(j,i) = 2/(dr(j)*(dr(j)+dr(j)))
         dw(j,i+1) = de(j,i)
  11  continue
  10  continue
      do 40 i=1,nz
      do 41 j=1,nr
         dc(j,i) = de(j,i) + dw(j,i) + dn(j,i) + ds(j,i)
   41 continue
   40 continue

c     setting up boundary conditions
c     top and bottom boundaries
      do 30 j=1,nr
         ds(j,1) = 0.d0
	 dn(j,nz) = 0.d0
c top of lower electrode
c         if (j.le.nra) dn(j,nz) = 0.d0
 30   continue
c-Vertical boundary conditions:
c--dphi/dr =0. at r= 0.
      do 91 i=1,nz
         dc(1,i) = dc(1,i) - dw(1,i)
         dw(1,i) = 0.d0
   91 continue
c--dphi/dr =0. at r= radius.
      do 20 i=1,nz
         dc(nr,i) = dc(nr,i) - de(nr,i)
         de(nr,i) = 0.d0
 20   continue

c-Bottom boundary conditions:
      do 90 j=1,nr
            dc(j,nz) = 1.0
            dn(j,nz) = 0.0
            ds(j,nz) = 0.0
            dw(j,nz) = 0.0
            de(j,nz) = 0.0
            df(j,nz) = 0.0
c	   call zerocoef(dc,df,de,dw,dn,ds,0.0,j,nz)
   90 continue
c Top boundary condition
      do 92 j=1,nr
            dc(j,1) = 1.0
            dn(j,1) = 0.0
            ds(j,1) = 0.0
            dw(j,1) = 0.0
            de(j,1) = 0.0
            df(j,1) = volt
   92 continue
c-Top electrode - voltage = volt 
      do 94 i=1,nz
      do 93 j=1,nr
          if (ica(j,i).eq.1) then
            dc(j,i) = 1.0
            dn(j,i) = 0.0
            ds(j,i) = 0.0
            dw(j,i) = 0.0
            de(j,i) = 0.0
            df(j,i) = volt
          endif
   93 continue
   94 continue
c-Bottom electrode - voltage = zero 
      do 95 i=1,nz
      do 96 j=1,nr
          if (ian(j,i).eq.1) then
            dc(j,i) = 1.0
            dn(j,i) = 0.0
            ds(j,i) = 0.0
            dw(j,i) = 0.0
            de(j,i) = 0.0
            df(j,i) = 0.0
          endif
   96 continue
   95 continue
c--dphi/dz =0. at z= bottom for radius > anode
c      do 90 j=nra+1,nr
c         dc(j,nz) = dc(j,nz) - dn(j,nz)
c         dn(1,i) = 0.d0
c   90 continue
      return
      end
c end of subroutine fill
c-------------------------------------------------------------------

c     subroutine elec
c     v.2 ls 4.11.98
c
c     this subroutine solves the electron continuity equation. the electron
c     density is determined by the ionisation, the three body recombination 
c     and the ambipolar diffusion. Patankar method and "fillde" subroutine not used.
c     the convection can also be added by removing the comments in the
c     end of the subroutine fillde
c-----------------------------------------------------------------------
      subroutine elec
      implicit double precision (a-h,o-z)
      logical rvpol
      parameter(mr=200, mz=100, nsh=21)
      dimension zdc(mr,mz),zdf(mr,mz)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz),
     r  cdz(0:mr,0:mz), cdr(0:mr,0:mz),cdez(0:mr,0:mz),cder(0:mr,0:mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/pol/
     r  rvpol
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr  
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi 

      data dpi /6.283185307/
c--------------------------------------------------------------------
 496  format('  i        rjnez     d2dz2       rjner       d2ndr2    
     r rion       att          dne       ne        dnp         np') 
 498  format('  j        rnez     d2dz2       rner       d2ndr2    
     r rion       att          dne       ne      dnp     np') 
 497  format(' temp(3,3) ',e12.3)
 499     format(i4,10E12.3)
      electron  = 1.60206d-19
c      write(*,496) 
      boltz     = 1.38044d-23      
      emass     = 9.1094d-31
      pi        = 3.1415925
      rn        = 2.5E19
c---calculate patankar coefficients  
c      call fillde
c---relaxation
 1072  format('  deltel =',2e12.3)
c      do 15 i=1,nz
         de(j,1) = 0.0
c      do 16 j=1,nr
c         zdc(j,i) =dc(j,i)
c         zdf(j,i)=df(j,i)
c         dc(j,i) = dc(j,i)/rlxne
c         df(j,i) = df(j,i) + (1.-rlxne)*dc(j,i)*rne(j,i)
c  16  continue
c  15  continue
c set anode and cathode electron density

c cathode electron density boundary
c ne = 4*jrich/(e * veth)
         de(j,1) = 0.0
 1071  format(/'  gammai = ',E12.3,'  gamma =',E12.3) 

c anode electron density boundary
c         if(ian(j,i).eq.1)then
c        call zerocoef(dc,df,de,dw,dn,ds,1.d0,j,i)
c	endif
c  14  continue   
c  17  continue   
c nr boundary
         de(j,1) = 0.0
      do 18 i=0,nzp
         rne(0,i)= rne(1,i) 
 18   continue
c beginning insert jetair8d.f
c diffusion from D/mu = kT/e ~ 1 eV; mu ~ 1 for ions and 100 for electrons
       rmob = 2.0
c       diff =100.0*rmob
        diff = 0.0
c diffusion omitted to avoid complexity
c Radial and axial fluxes with backward differences
       do 193 i=1,nz
c         de(j,1) = 0.0
       do 194 j=1,nr
          rion = 0.0
          att = 0.0
          rece = 0.0
          reci = 0.0
          wir = -rmob*er(j,i)
          wire = 100*wir
          if (wire.ge.1000000.) wire = 1000000.
          if (wire.le.-1000000.) wire = -1000000.
          if (wir.ge.100000.) wir = 100000.
          if (wir.le.-100000.) wir = -100000.
          if(er(j,i).le.0.0) then
              rjnir(j,i) = wir*rni(j,i)
              rjner(j,i)=wire*rne(j,i)-diff*(rne(j+1,i)-rne(j,i))*2/
     r	          (dr(j+1)+dr(j))
              rjnpr(j,i) = -wir*rnp(j+1,i)
          else
              rjnir(j,i) = wir*rni(j+1,i)
              rjner(j,i)=wire*rne(j+1,i)-diff*(rne(j+1,i)-rne(j,i))*2/
     r	          (dr(j+1)+dr(j))
              rjnpr(j,i) = -wir*rnp(j,i)
          endif
          if (j.le.1) rjner(j,i) = 0.0
c Backward differences modified for reversed flow of current after current zero.
c        wiz = -rmob*ez(j,i)
        wiz = rmob*ez(j,i)
        wize = 100*wiz
        if (wize.ge.10000000.) wize = 10000000.
        if (wize.le.-10000000.) wize = -10000000.
c          if (wiz.ge.100000.) wiz = 100000.
c          if (wiz.le.-100000.) wiz = -100000.
        if(ez(j,i).le.0.0) then
c           rjniz(j,i) = wiz*rni(j,i+1)
           rjniz(j,i) = wiz*rni(j,i-1)
c           rjnez(j,i) = wize*rne(j,i+1)
           rjnez(j,i) = wize*rne(j,i-1)
c     r	          -diff*(rne(j,i+1)-rne(j,i))/dz(i)
           rjnpz(j,i) = -wiz*rnp(j,i)
        else
           rjniz(j,i) = wiz*rni(j,i)
           rjnez(j,i) = wize*rne(j,i)
c     r	          -diff*(rne(j,i+1)-rne(j,i))/dz(i)
c           rjnpz(j,i) = -wiz*rnp(j,i+1)
           rjnpz(j,i) = -wiz*rnp(j,i-1)
      endif
 194    continue
 193    continue

c Boundary conditions for fluxes at centre and top and bottom electrodes
c flow symmetric around central axial axis at r = 0 when j = 0
      do 212 i = 1,nz
          rjnez(0,i) = rjnez(1,i)
          rjner(0,i) = rjner(1,i)
          rjniz(0,i) = rjniz(1,i)
          rjnir(0,i) = rjnir(1,i)
          rjnpz(0,i) = rjnpz(1,i)
          rjnpr(0,i) = rjnpr(1,i) 
 212  continue
c top electrode flow out or in at i=0 the same
      do 213 j = 1,nr
          rjnez(j,0) = rjnez(j,1)
          rjner(j,0) = rjner(j,1)
          rjniz(j,0) = rjniz(j,1) 
          rjnir(j,0) = rjnir(j,1)
          rjnpz(j,0) = rjnpz(j,1)
          rjnpr(j,0) = rjnpr(j,1) 
c bottom boundary does not emit electrons, neg ions or pos ions
          rjnez(j,nzp) = rjnez(j,1)
          rjner(j,nzp) = rjner(j,1)
          rjniz(j,nzp) = rjniz(j,1) 
          rjnir(j,nzp) = rjnir(j,1)
          rjnpz(j,nzp) = rjnpz(j,1)
          rjnpr(j,nzp) = rjnpr(j,1) 
          rjnez(j,nzp) = rjnez(j,nz)
          rjniz(j,nzp) = rjniz(j,nz) 
 213  continue  
c bottom electrode does not emit electrons, neg ions or pos ion= after currentzero
      do 214 j = 1,nr
       if (time.gt.timezero) then
       if (i.eq.70.and.j.le.nra) then
       if (ez(j,i).le.0.0) then
            rjnez(j,i) = 0.0
            rjniz(j,i) = 0.0
	    rjnpz(j,i) = rjnpz(j,i-1)
       else
            rjnez(j,i) = rjnez(j,i-1)
            rjniz(j,i) = rjniz(j,i-1)
	    rjnpz(j,i) = 0.0
       endif
       endif
       endif
 214  continue  
c bottom electrode  has particle densities set to zero
      do 216 i = 71,nzp
      do 215 j = 0,nra
            rne(j,i) = 0.0
            rni(j,i) = 0.0
	    rnp(j,i) = 0.0
            rjnez(j,i-1) = 0.0
            rjniz(j,i-1) = 0.0
 215  continue  
 216  continue  
c Increments of particle densities
      do 210 i=1,nz         
      do 211 j=1,nr
	  reci = gammai*rni(j,i)*rnp(j,i)
	  rece = gamma*rne(j,i)*rnp(j,i)
         eav = eabs(j,i)
             eon = eav/rn
             weav = 100*rmob*eav
          if (weav.ge.100000000.) weav = 100000000.
          if (weav.le.-100000000.) weav = -100000000.
c Tam representation of alpha/n(E,ne) from Petrova; eont = E/N(Tam)
        rn1 = rn - rne(j,i) -rnp(j,i) - rni(j,i)
	      if (rn1.le.0.0) rn1 = 1.0
	      rn2 = rne(j,i)
	      if (rne(j,i).le.0.0) rne(j,i) = 1.0
	      eont = eon/1.0E4
        tam = 7.66893E-19/(eont + 20*eont*log((rn1+rn2)/rn1)**0.25)
        alphaon = 2.42013*1.0E-16*exp(-tam)
          etaon = 0.078E-17 - 0.0029*eon
	  if (etaon.le.1.6E-19) etaon = 1.6E-19
c Ionization and attachment equal at 25 kV/cm; 100 Td.
          eta = etaon*rn
c adjustment to make attachment rate constant near zero electric field
          if (eon.le.3.0E-17) att = rne(j,i)*3.75E6
          alpha = alphaon*rn
          ttt = 1400.0
	  if (t1(j,i).le.ttt) then
            att = rne(j,i)*eta*weav
          else
            att = 0.0
          endif
c alpha for argon from A. von Engel "Ionized Gaes"  Oxford 1965;  p 18
c p is pressure in mm Hg
c      aaa = 12.0
c      bbb = 180.0
c      p = 760.0
c      alpha = aaa*p*exp(-bbb/(eav/p))
c      if (alpha.ge.1.0E2) alpha = 100.0
      if (alpha.ge.1.0E2) alpha = 100.0
      rion = rne(j,i)*alpha*weav
c	  phoi = 4.0E7*rne(j,i) photoionization term
	  ddrf = 1+dr(j)/(2*(r(j)+dr(j)/2))
	  ddrb = 1-dr(j)/(2*(r(j)+dr(j)/2))
c          dni(j,i) =-((rjnir(j,i)*ddrf-rjnir(j-1,i)*ddrb)/dr(j) -
c     r        (rjniz(j,i)-rjniz(j,i-1)))*deltel/dz(i)+att-reci 
          rner = -(rjner(j,i)*ddrf-rjner(j-1,i)*ddrb)/dr(j)
          rnez = - (rjnez(j,i)-rjnez(j,i-1))/dz(i)
          dne(j,i) =-((rjner(j,i)*ddrf-rjner(j-1,i)*ddrb)/dr(j) -
     r        (rjnez(j,i)-rjnez(j,i-1)))*deltel/dz(i)-att
c     r        +SS(j,i)+rion-rece 
c Omit thermal ionization and recombination, which are ~equal and produce little net charge  
c - SS(j,i) is for thermal ionization
	d2dr2= -(rjner(j,i)*ddrf-rjner(j-1,i)*ddrb)/dr(j)
          d2dz2= -(rjnez(j,i)-rjnez(j,i+1))/dz(i)
          if (rne(j,i).le.1E-3) rne(j,i) = 1.0E-3
          if (rnp(j,i).le.1E-3) rnp(j,i) = 1.0E-3
          if (rne(j,i).ge.rmaxne) rmaxne = rne(j,i)
          dnp(j,i) =-((rjnpr(j,i)*ddrf-rjnpr(j-1,i)*ddrb)/dr(j) -
     r        (rjnpz(j,i)-rjnpz(j,i-1)))*deltel/dz(i)
c     r        +SS(j,i)+rion-rece-reci 
c Omit thermal ionization and recombination, which are ~equal and produce little net charge  
c        if (i.eq.66)
c     r      write(*,499) j,rjnez(j,i),d2dz2,rjner(j,i),d2dr2,rion,
c     r    att,dne(j,i),rne(j,i),dnp(j,i),rnp(j,i) 
c        if (j.eq.3)
c     r      write(*,499) i,rjnez(j,i),d2dz2,rjner(j,i),d2dr2,rion,
c     r    att,dne(j,i),rne(j,i),dnp(j,i),rnp(j,i) 
 211  continue
 210  continue

c increment particle densities
       do 198 i=1,nz
       do 197 j=1,nr
          rne(j,i)=rne(j,i)+dne(j,i)
	  if (rne(j,i).le.1E-3) rne(j,i) = 1.0E-3
c	  rni(j,i)=rni(j,i)+dni(j,i)
          rnp(j,i)=rnp(j,i)+dnp(j,i)
          if (rnp(j,i).le.1E-3) rnp(j,i) = 1.0E-3
	  rnet(j,i) = rnp(j,i) - rne(j,i)
c          rnet(j,i) = rnp(j,i) - rne(j,i) - rni(j,i)
 197  continue
 198  continue
c        jj = 3
c     r      write(*,499) i,rjnez(jj,i),d2dz2,rjner(jj,i),d2dr2,rion,
c     r    att,dne(jj,i),rne(jj,i),dnp(jj,i),rnp(jj,i) 
c        i = 66
c        do 301 j=0,nr
c           write(*,499) j,rjnez(j,i),d2dz2,rjner(j,i),d2dr2,rion,
c     r    att,dne(j,i),rne(j,i),dnp(j,i),rnp(j,i) 
c  301  continue
       do 199 i=1,nz
          rne(0,i) = rne(1,i)
          rni(0,i) = rni(1,i)
          rnp(0,i) = rnp(1,i)
 199  continue
c---solve matrix only used in Patankar method, assumes no net flow
c        call solvene(1,nr,1,nz,de,dw,dn,ds,dc,df,rne,0,1)

c--calculate residual
c      resne = 0.d0
c      relne = 0.d0
c      do 40 i=2,nzm
c      do 41 j=1,nrm
c         if(ica(j,i).eq.0.and.ian(j,i).eq.0)then
c          dc(j,i)=zdc(j,i)
c          df(j,i)=zdf(j,i)
c	east= de(j,i)*rne(j+1,i)/dc(j,i)
c          west= dw(j,i)*rne(j-1,i)/dc(j,i)
c          rnorth= dn(j,i)*rne(j,i+1)/dc(j,i)
c          south= ds(j,i)*rne(j,i-1)/dc(j,i)
c          centre= df(j,i)/dc(j,i)
c          res = east+west+south+rnorth+centre-rne(j,i)
c          east = dabs(east)
c          west = dabs(west)
c          rnorth = dabs(rnorth)
c          south = dabs(south)
c          centre = dabs(centre)
c          azne = dabs(rne(j,i))
c          res = dabs(res)
c          rel = res/dmax1(east,west,south,rnorth,centre,azne)
c        if (j.eq.68.and.i.eq.47) write(*,*) res
c
c        ares = abs(res)
c        if (ares .gt. resne) then
c           resne = ares
c           ine = i
c           jne = j
c        end if
c        if (j.eq.68.and.i.eq.47) write(*,*) resne
c        if (rel .gt.relne) then
c            relne = rel
c            iine = i
c            jjne = j
c        endif
c----end of residu calculation
c       endif
c 41   continue
c 40   continue

cls----prevent negative electron densities 
      do 20 i=0,nzp
      do 21 j=0,nrp
         if (rne(j,i) .lt. 1.d0) rne(j,i) = 1.0 
         if (rnp(j,i) .lt. 1.d0) rnp(j,i) = 1.0 
  21  continue
  20  continue

      return
      end
c------------------------------------------------------------------
c end of subroutine elec
c-------------------------------------------------------------------

c-------------------------------------------------------------------
c     subroutine fillde
c     v.1 ls 18.9.98
c
c     fill calculates Patankar coefficients needed in the calculation 
c     of the electron continuity equation; subroutine not needed with new elec
c----------------------------------------------------------------------
      subroutine fillde
      implicit double precision (a-h,o-z)
      logical rvpol,actest
      parameter(mr=200, mz=100, nsh=21)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/pol/
     r  rvpol 
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
c----------------------------------------------------------------------
      e = 1.60206d-19
      epsilon = 8.85d-14
      do 5 i=1,nz
      do 6 j=1,nr
          dc(j,i) = 0.d0
          df(j,i) = 0.d0
          dw(j,i) = 0.d0
          de(j,i) = 0.d0
          dn(j,i) = 0.d0
          ds(j,i) = 0.d0
   6   continue
   5   continue

      do 10 i=1,nz
      do 11 j=0,nr
          temp = 2.*ae(j,i)/(dr(j) + dr(j+1))
          if (j.ne.0) de(j,i) = temp
          dw(j+1,i) = temp
  11  continue
  10  continue
      do 20 i=0,nz
      do 21 j=1,nr
          temp = 2.*an(j)/(dz(i)+dz(i+1))
          if (i .ne. 0) dn(j,i) = temp
          ds(j,i+1) = temp
c electric field from space charge effects
c if (i.ne.0) df(j,i) = (e/epsilon)*(rnp(j,i)-rne(j,i)-rni(j,i))*vol(j,i)
c          if (i.ne.0)  df(j,i) = 0.0
  21   continue
  20   continue
      do 82 i=1,nz
      do 83 j=1,nr
          rnp(j,i) = rne(j,i)
          rni(j,i) = 0.0
  83  continue
  82  continue
c Boundary conditions
c   dne/dr = 0.0 at r = 0 and r = R 
      do 90 i=1,nz
         dw(1,i) = 0.d0
         de(nr,i) = 0.d0
  90   continue
c   dne/dz = 0.0 at z = 0 and i = nz 
      do 92 j=1,nr
         ds(j,1) = 0.d0
         dn(j,nz) = 0.d0
  92   continue
      do 80 i=1,nz
      do 81 j=1,nr
          dc(j,i) = de(j,i) + dw(j,i) +dn(j,i) +ds(j,i)
c Thermal ionization and recombination 
          df(j,i) = (SS(j,i)-gamma*rne(j,i)*rnp(j,i))*vol(j,i)
     r       -gamma*rni(j,i)*rnp(j,i)
  81  continue
  80  continue
c--------------------------------------------------------------------
c -- source terme: sc = ionisation, sp*ne = recombination         
c  rec = alpha * ne = conventional gamma , sc = alpha * (ne^2/n0)e * ne * n0
c  sp  = -alpha * ne^2
c  alpha from hinnov for t<3200 and from hoffert for t>3200
c         if (t(j,i).lt.3200.) then
c           rec = 1.1d-8*rne(j,i)/t(j,i)**4.5
c         else
c      rec=1.29d-32*(135300./t(j,i)+2.)*dexp(47800./t(j,i))*rne(j,i)
c         endif
c         zran=(ran(j,i)+2.*(rneth(j,i)-rne(j,i)))/ran(j,i)
c         sc = rec*(rneth(j,i)*ppp(i))**2
c         sc is source term S for ionization
c         if (sc .lt. 1.d-10) then
c            sc=0.0
c         endif
c         apot = abs(phi(1,1))
c         if (apot.le.voltion)  then
c             sc1 =  0.0
c         else
c             sc1 = sc
c         endif
c         sp = rec*rne(j,i)
c---------------------------------------------------------------------
c Thermal ionization 
c          df(j,i) = (SS(j,i)-gamma*rne(j,i)*rnp(j,i))*vol(j,i)
c     r       -gamma*rni(j,i)*rnp(j,i)
c Electron diffusion
c         dife = 2.*ae(j,i)*diffar(j,i)/(dr(j)+dr(j+1))
c         if (i.gt.0) de(j,i) = dife
c         if (i.gt.0) dw(j+1,i) = dife
c         difn = 2.*an(j)*diffaz(j,i)/(dz(i)+dz(i+1))
c         if (i.gt.0) dn(j,i) = difn
c         ds(j,i+1) = difn
c attachment for edge where temp<tempattach
c	 if (t(j,i).le.tempattach) then
c	      attach1 = attach 
c         else
c`	      attach1 = 0.0
c        endif
c         if (i.gt.0) df(j,i) = sc1*vol(j,i)
c     r    - attach1*rne(j,i)*ran(j,i)*ppp(i)*vol(j,i)
c
c attachment for all neutrals ran(j,i)
c         if (i.gt.0) df(j,i) = SS(j,i)*vol(j,i)
c     r    - attach*rne(j,i)*ran(j,i)*ppp(i)*vol(j,i)

c convection term for plasma (neutral flow); not included
c--- convection r-direction
c            temp1= ae(j,i)*vr(j,i)
c            de(j,i) = de(j,i) + dmax1(0.d0,-temp1)
c            dw(j+1,i) = dw(j+1,i) + dmax1(0.d0,temp1)
c--- convection z-direction
c	  temp2= an(j)*vz(j,i)
c	  dn(j,i)= dn(j,i)+ dmax1(0.d0,-temp2)
c	  ds(j,i+1)= ds(j,i+1) +dmax1(0.d0,temp2)

c convection term  for electron convection not plasma convection
c---electron convection r-direction
cc           temp1= ae(j,i)*zmob(j,i)*(phi(j+1,i)-phi(j,i))/dr(j)
cc           de(j,i) = de(j,i) + dmax1(0.d0,-temp1)
cc           dw(j+1,i) = dw(j+1,i) + dmax1(0.d0,temp1)
c-- electron convection z-direction
cc           temp2= an(j)*zmob(j,i)*(phi(j,i+1)-phi(j,i))/dz(i)
cc	   dn(j,i)= dn(j,i)+ dmax1(0.d0,-temp2)
cc           ds(j,i+1)= ds(j,i+1) +dmax1(0.d0,temp2)
c recombination no longer on left of equation for stability    
c gamma is effective two term recombination coefficient for both 
c       electron and negative ion recombination with positive ions
c           rec = gamma*(rne(j,i)+rni(j,i))*rnp(j,i)
c           df(j,i) = -vol(j,i)*rec
c      if (i.gt.0) dc(j,i)=dn(j,i)+ds(j,i)+de(j,i)+dw(j,i)+sp*vol(j,i)
c        if (actest) df(j,i) = df(j,i)+rne1(j,i)*vol(j,i)/delt
c         if (actest) dc(j,i) = dc(j,i)+vol(j,i)/delt
c 11   continue
c 10   continue

      return
      end
c end of subroutine fillde
c-------------------------------------------------------------------
      subroutine filldi
c for negative ion indices in electron balance
      implicit double precision (a-h,o-z)
      logical rvpol,actest
      parameter(mr=200, mz=100, nsh=21)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/pol/
     r  rvpol 
      common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
c----------------------------------------------------------------------
c     setting up boundary values
      do 20 i=1,nz
         dw(1,i) = 0.d0
c         de(nr,i) = 0.d0
 20   continue
      do 30 j=1,nr
c         ds(j,1) = 0.d0
c         dn(j,nz) = 0.d0
 30   continue

      do 10 i=1,nz
      do 11 j=1,nr

c Thermal ionization 
c          df(j,i) = SS(j,i) - gamma*rne(j,i)*rnp(j,i)-gamma*rni(j,i)*rnp(j,i)
         if (i.gt.0) df(j,i) = sc1*vol(j,i)

c convection term  for electron convection not plasma convection
c---electron convection r-direction
           temp1= ae(j,i)*zmob(j,i)*(phi(j+1,i)-phi(j,i))/dr(j)
           de(j,i) = de(j,i) + dmax1(0.d0,-temp1)
           dw(j+1,i) = dw(j+1,i) + dmax1(0.d0,temp1)
c-- electron convection z-direction
           temp2= an(j)*zmob(j,i)*(phi(j,i+1)-phi(j,i))/dz(i)
	   dn(j,i)= dn(j,i)+ dmax1(0.d0,-temp2)
           ds(j,i+1)= ds(j,i+1) +dmax1(0.d0,temp2)
c gamma is effective two term recombination coefficient for both 
c       electron and negative ion recombination with positive ions
           rec = gamma*(rne(j,i)+rni(j,i))*rnp(j,i)
           df(j,i) = -vol(j,i)*rec
c      if (i.gt.0) dc(j,i)=dn(j,i)+ds(j,i)+de(j,i)+dw(j,i)+sp*vol(j,i)
c        if (actest) df(j,i) = df(j,i)+rne1(j,i)*vol(j,i)/delt
c         if (actest) dc(j,i) = dc(j,i)+vol(j,i)/delt
 11   continue
 10   continue

      return
      end
c end of subroutine filldi
c----------------------------------------------------------------------
      subroutine filldp
c for positive ion indices in electron balance
      implicit double precision (a-h,o-z)
      logical rvpol,actest
      parameter(mr=200, mz=100, nsh=21)
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz)
      common/elc/
     r  sig(0:mr,0:mz),rad(0:mr,0:mz),dgamdt(0:mr),
     r  sigr(0:mr,0:mz),sigz(0:mr,0:mz),zmob(0:mr,0:mz)
      common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
      common/diff/
     r  zet(0:mr,0:mz), zemr(0:mr,0:mz), zemz(0:mr,0:mz),
     r  zemr1(0:mr,0:mz), zemz1(0:mr,0:mz), schm(0:mr,0:mz),
     r  cp(0:mr,0:mz), cpz(0:mr,0:mz),  cpr(0:mr,0:mz),    
     r  diffa(0:mr,0:mz),  diffar(0:mr,0:mz),diffaz(0:mr,0:mz)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
c----------------------------------------------------------------------
c     setting up boundary values
      do 20 i=1,nz
         dw(1,i) = 0.d0
c         de(nr,i) = 0.d0
 20   continue
      do 30 j=1,nr
c         ds(j,1) = 0.d0
c         dn(j,nz) = 0.d0
 30   continue

      do 10 i=1,nz
      do 11 j=1,nr

c---------------------------------------------------------------------
c Thermal ionization and recombination 
c          df(j,i) = SS(j,i) - gamma*rne(j,i)*rnp(j,i)-gamma*rni(j,i)*rnp(j,i)
c Electron diffusion
         dife = 2.*ae(j,i)*diffar(j,i)/(dr(j)+dr(j+1))
         if (i.gt.0) de(j,i) = dife
         if (i.gt.0) dw(j+1,i) = dife
         difn = 2.*an(j)*diffaz(j,i)/(dz(i)+dz(i+1))
         if (i.gt.0) dn(j,i) = difn
         ds(j,i+1) = difn

c convection term for pos ion convection not plasma convection
c---pos ion convection r-direction
           temp1= ae(j,i)*zmob(j,i)*(phi(j+1,i)-phi(j,i))/dr(j)
           de(j,i) = de(j,i) + dmax1(0.d0,-temp1)
           dw(j+1,i) = dw(j+1,i) + dmax1(0.d0,temp1)
c-- pos ion convection z-direction
           temp2= an(j)*zmob(j,i)*(phi(j,i+1)-phi(j,i))/dz(i)
	   dn(j,i)= dn(j,i)+ dmax1(0.d0,-temp2)
           ds(j,i+1)= ds(j,i+1) +dmax1(0.d0,temp2)
c recombination no longer on left of equation for stability    
c gamma is effective two term recombination coefficient for both 
c       electron and negative ion recombination with positive ions
           rec = gamma*(rne(j,i)+rni(j,i))*rnp(j,i)
           df(j,i) = -vol(j,i)*rec
c      if (i.gt.0) dc(j,i)=dn(j,i)+ds(j,i)+de(j,i)+dw(j,i)+sp*vol(j,i)
c        if (actest) df(j,i) = df(j,i)+rne1(j,i)*vol(j,i)/delt
c         if (actest) dc(j,i) = dc(j,i)+vol(j,i)/delt
 11   continue
 10   continue

      return
      end
c end of subroutine filldp
 
c-----------------------------------------------------------------
c     subroutine solvene
c     v.1 ls 18.9.98
c
c     solvene is the solver for the electron density equation
c     solvene uses the line by line technique in conjonction with 
c     the block correction technique.
c-----------------------------------------------------------------

      subroutine solvene(m,nr,n,nz,de,dw,dn,ds,dc,df,phi,flag1,flag2)
      implicit double precision (a-h,o-z)
      logical actest
      parameter(mr=200, mz=100, nsh=21)
      dimension de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1   df(mr,mz),phi(0:mr,0:mz),v(mz),x(0:mr,0:mz)
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr  
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz) 
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest

      integer flag1,flag2
c------------------------------------------------------------------  
      nr1 = nr-1
      nz1 = nz-1
      do 5 i=0,mz
      do 6 j=0,mr
        x(j,i) = 0.d0
   6  continue
   5  continue
      do 10 i=0,nz+1
      do 11 j=0,nr+1
         x(j,i) = phi(j,i)
  11  continue
  10  continue
      rel = 1.85
      if (flag1 .eq. 0) then
c         gams = 0.00625d0
         gams = 0.00125d0
         res0 = 0.d0
         do 20 i=n,nz
         do 21 j=m,nr
c            if (x(j,i) .eq. 0.d0 .and. dc(j,i) .gt. 1.d6) go to 20
            res0 = res0 + (de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     1             dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) + df(j,i) -
     2             dc(j,i)*x(j,i))**2
  21     continue
  20     continue
      else
         gams = 0.d0
         iconv = 0
      end if
      do 200 iter=1,50
c
c     first step: block correction along lines of constant j
c

      do 100 j=m,nr1
         a(j) = 0.d0
         b(j) = 0.d0
         c(j) = 0.d0
         d(j) = 0.d0
 100   continue
      do 110 i=n,nz
      do 111 j=m,nr1
         a(j) = a(j) + dc(j,i) - dn(j,i) - ds(j,i)
         b(j) = b(j) + de(j,i)
         c(j) = c(j) + dw(j,i)
         d(j) = d(j) + dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) +
     1       de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     2       df(j,i) - dc(j,i)*x(j,i)
 111  continue
 110  continue
      call tdma(m,nr1,v)
      do 120 i=n,nz
      do 121 j=m,nr1
         x(j,i) = x(j,i) + v(j)
 121  continue
 120  continue
c
c     second step: block corrections along lines of constant i
c
      do 130 i=n,nz1
         a(i) = 0.d0
         b(i) = 0.d0
         c(i) = 0.d0
         d(i) = 0.d0
 130  continue
      do 140 i=n,nz1
      do 141 j=m,nr
         a(i) = a(i) + dc(j,i) - de(j,i) - dw(j,i)
         b(i) = b(i) + dn(j,i)
         c(i) = c(i) + ds(j,i)
         d(i) = d(i) + dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) +
     1       de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     2       df(j,i) - dc(j,i)*x(j,i)
 141  continue
 140  continue
      call tdma(n,nz1,v)
      do 150 i=n,nz1
      do 151 j=m,nr
         x(j,i) = x(j,i) + v(i)
 151  continue
 150  continue

c     third step: sweep from i=1 to nz
c
      do 30 i=n,nz
         do 35 j=m,nr
            a(j) = dc(j,i) - (rel-1.)*dn(j,i)
            b(j) = de(j,i)
            c(j) = dw(j,i)
            d(j) = df(j,i) + dn(j,i)*(x(j,i+1) - (rel-1.)*x(j,i))
     1              + ds(j,i)*x(j,i-1)
  35     continue
         call tdma(m,nr,v)
         do 40 j=m,nr
            x(j,i) = v(j)
  40     continue
  30  continue
c
c     fourth step: sweep from j=nr to 1
c
      do 70 j=nr,m,-1
         do 75 i=n,nz
            a(i) = dc(j,i) - (rel-1.)*dw(j,i)
            b(i) = dn(j,i)
            c(i) = ds(j,i)
            d(i) = df(j,i) + de(j,i)*x(j+1,i)
     1             + dw(j,i)*(x(j-1,i) - (rel-1.)*x(j,i))
  75     continue
         call tdma(n,nz,v)
         do 80 i=n,nz
            x(j,i) = v(i)
  80     continue
  70  continue
c
c     fifth step: sweep from i=nz to 1
c
      do 50 i=nz,n,-1
         do 55 j=m,nr
            a(j) = dc(j,i) - (rel-1.)*ds(j,i)
            b(j) = de(j,i)
            c(j) = dw(j,i)
            d(j) = df(j,i) + dn(j,i)*x(j,i+1) + ds(j,i)*(x(j,i-1) -
     1             (rel-1.)*x(j,i))
  55     continue
         call tdma(m,nr,v)
         do 60 j=m,nr
            x(j,i) = v(j)
  60     continue
  50  continue
c
c     sixth step: sweep from j=1 to nr
c
      do 90 j=m,nr
         do 95 i=n,nz
            a(i) = dc(j,i) - (rel-1.)*de(j,i)
            b(i) = dn(j,i)
            c(i) = ds(j,i)
            d(i) = df(j,i) + dw(j,i)*x(j-1,i) + de(j,i)*(x(j+1,i) -
     1             (rel-1.)*x(j,i))
  95     continue
         call tdma(n,nz,v)
         do 96 i=n,nz
            x(j,i) = v(i)
  96     continue
  90  continue
c
c     check convergence
c
      if (gams .ne. 0.d0) then
         res = 0.d0
         do 160 i=n,nz
         do 161 j=m,nr
            res = res + (de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     1             dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) + df(j,i) -
     2             dc(j,i)*x(j,i))**2
 161     continue
 160     continue
c         if (jscreen.eq.1) write(*,999) res,res0
cl	 if (mc.ge.mm1) write(5,999) res,res0
c  999     format(1x,'res =',1p,e16.8,'   res0 =',e16.8)
         if (res .le. gams*res0.and.iter.gt.2) then
            do 170 i=n,nz
            do 171 j=m,nr
               phi(j,i) = x(j,i)            
 171        continue
 170        continue
            return
         end if
      else
         do 260 i=n,nz
         do 261 j=m,nr
            if (x(j,i) .ne. 0.) then
               diff = (x(j,i)-phi(j,i))/x(j,i)
            else
               diff = -phi(j,i)
            end if
	    if (dabs(diff) .gt. 1.e-3) go to 265
  261    continue
  260    continue
         iconv = 1
  265    continue
c       if (jscreen.eq.1) write(*,998) iter,diff,i,j,x(j,i),phi(j,i)
  998    format(1x,'iter =',i4,'    diff =',1p,e16.8,
     1            '  at   (',i2,',',i2,')',2e12.3)
         do 270 i=n,nz
         do 271 j=m,nr
            phi(j,i) = x(j,i)
  271    continue
  270    continue
         if (iconv .eq. 1) return
      end if
 200  continue
      write(*,300)
c      write(5,300)
 300  format(1x,'convergence not attained in subroutine solvene')
      return
      end

c------------------------------------------------------------------
c end of subroutine solvene
c-------------------------------------------------------------------
 
c-----------------------------------------------------------------
c     subroutine solve
c     v.1 ls 18.9.98
c
c     solve is the solver for the poisson equation
c     solve uses the line by line technique in conjonction with 
c     the block correction technique.
c-----------------------------------------------------------------

      subroutine solve(m,nr,n,nz,de,dw,dn,ds,dc,df,phi,flag1,flag2)
      implicit double precision (a-h,o-z)
      logical actest
      parameter(mr=200, mz=100, nsh=21)
      dimension de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1   df(mr,mz),phi(0:mr,0:mz),v(mz),x(0:mr,0:mz)
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr  
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz) 
      common/ac/
     r   time,delt,t1(0:mr,0:mz),vr1(0:mr,0:mz),h1(0:mr,0:mz),
     1   vz1(0:mr,0:mz),rho1(0:mr,0:mz),rhmr1(0:mr,0:mz),deltel,
     2   rhmz1(0:mr,0:mz),rne1(0:mr,0:mz),rrrv,actest

      integer flag1,flag2
c------------------------------------------------------------------  
      nr1 = nr-1
      nz1 = nz-1
      do 5 i=0,mz
      do 6 j=0,mr
        x(j,i) = 0.d0
   6  continue
   5  continue
      do 10 i=0,nz+1
      do 11 j=0,nr+1
         x(j,i) = phi(j,i)
  11  continue
  10  continue
      rel = 1.85
      if (flag1 .eq. 0) then
         gams = 0.00625d0
         res0 = 0.d0
         do 20 i=n,nz
         do 21 j=m,nr
            if (x(j,i) .eq. 0.d0 .and. dc(j,i) .gt. 1.d6) go to 20
            res0 = res0 + (de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     1             dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) + df(j,i) -
     2             dc(j,i)*x(j,i))**2
  21     continue
  20     continue
      else
         gams = 0.d0
         iconv = 0
      end if
      do 200 iter=1,50
c
c     first step: block correction along lines of constant j
c
      do 100 j=m,nr1
         a(j) = 0.d0
         b(j) = 0.d0
         c(j) = 0.d0
	 d(j) = 0.d0
 100  continue
      do 110 i=n,nz
      do 111 j=m,nr1
         a(j) = a(j) + dc(j,i) - dn(j,i) - ds(j,i)
         b(j) = b(j) + de(j,i)
         c(j) = c(j) + dw(j,i)
         d(j) = d(j) + dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) +
     1       de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     2       df(j,i) - dc(j,i)*x(j,i)
 111  continue
 110  continue
      call tdma(m,nr1,v)
      do 120 i=n,nz
      do 121 j=m,nr1
         x(j,i) = x(j,i) + v(j)
 121  continue
 120  continue
c
c     second step: block corrections along lines of constant i
c
      do 130 i=n,nz1
         a(i) = 0.d0
         b(i) = 0.d0
         c(i) = 0.d0
         d(i) = 0.d0
 130  continue
      do 140 i=n,nz1
      do 141 j=m,nr
         a(i) = a(i) + dc(j,i) - de(j,i) - dw(j,i)
         b(i) = b(i) + dn(j,i)
         c(i) = c(i) + ds(j,i)
         d(i) = d(i) + dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) +
     1       de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     2       df(j,i) - dc(j,i)*x(j,i)
 141  continue
 140  continue
      call tdma(n,nz1,v)
      do 150 i=n,nz1
      do 151 j=m,nr
         x(j,i) = x(j,i) + v(i)
 151  continue
 150  continue
c
c     third step: sweep from i=1 to nz
c
      do 30 i=n,nz
         do 35 j=m,nr
            a(j) = dc(j,i) - (rel-1.)*dn(j,i)
            b(j) = de(j,i)
            c(j) = dw(j,i)
            d(j) = df(j,i) + dn(j,i)*(x(j,i+1) - (rel-1.)*x(j,i))
     1              + ds(j,i)*x(j,i-1)
  35     continue
         call tdma(m,nr,v)
         do 40 j=m,nr
            x(j,i) = v(j)
  40     continue
  30  continue
c
c     fourth step: sweep from j=nr to 1
c
      do 70 j=nr,m,-1
         do 75 i=n,nz
            a(i) = dc(j,i) - (rel-1.)*dw(j,i)
            b(i) = dn(j,i)
            c(i) = ds(j,i)
            d(i) = df(j,i) + de(j,i)*x(j+1,i)
     1             + dw(j,i)*(x(j-1,i) - (rel-1.)*x(j,i))
  75     continue
         call tdma(n,nz,v)
         do 80 i=n,nz
            x(j,i) = v(i)
  80     continue
  70  continue
c
c     fifth step: sweep from i=nz to 1
c
      do 50 i=nz,n,-1
         do 55 j=m,nr
            a(j) = dc(j,i) - (rel-1.)*ds(j,i)
            b(j) = de(j,i)
            c(j) = dw(j,i)
            d(j) = df(j,i) + dn(j,i)*x(j,i+1) + ds(j,i)*(x(j,i-1) -
     1             (rel-1.)*x(j,i))
  55     continue
         call tdma(m,nr,v)
         do 60 j=m,nr
            x(j,i) = v(j)
  60     continue
  50  continue
c
c     sixth step: sweep from j=1 to nr
c
      do 90 j=m,nr
         do 95 i=n,nz
            a(i) = dc(j,i) - (rel-1.)*de(j,i)
            b(i) = dn(j,i)
            c(i) = ds(j,i)
            d(i) = df(j,i) + dw(j,i)*x(j-1,i) + de(j,i)*(x(j+1,i) -
     1             (rel-1.)*x(j,i))
  95     continue
         call tdma(n,nz,v)
         do 96 i=n,nz
            x(j,i) = v(i)
  96     continue
  90  continue
c
c     check convergence
c
      if (gams .ne. 0.d0) then
         res = 0.d0
         do 160 i=n,nz
         do 161 j=m,nr
            res = res + (de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     1             dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) + df(j,i) -
     2             dc(j,i)*x(j,i))**2
 161     continue
 160     continue
c         if (jscreen.eq.1) write(*,999) res,res0
c	 if (mc.ge.mm1) write(5,999) res,res0
 999     format(1x,'res =',1p,e16.8,'   res0 =',e16.8)
c         if (res.le.gams*res0.or.res.lt.1.d-10) then
         if (res.le.gams*res0) then
            do 170 i=n,nz
            do 171 j=m,nr
               phi(j,i) = x(j,i)            
 171        continue
 170        continue
            return
         end if
      else
         do 260 i=n,nz
         do 261 j=m,nr
            if (x(j,i) .ne. 0.) then
               diff = (x(j,i)-phi(j,i))/x(j,i)
            else
               diff = -phi(j,i)
            end if
            if (dabs(diff) .gt. 1.e-3) go to 265
  261    continue
  260    continue
         iconv = 1
  265    continue
c       if (jscreen.eq.1) write(*,998) iter,diff,i,j,x(j,i),phi(j,i)
c	 if (mc.ge.mm1) write(5,998) iter,diff,i,j
  998    format(1x,'iter =',i4,'    diff =',1p,e16.8,
     1            '  at   (',i2,',',i2,')',2e12.3)
         do 270 i=n,nz
         do 271 j=m,nr
            phi(j,i) = x(j,i)
  271    continue
  270    continue
         if (iconv .eq. 1) return
      end if
 200  continue
      write(*,300)
c      write(5,300)
 300  format(1x,'convergence not attained in subroutine solve')
      return
      end

c------------------------------------------------------------------
c end of subroutine solve
c-------------------------------------------------------------------
 
c-------------------------------------------------------------------
c     subroutine solvpp
c     v.1 ls 18.9.98
c
c     solve is the solver for the pressure.
c     h p schmidt version
c     df cast to left if df/x < 0.0d0 for stability
c     ---------------------------------------------------------------

      subroutine solvpp(m,nr,n,nz,de,dw,dn,ds,dc,df,phi,flag1,flag2)
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      dimension de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz),phi(0:mr,0:mz),x(0:mr,0:mz),dcc(mr,mz),
     2  dfold(mr,mz)
      common/dummy/
     r  rlxp,rlxvz,rlxvr,rlxt,rlxne,rlxphi,resvz,relne,relvz,  
     r  resvr, relvr,resp,relp,rese,rele,resphi,relphi,resne,
     r  ivz,jvz, iivz,jjvz,iiph,jjph,ie,je,iie,jje,iph,jph,     
     r  ine,jne, iine,jjne,nnel, 
     r  mc,mm1,ifread,ip,jp,iip,jjp,ivr,jvr,iivr,jjvr  
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz)     

      integer flag1,flag2
c-----------------------------------------------------------------------

      gams = 0.0625d0 
      gamsi=gams*0.2d0
      iter2=1
      iterma=1
      iall= 20
      icountm=50
      
      do 5 i=0,mz
      do 6 j=0,mr
        x(j,i) = 0.0d0 
   6  continue  
   5  continue  

       
      do 7 i=1,mz
      do 8 j=1,mr          
         dfold(j,i)=df(j,i)
         dcc(j,i) = dc(j,i)
   8  continue
   7  continue
      
      do 10 i=0,nz+1
      do 11 j=0,nr+1
         x(j,i) = phi(j,i)
  11  continue
  10  continue
      
c     solve one time and check residiuum 
      call residu(m,nr,n,nz,de,dw,dn,ds,dcc,dfold,x,res0)
      call solvie(m,nr,n,nz,de,dw,dn,ds,dcc,dfold,x,iterma)
      call residu(m,nr,n,nz,de,dw,dn,ds,dcc,dfold,x,rest)
c      write(*,888) res,resold
      


      do 200 iter=1,iall
c          stop iteration if original matrix solved and converged----
          if (rest .le. gams*res0)then
                  goto 999
          endif
         
c         update coefficents; cast df to left for stability
           do 201 i=n,nz
           do 203 j=m,nr
               if ( dabs(x(j,i)) .gt. 1d-78) then
                 
                 dfdumy=dfold(j,i)/x(j,i)
                 if(dfdumy . lt. 0d0 ) then
                    dc(j,i)=dcc(j,i)-dfdumy
                    df(j,i)=0.0d0 
                 else
                   dc(j,i)=dcc(j,i)
                   df(j,i)=dfold(j,i)
                 endif
                 
               else
                 dc(j,i)=dcc(j,i)
                 df(j,i)=dfold(j,i)
               endif 
               
  203       continue  
  201       continue  
           
            
            call residu(m,nr,n,nz,de,dw,dn,ds,dc,df,x,resi0)
            resi=resi0
c           inner loop for one df(j,i)/x(j,i) to the left------------ 
            do 202 icount=1,icountm 
               if (resi.lt. gams*resi0) goto 202
c              stop if not converging- -------------------------------         
               if (resi.gt. resi0)then 
                 write(*,*) ' not converged for pp icount=',icount
                 goto 999 
               endif
               call solvie(m,nr,n,nz,de,dw,dn,ds,dc,df,x,iter2)
               call residu(m,nr,n,nz,de,dw,dn,ds,dc,df,x,resi)

c               write(*,887) resi,resi0                 

  202       continue                                       
c           calculate actual residuum -------------------------------
            call residu(m,nr,n,nz,de,dw,dn,ds,dcc,dfold,x,rest)
  200 continue               
c       
c
      

c 999       write(*,888) rest,res0
 999       continue
            do 170 i=n,nz
            do 171 j=m,nr 
               phi(j,i) = x(j,i)
               dc(j,i) = dcc(j,i)
               df(j,i) = dfold(j,i)
 171        continue
 170        continue
      return
      
 887  format(1x,'resi =',1p,e16.8,'   resi0=',e16.8)
 888  format(1x,'res =',1p,e16.8,'   res0 =',e16.8)      
      
      end

c----------------------------------------------------------------------
c end of subroutine solvpp
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c     subroutine solvie
c     v.1 ls 18.9.98
c
c     solvie is used in the subroutine solvpp
c     solvie solve a x= b matrix  
c     solvie does only one iteration  i.e solve from fbarc code
c
c     this subroutine uses the line by line technique in conjunction
c     with the block correction technique
c-----------------------------------------------------------------------

      subroutine solvie(m,nr,n,nz,de,dw,dn,ds,dc,df,x,imax)
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      dimension de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz),v(mz),x(0:mr,0:mz)
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz) 
c------------------------------------------------------------------------     

      nr1 = nr-1
      nz1 = nz-1
      rel = 1.85 
      gams = 0.1d0 
      
      do 200 iter=1,imax
c
c     first step: block correction along lines of constant j
c
      do 100 j=m,nr1
         a(j) = 0.d0
         b(j) = 0.d0
         c(j) = 0.d0
         d(j) = 0.d0
100   continue
      do 110 i=n,nz
      do 111 j=m,nr1
         a(j) = a(j) + dc(j,i) - dn(j,i) - ds(j,i)
         b(j) = b(j) + de(j,i)
         c(j) = c(j) + dw(j,i)
         d(j) = d(j) + dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) +
     1       de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     2       df(j,i) - dc(j,i)*x(j,i)
 111  continue
 110  continue
      call tdma(m,nr1,v)
      do 120 i=n,nz
      do 121 j=m,nr1
         x(j,i) = x(j,i) + v(j)
 121  continue
 120  continue
c
c     second step: block corrections along lines of constant i
c
      do 130 i=n,nz1
         a(i) = 0.d0
         b(i) = 0.d0
         c(i) = 0.d0
         d(i) = 0.d0
 130  continue
      do 140 i=n,nz1
      do 141 j=m,nr
         a(i) = a(i) + dc(j,i) - de(j,i) - dw(j,i)
         b(i) = b(i) + dn(j,i)
         c(i) = c(i) + ds(j,i)
         d(i) = d(i) + dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) +
     1       de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     2       df(j,i) - dc(j,i)*x(j,i)
 141  continue
 140  continue
      call tdma(n,nz1,v)
      do 150 i=n,nz1
      do 151 j=m,nr
         x(j,i) = x(j,i) + v(i)
 151  continue
 150  continue
c
c     third step: sweep from i=1 to nz
c
  25  do 30 i=n,nz
         do 35 j=m,nr
            a(j) = dc(j,i) - (rel-1.)*dn(j,i)
            b(j) = de(j,i)
            c(j) = dw(j,i)
            d(j) = df(j,i) + dn(j,i)*(x(j,i+1) - (rel-1.)*x(j,i))
     1              + ds(j,i)*x(j,i-1)
  35     continue
         call tdma(m,nr,v)
         do 40 j=m,nr
            x(j,i) = v(j)
  40     continue
  30  continue
c
c     fourth step: sweep from j=nr to 1
c
      do 70 j=nr,m,-1
         do 75 i=n,nz
            a(i) = dc(j,i) - (rel-1.)*dw(j,i)
            b(i) = dn(j,i)
            c(i) = ds(j,i)
            d(i) = df(j,i) + de(j,i)*x(j+1,i)
     1             + dw(j,i)*(x(j-1,i) - (rel-1.)*x(j,i))
  75     continue
         call tdma(n,nz,v)
         do 80 i=n,nz
            x(j,i) = v(i)
  80     continue
  70  continue
c
c     fifth step: sweep from i=nz to 1
c
      do 50 i=nz,n,-1
         do 55 j=m,nr
            a(j) = dc(j,i) - (rel-1.)*ds(j,i)
            b(j) = de(j,i)
            c(j) = dw(j,i)
            d(j) = df(j,i) + dn(j,i)*x(j,i+1) + ds(j,i)*(x(j,i-1) -
     1             (rel-1.)*x(j,i))
  55     continue
         call tdma(m,nr,v)
         do 60 j=m,nr
            x(j,i) = v(j)
  60     continue
  50  continue
c
c     sixth step: sweep from j=1 to nr
c
      do 90 j=m,nr
         do 95 i=n,nz
            a(i) = dc(j,i) - (rel-1.)*de(j,i)
            b(i) = dn(j,i)
            c(i) = ds(j,i)
            d(i) = df(j,i) + dw(j,i)*x(j-1,i) + de(j,i)*(x(j+1,i) -
     1             (rel-1.)*x(j,i))
  95     continue
         call tdma(n,nz,v)
         do 96 i=n,nz
            x(j,i) = v(i)
  96     continue
  90  continue
c

c    
 200  continue
      end

c-----------------------------------------------------------------------
c end of subroutine solvie
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     subroutine residu
c     v.1 ls 18.9.98
c
c    residu calcultes the residuum for x
c    residu is used in subroutine solvpp 
c-----------------------------------------------------------------------  
      
      subroutine residu(m,nr,n,nz,de,dw,dn,ds,dc,df,x,resid)
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      dimension de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1  df(mr,mz),x(0:mr,0:mz)
c------------------------------------------------------------------------
      resid=0.0d0
      do 10 i=n,nz
      do 11 j=m,nr
         if (x(j,i) .eq. 0.d0 .and. dc(j,i) .gt. 1.d+06) goto 10
         resid = resid + (de(j,i)*x(j+1,i) + dw(j,i)*x(j-1,i) +
     1          dn(j,i)*x(j,i+1) + ds(j,i)*x(j,i-1) + df(j,i) -
     2          dc(j,i)*x(j,i))**2
  11  continue              
  10  continue              
c              
      return
      end

c----------------------------------------------------------------------
c end of subroutine residu
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     subroutine tdma
c     v.1 ls 18.9.98
c
c     tdma solves equations of the form
c     a(i)*v(i) = b(i)*v(i+1) + c(i)*v(i-1) + d(i)  ,i=1,....,n
c     with c(1) = 0 and b(n) = 0
c-------------------------------------------------------------------------

      subroutine tdma(if,l,v)
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      dimension v(mz),p(mz),q(mz)
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz)  
c-------------------------------------------------------------------------

      l1 = l-1
      p(if) = b(if)/a(if)
      q(if) = d(if)/a(if)
      ifp1 = if + 1
      do 10 i=ifp1,l
         den = 1./(a(i) - c(i)*p(i-1))
         p(i) = b(i)*den
         q(i) = (d(i)+c(i)*q(i-1))*den
  10  continue
      v(l) = q(l)
      do 20 i=if,l1
         k = l1 -i + if
         v(k) = p(k)*v(k+1) + q(k)
  20  continue
      return
      end

c-----------------------------------------------------------
c end subroutine tdma
c-----------------------------------------------------------
c-----------------------------------------------------------------------
c       subroutine meshls
c       v.2 ls 16.04.99
c
c meshls calculates the mesh from the data of the file gridls.
c 
c     initial setting up of coordinate system.
c     there are nz axial grid points and nr radial grid points.
c     there are nz control volumes in the z direction and
c     nr control volumes in the r direction.
c     the j-th radial control volume contains the j-th
c        radial grid point and has a width of dr(j).
c     the i-th axial control volume contains the i-th
c        axial grid point and has a length of dz(i).
c     the first and the last control volumes in the z and r
c        directions have zero width or zero length.
c     r(0), r(nr), z(0), and z(nz) determine the domain
c        of integration.
c
c-----------------------------------------------------------------------
	subroutine meshls

      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      dimension  zb(0:mz),rb(0:mr)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
c------------------------------------------------------------------------
c read the data for the mesh generation in file gridls

        open(18,file='gridls')
        rewind 18
	read(18,*)
	read(18,*) cathgeom
        read(18,*)
        read(18,*)
        read(18,*) cathode, diam, wtip, angtip
	read(18,*)
	read(18,*) nntip, nzc, nuni
        read(18,*)
        read(18,*)
	read(18,*) anode, adiam, nanode,nra
        read(18,*)
        read(18,*)
        read(18,*) gap, ngap
        read(18,*)
        read(18,*) shca, shan
        read(18,*)
        read(18,*) radius, nrout
        read(18,*)
	read(18,*) dnozz
        read(18,*)
        read(18,*)
        read(18,*) rcathode,rdiam, nzcr, nrcr
	close(18)
c--------------------------------------------------------------------
c       some initialisation
	do 10 j = 0,mr
	   rb(j)=0.0
 10	continue   

	do 20 i = 0,mz
	   zb(i)=0.0
 20	continue
c---------------------------------------------
	if(cathgeom.eq.0) then
c       conical cathode

        cradius = diam/2.
        aradius = adiam/2.
c       ang is the half cathode tip angle
	ang = 3.14159*angtip/360.

	nr = nuni+nntip+nrout
	nrp = nr+1
	nrm = nr-1
	nz = nzc+ngap+nanode
	nzp = nz+1
	nzm = nz-1
	nza = nzc + ngap +1
	nrspc = 0

c----------------------------------------------------------------
c axial mesh calculation
c---------------------------------------------------------------
c
c calculation of the axial cell boundaries
c
c----------- axial cathode mesh -------------------------------
	zb(nzc) = cathode
c--------------------------------------------------------------
c 1) variable mesh for the conical tip with
c    shca for the first cell thickness from the tip
c    and with nntip points

	ipro = 1
	nn = nntip
	zz = (cradius-wtip)/tan(ang)
	db = shca
	call multcoef(ipro,zz,db,nn,a)	
	do 41 i= 1,nn
	   dz(i)= db*(a**float(i-1))
	   zb(nzc-i)=zb(nzc+1-i)-dz(i)
 41	continue
c-------------------------------------------------------------
c 2) variable mesh for the rest of the cathode with
c    the same thickness for the first point as for the last cell
c    of the conical tip and nzc-nntip points

	ipro = 1
	nn = nzc-nntip
	zz = cathode-(cradius-wtip)/tan(ang)
	db = dz(nntip)
	call multcoef(ipro,zz,db,nn,a)
	do 42 i = 1,nn
	   dz(i)= db*(a**float(i-1))
	   zb(nzc-nntip-i)=zb(nzc-nntip+1-i)-dz(i)
 42	continue
c----- end of axial cathode mesh -----------------------

c----- axial plasma gap mesh --------------------------

c--------------------------------------------------------
c 1) first half of the gap from the cathode tip.
c    varible mesh with shca for the first cell's thickness
c    and ngap/2 points

	ipro = 1
	nn = ngap/2
	zz = gap/2
	db = shca
	call multcoef(ipro,zz,db,nn,a)
	do 43 i= 1, ngap/2
	   dz(i) = db*(a**float(i-1))
	   zb(nzc+i)=zb(nzc-1+i)+dz(i)
 43	continue
c------------------------------------------------------
c 2) second half of the gap until the anode.
c    variable mesh with shan for the cell's thickness in 
c    front of the anode.
	zb(nza-1) = cathode + gap
	ipro = 1
	nn = ngap-ngap/2
	zz = gap/2
	db = shan
	call multcoef(ipro,zz,db,nn,a)
	do 44 i= 1, nn-1
	   dz(i) = db*(a**float(i-1))
	   zb(nza-1-i) = zb(nza-i)-dz(i)   
 44	continue
c------- end of the axial plasma gap mesh -------------
c
c------  anode axial mesh ----------------------------
c variable mesh with shan for the first cell's thickness
c and naonde point.

	ipro=1
	nn = nanode
	zz = anode
	db = shan
	call multcoef(ipro,zz,db,nn,a)
	do 45 i=0,nn-1
	   dz(i) = db*(a**float(i))
	   zb(nza+i) = zb(nza-1+i)+ dz(i)
 45	continue
c----- end of the anode axial mesh ---------------------
c
c axial calculation of the cell's thickness dz and
c of the coordinate of the middle point z of ech cell

	z(0) = 0.0
	dz(0) = 0.0
	zb(nzp)= zb(nz)

	do 46 i = 1,nzp
	   dz(i)=zb(i)-zb(i-1)
	   z(i) = zb(i-1)+0.5*dz(i)
 46	continue
c-----------------------------------------------
c end of the axial mesh calculation
c----------------------------------------------

c---------------------------------------------
c radial mesh calculation
c--------------------------------------------

c------radial cathode mesh---------------------------

c 1) mesh uniform over the flat tip radius wtip and
c    nuni poit 

	runi = wtip
	db= runi/(float(nuni)-0.5)
	rb(0)=0.
	rb(1)=db/2.
	do 30 j = 2, nuni
	   rb(j)= rb(j-1)+db
 30	continue
c----------------------------------------------
c 2) mesh until the end of the cathode radius
c    with the same number of point as for the axial 
c    conical part of the cathode (nntip).
c    the thickness of each cell is given by
c    dr = dz*tan(ang)
	do 31 j=1,nntip
	   dr(nuni+j)=dz(nzc+1-j)*tan(ang)
	   rb(nuni+j)=rb(nuni-1+j)+dr(nuni+j)
 31	continue
	nrc=nuni+nntip
c----- end radial cathode mesh ----------------------------
c
c---- radial mesh to the wall ----------------------------
c variable mesh until the radial wall with nrout point
c and the thickness of the first cell equal to the thickness
c of the last cell in the cathode. 
	ipro = 1
	nn = nrout
	zz = radius-rb(nrc)
	db= dr(nrc)
	call multcoef(ipro,zz,db,nn,a)
	do 32 j= nrc+1, nr
	   dr(j) = db*(a**float(j-nrc-1))
	   rb(j) = rb(j-1) + dr(j)
 32	continue
c------ end radial mesh to the wall ------------------------
c
c axial calculation of the cell's thickness dr and
c of the coordinate of the middle point r of ech cell

	r(0) = 0.0
	dr(0) = 0.0
	r(1) = 0.0
	dr(1) =2.* rb(1)
	rb(nrp)= rb(nr)

	do 33 j= 2,nrp
	   dr(j) = rb(j)-rb(j-1)
	   r(j) = rb(j-1)+0.5*dr(j)
 33	continue
c-------------------------------------------------
c end of the radial mesh calculation
c-------------------------------------------------
c
c cathode radial j boundary coordinate nrca irca
	if(ang.lt.3.10/2.) then
	   zza = -tan(ang)
	   zzb = wtip + cathode*tan(ang)
	else
	   zza = 0.0
	   zzb = radius
	endif

	do 51 i = 0, nzc
 	     zzr = zza *z(i)+zzb
	     irca(i)=0
	do 56 j = 0,nr
	      if (rb(j).lt.cradius-1.d-6.and.rb(j).lt.zzr) nrca(i)=j	      
 56	continue      
 51	continue      

c anode radial j boundary coordinate nran iran
	do 55 i = nza, nzp
	   iran(j)=0
        do 54 j = 0, nrp
	      if (r(j).lt.aradius) nran(i)=j
 54	continue
 55	continue
c-------------------------------------------------
	else
c       rectangular cathode

        cradius = rdiam/2.
        aradius = adiam/2.

	nr = nrcr+nrout
	nrp = nr+1
	nrm = nr-1
	nz = nzcr+ngap+nanode
	nzp = nz+1
	nzm = nz-1
	nza = nzcr + ngap +1
	nrspc = 0

c----------------------------------------------------------------
c axial mesh calculation
c---------------------------------------------------------------
c
c calculation of the axial cell boundaries
c
c----------- axial cathode mesh -------------------------------
	zb(nzcr) = rcathode
c--------------------------------------------------------------
c 1) variable mesh across the cathode
c    shca for the first cell thickness from the tip
c    and with nzcr points


	ipro = 1
	nn = nzcr
	zz = rcathode
	db = shca
	call multcoef(ipro,zz,db,nn,a)	
	do 141 i= 1,nn
	   dz(i)= db*(a**float(i-1))
	   zb(nzcr-i)=zb(nzcr+1-i)-dz(i)
 141	continue
c----- axial plasma gap mesh --------------------------

c--------------------------------------------------------
c 1) first half of the gap from the cathode tip.
c    varible mesh with shca for the first cell's thickness
c    and ngap/2 points

	ipro = 1
	nn = ngap/2
	zz = gap/2
	db = shca
	call multcoef(ipro,zz,db,nn,a)
	do 143 i= 1, ngap/2
	   dz(i) = db*(a**float(i-1))
	   zb(nzcr+i)=zb(nzcr-1+i)+dz(i)
 143	continue
c------------------------------------------------------
c 2) second half of the gap until the anode.
c    variable mesh with shan for the cell's thickness in 
c    front of the anode.
	zb(nza-1) = rcathode + gap
	ipro = 1
	nn = ngap-ngap/2
	zz = gap/2
	db = shan
	call multcoef(ipro,zz,db,nn,a)
	do 144 i= 1, nn-1
	   dz(i) = db*(a**float(i-1))
	   zb(nza-1-i) = zb(nza-i)-dz(i)   
 144	continue
c------- end of the axial plasma gap mesh -------------
c
c------  anode axial mesh ----------------------------
c variable mesh with shan for the first cell's thickness
c and naonde point.

	ipro=1
	nn = nanode
	zz = anode
	db = shan
	call multcoef(ipro,zz,db,nn,a)
	do 145 i=0,nn-1
	   dz(i) = db*(a**float(i))
	   zb(nza+i) = zb(nza-1+i)+ dz(i)
 145	continue
c----- end of the anode axial mesh ---------------------
c
c axial calculation of the cell's thickness dz and
c of the coordinate of the middle point z of ech cell

	z(0) = 0.0
	dz(0) = 0.0
	zb(nzp)= zb(nz)

	do 146 i = 1,nzp
	   dz(i)=zb(i)-zb(i-1)
	   z(i) = zb(i-1)+0.5*dz(i)
 146	continue
c-----------------------------------------------
c end of the axial mesh calculation
c----------------------------------------------

c---------------------------------------------
c radial mesh calculation
c--------------------------------------------

c------radial cathode mesh---------------------------

c 1) variable mesh across the cathode
c    shca for the first cell thickness from the tip
c    and with nzcr points
	rb(nrcr) = cradius

	ipro = 1
	nn = nrcr
	zz = cradius
	db = shca
	call multcoef(ipro,zz,db,nn,a)	
	do 241 j= 1,nn
	   dr(j)= db*(a**float(j-1))
	   rb(nrcr-j)=rb(nrcr+1-j)-dr(j)
 241	continue

c----- end radial cathode mesh ----------------------------
c
c---- radial mesh to the wall ----------------------------
c variable mesh until the radial wall with nrout point
c and the thickness of the first cell equal to the thickness
c of the last cell in the cathode. 
	ipro = 1
	nn = nrout
	zz = radius-rb(nrcr)
	db= shca
	call multcoef(ipro,zz,db,nn,a)
	do 232 j= nrcr+1, nr
	   dr(j) = db*(a**float(j-nrcr-1))
	   rb(j) = rb(j-1) + dr(j)
 232	continue
c------ end radial mesh to the wall ------------------------
c
c axial calculation of the cell's thickness dr and
c of the coordinate of the middle point r of ech cell

	r(0) = 0.0
	dr(0) = 0.0
	r(1) = 0.0
	dr(1) =2.* rb(1)
	rb(nrp)= rb(nr)

	do 233 j= 2,nrp
	   dr(j) = rb(j)-rb(j-1)
	   r(j) = rb(j-1)+0.5*dr(j)
 233	continue
c-------------------------------------------------
c end of the radial mesh calculation
c-------------------------------------------------
c
c cathode radial j boundary coordinate nrca irca

	do 251 i = 0, nzcr
	     irca(i)=0
	     nrca(i)=nrcr	      
 251	continue      

c anode radial j boundary coordinate nran iran
	do 252 i = nza, nzp
	   iran(j)=0
        do 253 j = 0, nrp
	      if (r(j).lt.aradius) nran(i)=j
 253	continue
 252	continue
	endif

c j coordinate for the nozzle nrspc
	do 53 j=1,nr
        if((rb(j).ge.dnozz/2.).and.(rb(j-1).lt.dnozz/2.))nrspc=j
 53	continue
        if(nrspc.le.1)nrspc=nr	
c----------------------------------------------------------
c calculate the volume and surfaces of each elementary cell
	call volume
c define the position of each cell,
	call cellpos

c----------------------------------------------------------
c save mesh data in file lsmesh
	open(19,file='lsmesh')
	write(19,*) 'j, r, rb, dr'
	do 779 j=0,nrp
	   write(19,*) j, r(j), rb(j), dr(j)
 779	continue

	write(19,*) 'i, z, zb, dz, nrca, nran'
	do 778 i=0,nzp
	   write(19,*) i, z(i), zb(i), dz(i), nrca(i), nran(i)
 778	continue   
	close(19)
c-----------------------------------------------------------
        return
      end

c-------------------------------------------------------
c end subroutine meshls       
c-------------------------------------------------------

c--------------------------------------------------------
c     subroutine multcoef
c     v.1 ls 15.9.98
c--------------------------------------------------------

	subroutine multcoef(ipro,z,fst,n,a)

      implicit double precision (a-h,o-z)
	alp=z/fst

	if(abs(z/fst-float(n)).lt.1.0e-4)then
	go to 12
	endif
	if(ipro.eq.1)then
	   x1=1.0001d0
	   y1=fun(x1,z,n,fst)
	   x2 = x1
  16	continue
		x2=x2+.1d0
		y2=fun(x2,z,n,fst)
		if(y1*y2.gt.0.0)go to 16
	else
	   x1=0.999
	   y1=fun(x1,z,n,fst)
	   x2=x1
  17    continue	
		x2=x2-0.1d0
		y2=fun(x2,z,n,fst)
		if(y1*y2.gt.0.0d0)go to 17
	endif

 10	x=x1-(x2-x1)*y1/(y2-y1)
	y=fun(x,z,n,fst)
	k=k+1
	if(abs(y).lt.2.0e-5)go to 12
	if(y*y1.gt.0.0d0)then
	x1=x
	y1=y
	go to 10
	else
	x2=x
	y2=y
	go to 10
	endif
  12	continue	
	a=x
	end

c---------------------------------------------------
c end subroutine mulcoef
c---------------------------------------------------

c----------------------------------------------------
c     function fun
c     v.1 ls 21.9.98
c----------------------------------------------------

	function fun(x,z,n,fst)	
      implicit double precision (a-h,o-z)
	fun=fst*((1.d0-x**float(n))/(1.d0-x))-z
	end

c----------------------------------------------------
c  end function fun 
c---------------------------------------------------
c-----------------------------------------------------
c     subroutine volume
c     v.1 ls 15.9.98
c 
c calculates the volume and the surfaces of each elementary
c cell. calculates also the ponderation factor used for the 
c geometric average calculation.
c an top surface of the cell
c ae external radial surface of the cell
c vol volume of the cell
c f radial ponderation factor
c g axial ponderation factor
c---------------------------------------------------- 

      subroutine volume
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)

      an(0) = 0.d0
      an(1) = 0.125*dr(1)*dr(1)
      do 30 j=2,nrp
         xa = r(j)*dr(j)
         an(j) = xa
  30  continue
      do 50 i=0,nz
      do 51 j=0,nr
         xa = (r(j)+0.5*dr(j))*dz(i)
         ae(j,i) = xa
  51  continue
  50  continue
      do 70 j=0,nr
         f(j) = dr(j)/(dr(j) + dr(j+1))
  70  continue
      do 80 i=0,nz
         g(i) = dz(i)/(dz(i) + dz(i+1))
  80  continue
      do 90 i=1,nz
      do 91 j=1,nr
        vol(j,i) = an(j)*dz(i)
  91  continue
  90  continue

      return
      end

c------------------------------------------------------
c end subroutine volume
c------------------------------------------------------

c------------------------------------------------------
c     subroutine cellpos
c     v.1 ls 21.9.98
c
c assigns number to cells in the calculation domain.
c ica=1 for cathode =0 elsewhere
c ian=1 for anode , =0 elsewhere
c-------------------------------------------------------

        subroutine cellpos

        implicit double precision (a-h,o-z)
        parameter(mr=200, mz=100, nsh=21)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
	
        do 9 i=0,mz
        do 10 j=0,mr
	   ica(j,i)=0
	   ian(j,i)=0
 10     continue
 9      continue
c--------------------------------------------------------------
c  ica=1 inside the cathode and 0 elsewhere
        do 121 i=0,nzc
	if(nrca(i).le.irca(i))go to 121
	ist=0
	ien=nrca(i)
	if(ien.eq.nr)ien=nrp
	if(irca(i).gt.1)ist=irca(i)
        do 21 j=ist,ien
	ica(j,i)=1
 21	continue
 121	continue

c  ian=1 inside the anode and 0 elsewhere
        do 122 i=nza,nzp
        do 22 j=0,nra
            ian(j,i)=1
 22	continue
 122	continue
c----------------------------------------------------------------
	return
        end

c-----------------------------------------------------------------
c end subroutine cellpos
c-----------------------------------------------------------------
c-------------------------------------------------------------------
c     subroutine readmatf
c     v.1 ls 16.9.98
c    material functions for plasma, solid cathode, solid anode and liquid weldpool
c-------------------------------------------------------------------

	subroutine readmatf
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      dimension akappaan(41),akappaca(41),akappaliqu(41),akappa(41)
      common/matf/
     r  ah(41), at(41), aet(41),azet(41),asig(41),arho(41),
     r  acp(41),au(41), aan(41),aneth(41), amue(41), aschm(41),
     r  casig(41),cazet(41),ansig(41),anzet(41),cacp(41),ancp(41),
     r  carho,anrho,attach,tempattach,gamma,gammai
      common/liqu/  asigliqu(41),tmla,rholiqu,
     1  rhocoeff,cpliqu,cpcoeff,npr 
      common/emiss/
     r  emissc,emissa
      common/rich/
     r   potion, cworkf, cworke, carich, aworkf, anrich, tcfi
      common/cath/
     r  dtdh, dsdt(0:mr, 0:mz),rc

      open(1,file='matf')
      rewind 1
c attach is attachment rate coef cm3s-1, tempattach is temperature above which attachment is zero
      read(1,*)
      read(1,*) attach,tempattach
      read(1,*)
      read(1,*) potion, cworkf, cworke, carich, aworkf, anrich
      read(1,*)	
c      read(1,*) (at(i),azet(i),asig(i),arho(i),ah(i),acp(i),aet(i),
c     1  au(i),aan(i),aneth(i),amue(i),i=1,41)
      read(1,*) (at(i),akappa(i),asig(i),arho(i),ah(i),acp(i),aet(i),
     1  au(i),aan(i),aneth(i),amue(i),i=1,41)
      read(1,*)
      read(1,*) (akappaca(i),i=1,41)
      read(1,*)
      read(1,*) (casig(i),i=1,41)
      read(1,*)
      read(1,*) (akappaan(i),i=1,41)
      read(1,*)
      read(1,*) (ansig(i),i=1,41)
      read(1,*)
      read(1,*)emissc ,carho, (cacp(i),i=1,41)
      read(1,*)
      read(1,*)emissa ,anrho, (ancp(i),i=1,41)
      close(1)
c       electron and negative ion recombination with positive ions
c recombination coefficients for positve ion  and electrons; gammai and gamma
      gammai = 1.0E-7
      gamma = 1.0E-7
c read weldpool liquid material functions;  jjl
      open(1,file='matfliquid')
      rewind 1
      read(1,*)
      read(1,*)
      read(1,*) (akappaliqu(i),i=1,41)
      read(1,*)
      read(1,*) (asigliqu(i),i=1,41)
      read(1,*)
      read(1,*)  tmla,rholiqu,rhocoeff
      read(1,*)
      read(1,*)  cpliqu,cpcoeff
      close(1)
c-----------------------------------------------------------------
	ah1=ah(1)
	do 10 i=1,41
	ah(i)=ah(i)-ah1
 10	continue
c------------------------------------------------------------------
c  input thermal conductivities, kappa, divided by cp for gas, 
c  use cp of gas even for electrodes
c  ener assumes zet is kappa/cp where cp is for gas; acp(1)=acp(2)
        dtdh = 1700.d0/(ah(3) - ah(1))
        do 30 j = 1, 41
           cazet(j) = akappaca(j)/acp(j)
           anzet(j) = akappaan(j)/acp(j)
	   azet(j) = akappa(j)/acp(j)
 30     continue
c-----------------------------------------------------------------------
cls calculation of the schmitz(T) function
	aschm(1)=0.d0
	do 20 i=2,41
	   zkappa=(azet(i)*acp(i)+azet(i-1)*acp(i-1))/2
	   aschm(i) = aschm(i-1)+zkappa*(at(i)-at(i-1))
 20     continue
c----------------------------------------------------------------------- 
      cazet(41) = cazet(40)
      anzet(41) = anzet(40)
	end
c end of subroutine readmatf
c-------------------------------------------------------------------

c-------------------------------------------------------------------
c     subroutine readarcin
c     v.1 ls 16.9.98
c-------------------------------------------------------------------
	
        subroutine readarcin(currm,flowin,rc,twall,changf,ifread,
     1  te,tc,rlxphi,rlxp,rlxvr,rlxvz,rlxt,rlxne,deltel,didt,delt,
     2  p0,jscreen,rlamt,tamb,rvpol,rrrv,rsdl,rldl,ntimesteps,mm,
     4  ielec,nra,nzci,istore,lte,nnn,npr,nnel,
     5  actest)
c	subroutine readarcin(currm,flowin,rc,twall,
c     1  te,tc,rlxphi,rlxp,rlxvr,rlxvz,rlxt,rlxne,mm,ifread,
c     2  p0,jscreen,rlamt,tamb,rvpol,ntimesteps,
c     3  npr, changf,actest,delt,didt,rrrv,rsdl,
c     4  rldl,ielec,nra,nzci,istore,lte,nnn,deltel)
      implicit double precision (a-h,o-z)
      logical rvpol,ltest,actest,test
      parameter(mr=200, mz=100, nsh=21)
      common/plot/ iplot(50),jplot(50),iplotn,jplotn,ib,jend
c rc=0.1 is source rel limit
c ltest, changf no longer used; jjl
c if jscreen = 1 writes faults to screen
	nzci=49
c	nra=44
      open(3,file='arcin')
      rewind 3
      read(3,*)	
c if lte = 1, lte used for sigma, not ne      
      read(3,*) currm, flowin,lte,nnn
      read(3,*)	
      read(3,*) rc,rlxphi,rlxp,rlxvr,rlxvz,rlxt,rlxne
      read(3,*)	
      read(3,*) twall,te,tc
      read(3,*)	
      read(3,*) mm, ifread, rsdl, istore
      read(3,*) ib,jend,iplotn,jplotn,p0,jscreen
c ib, jend first i and endj for contour plot;iplot,jplotn arrow numbers
      read(3,*) test,rlamt,tamb,ltest,rvpol,npr, changf
c if test true reads elements for vector plot; if j.ge.npr insulating layer 
      read(3,*) 
c voltion defines E from arc length. If below E, source term of ionisation, sc=S=0
      read(3,*) actest,delt,didt,ntimesteps,rrrv,deltel,nnel
c actest true, time dependent; delt time step; didt current change A/s; number time steps 
        close(3)
	end

c---------------------------------------------------------------------
c end of subroutine readarcin
c---------------------------------------------------------------------

c---------------------------------------------------------------
c     subroutine readfn
c     v.1 ls 16.9.98
c---------------------------------------------------------------


      subroutine readfn(ir,nfl,time,curr,js,je,is,ie)
      	implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
      common/var/
     r  dni(0:mr,0:mz),dnp(0:mr,0:mz),rnet(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  vr(0:mr,0:mz),vz(0:mr,0:mz),t(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),
     r  p(0:mr,0:mz),h(0:mr,0:mz),phi(0:mr,0:mz),rrne(0:mr,0:mz),
     r  rne(0:mr,0:mz),rneth(0:mr,0:mz),eohm(0:mr,0:mz),
     r  rni(0:mr,0:mz),rnp(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),
     r  ran(0:mr,0:mz),rtest(0:mr,0:mz),SS(0:mr,0:mz),dne(0:mr,0:mz)
            common/dens/
     r  rho(0:mr,0:mz), rhmr(0:mr,0:mz),rhmz(0:mr,0:mz)

        if(ir.eq.1)then
         read(nfl,496) time,curr
c         read(nfl,520) ((t(j,i),j=js,je),i=is,ie)
         read(nfl,520) ((h(j,i),j=js,je),i=is,ie)
         read(nfl,520) ((vz(j,i),j=js,je),i=is,ie)
         read(nfl,520) ((vr(j,i),j=js,je),i=is,ie)
         read(nfl,520) ((p(j,i),j=js,je),i=is,ie)
         read(nfl,520) ((phi(j,i),j=js,je),i=is,ie)
	 read(nfl,520) ((rne(j,i),j=js,je),i=is,ie)
	 read(nfl,520) ((rneth(j,i),j=js,je),i=is,ie)
        else
         write(nfl,496) time,curr
c         write(nfl,520) ((t(j,i),j=js,je),i=is,ie)
         write(nfl,520) ((h(j,i),j=js,je),i=is,ie)
         write(nfl,520) ((vz(j,i),j=js,je),i=is,ie)
         write(nfl,520) ((vr(j,i),j=js,je),i=is,ie)
         write(nfl,520) ((p(j,i),j=js,je),i=is,ie)
         write(nfl,520) ((phi(j,i),j=js,je),i=is,ie)
	 write(nfl,520) ((rne(j,i),j=js,je),i=is,ie)
	 write(nfl,520) ((rneth(j,i),j=js,je),i=is,ie)
        endif
 496     format(2e15.5)
 519     format (3e25.8)
 520     format (5e16.8)
        end

c----------------------------------------------------------
c end of subroutine readfn
c----------------------------------------------------------

c----------------------------------------------------------
c     subroutine cont
c     v.1 ls 18.9.98
c----------------------------------------------------------



c----------------------------------------------------------------
c     subroutine ppp
c     v.1 ls 16.9.98
c------------------------------------------------------------------
c     subroutine meshrec
c     v.1 ls 21.9.98
c------------------------------------------------------------------
	
	subroutine meshrec
      implicit double precision (a-h,o-z)
      parameter(mr=200, mz=100, nsh=21)
        common/coord/
     r  dr(0:mr),dz(0:mz),r(0:mr),z(0:mz),ppp(0:mz),
     r  an(0:mr),ae(0:mr,0:mz),vol(0:mr,0:mz),rr(0:mr),drr(0:mr),
     r  f(0:mr), g(0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,
     i  nzc, nrc, nza, nra, nrspc,nuni,nntip, 
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz), ian(0:mr,0:mz)
	
      open(4,file='coord')
        rewind 4
      write(4,30) (z(i),i=0,nzp)
      write(4,40) (r(j),j=0,nrp)
        close(4)
        open(901,file='mesh')
        rewind 901
c        write(901,902)nz,nr,nzc,nra,nrspc,nza
        do 904 i=0,nzc
          write(901,903)nrca(i)
 904  continue
        do 905 i=0,nz+1
            write(901,*)dz(i)
 905  continue
        do 906 i=0,nr+1
            write(901,*) dr(i) 
 906  continue
        close(901)
        open(901,file='zcoor')
        rewind 901
        do 907 i=0,nz+1
        write(901,*)i,nrca(i),dz(i),z(i)
 907    continue
        close(901)
        open(901,file='rcoor')
        rewind 901
        do 919 j=0,nr+1
        write(901,*)j,dr(j),r(j)
 919    continue
        close(901)
  30  format(1x,'nodal positions - axial direction'/(9f8.5))
  40  format(1x,'nodal positions - radial direction'/(9f8.5))
 902    format(6(1p,i3))
 903    format(i3)
	end

c--------------------------------------------------------------
c end subroutine meshrec
c-------------------------------------------------------------
