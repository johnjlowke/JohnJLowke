c jetair6 has initialzero subroutine 
c with backward differences
c line 447 to make attachment finite at zero field
c diffusion correction line 370; flux for j=0 given
c jetair3.f has zero charge emission inside pin
c arglow6 solves  poisson equation (dE/dx = (e/epsilon)= net charge) insted of current continuity.
c arglow4 & 5 removes backward differencing.
c arglow3.f uses backward differences
c max alpha increased from 100 to 1000
c arglow1.f series are for argon instead of air
c based on jetlow11a.f no attachment; no singlet delta excited state but rnm() not removed.
c also no negative ion species, rni retained but not used.
c jetlow11a.f has removal of backward difeerencing for convection. Instead considers flow into cell.
c jetlow9.f includes space charge effects
c jetlow8.f is same as letlow6a.f
c jetlow6.f has subroutine zerocoef  added to jetlow4.f Sets any value at j,i to be a fixed value dfval
c jet1.f has hail removed with dc voltage at the upper plane
c hail1 is copy of cyl2.f
c cyl has cylindrical cloud of negative ions to produce lightning field
c cyl2.f has upper boundary condition  voltage, potb(j) being boundary voltage
c contains correction to make attachment rate const near zero electric field
c cyl.f has whole electron clouud, not with d phi/dz = 0 at uppe boundary; json removed
c cylstep codes have d phi/dx = 0 at z = 0 to simulate z = 0 asplane through negative cloud centre
c cylstep4.f has printout of max ne for each iteration
c cylstep4.f creates json files for potential, neg and pos ions, net neg and pos ions, metasbles and electric field.
c cylstep4.f has alpha limited to 100
c cylstep1.f has metastable influence removed
c steplead6.f has dphi/dz = 0 at z = 0
c steplead5b.f removes limit on alpha
c steplead5a.f removes alpha limit of 500 to 1000
c steplead5.f has recombination for ions reduced to 1.0E-8
c steplead3.f has alpha limit of 500
c steplead2.f is for test of step leader resulting from metastable oxygen
c steplead1 is cloud8a.f
c cloud8a.f has low eta values
c cloud8.f has photoionization term added
c cloud6h.f has corrected printout for Abs field
c cloud6h.f has Petrova alphas but high eta values from 25 kV/cm sustaining arcin field.
c cloud6g.f has Tam representation of alphaon from Petrova
c cloud6fg.f has alpha limited to 10000.0
c cloud6f arcfn input not distorted by input cloud charges
c cloud6f corrected eta values for sustaining field of 5 kV/cm
c cloud6d has average electic field for production of ionization as average of field at the four faces
c cloud6c has gamma for electrons 5E-7 and gammai for ions 1E-7 for recombination
c cloud6a has removal of ion velocity cap and correction to "wire" lines 334, 337,338
c cloud6a has printout for ion fluxes
c cloud6 has increased array boundaries and includes streamer of uniform z grid size; L814 file 32 available
c cloud4a corrects for "wings" at top and bottom of cyinder
c cloud4 treats discharge development;  limits set on alpha and drift velocities to limit charge densities.
c cloud3 has voltage input at cloud plane for a sphere of negative charge, not just a plane
c cloud3 has average eabs determined by field at metal-atmosphere interface - only average is RMS Er and EZ
c cloud2 corrected for negative voltage and definitions, lines 124 692 693
c cloud1.f is leader5.f; ggg
c leader5 has net ionization zero at high E/N, line 370
c leader4 has electron diffusion from kT/e
c leader3 has tmemeta corrected; convective terms in radius corrected
c leader2 has ionization increased for metastables for time greater than timemeta
c leader2 reads timemeta in arcin for transition from 25 kV/cm to 5 kV/cm transition field
c leader1b has sig as relative permittivity
c leader1b has axial number of points defined by mz extended from 100 to 120
c leader1a has sig a multiplying factor for field calculation
c leader1.f is a combination of codes boeing5d.f and ballbirth6c.f
c boeing5d.f is code to calculate electric  fields around a cylinder capped
c    by hemispheres under a cloud for lightning calculations; 100MV over 2 km
c    boeing5d.f is in folders lightning, Loeb, Boeing, Boeing2
c ballbirth6c.f is code of Ball Lightning calculations  of J Geophys Res
c    ballbirth6c.f is in folders lightning, BallBirth, Short, Plot
c Voltage at top plane not at Volt(nz); see fill. Voltage at bottomn plane, not at Volt(1)
c sig e-3, e6
c boeing4.f has anode removed,nntop is top of object nnbot is bottom
c boeing3.f has hemisphere on cylinder
c boeing1.f has dv/dr = 0
c boeing0.f is loeb5a.f
c loeb4 high voltage on point
c loeb3 is for hemisphere on a rod
c Loeb2 has zero volts at z = 0; sig = 100000 and 1.0E-6 for metal and air
c Laplace equation solution between upper plane and lower sphere on rod.
c BernardiICLP Florence
c Aug 2011
c-----------------------------------------------------------------------
      program  brouwer
      implicit double precision (a-h,o-z)
c-----------------------------------------------------------------------
      character*40 stfile,fname
      parameter(mr=150, mz=220, nsh=21)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),df(mr,mz)
      common/elc/
     r  sig(0:mr,0:mz),sigr(0:mr,0:mz),sigz(0:mr,0:mz)
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz)
      common/var/
     r  phi(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rni(0:mr,0:mz),rne(0:mr,0:mz),rnp(0:mr,0:mz),
     r  rrne(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),rnm(0:mr,0:mz),
     r  rrnet(0:mr,0:mz),dne(0:mr,0:mz),dni(0:mr,0:mz),dnm(0:mr,0:mz),
     r  dnp(0:mr,0:mz),rnet(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),rlognetp(0:mr,0:mz),
     r  rlognetn(0:mr,0:mz),rlognm(0:mr,0:mz),rrnm(0:mr,0:mz),
     r  rphi(0:mr,0:mz)
      common/dummy/
     r  resphi,relphi,rlxphi,ie,je,iph,jph,ifread,mprint
      common/dummy2/m,mm2
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),rr(0:mr),drr(0:mr),
     r  an(0:mr), ae(0:mr,0:mz),f(0:mr),g(0:mz),vol(0:mr,0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,nri,nro,pinl,
     i  nzc, nrc, nnbot, nra,nnrad,ngapup,radius,volts,
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz)
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz)
      common/plot/iplot(50),jplot(50),iplotn,jplotn,ib,jend,iend
        stfile='arcst'
c gas number density at 1 bar
      rn = 2.5E19
c recombination coefficient, gamma for electrons, gammai for ions, for air
c        gamma = 5.0E-7
c        gammai = 1.0E-7
        gammai = 1.0E-8
c recombination coefficient, gamma for electrons and argon ions,taken as 5.0E-9
c    see von Engel "Ionizedf Gases"  1965 p163
        gamma = 5.0E-9
c ion mobility where wi = rmob*E in cm/s; E in V/cm
      rmob = 2.0
c quenching coefficient
      rkq = 2.22E-18
c detachment coefficient
      rkd = 2.0E-10
      epsilon = 8.85d-14
c-----------------------------------------------------------------
c initialise the variables to zero
        call initialise
c------------------------------------------------------------------
c reads the input parameters in file 'arcin'
        write(*,*)'started'
c        call readarcin(volts,rlxphi,mm1,ifread,iend,delt,
        call readarcin(volts,rlxphi,mm1,ifread,delt,
     r    etime,mm2,mprint,freq,meta)
c------------------------------------------------------------------
c generates the mesh from the data read in file 'gridbern'
        write(*,*)'read  arcin'
        call mesh
        write(*,*)'mesh generation completed'
        jend = nr
        iend = nz
        ib = 1
c-----------------------------------------------------------------------
       write(*,25) rlxphi,mm1
       if (ifread.eq.0) write(*,26)
       if (ifread.eq.1) write(*,27)
       if (ifread.eq.2) write(*,28)
c-----------------------------------------------------------------
c     initial conditions:
c if ifread = 0 start from scratch.
c if ifread = 1 start from file phiin
c if ifread = 2 start from arcfn
c----------------------------------------------------------------
      if (ifread .eq. 1) then
        open(2,file='phiin')
        read(2,526) ((phi(j,i),j=0,nrp),i=0,nzp)
        close(2)
      endif
      do 54 i=0,nz+1
      do 55 j=0,nr+1
          if (time.le.1.0e-11) then
	      rni(j,i) = 1.0
	      rne(j,i) = 1.0
	      rnp(j,i) = 1.0
	      rnm(j,i) = 1.0

	  endif
c	  rrr = r(j)
c	  rad = sqrt(zzz**2 + rrr**2)

  55  continue
  54  continue
      write(*,*) tlen, tlen10
      if (ifread .eq. 2) then
         open(304,file='arcfn')
         rewind 304
         read(304,526) time
         read(304,526) ((rni(j,i),j=1,nr),i=1,nz)
         read(304,526) ((rne(j,i),j=1,nr),i=1,nz)
         read(304,526) ((rnp(j,i),j=1,nr),i=1,nz)
         read(304,526) ((rnm(j,i),j=1,nr),i=1,nz)
         read(304,526) ((phi(j,i),j=1,nr),i=1,nz)
         close(304)
      endif
c-----------------------------------------------------------------
c     start  iterations for initial field
      do 200 m=1,mm1
         call poiss
         write(*,850) m, resphi,relphi
c         volts = phi(1,1)
 200  continue
c
c     end iterations for initial field
c
c-----------------------
   9  format(' Patankar code, solution of Laplace equn,sphere on rod.'/)
  10  format(1x,'Voltage =',f13.0,', nz=',i3,', nr=',i3,', nzc=',i3,
     1    ', nnbot=',i3)
  25  format(1x,'rlxn. factor =',f6.2,',',' no. of iterations =',i5)
  26  format(1x,'initial phi fields generated by program')
  27  format(1x,'initial phi fields from <phiin>')
  28  format(1x,'initial phi fields from <arcfn>')
  29  format(' mm1 =',i4,' volts =',e13.3)
  30  format(1x,'nodal positions-axial direction'/(10f11.2))
  34  format(2(e11.4))
  40  format(1x,'nodal positions-radial direction'/(10f11.3))
  50  format(1p6e13.5)
  53  format(5(e11.4))
 181  format(' time= ',e10.4,' dt = ',e10.4,'  mm2 = ',i6)
 463    format(2x,5(3x,f6.3),3x,i5)
 495  format(2e15.5)
 496     format(2e15.5)
 500     format(8f9.0)
 501     format(/)
 502     format(8f9.0)
 503  format(/'  rjnez')
 504  format(/'  rjner')
 505  format(/'  rjnpz')
 519  format(/'  rjnpr')
 527  format(/'  rjniz')
 528  format(/'  rjnir')
c 506  format('  hail')
 507  format('  potential')
 508  format('  electric field Ez')
 509  format('  electric field Er')
 510  format('  i =',i4,'   z =',f11.4)
 511  format(14f7.2)
 512  format(14f6.0)
 513  format(11f8.0)
 514  format(/'  Negative Ion Density')
 516  format(9f7.0)
 518  format (3e24.8)
 520  format (5e16.9)
 521  format (1pe10.1,1pe10.1,1pe10.1,1pe10.1,1pe10.1,1pe10.1,
     1       1pe10.1,1pe10.1,1pe10.1,1pe10.1)
 522  format(/'  Electron Density')
 523  format(/'  Positive Ion Density')
 524  format(/'  Net Charge Density')
 525  format(/'  Absolute Field')
 526  format (5e16.8)
 572  format(1p9e11.3)
 573  format(25i3)
 620  format(1i3,6e12.3)
 796  format(200e11.3)
 797  format(200f12.3)
 820  format(/'  Metastable Density')
 850  format(1x,'iteration no.',i6,
     * '   maximum residual =', 1pe9.2,'   max relerror =',e9.2)
 851  format(1x,'iteration no.',i6,
     * '   maximum residual =', 1pe9.2,'   max relerror =',e9.2,
     * ' time =',e9.2)
 852  format(1x,'iteration no.',i6,'  time =',e12.4,
     r '  max electron density =',1pe9.2)
 854  format(' max. res. phi   (0.1 %)    =',e12.3,' at node',2i4,e12.3)
 859  format('     res    je     ie      phi        rel    jje  ',
     r 'iie    phi')
 861  format(' resphi=',1pe12.2,2i4,1pe12.2)
 899     format(1x,'iteration no. ',i4)
 909    format(1x,i3,6(2pe12.4))
 910  format(/'  j    rion    alph   att   eta  rece  dne  rne '/)
 911  format(i4,9e12.3)
 978  format(i3.3)
 990  format(1x,'   j      rion        alpha        rec         dne  
     r       ne          rner        rnez')
 995  format(1p6e13.5)
 996  format(1pe10.2)

      fname='record'
      open(5,file=fname)
      rewind 5
      write(5,9)
      write(5,10) volts,nz,nr,nzc,nnbot
      write(5,25) rlxphi,mm1
      if (ifread .eq. 0) write(5,26)
      if (ifread .eq. 1) write(5,27)
      if (ifread .eq. 2) write(5,28)
        write(5,30) (z(i),i=0,nzp)
        write(5,40) (r(j),j=0,nrp)
        write(5,850) m, resphi,relphi
        write(5,861) resphi,jph,iph,phi(jph,iph)
 227  format(1x,' voltage   =',f13.0)
c       write(5,227) volts
         write(5,501)
      nzcp=nzc+nrca(nzc)-1
      nzc1=nzcp-nzc
 780  continue
         write(5,507)
      do 778 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (phi(j,i),j=0,nrp)
 778  continue
         write(5,501)
         write(5,508)
      do 792 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (ez(j,i),j=0,nrp)
 792  continue
         write(5,501)
         write(5,509)
      do 781 i=0,nzp
         write(5,510) i,z(i)
         write(5,572) (er(j,i),j=0,nrp)
 781  continue
c plot files for MATLAB or IDL
      open(20,file='r.dat')
      write(20,*) (r(j),j=1,jend)
      close(20)
       open(21,file='z.dat')
      write(21,520) (z(i),i=ib,nz)
      close(21)
c     time iterations start here
 189  if (ifread .eq. 0) time = 0.0
      do 187 i=0,nz
	  rjnir(0,i) = 0.0
	  rjner(0,i) = 0.0
	  rjnpr(0,i) = 0.0
 187  continue
      do 188 j=0,nr+1
	  rjniz(j,0) = 0.0
	  rjnez(j,0) = 0.0
	  rjnpz(j,0) = 0.0
 188  continue
      write(5,990)
c      mprint = 253
      mmm = 0
      do 201 m=1,mm2
c      if(m.eq.mm2) then
c          write(*,*) rne(2,33),rni(2,33),rnp(2,33),vol(2,33)
c          write(*,*) dn(2,33),ds(2,33),de(2,33),dw(2,33),dc(2,33),
c     r        df(2,33)
c          write(*,*) rne(2,32),rni(2,32),rnp(2,32),vol(2,32)
c          write(*,*) dn(2,32),ds(2,32),de(2,32),dw(2,32),dc(2,32),
c     r        df(2,32)
c      endif
      call poiss
c      call amf
      mm8 = (nzc+nnbot)/2
      do 791 i=1,nz
      do 784 j=1,nr+1
        ezz = (ez(j,i-1) + ez(j,i))/2
	err = (er(j-1,i) + er(j,i))/2
	eabs(j,i) = sqrt(ezz**2 + err**2)
  784  continue
  791  continue
c diffusion from D/mu = kT/e ~ 1 eV; mu ~ 1 for ions and 100 for electrons
       diff =100.0*rmob
c Radial and axial fluxes with backward differences
       do 193 i=1,nz
       do 194 j=1,nr
          rion = 0.0
          att = 0.0
          rece = 0.0
          reci = 0.0
          wir = -rmob*er(j,i)
          wire = 100*wir
c          if (wire.ge.10000000.) wire = 1000000.
          if (wire.ge.10000000.) wire = 10000000.
          if (wire.le.-10000000.) wire = -10000000.
c          if (wir.ge.100000.) wire = 100000.
c          if (wir.le.-100000.) wire = -100000.
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
c Axial current densities; for negative ions and electrons w = -mu E
c At pin surfaces charges can flow in but not out.
        wiz = -rmob*ez(j,i)
        wize = 100*wiz
        if (wize.ge.10000000.) wize = 10000000.
        if (wize.le.-10000000.) wize = -10000000.
c          if (wiz.ge.100000.) wiz = 100000.
c          if (wiz.le.-100000.) wiz = -100000.
        if(ez(j,i).le.0.0) then
           rjniz(j,i) = wiz*rni(j,i)
            rjnez(j,i) = wize*rne(j,i)
     r	          -diff*(rne(j,i+1)-rne(j,i))/dz(i)
             rjnpz(j,i) = -wiz*rnp(j,i+1)
        else
       rjniz(j,i) = wiz*rni(j,i+1)
       rjnez(j,i) = wize*rne(j,i+1)
     r	          -diff*(rne(j,i+1)-rne(j,i))/dz(i)
        rjnpz(j,i) = -wiz*rnp(j,i)
      endif
 194    continue
 193    continue
c Boundary conditions on fluxes
c At centre, j=1, gradients with radius are zero 
      do 195 i=1,nz
             rjnez(0,i) = rjnez(1,i)
             rjner(0,i) = rjner(1,i)
             rjniz(0,i) = rjniz(1,i)
             rjnir(0,i) = rjnir(1,i) 
             rjnpz(0,i) = rjnpz(1,i)
             rjnpr(0,i) = rjnpr(1,i) 
c At glass boundary no emission or absorption of particles.
             rjner(nr,i) = 0.0
             rjnir(nr,i) = 0.0
             rjnpr(nr,i) = 0.0
c At inside wall of pin ensure only fluxes into metal
             nrim = nri-1
          if (i.le.ngapup) then
             if (er(nrim,i).le.0.0) rjner(nrim,i) = 0.0
             if (er(nrim,i).le.0.0) rjnir(nrim,i) = 0.0
             if (er(nrim,i).ge.0.0) rjnpr(nrim,i) = 0.0
          endif
 195   continue
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
          alpha = alphaon*rn
	  att = rne(j,i)*eta*weav
c adjustment to make attachment rate constant neat zero electric field
          if (eon.le.3.0E-17) att = rne(j,i)*3.75E6
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
c Oxygen metastable coefficients
          if (eon.le.1.0E-17) rmetaon = 0.0
          if (eon.ge.1.0E-17.and.eon.le.7.0E-17) rmetaon =
     r        1.99E-18*(eon-1.0E-17)/9.0E-17
          if (eon.ge.7.0E-17.and.eon.le.70.0E-17) rmetaon = 2.0E-18
          if (eon.ge.70.0E-17.and.eon.le.200.0E-17) rmetaon=
     r        1.0E-17*(eon-7.0e-17)/130.0e-17
          if (eon.ge.200.0E-17) rmetaon= 1.0e-17
          rmeta = rmetaon*rn
          rmet = rne(j,i)*rmeta*weav
	  ddrf = 1+dr(j)/(2*r(j))
	  ddrb = 1-dr(j)/(2*r(j))
          dni(j,i) = -(rjnir(j,i)*ddrf-rjnir(j-1,i)*ddrb)/dr(j) -
     r        (rjniz(j,i)-rjniz(j,i-1))/dz(i)+att-reci -
     r        rkd*rni(j,i)*rnm(j,i)
          rner = -(rjner(j,i)*ddrf-rjner(j-1,i)*ddrb)/dr(j)
          rnez = - (rjnez(j,i)-rjnez(j,i-1))/dz(i)
          dne(j,i) = -(rjner(j,i)*ddrf-rjner(j-1,i)*ddrb)/dr(j) -
     r        (rjnez(j,i)-rjnez(j,i-1))/dz(i)+rion-att-rece +
     r        rkd*rni(j,i)*rnm(j,i)
          if (m.eq.mm2.and.i.eq.100) then
             write(5,911) j,rion,alpha,rece,dne(j,i),rne(j,i)
     r          ,rner,rnez
          endif
          d2dr2= -(rjner(j,i)*ddrf-rjner(j-1,i)*ddrb)/dr(j)
          d2dz2= -(rjnez(j,i)-rjnez(j,i-1))/dz(i)
          recom=   rkd*rni(j,i)*rnm(j,i)
          if (rne(j,i).le.1E-3) rne(j,i) = 1.0E-3
          if (rne(j,i).ge.rmaxne) rmaxne = rne(j,i)
          dnp(j,i) = -(rjnpr(j,i)*ddrf-rjnpr(j-1,i)*ddrb)/dr(j) -
     r        (rjnpz(j,i)-rjnpz(j,i-1))/dz(i)+rion-rece-reci
         if (meta.eq.1)
     r       dnm(j,i)=rmet-rkd*rnm(j,i)*rni(j,i)-rkq*rnm(j,i)*rn/5
 211  continue
 210  continue
c Keep charge densities in pin zero.
      do 190 i=1,ngapup
      do 191 j=nri,nro
          dne(j,i)= 0.0
          dni(j,i)= 0.0
          dnp(j,i)= 0.0
          dnm(j,i)= 0.0
 191  continue
 190  continue
c increment particle densities
      do 198 i=1,nz
      do 197 j=1,nr
          rne(j,i)=rne(j,i)+dne(j,i)*delt
          if (rne(j,i).le.1E-3) rne(j,i) = 1.0E-3
          rni(j,i)=rni(j,i)+dni(j,i)*delt
          rnp(j,i)=rnp(j,i)+dnp(j,i)*delt
          rnm(j,i)=rnm(j,i)+dnm(j,i)*delt
c	  rnet(j,i) = rnp(j,i) - rne(j,i)
	  rnet(j,i) = rnp(j,i) - rne(j,i) - rni(j,i)
 197  continue
 198  continue
      do 196 i=1,nz
          rni(nr+1,i) = rni(nr,i)
          rne(nr+1,i) = rne(nr,i)
          rne(nr+2,i) = rne(nr+1,i)
          rnp(nr+1,i) = rnp(nr,i)
          rnet(nr+1,i) = rnet(nr,i)
          rnm(nr+1,i) = rnm(nr,i)
          rni(0,i) = rni(1,i)
          rne(0,i) = rne(1,i)
          rnp(0,i) = rnp(1,i)
          rnet(0,i) = rnet(1,i)
          rnm(0,i) = rnm(1,i)
 196  continue
      do 199 j=1,nr
          rne(j,nz+1) = 0.0
 199  continue
      time = time + delt
      write(*,852) m,time,rmaxne
c prints every mprint times
      mmm = mmm+1
      if (mmm.eq.mprint) then
         mmm = 0
      endif
 201  continue
c     end internal iteration
      open(304,file='arcfnout')
        rewind 304
         write(304,526) time
         write(304,526) ((rni(j,i),j=1,nr),i=1,nz)
         write(304,526) ((rne(j,i),j=1,nr),i=1,nz)
         write(304,526) ((rnp(j,i),j=1,nr),i=1,nz)
	write(304,526) ((rnm(j,i),j=1,nr),i=1,nz)
	write(304,526) ((phi(j,i),j=1,nr),i=1,nz)
      close(304)
        write(5,181) time,delt,mm2
        write(5,859)
        write(5,40) (r(j),j=0,nrp)
         write(5,507)
      do 788 i=0,nz+1
         write(5,510) i,z(i)
         write(5,572) (phi(j,i),j=0,nr+1)
 788  continue
         write(5,501)
         write(5,514)
      do 785 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rni(j,i),j=0,nr+1)
 785  continue
         write(5,501)
         write(5,522)
      do 786 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rne(j,i),j=0,nr+1)
 786  continue
         write(5,501)
         write(5,523)
      do 787 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rnp(j,i),j=0,nr+1)
 787  continue
         write(5,501)
         write(5,820)
      do 798 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rnm(j,i),j=0,nr+1)
 798  continue
         write(5,501)
         write(5,524)
      do 795 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rnet(j,i),j=0,nr+1)
 795  continue
         write(5,501)
         write(5,508)
      do 779 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (ez(j,i),j=0,nr+1)
 779  continue
         write(5,501)
         write(5,509)
      do 789 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (er(j,i),j=0,nr+1)
 789  continue
         write(5,501)
         write(5,525)
	 mm8 = (nzc+nnbot)/2
      do 790 i=1,nz
         write(5,510) i,z(i)
         write(5,572) (eabs(j,i),j=1,nr+1)
 790  continue
         write(5,501)
         write(5,503)
      do 782 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rjnez(j,i),j=0,nr+1)
 782  continue
         write(5,501)
         write(5,504)
      do 783 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rjner(j,i),j=0,nr+1)
 783  continue
         write(5,501)
         write(5,505)
      do 816 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rjnpz(j,i),j=0,nr+1)
 816  continue
         write(5,501)
         write(5,519)
      do 817 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rjnpr(j,i),j=0,nr+1)
 817  continue
         write(5,501)
         write(5,527)
      do 819 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rjniz(j,i),j=0,nr+1)
 819  continue
         write(5,501)
         write(5,528)
      do 814 i=0,nz
         write(5,510) i,z(i)
         write(5,572) (rjnir(j,i),j=0,nr+1)
 814  continue
       close(5)
c files for MATLAB
      open(21,file='z.dat')
      write(21,797) (z(i),i=ib,iend)
      close(21)
c err and rr files for full plot from r = -R to r = +R
      i = ib
      jjend = 2*jend
      n = 0
        write(*,*)'got here1'
      rr(0) = 0.d0
      do 808 j=1,jjend
	  n = n + 1
          jj = jend-j+1
          if (jj.le.0) then
	      jj = j - jend
              rr(n) = r(jj)
              drr(n) = dr(jj)
          else
              rr(n) = -r(jj)
              drr(n) = dr(jj)
          endif
 808   continue
        write(*,*)'got here2'
       open(23,file='rr.dat')
       rewind 23
      write(23,*) (rr(j),j=1,jjend)
      close(23)
      open(28,file='rrni.dat')
      do 799 i=ib,iend
      do 800 j=1,jjend
	  jj = jend-j+1
          if (jj.le.0) jj = j - jend
          rrni(j,i) = rni(jj,i)
	  if (rrni(j,i).le.1.0E-7) rrni(n,i)=1.0E-7
          rlogni(j,i) = log10(rrni(j,i))
 800  continue
 799  continue
      do 793 i=ib,iend
          write(28,796) (rlogni(j,i),j=1,jjend)
 793  continue
      open(29,file='rrnp.dat')
      do 807 i=ib,iend
      do 809 j=1,jjend
	  jj = jend-j+1
          if (jj.le.0) jj = j - jend
          rrnp(j,i) = rnp(jj,i)
	  if (rrnp(j,i).le.1.0D-7) rrnp(j,i)=1.0D-7
          rlognp(j,i) = log10(rrnp(j,i))
 809  continue
 807  continue
      do 812 i=ib,iend
          write(29,796) (rlognp(j,i),j=1,jjend)
 812  continue
      close(29)
      open(36,file='rrnm.dat')
      do 38 i=ib,iend
      do 37 j=1,jjend
	  jj = jend-j+1
          if (jj.le.0) jj = j - jend
          rrnm(j,i) = rnm(jj,i)
	  if (rrnm(j,i).le.1.0D-7) rrnm(j,i)=1.0D-7
          rlognm(j,i) = log10(rrnm(j,i))
  37  continue
  38  continue
      do 39 i=ib,iend
          write(36,796) (rlognm(j,i),j=1,jjend)
  39  continue
      close(36)
      do 810 i=ib,iend
      do 811 j=1,jjend
	  jj = jend-j+1
          if (jj.le.0) jj = j - jend
          rrnet(j,i) = rnet(jj,i)
	  if (rrnet(j,i).le.1.0D-7.and.rrnet(j,i).ge.-1.0D-7)
     r	     rrnet(j,i) = 1.0D-7
          if (rrnet(j,i).gt.0.0) rlognetp(j,i) = log10(rrnet(j,i))
          if (rrnet(j,i).lt.0.0) then
              rrnet(j,i) = -rrnet(j,i)
              rlognetn(j,i) = log10(rrnet(j,i))
          endif
 811  continue
 810  continue
      open(31,file='netp.dat')
      do 813 i=ib,iend
          write(31,796) (rlognetp(j,i),j=1,jjend)
 813  continue
      close(31)
      open(33,file='netn.dat')
      do 815 i=ib,iend
          write(33,796) (rlognetn(j,i),j=1,jjend)
 815  continue
      close(33)
      open(30,file='rrne.dat')
      do 794 i=ib,iend
      do 801 j=1,jjend
          jj = jend - j + 1
          if (jj.le.0) jj = j - jend
	  rrne(j,i) = rne(jj,i)
	  if (rrne(j,i).le.1.0D-7) rrne(j,i)=1.0D-7
          rlogne(j,i) = log10(rrne(j,i))
 801  continue
 794  continue
      do 802 i=ib,iend
          write(30,796) (rlogne(j,i),j=1,jjend)
 802  continue
      close(30)
      open(27,file='rreabs.dat')
      do 803 i=ib,iend
      do 806 j=1,jjend
	  jj = jend-j+1
          if (jj.le.0) jj = j - jend
          rreabs(j,i) = eabs(jj,i)
 806  continue
 803  continue
      do 805 i=ib,iend
          write(27,796) (rreabs(j,i),j=1,jjend)
 805  continue
      close(27)
        write(*,*)'got here3'
      open(7,file='rphi.dat')
      do 700 i=ib,iend
      do 701 j=1,jjend
	  jj = jend-j+1
          if (jj.le.0) jj = j - jend
          rphi(j,i) = phi(jj,i)
 701  continue
 700  continue
        write(*,*)'got here3a'
        write(*,*) iend,jend,jjend,ib
      do 702 i=ib,iend
          write(7,796) (rphi(j,i),j=1,jjend)
 702  continue
      close(7)
c      open(6,file='rhail.dat')
c      do 703 i=ib,iend
c      do 704 j=1,jjend
c	  jj = jend-j+1
c          if (jj.le.0) jj = j - jend
c          rhail(j,i) = hail(jj,i)
c 704  continue
c 703  continue
c      do 705 i=ib,iend
c          write(6,796) (rhail(j,i),j=1,jjend)
c 705  continue
c      close(6)
c      open(6,file='phiout')
c         write(6,526) ((phi(j,i),j=1,jjend),i=ib,iend)
c      close(6)
c      call printjson
      end
c-------------------------------------------------
c  end main
c--------------------------------------------------------------------------
c  subroutine zerocoef(dc,df,de,dw,dn,ds,dfval,j,i)
c sets value at j,i to be dfval

        subroutine zerocoef
      implicit double precision (a-h,o-z)
      parameter(mr=150, mz=220, nsh=21)
      dimension de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1     df(mr,mz)
           dc(j,i) = 1.d0
           df(j,i) = dfval
           de(j,i) = 0.d0
           dw(j,i) = 0.d0
           dn(j,i) = 0.d0
           ds(j,i) = 0.d0
       end
c end of subroutine zerocoef






c----------------------------------------------------------------
c     subroutine initialise
c     v.1 ls 16.9.98
c----------------------------------------------------------------
c
        subroutine initialise
      implicit double precision (a-h,o-z)
      parameter(mr=150, mz=220, nsh=21)
      common/elc/
     r  sig(0:mr,0:mz),sigr(0:mr,0:mz),sigz(0:mr,0:mz)
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz)
      common/var/
     r  phi(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rni(0:mr,0:mz),rne(0:mr,0:mz),rnp(0:mr,0:mz),
     r  rrne(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),rnm(0:mr,0:mz),
     r  rrnet(0:mr,0:mz),dne(0:mr,0:mz),dni(0:mr,0:mz),dnm(0:mr,0:mz),
     r  dnp(0:mr,0:mz),rnet(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),rlognetp(0:mr,0:mz),
     r  rlognetn(0:mr,0:mz),rlognm(0:mr,0:mz),rrnm(0:mr,0:mz),
     r  rphi(0:mr,0:mz)
      common/dummy/
     r  resphi,relphi,rlxphi,ie,je,iph,jph,ifread,mprint
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),rr(0:mr),drr(0:mr),
     r  an(0:mr), ae(0:mr,0:mz),f(0:mr),g(0:mz),vol(0:mr,0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,nri,nro,pinl,
     i  nzc, nrc, nnbot, nra,nnrad,ngapup,radius,volts,
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz)
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz)

      do 605 i=0,mz
      do 606 j=0,mr
         phi(j,i) = -1000.d-01
         rne(j,i) = 0.0
         rnp(j,i) = 0.0
         rnm(j,i) = 0.0
         rni(j,i) = 0.0
         rjner(j,i) = 0.0
         rjnez(j,i) = 0.0
         rjnpr(j,i) = 0.0
         rjnpz(j,i) = 0.0
         rjnir(j,i) = 0.0
         rjniz(j,i) = 0.0
         dne(j,i) = 0.0
         dnp(j,i) = 0.0
         dni(j,i) = 0.0
         dnm(j,i) = 0.0
         eabs(j,i) = 0.0
         rreabs(j,i) = 0.0
         rrne(j,i) = 0.0
         rrnp(j,i) = 0.0
         rrni(j,i) = 0.0
         rrnm(j,i) = 0.0
  606 continue
  605 continue
      do 490 i=0,mz
      do 491 j=0,mr
          ica(j,i)=0.0
          ae(j,i)=0.0
          sig(j,i)=0.
          ez(j,i)=0.0
          er(j,i)=0.0
  491 continue
  490 continue
      do 494 i=1,mz
          a(i)=0.
          b(i)=0.
          c(i)=0.
          d(i)=0.
  494 continue
      do 2 j=0,mr
          dr(j)=0.
          f(j)=0.
          r(j)=0.
  2    continue
      do 3 i=0,mz
          dz(i)=0.
          g(i)=0.
          nrca(i)=0.
          iran(i)=0
          nran(i)=0
          irca(i)=0
          z(i)=0.
  3   continue
        end

c-------------------------------------------------
c end of subroutine initialise
c-------------------------------------------------

c---------------------------------------------------------------------
c     subroutine amf
c
c     this module is used to calculate the values of sigma to represent solid fuselage or cylinder.
c---------------------------------------------------------------------
      subroutine amf
      implicit double precision (a-h,o-z)
      parameter(mr=150, mz=220, nsh=21)
      common/dummy/
     r  resphi,relphi,rlxphi,ie,je,iph,jph,ifread,mprint
      common/elc/
     r  sig(0:mr,0:mz),sigr(0:mr,0:mz),sigz(0:mr,0:mz)
      common/var/
     r  phi(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rni(0:mr,0:mz),rne(0:mr,0:mz),rnp(0:mr,0:mz),
     r  rrne(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),rnm(0:mr,0:mz),
     r  rrnet(0:mr,0:mz),dne(0:mr,0:mz),dni(0:mr,0:mz),dnm(0:mr,0:mz),
     r  dnp(0:mr,0:mz),rnet(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),rlognetp(0:mr,0:mz),
     r  rlognetn(0:mr,0:mz),rlognm(0:mr,0:mz),rrnm(0:mr,0:mz) ,
     r  rphi(0:mr,0:mz)
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz)
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),rr(0:mr),drr(0:mr),
     r  an(0:mr), ae(0:mr,0:mz),f(0:mr),g(0:mz),vol(0:mr,0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,nri,nro,pinl,
     i  nzc, nrc, nnbot, nra,nnrad,ngapup,radius,volts,
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz)

c---cathode electrical conductivity
c   make sig high in metal when ica = 1 and very small elsewhere
c   make sigr and sigz very small at metal boundaries
      do 303 i=0,nzp
      do 304 j=0,nrp
              sig(j,i) = 1.0
c	  else
  304  continue
  303  continue
      do 306 i=0,nzp
      do 307 j=0,nrp
       sigr(j,i) = 1.0
       sigz(j,i) = 1.0
  307  continue
  306  continue
c
      return
      end

c  end subroutine amf
c---------------------------------------------------------------------
c     subroutine poiss(volts)
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
      subroutine poiss
      implicit double precision (a-h,o-z)
      parameter(mr=150, mz=220, nsh=21)
      dimension zdc(0:mr,0:mz),zdf(0:mr,0:mz)
      common/var/
     r  phi(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rni(0:mr,0:mz),rne(0:mr,0:mz),rnp(0:mr,0:mz),
     r  rrne(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),rnm(0:mr,0:mz),
     r  rrnet(0:mr,0:mz),dne(0:mr,0:mz),dni(0:mr,0:mz),dnm(0:mr,0:mz),
     r  dnp(0:mr,0:mz),rnet(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),rlognetp(0:mr,0:mz),
     r  rlognetn(0:mr,0:mz),rlognm(0:mr,0:mz),rrnm(0:mr,0:mz),
     r  rphi(0:mr,0:mz)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),df(mr,mz)
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),rr(0:mr),drr(0:mr),
     r  an(0:mr), ae(0:mr,0:mz),f(0:mr),g(0:mz),vol(0:mr,0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,nri,nro,pinl,
     i  nzc, nrc, nnbot, nra,nnrad,ngapup,radius,volts,
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz)
      common/elc/
     r  sig(0:mr,0:mz),sigr(0:mr,0:mz),sigz(0:mr,0:mz)
      common/cds/
     r  ez(0:mr,0:mz),er(0:mr,0:mz),eabs(0:mr,0:mz),rreabs(0:mr,0:mz)
      common/dummy/
     r  resphi,relphi,rlxphi,ie,je,iph,jph,ifread,mprint
      common/dummy2/m,mm2

c---relaxation
      epsilon = 8.85d-14
      call fill
      do 15 i=1,nz
      do 16 j=1,nr
         zdc(j,i) =dc(j,i)
         zdf(j,i)=df(j,i)
         dc(j,i) = dc(j,i)/rlxphi
         df(j,i) = df(j,i) + (1.d0-rlxphi)*dc(j,i)*phi(j,i)
  16  continue
  15  continue
c---solve matrix
c--set boundaries for voltage
      do 20 j=0,nrp
            phi(j,nzp) = phi(j,nz)
            phi(j,0) = phi(j,1)
  20  continue
      do 30 i=0,nz
         phi(nrp,i) = phi(nr,i)
         phi(0,i)=phi(1,i)
  30  continue
 572  format(1p9e11.3)
 510  format('  i =',i4,'   z =',f11.4)
c      do 776 i=0,nzp
c         write(*,510) i,z(i)
c         write(*,572) (phi(j,i),j=0,nrp)
c 776  continue
  17      call solve(1,nr,1,nz,de,dw,dn,ds,dc,df,phi,1,1)
      do 775 i=0,nzp
            phi(nrp,i) = phi(nr,i)
            phi(0,i) = phi(1,i)
 775  continue
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
c	  if (m.ne.mm2) go to 31
c	  if (i.eq.85.and.j.eq.3) write(*,*) east, west,south,rnorth
c     r	      centre
  31      ph = dabs(phi(j,i))
          res = dabs(res)
          rel = res/dmax1(east,west,south,rnorth,centre,ph)
          ares = abs(res)
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
      do 166 j=0,nr
c        if (j>=nri.and.j<=nro) then
c        er(j,i) = 0.d0
c        ez(j,i) = 0.d0
c        else
        er(j,i)=(phi(j,i)-phi(j+1,i))*2/(dr(j)+dr(j+1))
        ez(j,i)=(phi(j,i)-phi(j,i+1))*2/(dz(i)+dz(i+1))
c        endif
 166  continue
 165  continue
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
c     of the electric potential used by the subroutine poiss.
c----------------------------------------------------------------------
      subroutine fill
      implicit double precision (a-h,o-z)
      parameter(mr=150, mz=220, nsh=21)
      common/gcoef/
     r  de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),df(mr,mz)
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),rr(0:mr),drr(0:mr),
     r  an(0:mr), ae(0:mr,0:mz),f(0:mr),g(0:mz),vol(0:mr,0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,nri,nro,pinl,
     i  nzc, nrc, nnbot, nra,nnrad,ngapup,radius,volts,
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz)
      common/elc/
     r  sig(0:mr,0:mz),sigr(0:mr,0:mz),sigz(0:mr,0:mz)
      common/var/
     r  phi(0:mr,0:mz),rjniz(0:mr,0:mz),rjnez(0:mr,0:mz),
     r  rjnpz(0:mr,0:mz),rjnir(0:mr,0:mz),rjner(0:mr,0:mz),
     r  rjnpr(0:mr,0:mz),rni(0:mr,0:mz),rne(0:mr,0:mz),rnp(0:mr,0:mz),
     r  rrne(0:mr,0:mz),rrni(0:mr,0:mz),rrnp(0:mr,0:mz),rnm(0:mr,0:mz),
     r  rrnet(0:mr,0:mz),dne(0:mr,0:mz),dni(0:mr,0:mz),dnm(0:mr,0:mz),
     r  dnp(0:mr,0:mz),rnet(0:mr,0:mz),rlogne(0:mr,0:mz),
     r  rlogni(0:mr,0:mz),rlognp(0:mr,0:mz),rlognetp(0:mr,0:mz),
     r  rlognetn(0:mr,0:mz),rlognm(0:mr,0:mz),rrnm(0:mr,0:mz),
     r  rphi(0:mr,0:mz)
c----------------------------------------------------------------------
c     setting up boundary values
c      do 50 j=1,nr
c           call zerocoef(dc,df,de,dw,dn,ds,volts,j,1)
c 50   continue
      e = 1.60206d-19
      epsilon = 8.85d-14
      mu = 10
c mobility for argon electrons taken as mu=10
 510  format('  i =',i4,'   z =',f11.4)
 572  format(1p9e11.3)
      do 5 i=1,nz
      do 6 j=1,nr
	 dc(j,i) = 0.d0
	 de(j,i) = 0.d0
         dw(j,i) = 0.d0
         dn(j,i) = 0.d0
         ds(j,i) = 0.d0
	 df(j,i)=  0.d0
  6   continue
  5   continue
c      do 10 i=1,nz
c      do 11 j=1,nr
c j = sig E = ne e W = ne e W/E E; sig = ne e mu
c sig = (rne(j,i)+rne(j+1,i))*e*10/2
c         dife = 2.*ae(j,i)*sigr(j,i)/(dr(j)+dr(j+1))
c         dife = ae(j,i)*(rne(j+1,i)+rne(j,i))*e*mu/(dr(j)+dr(j+1))
c         de(j,i) = dife
c         dw(j+1,i) = dife
c         difn = an(j)*(rne(j,i+1)+rne(j,i))*e*mu/(dz(i)+dz(i+1))
c         dn(j,i) = difn
c         ds(j,i+1) = difn
c         temp = 2.*ae(j,i)*sigr(j,i)/(dr(j)+dr(j+1))
c         temp = 2.*ae(j,i)/(dr(j)+dr(j+1))
c         if (j .ne. 0) de(j,i) = temp
c         dw(j+1,i) = temp
c  11  continue
c  10  continue
c electric field from space charge - valid for Poisson equation, not current continuity
      do 20 i=1,nz
      do 21 j=1,nr
	 df(j,i)=(e/epsilon)*(rnp(j,i)-rne(j,i)-rni(j,i))
         dn(j,i) = 2/(dz(i)*(dz(i)+dz(i)))
         ds(j,i+1) = dn(j,i)
	 de(j,i) = 2/(dr(j)*(dr(j)+dr(j)))
         dw(j,i+1) = de(j,i)
  21  continue
  20  continue
      do 80 i=1,nz
      do 81 j=1,nr
         dc(j,i) = de(j,i) + dw(j,i) + dn(j,i) + ds(j,i)
   81 continue
   80 continue
c votages applied within pin:
c            PI = 3.1414
c         omg = 2* PI* freq
c--phi =volts for i up to ngapup within pin.
      do 91 i = 1, ngapup
      do 93 j = nri, nro
         ds(j,i) = 0.d0
         dn(j,i) = 0.d0
         de(j,i) = 0.0
         dw(j,i) = 0.0
         dc(j,i) = 1.0
         df(j,i) = volts
   93 continue
   91 continue
c      write(*,*) nri,nro
c-Vertical boundary conditions:
c--dphi/dr =0. at r= 0.
      do 90 i=1,nz
         dc(1,i) = dc(1,i) - dw(1,i)
         dw(1,i) = 0.d0
   90 continue
c--phi=0 at outer radius
      do 92 i=1,nz
         dc(nr,i) = dc(nr,i) - de(nr,i)
         de(nr,i) = 0.d0
         ds(nr,i) = 0.d0
         dn(nr,i) = 0.d0
         de(nr,i) = 0.0
         dw(nr,i) = 0.0
         dc(nr,i) = 1.0
         df(nr,i) = 0.0
   92 continue
c Horizontal boundary conditions
c         deltaV = (Eedge/(2*radius)*(radius**2 - r(j)**2))
c upper plane at z = 0, at voltage of volts to outer radius of pin nro.
c         dc(j,1) = dc(j,1) - dn(j,1)
c         dn(j,1) = 0.0
c      do 100 j=1,nro
c         ds(j,1) = 0.d0
c         dn(j,1) = 0.d0
c         de(j,1) = 0.0
c         dw(j,1) = 0.0
c         dc(j,1) = 1.0
c         df(j,1) = 0.0
c  100 continue
c upper plane at z = 0 from outer pin radius nro to wall nr, d phi/dz = 0
c         deltaV = (Eedge/(2*radius)*(radius**2 - r(j)**2))
      do 101 j=nro,nr
         dc(j,1) = dc(j,1) - dn(j,1)
         dn(j,1) = 0.0
         ds(j,1) = 0.d0
         dn(j,1) = 0.d0
         de(j,1) = 0.0
         dw(j,1) = 0.0
         dc(j,1) = 1.0
c         df(j,1) = 0.0
         df(j,1) = volts * (radius-r(j))/(radius-r(nro))
  101 continue
c lower plane at i = nz, phi set at ground voltage 0.0
      do 102 j = 1,nr
         ds(j,nz) = 0.d0
         dn(j,nz) = 0.d0
         de(j,nz) = 0.0
         dw(j,nz) = 0.0
         dc(j,nz) = 1.0
         df(j,nz) = 0.0
  102 continue

      return
      end
c end of subroutine fill

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
      parameter(mr=150, mz=220, nsh=21)
      dimension de(mr,mz),dw(mr,mz),dn(mr,mz),ds(mr,mz),dc(mr,mz),
     1   df(mr,mz),phi(0:mr,0:mz),v(mz),x(0:mr,0:mz)
      common/dummy/
     r  resphi,relphi,rlxphi,ie,je,iph,jph,ifread,mprint
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz)

      integer flag1,flag2
c------------------------------------------------------------------
      nr1 = nr-1
      nz1 = nz-1
      do 8 i=0,mz
      do 9 j=0,mr
         x(j,i) = 0.d0
  9   continue
  8   continue
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
      do 200 iter=1,300
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
c         if (flag1.eq.1) write(*,999) res,res0
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
         iiph = 0
         jjph = 0
         relphi = 0.0
         do 260 i=n,nz
         do 261 j=m,nr
            if (x(j,i) .ne. 0.) then
               diff = (x(j,i)-phi(j,i))/x(j,i)
            else
               diff = -phi(j,i)
            end if
	    adiff = abs(diff)
        if (adiff .gt.relphi) then
            relphi = adiff
            iiph = i
            jjph = j
        endif
            if (dabs(diff) .gt. 1.e-3) go to 265
  261    continue
  260    continue
         iconv = 1
c       if (flag1.eq.1) write(*,998) iter,relphi,iiph,jjph,
c     r     x(jjph,iiph),phi(jjph,iiph)
  265    continue
  998    format(1x,'iter =',i4,'    relphi =',1p,e16.8,
     1            '  at   (',i5,',',i5,')',2e12.3)
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
c       if (flag1.eq.1) write(*,998) iter,diff,i,j,x(j,i),phi(j,i)
 300  format(1x,'convergence not attained in subroutine solve')
      return
      end

c------------------------------------------------------------------
c end of subroutine solve
c-------------------------------------------------------------------

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
      parameter(mr=150, mz=220, nsh=21)
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
c       subroutine mesh
c       v.2 ls 16.04.99
c
c mesh calculates the mesh from the data of the file gridbern6.
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
c     gap is length non-uniform z mesh between cloud/ground and streamer tips
c
c-----------------------------------------------------------------------
	subroutine mesh

      implicit double precision (a-h,o-z)
      parameter(mr=150, mz=220, nsh=21)
      dimension  zb(0:mz),rb(0:mr)
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),rr(0:mr),drr(0:mr),
     r  an(0:mr), ae(0:mr,0:mz),f(0:mr),g(0:mz),vol(0:mr,0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,nri,nro,pinl,
     i  nzc, nrc, nnbot, nra,nnrad,ngapup,radius,volts,
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz)
c------------------------------------------------------------------------
c------------------------------------------------------------------------
c read the data for the mesh generation in file gridbern6

        open(18,file='gridbern')
        rewind 18
	read(18,*)
c radius, distance below pin, and pin length; cm
        read(18,*) radius,hight,pinl
	read(18,*)
c nri and nro are the j values for inner and outer radius of pin; nzl is i value for pin length
        read(18,*) nri, nro
        read(18,*)
c ngap up/down is number of points for nonuniform mesh axially from pin base; ngapradius is radial points
	read(18,*) ngapup,ngapdown,ngapradius
        read(18,*)
c shup and shdown are first grid  sizes for z at pin tip and shrad for r axes at pin tip
        read(18,*) shup, shdown, shrad
	close(18)
c nz is total axial points; must not exceed mz
       nz = ngapup + ngapdown
       nzp = nz+1
       nzm = nz-1
       open(20,file='mesh')
 12    format(' nz = ',i4,' ngap = ',i4,'  nnsh = ',i4,'  nnrad',i4,
     r    '  gap = ',E11.3)
 11    format(' nr = ',i4,' nzc = ',i4,'  nnbot = ',i4,' nstream = ',i4)
       if (nz.ge.mz) then
           write(*,*) 'nz =',nz,'> mz'
	   stop
       endif
c nr is total radial points; must not exceed mr
       nr = ngapradius
       nrp = nr+1
       nrm = nr-1
       if (nr.ge.mr) then
           write(*,*) 'nr =',nr,'> mr'
	   stop
       endif
c--------------------------------------------------------------------
c       initialisation
	do 10 j = 0,mr
	   rb(j)=0.0
 10	continue
	do 20 i = 0,mz
	   zb(i)=0.0
 20	continue
c---------------------------------------------
	write(20,*) 'mesh1'
c----------------------------------------------------------------
c axial mesh calculation - calculation of the axial cell boundaries
c---------------------------------------------------------------
c Gap from the cylinder top and bottom to streamers; nnsh points of nonuniform grid
c    shca for the first cell's thickness, dsh total sheath distance
	write(20,*) 'mesh3'
c	nn = nnsh
c	zz = dsh
c	db = shca
c	call multcoef(zz,db,nn,a)
c	do 47 i= 1, nnsh
c	   dz(i) = db*(a**float(i-1))
c	   zb(nzc+i)=zb(nzc-1+i)+dz(i)
c	   zb(nnbot-i-1) = zb(nnbot-i)-dz(i)
c 47	continue
c------------------------------------------------------
c Uniform grid of streamer
c	write(20,*) 'mesh4'
c	nn = nstream
c	zz = zb(nnbot+1)
c	do 43 i= 1, nn
c	   dz(i) = gstream
c	   zb(nzc+nnsh+i)=zb(nzc+nnsh-1+i)+dz(i)
c	   zb(nnbot-nnsh-i-1)=zb(nnbot-nnsh-i)-dz(i)
c 43	continue
c------------------------------------------------------
c Non uniform mesh from pin tip bottom to  z = 0.
 	write(20,*) 'mesh5'
	nn = ngapup
        z(ngapup) = pinl
        zz = z(ngapup)
	db = shup
	call multcoef(zz,db,nn,a)
        zb(nn) = zz
	do 48 i= 1, ngapup
	   dz(nn-i+1) = db*(a**float(i-1))
	   zb(nn-i) = zb(nn-i+1)-dz(nn-i+1)
 48	continue
c------------------------------------------------------
c Non uniform mesh from pin tip to ground, z = z(nz).
	write(20,*) 'mesh6'
	nn = ngapdown
        zz = hight
	db = shdown
	call multcoef(zz,db,nn,a)
	do 45 i= 1,ngapdown
	   dz(ngapup+i) = db*(a**float(i-1))
	   zb(ngapup+i)=zb(ngapup-1+i)+dz(ngapup+i)
 45	continue
c axial calculation of the cell's thickness dz and
c of the coordinate of the middle point z of ech cell
        zb(0) = 0.0
	z(0) = 0.0
	dz(0) = 0.0
	do 46 i = 1,nz
	   dz(i)=zb(i)-zb(i-1)
	   z(i) = zb(i-1)+0.5*dz(i)
 46	continue
	dz(nzp) = dz(nz)
	z(nzp) = z(nz)+dz(nzp)
	zb(nzp) = zb(nz)+dz(nzp)
c-----------------------------------------------
c end of the axial mesh calculation
c----------------------------------------------

c---------------------------------------------
c radial mesh calculation; calculation of radial cell boundaries
c---------------------------------------------------------
c radial mesh, cylinderto the wall; non-uniform grid; nrgap points
c and the thickness of the first cell equal to the thickness
c of the last cell in the cathode.
      	write(20,*) 'mesh7'
	nn = nr
	zz = radius
	db= shrad
	call multcoef(zz,db,nn,a)
	do 32 j= 1, nrp
	   dr(j) = db*(a**float(j-1))
	   rb(j) = rb(j-1) + dr(j)
 32	continue
c---------------------------------------------------------------------
c radial calculation of cell's thickness dr and mid points of ech cell
        rb(0) = 0.0
	r(0) = 0.0
	dr(0) = 0.0
        r(1) = dr(1)/2
	do 33 j= 2,nrp
	   dr(j) = rb(j)-rb(j-1)
           r(j) = rb(j-1)+0.5*dr(j)
 33	continue
	dr(nrp) = dr(nr)
	r(nrp) = r(nr)+dr(nrp)
	rb(nrp) = rb(nr)+dr(nrp)
c-------------------------------------------------
c end of the radial mesh calculation
c-------------------------------------------------
c calculate the volume and surfaces of each elementary cell
	call volume
c define the position of each cell,
c	call cellpos

c----------------------------------------------------------
c save mesh data in file mesh
	write(20,*) 'j, r, rb, dr'
	do 779 j=0,nrp
	   write(20,*) j, r(j), rb(j), dr(j)
 779	continue
	write(20,*) 'i, z, zb, dz'
	do 378 i=0,nzp
	   write(20,*) i, z(i), zb(i), dz(i)
 378	continue
	close(20)
c-----------------------------------------------------------
        return
      end

c-------------------------------------------------------
c end subroutine mesh
c-------------------------------------------------------
c-----------------------------------------------------
c     subroutine volume
c
c calculates the volume and the surfaces of each elementary
c cell. calculates also the ponderation factor used for the
c geometric average calculation.
c an top surface of the cell
c ae external radial surface of the cell
c vol volume of the cell
c----------------------------------------------------

      subroutine volume
      implicit double precision (a-h,o-z)
      parameter(mr=150, mz=220, nsh=21)
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),rr(0:mr),drr(0:mr),
     r  an(0:mr), ae(0:mr,0:mz),f(0:mr),g(0:mz),vol(0:mr,0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,nri,nro,pinl,
     i  nzc, nrc, nnbot, nra,nnrad,ngapup,radius,volts,
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz)
      common/tri/
     r    a(mz),b(mz),c(mz),d(mz)
      an(0) = 0.d0
      an(1) = 0.5*dr(1)*dr(1)
      do 30 j=2,nr
         xa = r(j)*dr(j)
         an(j) = xa
  30  continue
      do 50 i=0,nz
      do 51 j=0,nr
         xa = (r(j)+0.5*dr(j))*dz(i)
         ae(j,i) = xa
  51  continue
  50  continue
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

c     subroutine multcoef
c     v.1 ls 15.9.98
c--------------------------------------------------------

	subroutine multcoef(z,fst,n,a)

      implicit double precision (a-h,o-z)
	alp=abs(z/fst)

	if(abs(z/fst-float(n)).lt.1.0e-4)then
	a = fst
	return
	endif
        if(alp.gt.n) then
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
	return
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

c------------------------------------------------------
c end subroutine volume
c------------------------------------------------------

c------------------------------------------------------
c     subroutine cellpos
c     v.1 ls 21.9.98
c
c assigns number to cells in the calculation domain.
c ica=1 for cathode =0 elsewhere
c-------------------------------------------------------

        subroutine cellpos

        implicit double precision (a-h,o-z)
      parameter(mr=150, mz=220, nsh=21)
      common/coord/
     r  dr(0:mr), dz(0:mz), r(0:mr), z(0:mz),rr(0:mr),drr(0:mr),
     r  an(0:mr), ae(0:mr,0:mz),f(0:mr),g(0:mz),vol(0:mr,0:mz),
     i  nr, nrp, nrm, nz, nzp, nzm,nri,nro,pinl,
     i  nzc, nrc, nnbot, nra,nnrad,ngapup,radius,volts,
     i  irca(0:mz), nrca(0:mz), ica(0:mr,0:mz),
     i  iran(0:mz), nran(0:mz)
        do 9 i=0,mz
        do 10 j=0,mr
	ica(j,i)=0
 10      continue
 9       continue
c--------------------------------------------------------------
c  ica=1 inside the cathode and 0 elsewhere
c        do 121 i=0,nzc
	do 121 i=nnbot,nzc
	if(nrca(i).le.irca(i))go to 121
	ist=0
	ien=nrca(i)
	if(ien.eq.nr)ien=nrp
	if(irca(i).gt.1)ist=irca(i)
 121	continue

        do 122 i=nnbot,nzp
	if(nran(i).le.iran(i))go to 122
	ist=0
	ien=nran(i)
	if(ien.eq.nr)ien=nrp
	if(iran(i).gt.1)ist=iran(i)
 122	continue
c----------------------------------------------------------------
	return
        end

c-----------------------------------------------------------------
c end subroutine cellpos
c-----------------------------------------------------------------

c-------------------------------------------------------------------
c     subroutine readarcin
c     v.1 ls 16.9.98
c-------------------------------------------------------------------
      subroutine readarcin(volts,rlxphi,mm1,ifread,delt,
     r    etime,mm2,mprint,freq,meta)
      implicit double precision (a-h,o-z)
      logical test
      parameter(mr=150, mz=220, nsh=21)
      common/plot/iplot(50),jplot(50),iplotn,jplotn,ib,jend,iend
c rc=0.1 is source rel limit
      open(3,file='arcin')
      rewind 3
      read(3,*)
      read(3,*) volts,delt,freq,meta
      read(3,*)
      read(3,*) rc,rlxphi
      read(3,*)
      read(3,*) mm1, ifread, mm2,mprint
c      read(3,*) ib,jend,iplotn,jplotn,iend,jscreen
c ib, jend first i and endj for contour plot;iplot,jplotn arrow numbers
c      read(3,*) test
c if test true reads elements for vector plot;
c      if (test) then
c          read(3,*) (iplot(i), i=1,iplotn)
c          read(3,*) (jplot(j), j=1,jplotn)
c      endif
        close(3)
	end

c---------------------------------------------------------------------
c end of subroutine readarcin
c---------------------------------------------------------------------

c----------------------------------------------------------------
c     subroutine meshrec
c     v.1 ls 21.9.98
c------------------------------------------------------------------
