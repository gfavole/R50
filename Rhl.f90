  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! This code takes the output table of SAG_inputs_radii.py
  ! and computes the half-light radii
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  implicit   none
  CHARACTER(450), PARAMETER   :: cat='/Users/gfavole/Documents/work/test_radii/SAG_0.9_basecat.txt'
  CHARACTER(450), PARAMETER   :: outall='../outputs_last/SAG_0.9_outputs.txt'
  real*8, parameter  :: pi=3.141592653589793
  real       r,r0,I0,re,Ie,logm,hsb
  common     /disk/r0,I0
  common     /bulge/re,Ie
  real       i_d,i_b,l_d,l_b,t_l_d,t_l_b,Qp,Rp,Fp,Rf
  external   i_d,i_b,l_d,l_b,t_l_d,t_l_b,Qp,Rp,Fp,Rf
  real       mabs_t,mabs_d,mabs_b,r50,r90,c,fb,mu0,mue
  !........................................................................................


REAL*8, dimension(1:2e6) :: rad0,rade ,int0,inte,logmass,rhalfmassarcsec,magr,x,y,z,logMh,Vpeak,haloID,sfr
real*8   :: temp_rad0,temp_rade,temp_int0,temp_inte,temp_logmass,temp_rhalfmassarcsec,temp_magr,temp_x,temp_y,temp_z,temp_logMh,temp_Vpeak,temp_haloID,temp_sfr

integer   :: io,i, Ndata,p

OPEN(1,FILE=cat,STATUS='OLD',action='read')
p =1
DO WHILE (io >= 0)   
   READ(1,*,IOSTAT=io) temp_rad0,temp_int0,temp_rade,temp_inte,temp_rhalfmassarcsec,temp_magr,temp_x,temp_y,temp_z,temp_logMh,temp_Vpeak,temp_haloID,temp_sfr
   rad0(p) = temp_rad0 !in kpc
   int0(p) = temp_int0
   rade(p) = temp_rade !in kpc
   inte(p) = temp_inte
   rhalfmassarcsec(p) = temp_rhalfmassarcsec
   magr(p) = temp_magr
   x(p)=temp_x
   y(p)=temp_y
   z(p)=temp_z
   logMh(p)=temp_logMh
   Vpeak(p)=temp_Vpeak
   haloID(p)=temp_haloID
   sfr(p)=temp_sfr
   p = p+1
ENDDO
Ndata = p-2
write(*,*) Ndata
close(1)



open(13,file=outall)
do i=1,Ndata
   r0=rad0(i)
   re=rade(i)
   I0=10.**(-int0(i))
   Ie=10.**(-inte(i))

   hsb=magr(i)+2.5*log10(2.*pi*Rf(r0,I0,re,Ie,0.5)*Rf(r0,I0,re,Ie,0.5))
   write(*,*) i
     if (hsb .LE. 24.5) then
        write(13,'(f15.5,x,f15.5,x,f15.5,x,f15.5,x,f15.5,x,f15.5,x,f15.5,x,i13)') x(i), &
        y(i), z(i), logMh(i), Vpeak(i), Rf(r0,I0,re,Ie,0.5), sfr(i), haloID(i)
      endif
   endif
enddo
close(13)

end program








!****************************************************************************************
!       Calcula el radio que engloba una fraccion f del flujo de Petrosian
real     function Rf(r0,I0,re,Ie,f)
  implicit none
  real     r0,I0,re,Ie,f,flux
  real     Fp,l_d,l_b
  external Fp,l_d,l_b
  real     rmax,rmin,r1,r2,r3,f2
  flux=f*Fp(r0,I0,re,Ie)
  flux=-2.5*log10(flux)
  rmax=1000.  !1 Mpc in kpc
  rmin=0.01   !10 pc in kpc
  r1=rmin
  r3=rmax
  f2=1000.
  do while (abs(f2).gt.1.e-4)
20   r2=0.5*(r1+r3)
     f2=-2.5*log10(l_d(r0,I0,r2)+l_b(re,Ie,r2))-flux
     if(abs(f2).lt.1.e-4) go to 10
     if(f2.gt.0.) then
        r1=r2
        go to 20
     else
        r3=r2
        go to 20
     end if
  end do
10 Rf=r2
  return
end function Rf
!****************************************************************************************
!       Calcula el flujo de Petrosian, i.e., el flujo dentro de 2 radios
!       de Petrosian de acuardo a SDSS
real     function Fp(r0,I0,re,Ie)
  implicit none
  real     r0,I0,re,Ie,r_pet2
  real     l_d,l_b,Rp
  external l_d,l_b,Rp
  r_pet2=2.*Rp(r0,I0,re,Ie)
  Fp=l_d(r0,I0,r_pet2)+l_b(re,Ie,r_pet2)
  return
end function Fp
!****************************************************************************************
!       Calcula el radio de Petrosian, el cual, SDSS definen como el
!       radio para el cual el cociente de Petrosian es 0.2
real     function Rp(r0,I0,re,Ie)
  implicit none
  real     r0,I0,re,Ie,r
  real     Qp,rmax,rmin,r1,r2,r3,f2
  external Qp
  rmax=1000.  !1 Mpc
  rmin=0.01   !10 pc
  r1=rmin
  r3=rmax
  f2=1000.
  do while (abs(f2).gt.1.e-4)
20   r2=0.5*(r1+r3)
     f2=Qp(r0,I0,re,Ie,r2)-0.2
     if(abs(f2).lt.1.e-4) go to 10
     if(f2.gt.0.) then
        r1=r2
        go to 20
     else
        r3=r2
        go to 20
     end if
  end do
10 Rp=r2
  return
end function Rp
!****************************************************************************************
!       Calcula el cociente de Petrosian a un radio r, definido como el cociente
!       entre el brillo superficial medio en un anillo comprendido entre
!       0.8r y 1.25r y el brillo superficial medio hasta r
real     function Qp(r0,I0,re,Ie,r)
  implicit none
  real     r0,I0,re,Ie,r,up,down
  real     l_d,l_b
  external l_d,l_b
  up=l_d(r0,I0,1.25*r)+l_b(re,Ie,1.25*r)-l_d(r0,I0,0.80*r)-l_b(re,Ie,0.80*r)
  down=l_d(r0,I0,r)+l_b(re,Ie,r)
  Qp=up/down/0.92250
  return
end function Qp
!****************************************************************************************
!       Perfil de intensidad de un Disco Exponencial
!       con radio de escala r0 e intensidad central I0
real     function i_d(r0,I0,r)
  implicit none
  real     r,r0,I0
  i_d=I0*exp(-r/r0)
  return
end function i_d
!****************************************************************************************
!       Perfil de intensidad de un Bulge que sigue la ley r^1/4
!       con radio de escala re e intensidad central Ie
real     function i_b(re,Ie,r)
  implicit none
  real     r,re,Ie
  i_b=Ie*exp(-7.688*((r/re)**0.25-1.))
  return
end function i_b
!****************************************************************************************
!       Luminosidad integrada de un Disco Exponencial dentro de un radio r
!       con radio de escala r0 e intensidad central I0
real     function l_d(r0,I0,r)
  implicit none
  real     r,r0,I0,x
  x=r/r0
  l_d=8.*atan(1.)*I0*r0**2*(1.-(1.+x)*exp(-x))
  return
end function l_d
!****************************************************************************************
!       Luminosidad integrada de un Bulge r^1/4 dentro de un radio r
!       con radio de escala re e intensidad central Ie
real     function l_b(re,Ie,r)
  implicit none
  real     r,re,Ie,x,poly,b
  b=7.688
  x=r/re
  poly=4./b*x**1.75+28./b**2*x**1.5+168./b**3*x**1.25+840./b**4*x+ &
  3360./b**5*x**0.75+10080./b**6*x**0.5+20160./b**7*x**0.25+20160/b**8
  l_b=8.*atan(1.)*Ie*re**2*exp(b)*(20160./b**8-exp(-b*x**0.25)*poly)
  return
end function l_b
!****************************************************************************************
!       Luminosidad total de un Disco Exponencial
!       con radio de escala r0 e intensidad central I0
real     function t_l_d(r0,I0)
  implicit none
  real     r0,I0
  t_l_d=8.*atan(1.)*I0*r0**2
  return
end function t_l_d
!****************************************************************************************
!       Luminosidad total de un Bulge r^1/4
!       con radio de escala re e intensidad central Ie
real     function t_l_b(re,Ie)
  implicit none
  real     re,Ie
  t_l_b=28.836*atan(1.)*Ie*re**2
  return
end function t_l_b
!****************************************************************************************
