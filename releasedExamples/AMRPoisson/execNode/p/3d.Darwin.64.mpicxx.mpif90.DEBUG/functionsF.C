#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      real_t  function getphirzfunc(radius)
      implicit none
      real_t radius
      getphirzfunc = radius*radius
      return
      end
      real_t  function getgradphirzfunc(radius)
      implicit none
      real_t radius
      getgradphirzfunc = two*radius
      return
      end
      real_t  function getlaplphirzfunc(radius)
      implicit none
      real_t radius
      getlaplphirzfunc = four
      return
      end
        subroutine GETPHI(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2)
      REAL_T freq(0:2)
      REAL_T dx(0:2)
      REAL_T problo(0:2)
      REAL_T probhi(0:2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
        integer i,j,k
        real_t x(0:CH_SPACEDIM-1)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          x(2) = (k+half)*dx(2) + problo(2)
          call getphipoint(phi(i,j,k),freq,x)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine GETMAGRESIST(
     &           mag
     &           ,imaglo0,imaglo1,imaglo2
     &           ,imaghi0,imaghi1,imaghi2
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,icomp
     &           ,whichmag
     &           )

      implicit none
      integer imaglo0,imaglo1,imaglo2
      integer imaghi0,imaghi1,imaghi2
      REAL_T mag(
     &           imaglo0:imaghi0,
     &           imaglo1:imaghi1,
     &           imaglo2:imaghi2)
      REAL_T freq(0:2)
      REAL_T dx(0:2)
      REAL_T problo(0:2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer icomp
      integer whichmag
        integer i,j,k
        real_t x(0:CH_SPACEDIM-1)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          x(2) = (k+half)*dx(2) + problo(2)
          call getmagpointresist(mag(i,j,k),freq,x,
     $         icomp, whichmag)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine GETMAGPOINTRESIST(
     &           mag
     &           ,freq
     &           ,xval
     &           ,icomp
     &           ,whichmag
     &           )

      implicit none
      REAL_T mag
      REAL_T freq(0:2)
      REAL_T xval(0:2)
      integer icomp
      integer whichmag
        REAL_T  x,y,z
        integer i,j,k
#if CH_SPACEDIM==2
        if(icomp.eq.2) then
           mag = zero
           return
        endif
#endif
        
        i = icomp
        j = max(1-icomp, 0)
        k = 3-icomp-j
        if(whichmag.eq. 2) then
           
           x = freq(i)*xval(i)
           y = freq(j)*xval(j)
           z = freq(k)*xval(k)
#if CH_SPACEDIM==2
           mag = sin(y)
#elif CH_SPACEDIM==3
           mag = sin(y) + sin(z)
#else
           mag = x
#endif
        elseif(whichmag.eq. 1) then
           x = freq(icomp)*xval(icomp)
           mag = sin(x)
        elseif(whichmag.eq.0) then
           x = xval(icomp)
           mag = x*x
        elseif(whichmag.eq.4) then
           x = xval(icomp)
           if(icomp .eq. 0) then
              mag = x
           else
              mag = zero
           endif
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETETARESIST(
     &           eta
     &           ,ietalo0,ietalo1,ietalo2
     &           ,ietahi0,ietahi1,ietahi2
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,idir
     &           ,eps
     &           ,whicheta
     &           )

      implicit none
      integer ietalo0,ietalo1,ietalo2
      integer ietahi0,ietahi1,ietahi2
      REAL_T eta(
     &           ietalo0:ietahi0,
     &           ietalo1:ietahi1,
     &           ietalo2:ietahi2)
      REAL_T freq(0:2)
      REAL_T dx(0:2)
      REAL_T problo(0:2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer idir
      REAL_T eps
      integer whicheta
        integer i,j,k,jdir
        real_t x(0:CH_SPACEDIM-1)
        integer iv(0:CH_SPACEDIM-1)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

           
           iv(0) = i
           iv(1) = j
           iv(2) = k
           do jdir = 0, CH_SPACEDIM-1
              if(idir .eq. jdir) then
                 x(jdir) = iv(jdir)*dx(jdir) + problo(jdir)
              else
                 x(jdir) = (iv(jdir)+half)*dx(jdir) + problo(jdir)
              endif
           enddo
          call getetapointresist(eta(i,j,k),freq,x,
     $         idir, eps, whicheta)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine GETETAPOINTRESIST(
     &           eta
     &           ,freq
     &           ,xval
     &           ,idir
     &           ,eps
     &           ,whicheta
     &           )

      implicit none
      REAL_T eta
      REAL_T freq(0:2)
      REAL_T xval(0:2)
      integer idir
      REAL_T eps
      integer whicheta
        REAL_T x, y, z
        if(whicheta.eq. 1) then
           
           x = freq(0)*xval(0)
           y = freq(1)*xval(1)
           z = freq(2)*xval(2)
           eta  = one + eps*(sin(x) + sin(y) +sin(z))
        elseif(whicheta.eq.0) then
           eta = one
        elseif(whicheta.eq.3) then
           eta = half
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETBETAVISCOUS(
     &           beta
     &           ,ibetalo0,ibetalo1,ibetalo2
     &           ,ibetahi0,ibetahi1,ibetahi2
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,eps
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,whichbeta
     &           )

      implicit none
      integer ibetalo0,ibetalo1,ibetalo2
      integer ibetahi0,ibetahi1,ibetahi2
      REAL_T beta(
     &           ibetalo0:ibetahi0,
     &           ibetalo1:ibetahi1,
     &           ibetalo2:ibetahi2)
      REAL_T freq(0:2)
      REAL_T dx(0:2)
      REAL_T problo(0:2)
      REAL_T eps
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer whichbeta
        integer i,j,k,jdir
        real_t x(0:CH_SPACEDIM-1)
        integer iv(0:CH_SPACEDIM-1)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

           
           iv(0) = i
           iv(1) = j
           iv(2) = k
           do jdir = 0, CH_SPACEDIM-1
              x(jdir) = (iv(jdir)+half)*dx(jdir) + problo(jdir)
           enddo
          call getbetapointviscous(beta(i,j,k),freq,x, eps, whichbeta)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine GETBETAPOINTVISCOUS(
     &           beta
     &           ,freq
     &           ,xval
     &           ,eps
     &           ,whichbeta
     &           )

      implicit none
      REAL_T beta
      REAL_T freq(0:2)
      REAL_T xval(0:2)
      REAL_T eps
      integer whichbeta
        REAL_T x, y, z
        if(whichbeta.eq. 1) then
           
           x = freq(0)*xval(0)
           y = freq(1)*xval(1)
           z = freq(2)*xval(2)
           beta  = one + eps*(sin(x) + sin(y) +sin(z))
        elseif(whichbeta.eq.0) then
           beta = one
        else
           call MayDay_Error()
        endif
        return
        end
        subroutine GETKLBRESIST(
     &           klb
     &           ,iklblo0,iklblo1,iklblo2
     &           ,iklbhi0,iklbhi1,iklbhi2
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,alpha
     &           ,beta
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,icomp
     &           ,eps
     &           ,whichmag
     &           ,whicheta
     &           )

      implicit none
      integer iklblo0,iklblo1,iklblo2
      integer iklbhi0,iklbhi1,iklbhi2
      REAL_T klb(
     &           iklblo0:iklbhi0,
     &           iklblo1:iklbhi1,
     &           iklblo2:iklbhi2)
      REAL_T freq(0:2)
      REAL_T dx(0:2)
      REAL_T problo(0:2)
      REAL_T alpha
      REAL_T beta
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer icomp
      REAL_T eps
      integer whichmag
      integer whicheta
        integer i,j,k
        real_t x(0:CH_SPACEDIM-1)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          x(2) = (k+half)*dx(2) + problo(2)
          call getklbpointresist(klb(i,j,k),freq,x,
     $         alpha, beta, icomp, eps, whichmag, whicheta)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine GETKLBPOINTRESIST(
     &           klb
     &           ,freq
     &           ,xvec
     &           ,alpha
     &           ,beta
     &           ,icomp
     &           ,eps
     &           ,whichmag
     &           ,whicheta
     &           )

      implicit none
      REAL_T klb
      REAL_T freq(0:2)
      REAL_T xvec(0:2)
      REAL_T alpha
      REAL_T beta
      integer icomp
      REAL_T eps
      integer whichmag
      integer whicheta
        REAL_T  fx,fy,fz
        REAL_T  x,y,z,  termone
        REAL_T  freqx,freqy,freqz, mag, divf, eta
        integer i,j,k
#if CH_SPACEDIM==2
        if(icomp.eq.2) then
           klb = zero
           return
        endif
#endif
        
        i = icomp
        j = max(1-icomp, 0)
        k = 3-icomp-j
        call getetapointresist(eta,freq,xvec, icomp, eps, whicheta)
        call getmagpointresist(mag,freq,xvec, icomp, whichmag)
        
        freqx = freq(i)
        freqy = freq(j)
        freqz = freq(k)
        
        x = xvec(i)
        y = xvec(j)
        z = xvec(k)
        if((whichmag.eq. 2).and.(whicheta.eq.0)) then
#if(CH_SPACEDIM==1)
           divf = one
#else
           divf = -(freqy*freqy*sin(freqy*y))
#if CH_SPACEDIM==3
           divf = divf - (freqz*freqz*sin(freqz*z))
#endif
#endif
        elseif((whichmag.eq. 3).and.(whicheta.eq.0)) then
           call getlofphipoint(klb, freq, xvec, alpha, beta)
           goto  123
        elseif((whichmag.eq. 2).and.(whicheta.eq.1)) then
           
           fx = freqx*x
           fy = freqy*y
           fz = freqz*z
#if(CH_SPACEDIM==1)
           divf = one
#else
           divf = (freqy*cos(fy) - freqx*cos(fx))*(eps*freqy*cos(fy)) - freqy*freqy*eta*sin(fy)
#if CH_SPACEDIM==3
           divf = eps*freqy*cos(fy)*(freqy*cos(fy) - freqx*cos(fx))  - eta*freqy*freqy*sin(fy)
     $          + eps*freqz*cos(fz)*(freqz*cos(fz) - freqx*cos(fx))  - eta*freqz*freqz*sin(fz)
#endif
#endif
        elseif((whichmag.eq. 1).and.(whicheta.eq.1)) then
           termone =  
     $          freqx*cos(freqx*x) +
     $          freqy*cos(freqy*y) +
     $          freqz*cos(freqz*z)
          divf = eps*freqx*cos(freqx*x)*termone
     $         -freqx*freqx*sin(freqx*x)*eta
        elseif((whichmag.eq.0).and.(whicheta.eq.0)) then
           divf = two
        elseif((whichmag.eq.4).and.(whicheta.eq.0)) then
           divf = zero
        elseif((whichmag.eq.1).and.(whicheta.eq.0)) then
           divf = -freqx*freqx*sin(freqx*x)
        else
           call MayDay_Error()
        endif
        klb = alpha*mag  + beta*divf
  123   continue
        return
        end
        subroutine GETKLVVISCOUS(
     &           klb
     &           ,iklblo0,iklblo1,iklblo2
     &           ,iklbhi0,iklbhi1,iklbhi2
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,alpha
     &           ,beta
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,icomp
     &           ,eps
     &           ,whichvel
     &           ,whicheta
     &           ,whichlambda
     &           ,lambdafactor
     &           )

      implicit none
      integer iklblo0,iklblo1,iklblo2
      integer iklbhi0,iklbhi1,iklbhi2
      REAL_T klb(
     &           iklblo0:iklbhi0,
     &           iklblo1:iklbhi1,
     &           iklblo2:iklbhi2)
      REAL_T freq(0:2)
      REAL_T dx(0:2)
      REAL_T problo(0:2)
      REAL_T alpha
      REAL_T beta
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer icomp
      REAL_T eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL_T lambdafactor
        integer i,j,k
        real_t x(0:CH_SPACEDIM-1)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          x(2) = (k+half)*dx(2) + problo(2)
          call getklvpointviscous(klb(i,j,k),
     $         freq,x, alpha, beta, icomp, eps,
     $         whichvel, whicheta, whichlambda, lambdafactor)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine GETKLVPOINTVISCOUS(
     &           klv
     &           ,freq
     &           ,xvec
     &           ,alpha
     &           ,beta
     &           ,icomp
     &           ,eps
     &           ,whichvel
     &           ,whicheta
     &           ,whichlambda
     &           ,lambdafactor
     &           )

      implicit none
      REAL_T klv
      REAL_T freq(0:2)
      REAL_T xvec(0:2)
      REAL_T alpha
      REAL_T beta
      integer icomp
      REAL_T eps
      integer whichvel
      integer whicheta
      integer whichlambda
      REAL_T lambdafactor
        REAL_T  x,y,z,    lambda
        REAL_T  freqx,freqy,freqz, vel, divf, eta
        real_t fx,fy,fz
        integer i,j,k
        
        i = icomp
        j = max(1-icomp, 0)
        k = 3-icomp-j
        call getetapointresist(   eta,freq,xvec, icomp, eps, whicheta)
        if(whichlambda .eq. 2) then
           lambda = -lambdafactor*eta
        else
           call getetapointresist(lambda,freq,xvec, icomp, eps, whichlambda)
        endif
        call getmagpointresist(   vel,freq,xvec, icomp, whichvel)
        
        freqx = freq(i)
        freqy = freq(j)
        freqz = freq(k)
        
        x = xvec(i)
        y = xvec(j)
        z = xvec(k)
        
        fx = freqx*x
        fy = freqy*y
        fz = freqz*z
        if((whichvel.eq.2).and.(whicheta.eq.0)) then
#if(CH_SPACEDIM==1)
           divf = one
#else
           divf = -freqy*freqy*sin(fy)
#if CH_SPACEDIM==3
           divf = divf -freqz*freqz*sin(fz)
#endif
#endif
        else if((whichvel.eq.1).and.(whicheta.eq.3).and.(whichlambda.eq.3)) then
           divf = -three*half*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.0)) then
           divf = -three*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.0).and.(whichlambda.eq.2)) then
           divf = -(two - lambdafactor)*freqx*freqx*sin(fx)
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.1)) then
#if CH_SPACEDIM==1
           divf = one
#else
           divf =       -three*eta*freqx*freqx*sin(fx)
           divf = divf + three*eps*freqx*freqx*cos(fx)*cos(fx)
           divf = divf +       eps*freqx*freqy*cos(fx)*cos(fy)
#if CH_SPACEDIM==3
           divf = divf +       eps*freqx*freqz*cos(fx)*cos(fz)
#endif
#endif
        else if((whichvel.eq.1).and.(whicheta.eq.1).and.(whichlambda.eq.2)) then
#if CH_SPACEDIM==1
           divf = one
#else
           divf =       -(two - lambdafactor)*eta*freqx*freqx*sin(fx)
           divf = divf + (two - lambdafactor)*eps*freqx*freqx*cos(fx)*cos(fx)
           divf = divf -       (lambdafactor)*eps*freqx*freqy*cos(fx)*cos(fy)
#if CH_SPACEDIM==3
           divf = divf -       (lambdafactor)*eps*freqx*freqz*cos(fx)*cos(fz)
#endif
#endif
        else if((whichvel.eq.2).and.(whicheta.eq.1)) then
#if  CH_SPACEDIM==1
           divf  = one
#else
           divf =       - eta*freqy*freqy*sin(fy)
#if CH_SPACEDIM==3
           divf = divf  - eta*freqz*freqz*sin(fz)
#endif
#endif
#if  CH_SPACEDIM==1
           divf = one
#else
           divf = divf + eps*freqy*cos(fy)*(freqy*cos(fy) + freqx*cos(fx))
#if CH_SPACEDIM==3
           divf = divf + eps*freqz*cos(fz)*(freqz*cos(fz) + freqx*cos(fx))
#endif
#endif
        else if(whichvel.eq.4) then
           divf = zero
        else
           call MayDay_Error()
        endif
        klv = alpha*vel  + beta*divf
        return
        end
        subroutine GETPHIPOINT(
     &           phi
     &           ,freq
     &           ,x
     &           )

      implicit none
      REAL_T phi
      REAL_T freq(0:2)
      REAL_T x(0:2)
        phi = sin(freq(0)*x(0))
     &                * sin(freq(1)*x(1))
     &                * sin(freq(2)*x(2))
        return
        end
        subroutine GETLOFPHIZPOLY(
     &           lofphi
     &           ,x
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      REAL_T lofphi
      REAL_T x(0:2)
      REAL_T alpha
      REAL_T beta
        real_t phi, laplphi
        real_t dist
        external getlaplphirzfunc
        real_t getlaplphirzfunc
        external getphirzfunc
        real_t getphirzfunc
        dist = abs(x(0))
        phi = getphirzfunc(dist)
        laplphi = getlaplphirzfunc(dist)
        lofphi = alpha*phi + beta*laplphi
        return
        end
        subroutine GETPHIRZPOLY(
     &           phi
     &           ,x
     &           )

      implicit none
      REAL_T phi
      REAL_T x(0:2)
        real_t dist
        external getphirzfunc
        real_t getphirzfunc
        dist =abs(x(0))
        phi = getphirzfunc(dist)
        return
        end
      subroutine GETGRADPHIRZPOLY(
     &           gradphi
     &           ,x
     &           )

      implicit none
      REAL_T gradphi(0:2)
      REAL_T x(0:2)
        real_t dist
        external getgradphirzfunc
        real_t getgradphirzfunc
        dist = abs(x(0))
        
        gradphi(0) = getgradphirzfunc(dist)
        gradphi(1) = zero
        gradphi(2) = zero
        return
        end
        subroutine GETGRADPHIPOINT(
     &           gradphi
     &           ,freq
     &           ,x
     &           )

      implicit none
      REAL_T gradphi(0:2)
      REAL_T freq(0:2)
      REAL_T x(0:2)
        
        gradphi(0) = freq(0) * cos(freq(0)*x(0)) * sin(freq(1)*x(1)) * sin(freq(2)*x(2))
        gradphi(1) = freq(1) * sin(freq(0)*x(0)) * cos(freq(1)*x(1)) * sin(freq(2)*x(2))
        gradphi(2) = freq(2) * sin(freq(0)*x(0)) * sin(freq(1)*x(1)) * cos(freq(2)*x(2))
        return
        end
        subroutine GETLOFPHI(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1,ilofphilo2
     &           ,ilofphihi0,ilofphihi1,ilofphihi2
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,alpha
     &           ,beta
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer ilofphilo0,ilofphilo1,ilofphilo2
      integer ilofphihi0,ilofphihi1,ilofphihi2
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1,
     &           ilofphilo2:ilofphihi2)
      REAL_T freq(0:2)
      REAL_T dx(0:2)
      REAL_T problo(0:2)
      REAL_T probhi(0:2)
      REAL_T alpha
      REAL_T beta
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
        integer i,j,k
        real_t x(0:CH_SPACEDIM-1)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          x(2) = (k+half)*dx(2) + problo(2)
          call getlofphipoint(lofphi(i,j,k),freq,x,alpha,beta)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine GETLOFPHIPOINT(
     &           lofphi
     &           ,freq
     &           ,x
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      REAL_T lofphi
      REAL_T freq(0:2)
      REAL_T x(0:2)
      REAL_T alpha
      REAL_T beta
        real_t fac,phi
        fac = -(freq(0)**2
     &                  + freq(1)**2
     &                  + freq(2)**2)
        phi = (sin(freq(0)*x(0))
     &                 * sin(freq(1)*x(1))
     &                 * sin(freq(2)*x(2)))
        lofphi = fac*phi
        lofphi = alpha*phi + beta*lofphi
        return
        end
        subroutine GETDBGPHI(
     &           dbgphi
     &           ,idbgphilo0,idbgphilo1,idbgphilo2
     &           ,idbgphihi0,idbgphihi1,idbgphihi2
     &           ,beta
     &           ,ibetalo0,ibetalo1,ibetalo2
     &           ,ibetahi0,ibetahi1,ibetahi2
     &           ,freq
     &           ,dx
     &           ,problo
     &           ,probhi
     &           ,alpha
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer idbgphilo0,idbgphilo1,idbgphilo2
      integer idbgphihi0,idbgphihi1,idbgphihi2
      REAL_T dbgphi(
     &           idbgphilo0:idbgphihi0,
     &           idbgphilo1:idbgphihi1,
     &           idbgphilo2:idbgphihi2)
      integer ibetalo0,ibetalo1,ibetalo2
      integer ibetahi0,ibetahi1,ibetahi2
      REAL_T beta(
     &           ibetalo0:ibetahi0,
     &           ibetalo1:ibetahi1,
     &           ibetalo2:ibetahi2)
      REAL_T freq(0:2)
      REAL_T dx(0:2)
      REAL_T problo(0:2)
      REAL_T probhi(0:2)
      REAL_T alpha
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
        integer i,j,k
        real_t x(0:CH_SPACEDIM-1)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          x(2) = (k+half)*dx(2) + problo(2)
          call getdbgphipoint(dbgphi(i,j,k),
     &        beta(i,j,k),freq,x,alpha)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine GETDBGPHIPOINT(
     &           dbgphi
     &           ,beta
     &           ,freq
     &           ,x
     &           ,alpha
     &           )

      implicit none
      REAL_T dbgphi
      REAL_T beta
      REAL_T freq(0:2)
      REAL_T x(0:2)
      REAL_T alpha
        real_t gradphi(0:CH_SPACEDIM-1),gradbeta(0:CH_SPACEDIM-1)
        real_t alphaphiplusbetalapphi,gradbetadotgradphi
        call getbetapoint(beta,freq,x)
        call getlofphipoint(alphaphiplusbetalapphi,freq,x,alpha,beta)
        call getgradbetapoint(gradbeta,freq,x)
        call getgradphipoint(gradphi,freq,x)
        gradbetadotgradphi = gradbeta(0)*gradphi(0)
     &                               + gradbeta(1)*gradphi(1)
     &                               + gradbeta(2)*gradphi(2)
        dbgphi = alphaphiplusbetalapphi
        dbgphi = dbgphi + gradbetadotgradphi
        return
        end
        subroutine GETBETAPOINT(
     &           beta
     &           ,freq
     &           ,x
     &           )

      implicit none
      REAL_T beta
      REAL_T freq(0:2)
      REAL_T x(0:2)
        beta = x(0)*x(0)
     &                 + x(1)*x(1)
     &                 + x(2)*x(2)
        return
        end
        subroutine GETGRADBETAPOINT(
     &           gradbeta
     &           ,freq
     &           ,x
     &           )

      implicit none
      REAL_T gradbeta(0:2)
      REAL_T freq(0:2)
      REAL_T x(0:2)
        integer idir
        do idir = 0, CH_SPACEDIM-1
            gradbeta(idir) = two*x(idir)
        enddo
        return
        end
        subroutine GETBETAGRADPHIPOINT(
     &           gradphi
     &           ,freq
     &           ,x
     &           )

      implicit none
      REAL_T gradphi(0:2)
      REAL_T freq(0:2)
      REAL_T x(0:2)
        integer idir
        real_t beta
        call getbetapoint(beta,freq,x)
        call getgradphipoint(gradphi,freq,x)
        do idir = 0, CH_SPACEDIM-1
           gradphi(idir) = gradphi(idir)*beta
        enddo
        return
        end
        subroutine GETSRC(
     &           src
     &           ,isrclo0,isrclo1,isrclo2
     &           ,isrchi0,isrchi1,isrchi2
     &           ,freq
     &           ,dx
     &           ,diffconst
     &           ,problo
     &           ,probhi
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer isrclo0,isrclo1,isrclo2
      integer isrchi0,isrchi1,isrchi2
      REAL_T src(
     &           isrclo0:isrchi0,
     &           isrclo1:isrchi1,
     &           isrclo2:isrchi2)
      REAL_T freq(0:2)
      REAL_T dx(0:2)
      REAL_T diffconst
      REAL_T problo(0:2)
      REAL_T probhi(0:2)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
        integer i,j,k
        real_t x(0:CH_SPACEDIM-1)
        
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          
          x(0) = (i+half)*dx(0) + problo(0)
          x(1) = (j+half)*dx(1) + problo(1)
          x(2) = (k+half)*dx(2) + problo(2)
          call getsrcpoint(src(i,j,k),freq,x,diffconst)
        
      enddo
      enddo
      enddo
        return
        end
        subroutine GETSRCPOINT(
     &           src
     &           ,freq
     &           ,x
     &           ,diffconst
     &           )

      implicit none
      REAL_T src
      REAL_T freq(0:2)
      REAL_T x(0:2)
      REAL_T diffconst
        real_t fac,phi
        fac = -(freq(0)**2
     &                  + freq(1)**2
     &                  + freq(2)**2)
        phi = (sin(freq(0)*x(0))
     &                 * sin(freq(1)*x(1))
     &                 * sin(freq(2)*x(2)))
        src = (-fac*diffconst)*phi
        return
        end
