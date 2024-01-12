      subroutine LTDget(p1,p21,p32,p3,p2,p13,m02,m12,m22,m32)
      implicit none
#include "../../util/DMNLO_LoopIntegrals.h"
#include "looptools.h"
            
      memindex XDget
      integer epsi
      double precision p1,p21,p32,p3,p2,p13,m02,m12,m22,m32,fac
      double complex xd0Eps, xd1Eps,xd2Eps,xd3Eps,xd00Eps ,xd11Eps, xd12Eps,xd13Eps ,xd22Eps,xd23Eps 
      double complex xd33Eps ,xd001Eps,xd002Eps ,xd003Eps,xd111Eps,xd112Eps,xd122Eps,xd113Eps,xd123Eps,xd133Eps
      double complex xd222Eps,xd223Eps,xd233Eps ,xd333Eps 
 
      if(EpsPole.le.0d0) then
 
      call setuvdiv(0d0)
      call setlambda(EpsPole)

      epsi = getepsi()
 
      XDget = Dget(p1,p21,p32,p3,p2,p13,m02,m12,m22,m32)

      xd0 = Dval(dd0+epsi,XDget)

      xd1 = Dval(dd1+epsi,XDget)	
      xd2 = Dval(dd2+epsi,XDget)
      xd3 = Dval(dd3+epsi,XDget)

      xd00 = Dval(dd00+epsi,XDget)	
      xd11 = Dval(dd11+epsi,XDget)		
      xd12 = Dval(dd12+epsi,XDget)	
      xd13 = Dval(dd13+epsi,XDget)		
      xd22 = Dval(dd22+epsi,XDget)
      xd23 = Dval(dd23+epsi,XDget)	
      xd33 = Dval(dd33+epsi,XDget)	 	

      xd001 = Dval(dd001+epsi,XDget)	 
      xd002 = Dval(dd002+epsi,XDget)	  	
      xd003 = Dval(dd003+epsi,XDget)	  	
      xd111 = Dval(dd111+epsi,XDget)	 
      xd112 = Dval(dd112+epsi,XDget)
      xd122 = Dval(dd122+epsi,XDget)	 
      xd113 = Dval(dd113+epsi,XDget)	 
      xd123 = Dval(dd123+epsi,XDget)	 
      xd133 = Dval(dd133+epsi,XDget)	 
      xd222 = Dval(dd222+epsi,XDget)	 
      xd223 = Dval(dd223+epsi,XDget)	 
      xd233 = Dval(dd233+epsi,XDget)	 
      xd333 = Dval(dd333+epsi,XDget)	 

      if(EpsPole.eq.0d0) then
 
      call setuvdiv(0d0)

! Compute coefficients of double poles
      call setlambda(-2d0)

      epsi = getepsi()
 
      XDget = Dget(p1,p21,p32,p3,p2,p13,m02,m12,m22,m32)

      xd0Eps = Dval(dd0+epsi,XDget)

      xd1Eps = Dval(dd1+epsi,XDget)	
      xd2Eps = Dval(dd2+epsi,XDget)
      xd3Eps = Dval(dd3+epsi,XDget)

      xd00Eps = Dval(dd00+epsi,XDget)	
      xd11Eps = Dval(dd11+epsi,XDget)		
      xd12Eps = Dval(dd12+epsi,XDget)	
      xd13Eps = Dval(dd13+epsi,XDget)		
      xd22Eps = Dval(dd22+epsi,XDget)
      xd23Eps = Dval(dd23+epsi,XDget)	
      xd33Eps = Dval(dd33+epsi,XDget)	 	

      xd001Eps = Dval(dd001+epsi,XDget)	 
      xd002Eps = Dval(dd002+epsi,XDget)	  	
      xd003Eps = Dval(dd003+epsi,XDget)	  	
      xd111Eps = Dval(dd111+epsi,XDget)	 
      xd112Eps = Dval(dd112+epsi,XDget)	
      xd122Eps = Dval(dd122+epsi,XDget)	
      xd113Eps = Dval(dd113+epsi,XDget)	 
      xd123Eps = Dval(dd123+epsi,XDget)	 
      xd133Eps = Dval(dd133+epsi,XDget)	 
      xd222Eps = Dval(dd222+epsi,XDget)	 
      xd223Eps = Dval(dd223+epsi,XDget)	 
      xd233Eps = Dval(dd233+epsi,XDget)	 
      xd333Eps = Dval(dd333+epsi,XDget)	

! transfer conventions

      fac = LPi**2/6d0

      xd0 = xd0 - fac*xd0Eps

      xd1 = xd1 - fac*xd1Eps
      xd2 = xd2 - fac*xd2Eps
      xd3 = xd3 - fac*xd3Eps

      xd00 = xd00 - fac*xd00Eps
      xd11 = xd11 - fac*xd11Eps
      xd12 = xd12 - fac*xd12Eps
      xd13 = xd13 - fac*xd13Eps	
      xd22 = xd22 - fac*xd22Eps
      xd23 = xd23 - fac*xd23Eps
      xd33 = xd33 - fac*xd33Eps

      xd001 = xd001 - fac*xd001Eps
      xd002 = xd002 - fac*xd002Eps	  	
      xd003 = xd003 - fac*xd003Eps	  	
      xd111 = xd111 - fac*xd111Eps	 
      xd112 = xd112 - fac*xd112Eps	 
      xd122 = xd122 - fac*xd112Eps	 
      xd113 = xd113 - fac*xd113Eps	 
      xd123 = xd123 - fac*xd123Eps	 
      xd133 = xd133 - fac*xd133Eps	 
      xd222 = xd222 - fac*xd222Eps	 
      xd223 = xd223 - fac*xd223Eps	 
      xd233 = xd233 - fac*xd233Eps	 
      xd333 = xd333 - fac*xd333Eps	 

      endif

      else if (EpsPole.eq.1d0) then

! Sets all coefficents to zero since there are no UV divergences. 
      call myDget(p1,p21,p32,p3,p2,p13,m02,m12,m22,m32) 
            
      else

      write(*,*) "in LTDget: value vor EpsPole is not valid!"

      endif

      
      end



!      subroutine LTCget(p1,p21,p2,m02,m12,m22)
!      implicit none
!#include "../../util/DMNLO_LoopIntegrals.h"
!#include "looptools.h"
!            
!      memindex XCget
!      integer epsi
!      double precision p1,p21,p2,m02,m12,m22,fac
!      double complex xc0Eps,xc1Eps,xc2Eps,xc00Eps,xc11Eps,xc12Eps
!      double complex xc22Eps,xc001Eps,xc002Eps,xc111Eps,xc122Eps,xc112Eps,xc222Eps
!
! 
!      if(EpsPole.le.0d0) then
! 
!      call ltini
! 
!      call setuvdiv(0d0)
!      call setmudim(muSc**2)
!      call setlambda(EpsPole)
!
!      epsi = getepsi()
! 
!      XCget = Cget(p1,p21,p2,m02,m12,m22)
!
!      xc0 = Cval(cc0+epsi,XCget)
!
!      xc1 = Cval(cc1+epsi,XCget)	
!      xc2 = Cval(cc2+epsi,XCget)
!      
!      xc00 = Cval(cc00+epsi,XCget)	
!      xc11 = Cval(cc11+epsi,XCget)		
!      xc12 = Cval(cc12+epsi,XCget)		
!      xc22 = Cval(cc22+epsi,XCget)
!             
!      xc001 = Cval(cc001+epsi,XCget)	 
!      xc002 = Cval(cc002+epsi,XCget)	  		  	
!      xc111 = Cval(cc111+epsi,XCget)	
!      xc112 = Cval(cc112+epsi,XCget)
!      xc122 = Cval(cc122+epsi,XCget)	 
!      xc222 = Cval(cc222+epsi,XCget)	 
!             
!      call ltexi
!
!      if(EpsPole.eq.0d0) then
! 
!      call ltini
! 
!      call setuvdiv(0d0)
!      call setmudim(muSc**2)
!
!! Compute coefficients of double poles
!      call setlambda(-2d0)
!
!      epsi = getepsi()
! 
!      XCget = Cget(p1,p21,p2,m02,m12,m22)
!
!      xc0Eps = Cval(cc0+epsi,XCget)
!
!      xc1Eps = Cval(cc1+epsi,XCget)	
!      xc2Eps = Cval(cc2+epsi,XCget)
!      
!      xc00Eps = Cval(cc00+epsi,XCget)	
!      xc11Eps = Cval(cc11+epsi,XCget)		
!      xc12Eps = Cval(cc12+epsi,XCget)		
!      xc22Eps = Cval(cc22+epsi,XCget)
!             
!      xc001Eps = Cval(cc001+epsi,XCget)	 
!      xc002Eps = Cval(cc002+epsi,XCget)	  		  	
!      xc111Eps = Cval(cc111+epsi,XCget)	 
!      xc112Eps = Cval(cc112+epsi,XCget)	
!      xc122Eps = Cval(cc122+epsi,XCget) 
!      xc222Eps = Cval(cc222+epsi,XCget)	
!
!      call ltexi
!
!! transfer conventions
!
!      fac = LPi**2/6d0
!
!      xc0 = xc0 -  fac*xc0Eps
!
!      xc1 = xc1 - fac*xc1Eps
!      xc2 = xc2 - fac*xc2Eps
!      
!      xc00 = xc00 - fac*xc00Eps
!      xc11 = xc11 - fac*xc11Eps		
!      xc12 = xc12 - fac*xc12Eps		
!      xc22 = xc22 - fac*xc22Eps
!             
!      xc001 = xc001 - fac*xc001Eps	 
!      xc002 = xc002 - fac*xc002Eps	  		  	
!      xc111 = xc111 - fac*xc111Eps	 
!      xc112 = xc112 - fac*xc112Eps	 
!      xc122 = xc122 - fac*xc122Eps	 
!      xc222 = xc222 - fac*xc222Eps
!
!      endif
!
!      else if (EpsPole.eq.1d0) then
!
!      call myCget(p1,p21,p2,m02,m12,m22)
!            
!      else
!
!      write(*,*) "in myCget: value vor EpsPole is not valid!"
!
!      endif
!
!      
!      end