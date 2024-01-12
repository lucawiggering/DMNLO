

      subroutine QQ2xx_M2AmpDigestHH(argxh,argyh)

      implicit none

      integer io_error,argxh,argyh
      double precision HH_QX(6),HH_SX(6),HH_VX(6),HH_TX(6),HH_UX(6)


!initialize O(alpha_s) corrections
      call QQ2xx_CorrQXS(0,argxh,argyh)
      call QQ2xx_CorrTXSA(0,argxh)
      call QQ2xx_CorrUXSB(0,argxh)                  
      call QQ2xx_CorrTXSB(0,argyh)            
      call QQ2xx_CorrUXSA(0,argyh)

      open(unit=1,file='AmpDigest.txt',status='old',action='write',iostat=io_error,position='append')

      if (io_error.eq.0) then
        write(1,*) '*****************************************************************************************************************'
        write(1,*) '*****************************************************************************************************************'

        write(1,*) 'argxh= ',argxh,'argyh= ',argyh

        write(1,*) '*****************************************************************************************************************'

        call HH_QX_(1,argxh,argyh,HH_QX)    
      
        write(1,1) 'Tree:                   ',' QQ= ',HH_QX(1),' QS= ',HH_QX(2),' QV= ',HH_QX(3),' QT= ',HH_QX(4),' QU= ',HH_QX(5)

        call HH_QX_(2,argxh,argyh,HH_QX)

        write(1,1) 'Vertex_Gluon:           ',' QQ= ',HH_QX(1),' QS= ',HH_QX(2),' QV= ',HH_QX(3),' QT= ',HH_QX(4),' QU= ',HH_QX(5)

        call HH_QX_(3,argxh,argyh,HH_QX)

        write(1,1) 'Vertex_ScalarLoop:      ',' QQ= ',HH_QX(1),' QS= ',HH_QX(2),' QV= ',HH_QX(3),' QT= ',HH_QX(4),' QU= ',HH_QX(5)

        call HH_QX_(4,argxh,argyh,HH_QX)

        write(1,1) 'Vertex_Counterterm:     ',' QQ= ',HH_QX(1),' QS= ',HH_QX(2),' QV= ',HH_QX(3),' QT= ',HH_QX(4),' QU= ',HH_QX(5)

        write(1,*) '*****************************************************************************************************************'

        call HH_SX_(1,argxh,argyh,HH_SX)    
      
        write(1,1) 'Tree:                   ',' SQ= ',HH_SX(1),' SS= ',HH_SX(2),' SV= ',HH_SX(3),' ST= ',HH_SX(4),' SU= ',HH_SX(5)

        call HH_SX_(2,argxh,argyh,HH_SX)

        write(1,1) 'Vertex_Gluon:           ',' SQ= ',HH_SX(1),' SS= ',HH_SX(2),' SV= ',HH_SX(3),' ST= ',HH_SX(4),' SU= ',HH_SX(5)

        call HH_SX_(3,argxh,argyh,HH_SX)

        write(1,1) 'Vertex_Gluino:          ',' SQ= ',HH_SX(1),' SS= ',HH_SX(2),' SV= ',HH_SX(3),' ST= ',HH_SX(4),' SU= ',HH_SX(5)

        call HH_SX_(4,argxh,argyh,HH_SX)

        write(1,1) 'Vertex_ScalarLoop:      ',' SQ= ',HH_SX(1),' SS= ',HH_SX(2),' SV= ',HH_SX(3),' ST= ',HH_SX(4),' SU= ',HH_SX(5)

        call HH_SX_(5,argxh,argyh,HH_SX)

        write(1,1) 'Vertex_Counterterm:     ',' SQ= ',HH_SX(1),' SS= ',HH_SX(2),' SV= ',HH_SX(3),' ST= ',HH_SX(4),' SU= ',HH_SX(5)

        write(1,*) '*****************************************************************************************************************'

        call HH_VX_(1,argxh,argyh,HH_VX)    
      
        write(1,1) 'Tree:                   ',' VQ= ',HH_VX(1),' VS= ',HH_VX(2),' VV= ',HH_VX(3),' VT= ',HH_VX(4),' VU= ',HH_VX(5)

        call HH_VX_(2,argxh,argyh,HH_VX)

        write(1,1) 'Vertex_Gluon:           ',' VQ= ',HH_VX(1),' VS= ',HH_VX(2),' VV= ',HH_VX(3),' VT= ',HH_VX(4),' VU= ',HH_VX(5)

        call HH_VX_(3,argxh,argyh,HH_VX)

        write(1,1) 'Vertex_Gluino:          ',' VQ= ',HH_VX(1),' VS= ',HH_VX(2),' VV= ',HH_VX(3),' VT= ',HH_VX(4),' VU= ',HH_VX(5)

        call HH_VX_(4,argxh,argyh,HH_VX)

        write(1,1) 'Vertex_SSSSCoupl:       ',' VQ= ',HH_VX(1),' VS= ',HH_VX(2),' VV= ',HH_VX(3),' VT= ',HH_VX(4),' VU= ',HH_VX(5)

        call HH_VX_(5,argxh,argyh,HH_VX)

        write(1,1) 'Vertex_SSVVCoupl_p1:    ',' VQ= ',HH_VX(1),' VS= ',HH_VX(2),' VV= ',HH_VX(3),' VT= ',HH_VX(4),' VU= ',HH_VX(5)

        call HH_VX_(6,argxh,argyh,HH_VX)

        write(1,1) 'Vertex_SSVVCoupl_p2:    ',' VQ= ',HH_VX(1),' VS= ',HH_VX(2),' VV= ',HH_VX(3),' VT= ',HH_VX(4),' VU= ',HH_VX(5)

        call HH_VX_(7,argxh,argyh,HH_VX)

        write(1,1) 'Vertex_Counterterm:     ',' VQ= ',HH_VX(1),' VS= ',HH_VX(2),' VV= ',HH_VX(3),' VT= ',HH_VX(4),' VU= ',HH_VX(5)

        write(1,*) '*****************************************************************************************************************'

        call HH_TX_(1,1,1,argxh,argyh,HH_TX)    
      
        write(1,1) 'Tree:                   ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(1,1,2,argxh,argyh,HH_TX)

        write(1,1) 'Prop_Gluon:             ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)
  
        call HH_TX_(1,1,3,argxh,argyh,HH_TX)

        write(1,1) 'Prop_Gluino:            ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(1,1,4,argxh,argyh,HH_TX)

        write(1,1) 'Prop_Scalar:            ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(1,1,5,argxh,argyh,HH_TX)

        write(1,1) 'Prop_Vector:            ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(1,1,6,argxh,argyh,HH_TX)

        write(1,1) 'Prop_Counter:           ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(2,1,1,argxh,argyh,HH_TX)

        write(1,1) 'Vertex_Gluon_Atype:     ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(1,2,1,argxh,argyh,HH_TX)

        write(1,1) 'Vertex_Gluon_Btype:     ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(3,1,1,argxh,argyh,HH_TX)

        write(1,1) 'Vertex_Gluino_Atype:    ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(1,3,1,argxh,argyh,HH_TX)

        write(1,1) 'Vertex_Gluino_Btype:    ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(4,1,1,argxh,argyh,HH_TX)

        write(1,1) 'Vertex_SSSS_Atype:      ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(1,4,1,argxh,argyh,HH_TX)

        write(1,1) 'Vertex_SSSS_Btype:      ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(5,1,1,argxh,argyh,HH_TX)

        write(1,1) 'Vertex_Counter_Atype:   ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call HH_TX_(1,5,1,argxh,argyh,HH_TX)

        write(1,1) 'Vertex_Counter_Btype:   ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call QQ2xx_M2GluonBoxHH(argxh,argyh,HH_TX,HH_UX)

        write(1,1) 'Box_Gluon:              ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call QQ2xx_M2GluinoBoxHH(argxh,argyh,HH_TX,HH_UX)

        write(1,1) 'Box_Gluino:             ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        call QQ2xx_M2TriangleBoxHH(argxh,argyh,HH_TX,HH_UX)

        write(1,1) 'Box_Triangle:           ',' TQ= ',HH_TX(1),' TS= ',HH_TX(2),' TV= ',HH_TX(3),' TT= ',HH_TX(4),' TU= ',HH_TX(5)

        write(1,*) '*****************************************************************************************************************'

        call HH_UX_(1,1,1,argxh,argyh,HH_UX)    
      
        write(1,1) 'Tree:                   ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(1,1,2,argxh,argyh,HH_UX)

        write(1,1) 'Prop_Gluon:             ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)
  
        call HH_UX_(1,1,3,argxh,argyh,HH_UX)

        write(1,1) 'Prop_Gluino:            ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(1,1,4,argxh,argyh,HH_UX)

        write(1,1) 'Prop_Scalar:            ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(1,1,5,argxh,argyh,HH_UX)

        write(1,1) 'Prop_Vector:            ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(1,1,6,argxh,argyh,HH_UX)

        write(1,1) 'Prop_Counter:           ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(2,1,1,argxh,argyh,HH_UX)

        write(1,1) 'Vertex_Gluon_Atype:     ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(1,2,1,argxh,argyh,HH_UX)

        write(1,1) 'Vertex_Gluon_Btype:     ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(3,1,1,argxh,argyh,HH_UX)

        write(1,1) 'Vertex_Gluino_Atype:    ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(1,3,1,argxh,argyh,HH_UX)

        write(1,1) 'Vertex_Gluino_Btype:    ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(4,1,1,argxh,argyh,HH_UX)

        write(1,1) 'Vertex_SSSS_Atype:      ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(1,4,1,argxh,argyh,HH_UX)

        write(1,1) 'Vertex_SSSS_Btype:      ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(5,1,1,argxh,argyh,HH_UX)

        write(1,1) 'Vertex_Counter_Atype:   ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call HH_UX_(1,5,1,argxh,argyh,HH_UX)

        write(1,1) 'Vertex_Counter_Btype:   ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call QQ2xx_M2GluonBoxHH(argxh,argyh,HH_TX,HH_UX)

        write(1,1) 'Box_Gluon:              ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call QQ2xx_M2GluinoBoxHH(argxh,argyh,HH_TX,HH_UX)

        write(1,1) 'Box_Gluino:             ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

        call QQ2xx_M2TriangleBoxHH(argxh,argyh,HH_TX,HH_UX)

        write(1,1) 'Box_Triangle:           ',' UQ= ',HH_UX(1),' US= ',HH_UX(2),' UV= ',HH_UX(3),' UT= ',HH_UX(4),' UU= ',HH_UX(5)

  1   Format( A24,  5(A5,E13.4))

      else

           write(*,*)  'Error during writing into AmpDigest.txt'

      end if

      close(unit=1)

      end





!**********************************************************************************************************************************
!**********************************************************************************************************************************





      subroutine QQ2xx_M2AmpDigestVH(argx,argyh)

      implicit none

      integer io_error,argx,argyh
      double precision VH_SX(6),VH_VX(6),VH_TX(6),VH_UX(6)


!initialize further argx/argyh dependent O(alpha_s) corrections
      call QQ2xx_CorrTXSB(0,argyh)
      call QQ2xx_CorrUXSA(0,argyh) 
      call QQ2xx_CorrTXVA(0,argx)
      call QQ2xx_CorrUXVB(0,argx)  

      open(unit=1,file='AmpDigest.txt',status='old',action='write',iostat=io_error,position='append')

      if (io_error.eq.0) then
        write(1,*) '*****************************************************************************************************************'
        write(1,*) '*****************************************************************************************************************'

        write(1,*) 'argx= ',argx,'argyh= ',argyh

        write(1,*) '*****************************************************************************************************************'

        call VH_SX_(1,argx,argyh,VH_SX)    
      
        write(1,1) 'Tree:                   ',' SQ= ',VH_SX(1),' SS= ',VH_SX(2),' SV= ',VH_SX(3),' ST= ',VH_SX(4),' SU= ',VH_SX(5)

        call VH_SX_(2,argx,argyh,VH_SX)

        write(1,1) 'Vertex_Gluon:           ',' SQ= ',VH_SX(1),' SS= ',VH_SX(2),' SV= ',VH_SX(3),' ST= ',VH_SX(4),' SU= ',VH_SX(5)

        call VH_SX_(3,argx,argyh,VH_SX)

        write(1,1) 'Vertex_Gluino:          ',' SQ= ',VH_SX(1),' SS= ',VH_SX(2),' SV= ',VH_SX(3),' ST= ',VH_SX(4),' SU= ',VH_SX(5)

        call VH_SX_(4,argx,argyh,VH_SX)

        write(1,1) 'Vertex_ScalarLoop:      ',' SQ= ',VH_SX(1),' SS= ',VH_SX(2),' SV= ',VH_SX(3),' ST= ',VH_SX(4),' SU= ',VH_SX(5)

        call VH_SX_(5,argx,argyh,VH_SX)

        write(1,1) 'Vertex_Counterterm:     ',' SQ= ',VH_SX(1),' SS= ',VH_SX(2),' SV= ',VH_SX(3),' ST= ',VH_SX(4),' SU= ',VH_SX(5)

        write(1,*) '*****************************************************************************************************************'

        call VH_VX_(1,argx,argyh,VH_VX)    
      
        write(1,1) 'Tree:                   ',' VQ= ',VH_VX(1),' VS= ',VH_VX(2),' VV= ',VH_VX(3),' VT= ',VH_VX(4),' VU= ',VH_VX(5)

        call VH_VX_(2,argx,argyh,VH_VX)

        write(1,1) 'Vertex_Gluon:           ',' VQ= ',VH_VX(1),' VS= ',VH_VX(2),' VV= ',VH_VX(3),' VT= ',VH_VX(4),' VU= ',VH_VX(5)

        call VH_VX_(3,argx,argyh,VH_VX)

        write(1,1) 'Vertex_Gluino:          ',' VQ= ',VH_VX(1),' VS= ',VH_VX(2),' VV= ',VH_VX(3),' VT= ',VH_VX(4),' VU= ',VH_VX(5)

        call VH_VX_(4,argx,argyh,VH_VX)

        write(1,1) 'Vertex_SSSSCoupl:       ',' VQ= ',VH_VX(1),' VS= ',VH_VX(2),' VV= ',VH_VX(3),' VT= ',VH_VX(4),' VU= ',VH_VX(5)

        call VH_VX_(5,argx,argyh,VH_VX)

        write(1,1) 'Vertex_SSVVCoupl_p1:    ',' VQ= ',VH_VX(1),' VS= ',VH_VX(2),' VV= ',VH_VX(3),' VT= ',VH_VX(4),' VU= ',VH_VX(5)

        call VH_VX_(6,argx,argyh,VH_VX)

        write(1,1) 'Vertex_SSVVCoupl_p2:    ',' VQ= ',VH_VX(1),' VS= ',VH_VX(2),' VV= ',VH_VX(3),' VT= ',VH_VX(4),' VU= ',VH_VX(5)


        call VH_VX_(7,argx,argyh,VH_VX)

        write(1,1) 'Vertex_Counterterm:     ',' VQ= ',VH_VX(1),' VS= ',VH_VX(2),' VV= ',VH_VX(3),' VT= ',VH_VX(4),' VU= ',VH_VX(5)

        write(1,*) '*****************************************************************************************************************'

        call VH_TX_(1,1,1,argx,argyh,VH_TX)    
      
        write(1,1) 'Tree:                   ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(1,1,2,argx,argyh,VH_TX)

        write(1,1) 'Prop_Gluon:             ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)
  
        call VH_TX_(1,1,3,argx,argyh,VH_TX)

        write(1,1) 'Prop_Gluino:            ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(1,1,4,argx,argyh,VH_TX)

        write(1,1) 'Prop_Scalar:            ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(1,1,5,argx,argyh,VH_TX)

        write(1,1) 'Prop_Vector:            ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(1,1,6,argx,argyh,VH_TX)

        write(1,1) 'Prop_Counter:           ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5) 

        call VH_TX_(2,1,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_Gluon_Atype:     ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(1,2,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_Gluon_Btype:     ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(3,1,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_Gluino_Atype:    ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(1,3,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_Gluino_Btype:    ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(4,1,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_SSSS_Atype:      ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(1,4,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_SSSS_Btype:      ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(5,1,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_SSVVp1_Atype:    ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(6,1,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_SSVVprop_Atype:  ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(7,1,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_Counter_Atype:   ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call VH_TX_(1,5,1,argx,argyh,VH_TX)

        write(1,1) 'Vertex_Counter_Btype:   ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call QQ2xx_M2GluonBoxVH(argx,argyh,VH_TX,VH_UX)

        write(1,1) 'Box_Gluon:              ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call QQ2xx_M2GluinoBoxVH(argx,argyh,VH_TX,VH_UX)

        write(1,1) 'Box_Gluino:             ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call QQ2xx_M2GluonTriangleBoxVH(argx,argyh,VH_TX,VH_UX)

        write(1,1) 'Box_GluonTriangle:      ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        call QQ2xx_M2sQuarkTriangleBoxVH(argx,argyh,VH_TX,VH_UX)

        write(1,1) 'Box_sQuarkTriangle:     ',' TQ= ',VH_TX(1),' TS= ',VH_TX(2),' TV= ',VH_TX(3),' TT= ',VH_TX(4),' TU= ',VH_TX(5)

        write(1,*) '*****************************************************************************************************************'

        call VH_UX_(1,1,1,argx,argyh,VH_UX)    
      
        write(1,1) 'Tree:                   ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,1,2,argx,argyh,VH_UX)

        write(1,1) 'Prop_Gluon:             ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)
  
        call VH_UX_(1,1,3,argx,argyh,VH_UX)

        write(1,1) 'Prop_Gluino:            ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,1,4,argx,argyh,VH_UX)

        write(1,1) 'Prop_Scalar:            ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,1,5,argx,argyh,VH_UX)

        write(1,1) 'Prop_Vector:            ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,1,6,argx,argyh,VH_UX)

        write(1,1) 'Prop_Counter:           ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(2,1,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_Gluon_Atype:     ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,2,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_Gluon_Btype:     ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5) 

        call VH_UX_(3,1,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_Gluino_Atype:    ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,3,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_Gluino_Btype:    ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(4,1,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_SSSS_Atype:      ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,4,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_SSSS_Btype:      ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,5,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_SSVVprop_Btype:  ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,6,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_SSVVp2_Btype:    ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(1,7,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_Counter_Btype:   ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call VH_UX_(5,1,1,argx,argyh,VH_UX)

        write(1,1) 'Vertex_Counter_Btype:   ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call QQ2xx_M2GluonBoxVH(argx,argyh,VH_TX,VH_UX)

        write(1,1) 'Box_Gluon:              ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call QQ2xx_M2GluinoBoxVH(argx,argyh,VH_TX,VH_UX)

        write(1,1) 'Box_Gluino:             ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call QQ2xx_M2GluonTriangleBoxVH(argx,argyh,VH_TX,VH_UX)

        write(1,1) 'Box_GluonTriangle:      ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

        call QQ2xx_M2sQuarkTriangleBoxVH(argx,argyh,VH_TX,VH_UX)

        write(1,1) 'Box_sQuarkTriangle:     ',' UQ= ',VH_UX(1),' US= ',VH_UX(2),' UV= ',VH_UX(3),' UT= ',VH_UX(4),' UU= ',VH_UX(5)

  1   Format( A24,  5(A5,E13.4))

      else

           write(*,*)  'Error during writing into AmpDigest.txt'

      end if

      close(unit=1)

      end




!**********************************************************************************************************************************
!**********************************************************************************************************************************





      subroutine QQ2xx_M2AmpDigestHV(argxh,argy)

      implicit none

      integer io_error,argxh,argy
      double precision HV_SX(6),HV_VX(6),HV_TX(6),HV_UX(6)

!initialize further argxh/argy dependent O(alpha_s) corrections
      call QQ2xx_CorrTXSA(0,argxh)
      call QQ2xx_CorrUXSB(0,argxh) 
      call QQ2xx_CorrTXVB(0,argy)
      call QQ2xx_CorrUXVA(0,argy)  


      open(unit=1,file='AmpDigest.txt',status='old',action='write',iostat=io_error,position='append')

      if (io_error.eq.0) then
        write(1,*) '*****************************************************************************************************************'
        write(1,*) '*****************************************************************************************************************'

        write(1,*) 'argxh= ',argxh,'argy= ',argy

        write(1,*) '*****************************************************************************************************************'

        call HV_SX_(1,argxh,argy,HV_SX)    
      
        write(1,1) 'Tree:                   ',' SQ= ',HV_SX(1),' SS= ',HV_SX(2),' SV= ',HV_SX(3),' ST= ',HV_SX(4),' SU= ',HV_SX(5)

        call HV_SX_(2,argxh,argy,HV_SX)

        write(1,1) 'Vertex_Gluon:           ',' SQ= ',HV_SX(1),' SS= ',HV_SX(2),' SV= ',HV_SX(3),' ST= ',HV_SX(4),' SU= ',HV_SX(5)

        call HV_SX_(3,argxh,argy,HV_SX)

        write(1,1) 'Vertex_Gluino:          ',' SQ= ',HV_SX(1),' SS= ',HV_SX(2),' SV= ',HV_SX(3),' ST= ',HV_SX(4),' SU= ',HV_SX(5)

        call HV_SX_(4,argxh,argy,HV_SX)

        write(1,1) 'Vertex_ScalarLoop:      ',' SQ= ',HV_SX(1),' SS= ',HV_SX(2),' SV= ',HV_SX(3),' ST= ',HV_SX(4),' SU= ',HV_SX(5)

        call HV_SX_(5,argxh,argy,HV_SX)

        write(1,1) 'Vertex_Counterterm:     ',' SQ= ',HV_SX(1),' SS= ',HV_SX(2),' SV= ',HV_SX(3),' ST= ',HV_SX(4),' SU= ',HV_SX(5)

        write(1,*) '*****************************************************************************************************************'

        call HV_VX_(1,argxh,argy,HV_VX)    
      
        write(1,1) 'Tree:                   ',' VQ= ',HV_VX(1),' VS= ',HV_VX(2),' VV= ',HV_VX(3),' VT= ',HV_VX(4),' VU= ',HV_VX(5)

        call HV_VX_(2,argxh,argy,HV_VX)

        write(1,1) 'Vertex_Gluon:           ',' VQ= ',HV_VX(1),' VS= ',HV_VX(2),' VV= ',HV_VX(3),' VT= ',HV_VX(4),' VU= ',HV_VX(5)

        call HV_VX_(3,argxh,argy,HV_VX)

        write(1,1) 'Vertex_Gluino:          ',' VQ= ',HV_VX(1),' VS= ',HV_VX(2),' VV= ',HV_VX(3),' VT= ',HV_VX(4),' VU= ',HV_VX(5)

        call HV_VX_(4,argxh,argy,HV_VX)

        write(1,1) 'Vertex_SSSSCoupl:       ',' VQ= ',HV_VX(1),' VS= ',HV_VX(2),' VV= ',HV_VX(3),' VT= ',HV_VX(4),' VU= ',HV_VX(5)

        call HV_VX_(5,argxh,argy,HV_VX)

        write(1,1) 'Vertex_SSVVCoupl_p1:    ',' VQ= ',HV_VX(1),' VS= ',HV_VX(2),' VV= ',HV_VX(3),' VT= ',HV_VX(4),' VU= ',HV_VX(5)

        call HV_VX_(6,argxh,argy,HV_VX)

        write(1,1) 'Vertex_SSVVCoupl_p2:    ',' VQ= ',HV_VX(1),' VS= ',HV_VX(2),' VV= ',HV_VX(3),' VT= ',HV_VX(4),' VU= ',HV_VX(5)

        call HV_VX_(7,argxh,argy,HV_VX)

        write(1,1) 'Vertex_Counterterm:     ',' VQ= ',HV_VX(1),' VS= ',HV_VX(2),' VV= ',HV_VX(3),' VT= ',HV_VX(4),' VU= ',HV_VX(5)

        write(1,*) '*****************************************************************************************************************'

        call HV_TX_(1,1,1,argxh,argy,HV_TX)    
      
        write(1,1) 'Tree:                   ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,1,2,argxh,argy,HV_TX)

        write(1,1) 'Prop_Gluon:             ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)
  
        call HV_TX_(1,1,3,argxh,argy,HV_TX)

        write(1,1) 'Prop_Gluino:            ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,1,4,argxh,argy,HV_TX)

        write(1,1) 'Prop_Scalar:            ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,1,5,argxh,argy,HV_TX)

        write(1,1) 'Prop_Vector:            ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,1,6,argxh,argy,HV_TX)

        write(1,1) 'Prop_Counter:           ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(2,1,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_Gluon_Atype:     ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,2,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_Gluon_Btype:     ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(3,1,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_Gluino_Atype:    ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,3,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_Gluino_Btype:    ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(4,1,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_SSSS_Atype:      ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,4,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_SSSS_Btype:      ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,5,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_SSVVprop_Btype:  ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,6,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_SSVVp2_Btype:    ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(7,1,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_Counter_Atype:   ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call HV_TX_(1,5,1,argxh,argy,HV_TX)

        write(1,1) 'Vertex_Counter_Btype:   ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call QQ2xx_M2GluonBoxHV(argxh,argy,HV_TX,HV_UX)

        write(1,1) 'Box_Gluon:              ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call QQ2xx_M2GluinoBoxHV(argxh,argy,HV_TX,HV_UX)

        write(1,1) 'Box_Gluino:             ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call QQ2xx_M2GluonTriangleBoxHV(argxh,argy,HV_TX,HV_UX)

        write(1,1) 'Box_GluonTriangle:      ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        call QQ2xx_M2sQuarkTriangleBoxHV(argxh,argy,HV_TX,HV_UX)

        write(1,1) 'Box_sQuarkTriangle:     ',' TQ= ',HV_TX(1),' TS= ',HV_TX(2),' TV= ',HV_TX(3),' TT= ',HV_TX(4),' TU= ',HV_TX(5)

        write(1,*) '*****************************************************************************************************************'

        call HV_UX_(1,1,1,argxh,argy,HV_UX)    
      
        write(1,1) 'Tree:                   ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(1,1,2,argxh,argy,HV_UX)

        write(1,1) 'Prop_Gluon:             ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)
  
        call HV_UX_(1,1,3,argxh,argy,HV_UX)

        write(1,1) 'Prop_Gluino:            ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(1,1,4,argxh,argy,HV_UX)

        write(1,1) 'Prop_Scalar:            ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(1,1,5,argxh,argy,HV_UX)

        write(1,1) 'Prop_Vector:            ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(1,1,6,argxh,argy,HV_UX)

        write(1,1) 'Prop_Counter:           ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(2,1,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_Gluon_Atype:     ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(1,2,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_Gluon_Btype:     ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(3,1,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_Gluino_Atype:    ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(1,3,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_Gluino_Btype:    ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(4,1,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_SSSS_Atype:      ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(1,4,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_SSSS_Btype:      ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(5,1,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_SSVVp1_Atype:    ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(6,1,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_SSVVprop_Atype:  ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(1,7,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_Counter_Atype:   ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call HV_UX_(5,1,1,argxh,argy,HV_UX)

        write(1,1) 'Vertex_Counter_Btype:   ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call QQ2xx_M2GluonBoxHV(argxh,argy,HV_TX,HV_UX)

        write(1,1) 'Box_Gluon:              ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call QQ2xx_M2GluinoBoxHV(argxh,argy,HV_TX,HV_UX)

        write(1,1) 'Box_Gluino:             ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call QQ2xx_M2GluonTriangleBoxHV(argxh,argy,HV_TX,HV_UX)

        write(1,1) 'Box_GluonTriangle:      ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

        call QQ2xx_M2sQuarkTriangleBoxHV(argxh,argy,HV_TX,HV_UX)

        write(1,1) 'Box_sQuarkTriangle:     ',' UQ= ',HV_UX(1),' US= ',HV_UX(2),' UV= ',HV_UX(3),' UT= ',HV_UX(4),' UU= ',HV_UX(5)

  1   Format( A24,  5(A5,E13.4))

      else

           write(*,*)  'Error during writing into AmpDigest.txt'

      end if

      close(unit=1)

      end



!**********************************************************************************************************************************
!**********************************************************************************************************************************


      subroutine QQ2xx_M2AmpDigestVV(argx,argy)

      implicit none

      integer io_error,argx,argy
      double precision VV_QX(6),VV_SX(6),VV_VX(6),VV_TX(6),VV_UX(6)


      open(unit=1,file='AmpDigest.txt',status='old',action='write',iostat=io_error,position='append')

!initialize O(alpha_s) corrections
      call QQ2xx_CorrQXV(0,argx,argy)
      call QQ2xx_CorrTXVA(0,argx)
      call QQ2xx_CorrUXVB(0,argx)                  
      call QQ2xx_CorrTXVB(0,argy)            
      call QQ2xx_CorrUXVA(0,argy)

      if (io_error.eq.0) then
        write(1,*) '*****************************************************************************************************************'
        write(1,*) '*****************************************************************************************************************'

        write(1,*) 'argx= ',argx,'argy= ',argy

        write(1,*) '*****************************************************************************************************************'

        call VV_QX_(1,argx,argy,VV_QX)    
      
        write(1,1) 'Tree:                   ',' QQ= ',VV_QX(1),' QS= ',VV_QX(2),' QV= ',VV_QX(3),' QT= ',VV_QX(4),' QU= ',VV_QX(5)

        call VV_QX_(2,argx,argy,VV_QX)

        write(1,1) 'Vertex_Gluon:           ',' QQ= ',VV_QX(1),' QS= ',VV_QX(2),' QV= ',VV_QX(3),' QT= ',VV_QX(4),' QU= ',VV_QX(5)

        call VV_QX_(3,argx,argy,VV_QX)

        write(1,1) 'Vertex_ScalarLoop:      ',' QQ= ',VV_QX(1),' QS= ',VV_QX(2),' QV= ',VV_QX(3),' QT= ',VV_QX(4),' QU= ',VV_QX(5)

        call VV_QX_(4,argx,argy,VV_QX)

        write(1,1) 'Vertex_Counterterm:     ',' QQ= ',VV_QX(1),' QS= ',VV_QX(2),' QV= ',VV_QX(3),' QT= ',VV_QX(4),' QU= ',VV_QX(5)

        write(1,*) '*****************************************************************************************************************'

        call VV_SX_(1,argx,argy,VV_SX)    
      
        write(1,1) 'Tree:                   ',' SQ= ',VV_SX(1),' SS= ',VV_SX(2),' SV= ',VV_SX(3),' ST= ',VV_SX(4),' SU= ',VV_SX(5)

        call VV_SX_(2,argx,argy,VV_SX)

        write(1,1) 'Vertex_Gluon:           ',' SQ= ',VV_SX(1),' SS= ',VV_SX(2),' SV= ',VV_SX(3),' ST= ',VV_SX(4),' SU= ',VV_SX(5)

        call VV_SX_(3,argx,argy,VV_SX)

        write(1,1) 'Vertex_Gluino:          ',' SQ= ',VV_SX(1),' SS= ',VV_SX(2),' SV= ',VV_SX(3),' ST= ',VV_SX(4),' SU= ',VV_SX(5)

        call VV_SX_(4,argx,argy,VV_SX)

        write(1,1) 'Vertex_ScalarLoop:      ',' SQ= ',VV_SX(1),' SS= ',VV_SX(2),' SV= ',VV_SX(3),' ST= ',VV_SX(4),' SU= ',VV_SX(5)

        call VV_SX_(5,argx,argy,VV_SX)

        write(1,1) 'Vertex_Counterterm:     ',' SQ= ',VV_SX(1),' SS= ',VV_SX(2),' SV= ',VV_SX(3),' ST= ',VV_SX(4),' SU= ',VV_SX(5)

        write(1,*) '*****************************************************************************************************************'

        call VV_VX_(1,argx,argy,VV_VX)    
      
        write(1,1) 'Tree:                   ',' VQ= ',VV_VX(1),' VS= ',VV_VX(2),' VV= ',VV_VX(3),' VT= ',VV_VX(4),' VU= ',VV_VX(5)

        call VV_VX_(2,argx,argy,VV_VX)

        write(1,1) 'Vertex_Gluon:           ',' VQ= ',VV_VX(1),' VS= ',VV_VX(2),' VV= ',VV_VX(3),' VT= ',VV_VX(4),' VU= ',VV_VX(5)

        call VV_VX_(3,argx,argy,VV_VX)

        write(1,1) 'Vertex_Gluino:          ',' VQ= ',VV_VX(1),' VS= ',VV_VX(2),' VV= ',VV_VX(3),' VT= ',VV_VX(4),' VU= ',VV_VX(5)

        call VV_VX_(4,argx,argy,VV_VX)

        write(1,1) 'Vertex_SSSSCoupl:       ',' VQ= ',VV_VX(1),' VS= ',VV_VX(2),' VV= ',VV_VX(3),' VT= ',VV_VX(4),' VU= ',VV_VX(5)

        call VV_VX_(5,argx,argy,VV_VX)

        write(1,1) 'Vertex_SSVVCoupl_p1:    ',' VQ= ',VV_VX(1),' VS= ',VV_VX(2),' VV= ',VV_VX(3),' VT= ',VV_VX(4),' VU= ',VV_VX(5)

        call VV_VX_(6,argx,argy,VV_VX)

        write(1,1) 'Vertex_SSVVCoupl_p2:    ',' VQ= ',VV_VX(1),' VS= ',VV_VX(2),' VV= ',VV_VX(3),' VT= ',VV_VX(4),' VU= ',VV_VX(5)


        call VV_VX_(7,argx,argy,VV_VX)

        write(1,1) 'Vertex_Counterterm:     ',' VQ= ',VV_VX(1),' VS= ',VV_VX(2),' VV= ',VV_VX(3),' VT= ',VV_VX(4),' VU= ',VV_VX(5)

        write(1,*) '*****************************************************************************************************************'

        call VV_TX_(1,1,1,argx,argy,VV_TX)    
      
        write(1,1) 'Tree:                   ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,1,2,argx,argy,VV_TX)

        write(1,1) 'Prop_Gluon:             ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)
  
        call VV_TX_(1,1,3,argx,argy,VV_TX)

        write(1,1) 'Prop_Gluino:            ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,1,4,argx,argy,VV_TX)

        write(1,1) 'Prop_Scalar:            ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,1,5,argx,argy,VV_TX)

        write(1,1) 'Prop_Vector:            ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,1,6,argx,argy,VV_TX)

        write(1,1) 'Prop_Counter:           ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(2,1,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_Gluon_Atype:     ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,2,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_Gluon_Btype:     ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(3,1,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_Gluino_Atype:    ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,3,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_Gluino_Btype:    ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(4,1,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_SSSS_Atype:      ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,4,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_SSSS_Btype:      ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(5,1,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_SSVVp1_Atype:    ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,5,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_SSVVprop_Btype:  ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(6,1,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_SSVVprop_Atype:  ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,6,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_SSVVp2_Btype:    ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(7,1,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_Counter_Atype:   ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call VV_TX_(1,7,1,argx,argy,VV_TX)

        write(1,1) 'Vertex_Counter_Btype:   ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call QQ2xx_M2GluonBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_Gluon:              ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call QQ2xx_M2GluinoBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_Gluino:             ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call QQ2xx_M2GluonTriangleBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_GluonTriangle:      ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call QQ2xx_M2GluonsQuarkLoopBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_GluonsQuarkLoopBox: ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        call QQ2xx_M2sQuarkTriangleBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_sQuarkTriangle:     ',' TQ= ',VV_TX(1),' TS= ',VV_TX(2),' TV= ',VV_TX(3),' TT= ',VV_TX(4),' TU= ',VV_TX(5)

        write(1,*) '*****************************************************************************************************************'

        call VV_UX_(1,1,1,argx,argy,VV_UX)    
      
        write(1,1) 'Tree:                   ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(1,1,2,argx,argy,VV_UX)

        write(1,1) 'Prop_Gluon:             ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)
  
        call VV_UX_(1,1,3,argx,argy,VV_UX)

        write(1,1) 'Prop_Gluino:            ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(1,1,4,argx,argy,VV_UX)

        write(1,1) 'Prop_Scalar:            ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(1,1,5,argx,argy,VV_UX)

        write(1,1) 'Prop_Vector:            ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(1,1,6,argx,argy,VV_UX)

        write(1,1) 'Prop_Counter:           ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(2,1,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_Gluon_Atype:     ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(1,2,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_Gluon_Btype:     ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(3,1,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_Gluino_Atype:    ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(1,3,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_Gluino_Btype:    ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(4,1,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_SSSS_Atype:      ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(1,4,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_SSSS_Btype:      ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(5,1,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_SSVVp1_Atype:    ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' TU= ',VV_UX(5)

        call VV_UX_(1,5,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_SSVVprop_Btype:  ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' TU= ',VV_UX(5)

        call VV_UX_(6,1,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_SSVVprop_Atype:  ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' TU= ',VV_UX(5)

        call VV_UX_(1,6,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_SSVVp2_Btype:    ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' TU= ',VV_UX(5)

        call VV_UX_(7,1,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_Counter_Atype:   ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call VV_UX_(1,7,1,argx,argy,VV_UX)

        write(1,1) 'Vertex_Counter_Btype:   ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call QQ2xx_M2GluonBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_Gluon:              ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call QQ2xx_M2GluinoBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_Gluino:             ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call QQ2xx_M2GluonTriangleBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_Triangle:           ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call QQ2xx_M2GluonsQuarkLoopBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_GluonsQuarkLoopBox: ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

        call QQ2xx_M2sQuarkTriangleBoxVV(argx,argy,VV_TX,VV_UX)

        write(1,1) 'Box_sQuarkTriangle:     ',' UQ= ',VV_UX(1),' US= ',VV_UX(2),' UV= ',VV_UX(3),' UT= ',VV_UX(4),' UU= ',VV_UX(5)

  1   Format( A24,  5(A5,E13.4))

      else

           write(*,*)  'Error during writing into AmpDigest.txt'

      end if

      close(unit=1)

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_M2AmpDigestFFbar(f3_type,f3_gen,f4_type,f4_gen)

      implicit none

      integer io_error,f3_type,f3_gen,f4_type,f4_gen
      double precision FFbar_SX(6),FFbar_VX(6),FFbar_TX(6)


      open(unit=1,file='AmpDigest.txt',status='old',action='write',iostat=io_error,position='append')

      if (io_error.eq.0) then
        write(1,*) '*****************************************************************************************************************'
        write(1,*) '*****************************************************************************************************************'

        write(1,*) 'f3type= ',f3_type,'f3gen= ',f3_gen,'f4type= ',f4_type,'f4gen= ',f4_gen

        write(1,*) '*****************************************************************************************************************'

        call FFbar_SX_(1,f3_type,f3_gen,f4_type,f4_gen,FFbar_SX)    
      
        write(1,1) 'Tree:                   ',' SQ= ',FFbar_SX(1),' SS= ',FFbar_SX(2),' SV= ',FFbar_SX(3),' ST= ',FFbar_SX(4),' SU= ',FFbar_SX(5)

        call FFbar_SX_(2,f3_type,f3_gen,f4_type,f4_gen,FFbar_SX)

        write(1,1) 'Vertex_Gluon:           ',' SQ= ',FFbar_SX(1),' SS= ',FFbar_SX(2),' SV= ',FFbar_SX(3),' ST= ',FFbar_SX(4),' SU= ',FFbar_SX(5)

        call FFbar_SX_(3,f3_type,f3_gen,f4_type,f4_gen,FFbar_SX)

        write(1,1) 'Vertex_Gluino:          ',' SQ= ',FFbar_SX(1),' SS= ',FFbar_SX(2),' SV= ',FFbar_SX(3),' ST= ',FFbar_SX(4),' SU= ',FFbar_SX(5)

        call FFbar_SX_(4,f3_type,f3_gen,f4_type,f4_gen,FFbar_SX)

        write(1,1) 'Vertex_ScalarLoop:      ',' SQ= ',FFbar_SX(1),' SS= ',FFbar_SX(2),' SV= ',FFbar_SX(3),' ST= ',FFbar_SX(4),' SU= ',FFbar_SX(5)

        call FFbar_SX_(5,f3_type,f3_gen,f4_type,f4_gen,FFbar_SX)

        write(1,1) 'Vertex_Counterterm:     ',' SQ= ',FFbar_SX(1),' SS= ',FFbar_SX(2),' SV= ',FFbar_SX(3),' ST= ',FFbar_SX(4),' SU= ',FFbar_SX(5)

        write(1,*) '*****************************************************************************************************************'

        call FFbar_VX_(1,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)    
      
        write(1,1) 'Tree:                   ',' VQ= ',FFbar_VX(1),' VS= ',FFbar_VX(2),' VV= ',FFbar_VX(3),' VT= ',FFbar_VX(4),' VU= ',FFbar_VX(5)

        call FFbar_VX_(2,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)

        write(1,1) 'Vertex_Gluon:           ',' VQ= ',FFbar_VX(1),' VS= ',FFbar_VX(2),' VV= ',FFbar_VX(3),' VT= ',FFbar_VX(4),' VU= ',FFbar_VX(5)

        call FFbar_VX_(3,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)

        write(1,1) 'Vertex_Gluino:          ',' VQ= ',FFbar_VX(1),' VS= ',FFbar_VX(2),' VV= ',FFbar_VX(3),' VT= ',FFbar_VX(4),' VU= ',FFbar_VX(5)

        call FFbar_VX_(4,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)

        write(1,1) 'Vertex_SSSSCoupl:       ',' VQ= ',FFbar_VX(1),' VS= ',FFbar_VX(2),' VV= ',FFbar_VX(3),' VT= ',FFbar_VX(4),' VU= ',FFbar_VX(5)

        call FFbar_VX_(5,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)

        write(1,1) 'Vertex_SSVVCoupl_p1:    ',' VQ= ',FFbar_VX(1),' VS= ',FFbar_VX(2),' VV= ',FFbar_VX(3),' VT= ',FFbar_VX(4),' VU= ',FFbar_VX(5)

        call FFbar_VX_(6,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)

        write(1,1) 'Vertex_SSVVCoupl_p2:    ',' VQ= ',FFbar_VX(1),' VS= ',FFbar_VX(2),' VV= ',FFbar_VX(3),' VT= ',FFbar_VX(4),' VU= ',FFbar_VX(5)

        call FFbar_VX_(7,f3_type,f3_gen,f4_type,f4_gen,FFbar_VX)

        write(1,1) 'Vertex_Counterterm:     ',' VQ= ',FFbar_VX(1),' VS= ',FFbar_VX(2),' VV= ',FFbar_VX(3),' VT= ',FFbar_VX(4),' VU= ',FFbar_VX(5)

        write(1,*) '*****************************************************************************************************************'

        call FFbar_TX_(f3_type,f3_gen,f4_type,f4_gen,FFbar_TX)    
      
        write(1,1) 'Tree:                   ',' TQ= ',FFbar_TX(1),' TS= ',FFbar_TX(2),' TV= ',FFbar_TX(3),' TT= ',FFbar_TX(4),' TU= ',FFbar_TX(5)
        
  1   Format( A24,  5(A5,E13.4))

      else

           write(*,*)  'Error during writing into AmpDigest.txt'

      end if

      close(unit=1)

      end


!**********************************************************************************************************************************
!**********************************************************************************************************************************



      subroutine QQ2xx_M2AmpDigestFF(f3_type,f3_gen,f4_type,f4_gen)

      implicit none

      integer io_error,f3_type,f3_gen,f4_type,f4_gen
      double precision FF_TX(6),FF_UX(6)


      open(unit=1,file='AmpDigest.txt',status='old',action='write',iostat=io_error,position='append')

      if (io_error.eq.0) then
        write(1,*) '*****************************************************************************************************************'
        write(1,*) '*****************************************************************************************************************'

        call FF_TX_(f3_type,f3_gen,f4_type,f4_gen,FF_TX)    
      
        write(1,1) 'Tree:                   ',' TQ= ',FF_TX(1),' TS= ',FF_TX(2),' TV= ',FF_TX(3),' TT= ',FF_TX(4),' TU= ',FF_TX(5)

        call FF_UX_(f3_type,f3_gen,f4_type,f4_gen,FF_UX)    
      
        write(1,1) 'Tree:                   ',' UQ= ',FF_UX(1),' US= ',FF_UX(2),' UV= ',FF_UX(3),' UT= ',FF_UX(4),' UU= ',FF_UX(5)
        
  1   Format( A24,  5(A5,E13.4))

      else

           write(*,*)  'Error during writing into AmpDigest.txt'

      end if

      close(unit=1)

      end