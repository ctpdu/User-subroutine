      subroutine uanisohyper_inv (ainv, ua, zeta, nfibers, ninv, 
     &                            ui1, ui2, ui3, temp, noel, cmname, 
     &                            incmpflag, ihybflag, numstatev, 
     &                            statev, numfieldv, fieldv, 
     &                            fieldvinc, numprops, props)
C
      include 'aba_param.inc'
C
      character*80 cmname
      dimension ua(2), ainv(ninv), ui1(ninv),
     &          ui2(ninv*(ninv+1)/2), ui3(ninv*(ninv+1)/2),
     &          statev(numstatev), fieldv(numfieldv),
     &          fieldvinc(numfieldv), props(numprops)
C
C     ainv  : invariants
C     ua    : energies
C     ua(1) : utot, 
C     ua(2) : dev
C     ui1   : dU/dI
C     ui2   : d2U/dIdJ
C
      parameter ( zero  = 0.d0, 
     &            one   = 1.d0, 
     &            two   = 2.d0, 
     &            three = 3.d0,
     &            four  = 4.d0, 
     &            five  = 5.d0, 
     &            six   = 6.d0, 
C
     &            index_I1 = 1, 
     &            index_J  = 3 )
C
C--------------------neo-Hookean model-----------------------
C
      C10 = props(1)
      D1  = props(2)
C
C     Deviatoric energy
C
      ua(2) = zero
      ua(2) = ua(2) + C10 * (ainv(index_i1) - three)
C
C     Compute derivatives
C
      ui1(index_i1) =  C10
      ui2(indx(index_I1,index_I1))= zero
C     
C     Consider compressibility
C
      if(D1.gt.zero) then
         Dinv = one / D1
         det = ainv(index_J)
         DinJ = det - one
         ua(1) = ua(2) + Dinv * DinJ * DinJ
         ui1(index_J) = two * Dinv * DinJ
         ui2(indx(index_J,index_J))= two * Dinv
         if (hybflag.eq.1) then
           ui3(indx(index_J,index_J))= zero
         end if
      end if
C
      return
      end
C-------------------------------------------------------------
C     Function to map index from Square to Triangular storage 
C 		 of symmetric matrix
C
      integer function indx( i, j )
      include 'aba_param.inc'
      ii = min(i,j)
      jj = max(i,j)
      indx = ii + jj*(jj-1)/2
      return
      end