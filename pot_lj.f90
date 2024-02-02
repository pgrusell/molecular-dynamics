






    subroutine pot_lj(rx,ry,rz,ax,ay,az,epot,depotr,ddepotrr)

        use variables_comunes

        use def_precision

        implicit none

        real (kind=doblep), dimension(N):: rx,ry,rz,ax,ay,az
        integer (kind=entero):: i, j
        real(kind=doblep):: rxi,ryi,rzi, rxij, ryij, rzij, d2, d2_inv, d6_inv, d12_inv, epot, fmod, depotr, ddepotrr, factor

!       ===============================================================================================================
!       Puesta a 0 de los contadores
!       ===============================================================================================================

        rxi=0.d00
        ryi=0.d00
        rzi=0.d00
        rxij=0.d00
        ryij=0.d00
        rzij=0.d00
        epot=0.d00
        depotr=0.d00
        ddepotrr=0.d00
        ax=0d00
        ay=0d00
        az=0d00

!       ===============================================================================================================
!       Calculamos la energía potencial y las fuerzas
!       ===============================================================================================================
        do i=1,N-1

            rxi=rx(i)
            ryi=ry(i)
            rzi=rz(i)

            do j=i+1,N

                rxij=rxi-rx(j)
                ryij=ryi-ry(j)
                rzij=rzi-rz(j)

!               Calculamos rij como la distancia entre la partícula i y la imagen más próxima de la partícula j

                rxij=rxij-L*dnint(rxij*L_inv)
                ryij=ryij-L*dnint(ryij*L_inv)
                rzij=rzij-L*dnint(rzij*L_inv)

                d2=rxij*rxij+ryij*ryij+rzij*rzij

                d2_inv=1/d2
                d6_inv=d2_inv*d2_inv*d2_inv
                d12_inv=d6_inv*d6_inv

!               Si la imagen de la partícula j más próxima está dentro del radio de corte estas interactúan y calculamos
!               la energía potencial y la fuerza (aceleración, m=1) de la partícula i.

                if (d2<=rc2) then
!
!                       Sumamos a la energía potencial, la de la partícula

                        factor=-2*d12_inv+d6_inv
                        epot=epot+d12_inv-d6_inv
                        depotr=depotr+factor
                        ddepotrr=26*d12_inv-7*d6_inv
                        fmod=(-factor)*d2_inv

!                       Calculamos la fuerza

                        ax(i)=ax(i)+fmod*rxij
                        ax(j)=ax(j)-fmod*rxij
!
                        ay(i)=ay(i)+fmod*ryij
                        ay(j)=ay(j)-fmod*ryij
!
                        az(i)=az(i)+fmod*rzij
                        az(j)=az(j)-fmod*rzij
!
                end if
            end do
        end do

        epot=4*epot
        depotr=24*depotr
        ddepotrr=24*ddepotrr
        ax=24*ax
        ay=24*ay
        az=24*az

!       ===============================================================================================================
!       Añadimos la corrección a la energía
!       ===============================================================================================================

        epot=epot+ecorr
        depotr=depotr+depotr_corr
        ddepotrr=ddepotrr+ddepotrr_corr

        return

    end subroutine
