






    subroutine velocity_verlet(rx,ry,rz,vx,vy,vz,ax,ay,az,ecin,epot,etot,depotr,ddepotrr)

        use variables_comunes

        use def_precision

        implicit none

        real (kind=doblep), dimension (N), intent(inout):: rx, ry, rz, vx, vy, vz, ax, ay, az
        real (kind=doblep)::  deltat_cm, deltat_m, ecin, epot, etot, depotr, ddepotrr

        deltat_cm=0.5*deltat*deltat
        deltat_m=0.5*deltat

        rx=rx+vx*deltat+ax*deltat_cm
        ry=ry+vy*deltat+ay*deltat_cm
        rz=rz+vz*deltat+az*deltat_cm

        vx=vx+deltat_m*ax
        vy=vy+deltat_m*ay
        vz=vz+deltat_m*az


        call pot_lj(rx,ry,rz,ax,ay,az,epot,depotr,ddepotrr)

        vx=vx+deltat_m*ax
        vy=vy+deltat_m*ay
        vz=vz+deltat_m*az

        ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

        etot=ecin+epot

        return

    end subroutine
