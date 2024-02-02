





    program crea_red

        use def_precision

	use variables_comunes

        implicit none

        real(kind=doblep):: mi_random, b, epot, ecin, etot, depotr, ddepotrr
        real(kind=doblep):: ecin_, factor, sx, sy, sz, tiempo, pref
        character(40):: date, time, zone, fname1, fname2,fname3
        integer(kind=entero):: i, j, m, p1, p2, iduma, Npasos
        real (kind=doblep), dimension(N):: rx, ry, rz, vx, vy, vz, ax, ay, az
        integer, dimension(8)::values

	call configura()

!       ===================================================================================================================
!       GLOSARIO (variables del archivo principal)
!       ===================================================================================================================
!       mi_random: Funcion que va a devolver, a partir de un número entero, un numero real en [0,1)
!       b: Coeficiente que controla cuanto vamos a descolocar las particulas.
!       epot: Energia potencial del conjunto de partículas
!       ecin: Energia cinetica de las particulas (calculada a partir de las velocidades)
!       etot: Energia total del sistema (calculada como etot=ecin+epot)
!       ecin_: Energia cinetica necesaria para que se conserve la energia (diferencia entre epot y etot_, ver en el glosario
!              de variables_comunes).
!       factor: Raiz cuadrada del cociente entre ecin_ y ecin (sirve para normalizar la energia cinetica al valor deseado).
!       sx, sy, sz: Contadores de la suma de sendas componentes del momento lineal (sirve para poder fijar a cero el momento
!               lineal total).
!       date, time, zone: argumentos necesarios para poder obtener la hora del ordenador.
!       fname1, 2 y 3: Variable que almacena el nombre de los archivos.
!       i, j: Contadores genericos
!       m: Combinacion lineal de i y j que nos evita repetir dicha expresión varias veces.
!       p1 y p2: Coeficientes necesarios para la colocación de las posiciones en la fcc, (discutidos en el PDF).
!       iduma: Numero entero que usaremos para inicializar mi_random.
!       Npasos: Numero de pasos que vamos a dejar correr la simulacion, en esta su valor es de 5000.
!       tiempo: Producto de deltat*i (donde i es el número del paso en el que esta).
!       rx, ry, rz: Vectores que contienen las posiciones de las particulas. Inicialmente les vamos a asignar los valores
!               de una red fcc, posteriormente los descolocaremos y finalmente se dejara que velocity_verlet los haga
!               evolucionar.
!       vx, vy, vz: Velocidades de las particulas. Inicialemnte se asignan aleatoriamente y se corregiran para imponer con-
!                   servacion del momento lineal y normalizarán para imponer conservacion de la energia.
!       ax, ay, az: Aceleraciones de las particulas. Son calculadas por la subrutina pot_lj a partir de las posiciones de las
!                   particulas, suponiendo que están sometidas a un potencial LJ truncado en sigma=5.
!       values: vector de ocho componentes que contienen: (1) año, (2) mes, (3) dia, (4) diferencia en minutos respecto a UTC,
!               (5) hora del dia, (6) minutos, (7) segundos, (8) milisegundos. El motivo por el que queremos conocer estos valo-
!               res es, como son relativamente aleatorios, usar su suma, como inicializador de mi_random(). No es necesario,
!               podemos introducir nosotros uno a mano, pero es más automatico asi y genera una configuración inicial aleatoria.
!
!       ===================================================================================================================

!       ===================================================================================================================
!       Introducimos los valores de las variables (algunas de ellas definidas en variables_comunes)
!       ===================================================================================================================

        k=ceiling((N/4.d00)**(1.d00/3.d00),entero)
        a=L/k
        rc=L/2
        rc2=rc**2
        rc3=rc2*rc
        rc6=rc3*rc3
        b=0.25*a
        Npasos=5000
        deltat=1.d-04
        pref=N*N/V*16*pi/(rc3)

        ecorr=pref/6*(1/(3*rc6)-1)
        depotr_corr=pref*(1-2/(3*rc6))
        ddepotrr_corr=pref*(26/(3*rc6)-7)

!       Llamamos a la subrutina date_and_time que nos va a devolver los datos temporales del ordenador.

        call date_and_time(date, time, zone, values)

!       Definimos iduma como la suma de las componentes del vector values

        iduma=sum(values)

!       Y ahora los nombres de los archivos

        fname1='../resultados/datos_DM1.txt'
        fname2='../resultados/datos_DM2.dat'
        fname3='../resultados/datos_DM3.txt'

!       ===================================================================================================================
!       Ponemos a cero todos los contadores
!       ===================================================================================================================
        sx=0.d00
        sy=0.d00
        sz=0.d00

        rx=0d00
        ry=0d00
        rz=0d00

        vx=0d00
        vy=0d00
        vz=0d00

        ax=0d00
        ay=0d00
        az=0d00

        i=0
        j=0
        m=0

!       ===================================================================================================================
!       Creamos la celda unitaria de la fcc
!       ===================================================================================================================

        rx(1)=0
        rx(2)=a/2
        rx(3)=a/2
        rx(4)=0

        ry(1)=0
        ry(2)=a/2
        ry(3)=0
        ry(4)=a/2

        rz(1)=0
        rz(2)=0
        rz(3)=a/2
        rz(4)=a/2


!       Repetimos k-1 veces en cada dirección

        do i=1,(k-1)
            do j=1,4
                m=4*i+j
                rx(m)=rx(j)+a*i
                rz(m)=rz(j)
                ry(m)=ry(j)
            end do
        end do

        p1=m
        do i=1,(k-1)
            do j=1,p1
                m=p1*i+j
                rx(m)=rx(j)
                ry(m)=ry(j)+a*i
                rz(m)=rz(j)
            end do
        end do

        p2=m
        do i=1,(k-1)
            do j=1,p2
                m=p2*i+j
                rx(m)=rx(j)
                ry(m)=ry(j)
                rz(m)=rz(j)+a*i
            end do
        end do

!       Desviamos una cantidad entre -a/4 y a/4 a cada una de las particulas de su posicion de equilibrio
!       Ademas, sumamos a/4 en todas las direcciones para centrar las particulas.

        open(10, file=fname2)

        do i=1,N
            rx(i)=rx(i)+b*(0.5-mi_random(iduma))+a/4
            ry(i)=ry(i)+b*(0.5-mi_random(iduma))+a/4
            rz(i)=rz(i)+b*(0.5-mi_random(iduma))+a/4

            write(10,*) rx(i), ry(i), rz(i)
        end do

        close(10)


!       ===================================================================================================================
!       Calculo de la energia potencial: Llamada a pot_lj
!       ===================================================================================================================

        call pot_lj(rx,ry,rz,ax,ay,az,epot,depotr,ddepotrr)

!       ===================================================================================================================
!       Asignacion de velocidades
!       ===================================================================================================================

!       Asignamos velocidades aleatorias entre -0.5 y 0.5 aleatoriamente

        do i=1,N

            vx(i)=0.5-mi_random(iduma)
            vy(i)=0.5-mi_random(iduma)
            vz(i)=0.5-mi_random(iduma)

        end do

!       Debemos imponer que el momento lineal sea nulo, primer sumamos el momento en cada componente

        sx=sum(vx)
        sy=sum(vy)
        sz=sum(vz)

!       Corregimos el momento lineal de cada particula para que el momento lineal sea 0

        vx=vx-sx/N
        vy=vy-sy/N
        vz=vz-sz/N

!       Calculamos la nueva energia cinetica

        ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

!       Calculamos la energia cinetica que necesitamos para que etot sea la que queremos y obtenemos la relacion entre
!       la energia cinetica que queremos y la que tenemos.

        ecin_=etot_-epot
        factor=dsqrt(ecin_/ecin)

!       Multiplicamos cada componente por ese factor

        vx=vx*factor
        vy=vy*factor
        vz=vz*factor

!       Y ya tenemos la energia cinetica

        ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))


!       ===================================================================================================================
!       Con la energia potenicial y la cinetica calculamos la energía total
!       ===================================================================================================================

        etot=ecin+epot

!       ===================================================================================================================
!       Evolucion de la simulacion: Npasos(5000)*deltat(10**-4).
!       ===================================================================================================================

        open(10, file=fname1)

        do i=1, Npasos

!           Llamamos a velocity_verlet que calcula, a partir de r's, v's y a's antiguas, las nuevas, ademas de calcular las
!           nuevas energias, etot, ecin y epot.

            call velocity_verlet(rx,ry,rz,vx,vy,vz,ax,ay,az,ecin,epot,etot,depotr,ddepotrr)

!           Grabamos en un archivo cada 100 pasos, el valor del tiempo, el paso y las energias.

            if (mod(i,1)==0) then

                tiempo=dble(i)*dble(deltat)

                print*, i, tiempo, ecin, etot, epot

                write(10,'(e13.6, 2x, e13.6, 2x, e13.6, 2x,e13.6)') tiempo, ecin, etot, epot

            end if

        end do

        close(10)

        open(20, file=fname2, form='unformatted')

            write(20) N, L, V, rc, ecin, etot, epot, rx, ry, rz, vx, vy, vz, ax, ay, az
            write(20) Npasos, deltat, ecorr, depotr_corr, ddepotrr_corr

        close(20)

        open(30, file=fname3)

            write(30, '(I4,2x,I7)') N, Npasos
            write(30,'(9(2x,e13.6))') rx, ry, rz, vx, vy, vz, ax, ay, az
            write(30,'(7(2x,e13.6))') L, V, rc, ecorr, depotr_corr, ddepotrr_corr, deltat
            write(30,'(5(2x,e13.6))') ecin, etot, epot, depotr, ddepotrr

        close(30)

        print*, 'FIN'

    end program


