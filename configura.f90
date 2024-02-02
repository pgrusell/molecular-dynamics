subroutine configura()

    use variables_comunes

    implicit none

!   Cargamos los valores de las variables del archivo
    open(10, file='../src/configuracion.txt', status='old', action='read')
	read(10, *) L, etot_
    close(10)

!   Calculamos el resto de variables
    V = L*L*L
    L_inv = 1.d00/L

end subroutine