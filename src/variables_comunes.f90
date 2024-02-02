module variables_comunes

    use def_precision

    implicit none

    integer (kind=entero),parameter:: N=500
    real (kind=doblep):: L, V, L_inv, etot_
    real (kind=doblep), parameter:: pi=3.14159265358970
    integer (kind=entero):: k, Npasos
    real (kind=doblep):: a, rc, rc2, rc3, rc6, deltat, ecorr, depotr_corr, ddepotrr_corr

!   ===================================================================================================================
!   GLOSARIO (variables comunes)
!   ===================================================================================================================

!   L: Lado de la caja. En nuestro caso 10 en unidades reducidas.
!   V: Volumen del cubo.
!   L_inv: 1/L, se utiliza para calculos intermedios.
!   etot_: Energia total que queremos que tenga nuestro sistema. En general es distinta de etot=ecin+epot.
!   pi: Valor tipico de pi, empleado para calculos intermedios.
!   N: Numero total de particulas. En nuestro caso 500 particulas.
!   k: Numero de repeticiones de la celda unidad para construir la caja completa. Se calcula como
!   ceiling((N/4.d00)**(1.d00/3.d00),entero). Que sea entero impone que N tenga que valer, 4, 32, 108, 256, 500, 864...
!   a: Longitud de la celda unidad, calculada como L/k.
!   rc: Radio de corte de nuestro potencial, en nuestro caso vale L/2, este es el valor maximo posible, ya que, en caso
!   contrario podrían interactuar particulas de celdas virtuales.
!   rc2, rc3, rc6: Valores intermedios para los calculos.
!   deltat: Distancia temporal entre las dos configuraciones mas cercanas, en nuestro caso 10**-4
!   ecorr, depotr_corr, ddepotrr_corr: Correcciones a la energia potencial y sus derivadas respecto del volumen.
!   ===================================================================================================================

end module variables_comunes