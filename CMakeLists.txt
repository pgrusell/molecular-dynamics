cmake_minimum_required(VERSION 3.12)
project(MolecularDinamics LANGUAGES Fortran)

# Establecemos la configuración de Fortran
enable_language(Fortran)

# Agregamos los archivos fuente al ejecutable
add_executable(crea_red
	src/decl_precision.f90
	src/main.f90
	src/configura.f90
	src/mi_random.f90
	src/pot_lj.f90
	src/variables_comunes.f90
	src/velocity_verlet.f90
)

# Opciones de compilación
set(CMAKE_Fortran_FLAGS "-O2 -Wall")

# Especificamos el directorio de salida del ejecutable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)
