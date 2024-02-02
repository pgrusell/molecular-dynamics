# Creamos el directorio build
mkdir build

# Cambiamos al directorio build
cd build

# Ejecutamos cmake para configurar el proyecto
cmake ..

# Ejecutamos make para compilar el proyecto
make

# Creamos el directorio donde almacenaremos los resultados
cd .. 
mkdir resultados

# Abrimos el ejecutable
cd build
./crea_red


