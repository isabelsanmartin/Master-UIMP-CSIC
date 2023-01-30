### Construcción de mapas básicos de distribución en R
### Script v22  (created by R. Riina 29 Jan 2022; updated 30 Jan 2023)

### Vamos a explorar algunas de las opciones para hacer mapas disponibles en R
### Para ello, necesitamos instalar algunos paquetes y sus dependencias. 

### Instalar paquetes más importantes
### y sud dependencias

install.packages(c("rworldmap", "sp", "dismo", "maps", "prettymapr"), dependencies=TRUE) 


### Usaremos estas "libraries" y otras más indicadas más adelante

library(sp)  # classes for spatial data
library(rworldmap) # simple maps
library(maps) # para añadir la escala
library(prettymapr) # para añadir detalles (e.g., el norte)
library(dismo)


############################################################################
####  Ejercicio 1:  Crear un mapa de distribución de Euphorbia terracina ###
####  con datos de GBIF                                                  ###
############################################################################

### Comenzamos por crear un mapa general del mundo

mymap <- getMap(resolution = "low")
plot(mymap)

### opciones de resolución: "coarse","low","less islands","li","high". 
### para usar "high" necesitariamos instalar el paquete rworldxtra

### Ahora vamos a obtener los datos de la especie de interés del portal GBIF 
### y almacenarlos en la tabla o dataframe que llamaremos 'myspecies'

myspecies <- gbif("Euphorbia", "terracina")

### Usamos "head" para ver las primeras 6 filas y las columnas de nuestra  
### tabla 'myspecies'; veremos que hay muchas columnas que no nos interesan

head(myspecies)  

### Creamos un nuevo dataframe 'localidades' que contendrá solo las columnas  
### “lat”, “lon”, y “country” 

localidades <- subset(myspecies, select=c("country", "lat", "lon"))

head(localidades)  ### para ver las primeras filas de la nueva tabla 'localidades'

str(localidades)   ### para ver su estructura  

### Como hay algunas filas que contienen valores nulos, las eliminamos   
### y creamos una nueva tabla que llamaremos 'locs'
### la función complete.cases la usamos para excluir las filas con datos nulos

locs <- localidades[complete.cases(localidades), ]

head(locs)  ### para ver las primeras filas de la nueva tabla 'locs'

str(locs)  ### para ver su estructura 

### Para que se reconozcan como coordenadas y no solo como variables, 
### indicamos qué columnas ("lon", "lat") corresponden a la longitud  
### y latitud de nuestros puntos

coordinates(locs) <- c("lon", "lat")

### Añadimos la proyección espacial a nuestros datos
### 'crs' stands for 'coordinate reference system'
### El sistema más común es el WGS84 (World Geodesic System 1984) 

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  

### Definimos el sistema de proyección de datos 

proj4string(locs) <- crs.geo  

### Vemos un resumen de la tabla 'locs'

summary(locs)

### Dibujamos los puntos

plot(locs)

### Ponemos juntos puntos y mapa y modificamos atributos: tamaño de punto (cex), 
### color (col), añadimos título (main), etc.
### ver símbolos (pch) alternativos en
### http://www.statmethods.net/advgraphs/images/points.png
### Podemos cambiar el color (col) de los puntos indicando el nombre (blue, dark green)
### cex define el tamaño de los puntos en el mapa

plot(locs, pch=20, cex = 1, col="dark green", main="Localidades - Euphorbia terracina")

### Dibujamos mapa/bordes

plot(coastsCoarse, add =TRUE)

### Observa el mapa obtenido y trata de indentificar posibles errores. 
### Info: Euphorbia terracina es una especie nativa de la Región Mediterránea, 
### pero se considera invasora fuera de su área nativa 
### (pista: invade otras regiones mediterráneas del globo)
### ¿Cuáles puntos crees deberían ser investigados y verificados?


###########################################################################
### Ejercício 1.2: Construir mapa de la especie en un país o área de    ###
### de interés (España)                                                 ###
###########################################################################

### Usamos la función 'table'
###  '$' permite extraer elementos por nombre de una lista

table(locs$country)

locs.spain <- subset(locs, locs$country == "Spain")  

plot(locs.spain, pch = 20, cex = 1, col = "red")

title("Euphorbia terracina en España")

plot(countriesLow, add = T)

### En los próximos ejercícios haremos los mapas con todos los atributos 
### necesarios para publicar (e.g., artículo o website).

############################################################################
### Ejercício 1.3: Crear un mapa de una parte del país (Islas Baleares) ####
############################################################################

### Definimos los límites espaciales de nuestro mapa usando xlim, ylim

library(maps)    

plot(locs, pch = 20, col = "steelblue", xlim = c(1, 4.7), ylim = c(38.5,40))

title("Euphorbia terracina, Islas Baleares")

plot(countriesLow, add = T)

### añadimos los ejes de coordenadas y borde del mapa

map.axes(cex.axis=1)

### añadimos el norte y ecala en Km

library(prettymapr)

addnortharrow(pos="bottomright", scale=0.8, padin = c(0.5,0.3))

map.scale(x=3, y=38.6, metric = TRUE, ratio=FALSE, relwidth=0.15, cex=0.8) 

### “pos” posición del símbolo del Norte en el mapa
### “scale” define el tamaño del símbolo del Norte (prueba con 1, 0.7, 3, etc.) 

### Observa el mapa obtenido y trata de indentificar posibles puntos erróneos. 


#####################################################################################
###  Ejercicio 2: Crear un mapa con datos de un muestreo                          ###
###               realizado en 2018 de la especie Ricinus communis                ###
#####################################################################################
### En este caso, no usaremos datos de BBIF sino datos de nuestro propio muestreo ###  
### de campo, que tenemos en el archivo “Recolecciones_Ricinus.txt”               ###

library(dismo)
library(rworldmap)
library(prettymapr)
library(maps)

### Establecer directorio (el de este script y "Recolecciones_Ricinus.txt”)
### En mi caso tengo una carpeta en mi escritorio ("Ejercicio_Mapas")

setwd("~/Desktop/Ejercicio_Mapas")

tabla<-read.table("Recolecciones_Ricinus.txt", sep="\t", header=TRUE)

### Ver el contenido de 'tabla'
tabla 

### Crear un nuevo dataframe; un subset de 'tabla'

locs<-subset(tabla, select=c("lat_dec", "long_dec", "BOLSA_CODE"))

coordinates(locs) <- c("lat_dec", "long_dec")

### Proyección geográfica, datum WGS84

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 

### Define el sistema de proyección 

proj4string(locs) <- crs.geo                           

summary(locs)

plot(locs)

newmap <- getMap(resolution = "low") 

### En los ejemplos arriba usamos "low"
### para usar "high" necesitariamos instalar el paquete rworldxtra

plot(newmap, xlim=c(-13, 5), ylim = c(37, 39), asp = 1) 

title("Muestreo de Ricinus communis en 2018")

points(tabla$long_dec, tabla$lat_dec, col = "blue", cex = 1.5, pch=20)

map.axes(cex.axis=0.5) # añadimos las coordenadas al mapa

### añadimos símbolo Norte

addnortharrow(pos="bottomright", scale=0.5, padin = c(0.1,0.5)) 

### añadimos escala km

map.scale(x=-0, y=34, metric = TRUE, ratio=FALSE, relwidth=0.15, cex=0.5) 

### Puedes exportar el mapa (mejor como pdf, seleccionar "portrait")  
### usando "Export/Exportar" desde la ventana "Plots" de R Studio.  
### El pdf se puede editar y mejorar en programas 
### gráficos como Inkscape (añadir fotos de organismo, etc.)