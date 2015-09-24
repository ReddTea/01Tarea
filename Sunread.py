##########################
#Importar paquetes a usar#
##########################
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as unit

##############################################################################
#Arreglo de los datos, nótese que se pasa la longitud de onda a micrones y el#
#flujo a cgs, que según wikipedia en radiación es (erg·cm−2·s−1·um-1)             #
##############################################################################
datos=np.loadtxt("sun_AM0.dat")
wavelength=(datos[:,0]*unit.nm).to('um')
flujo=(datos[:,1]*unit.W*unit.m**(-2)*unit.nm**(-1)).to('erg/(s cm2 um)')

#######################
#se plotea la cuestión#
#######################

plt.plot(wavelength,flujo)
plt.xlabel('Longitud de onda [$um$]',fontsize=18)
plt.xlim(0,8) #se agrega límite
plt.ylabel('Flujo [$\\frac{erg}{s*cm^2*um^-1}$]',fontsize=18)
plt.title('Radiacion de Cuerpo Negro del Sol',fontsize=22)
plt.savefig('Spectre.png')
plt.show()

###################################################################
###################################################################
########################### FIN PARTE 1 ###########################
###################################################################
###################################################################

##########################################################################
#Nótese que con el método de los trapecios, si el ancho fuera h, la parte#
#cuadrada es y(i)*h y la triangular h*(y(i+1)-y(i))/2, sumando se tiene  #
#que cuadrado+triangulo=h*(y(i)+y(i+1))/2, como h= x(i+1)-x(i)           #
#queda que el trapecio es (y(i)+y(i+1))*(x(i+1)-x(i))/2                  #
##########################################################################



suma=0 #this! Me demoré demasiado en notar que se pone afuera del for u.u
for i in range(len(wavelength)-1):
    trapecitos= (wavelength[i+1]-wavelength[i])*(flujo[i]+flujo[i+1])/2.0
    suma+=trapecitos
print suma

###################################################################
###################################################################
########################### FIN PARTE 2 ###########################
###################################################################
###################################################################
