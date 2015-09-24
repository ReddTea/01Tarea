##########################
#Importar paquetes a usar#
##########################
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as cte # se agrega esto
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
'''
plt.plot(wavelength,flujo)
plt.xlabel('Longitud de onda [$um$]',fontsize=18)
plt.xlim(0,8) #se agrega límite
plt.ylabel('Flujo [$\\frac{erg}{s*cm^2*um^-1}$]',fontsize=18)
plt.title('Radiacion de Cuerpo Negro del Sol',fontsize=22)
plt.savefig('Spectre.png')
plt.show()
'''
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


#'''
energia_total=0 #this! Me demoré demasiado en notar que se pone afuera del for u.u
for i in range(len(wavelength)-1):
    trapecitos= (wavelength[i+1]-wavelength[i])*(flujo[i]+flujo[i+1])/2.0
    energia_total+=trapecitos
print 'energia total='+str(energia_total)

#la luminosidad total corresponde a la superficie de la esfera x su flujo

luminosidad_total= 4*np.pi*energia_total*(cte.au).to('cm')**2
#print 'luminosidad_total'
#print(luminosidad_total)
estabien=luminosidad_total.to('W')-cte.L_sun #nótese que la diferencia es de 3 órdenes de magnitud menor...
#'''
###################################################################
###################################################################
########################### FIN PARTE 2 ###########################
###################################################################
###################################################################

##########################################################################
#Calculamos la integral que nos dan en la tarea,con el cambio de variable#
#  y=arctan(x) y usamos el método de Simpson                             #
#parametros iniciales y vectores sobre los que integrar                  #
##########################################################################

a=0.01
b=np.pi/2.0

planx=np.arange(a,b,a)
blackbody_function=(np.tan(planx)**3)/((np.cos(planx)**2)*(np.exp(np.tan(planx))-1))
#|#|#|#|#|#| ERROR DE #$&%€ ACA!!!!

#Integramos usando el método del trapecio
blackbody_sum=0

for i in range(len(planx)-1):
    x2=planx[i+1]
    x1=planx[i]
    y2=blackbody_function[i+1]
    y1=blackbody_function[i]
    kkk=((x2-x1)/2)*(y2+y1)
    blackbody_sum+=kkk

###################VERIFICANDO SI ESTÁ BIEN HASTA AHORA########################
#print 'blackbody_sum= '+str(print blackbody_sum)                             #
#estabien2= np.pi**4/15.0 -blackbody_sum #el error es de 7 magnitudes menor!!#
##################################LO ESTÁ #####################################

#ahora la integral falta multiplicarla por las constantes omitidas y convertirla
#a unidades adecuadas

megaconstante=2*np.pi*((cte.k_B*5778*unit.K)**4)/((cte.c**2)*(cte.h**3))
bbrcgs=(blackbody_sum*megaconstante).to('erg / cm2 s') #radiacion en cgs

#cálculo del radio solar
radio_solar=((energia_total/bbrcgs)*(cte.au).to('cm')**2)**0.5
print 'radio_solar en cm= '+ str(radio_solar)

estabien3=radio_solar-cte.R_sun #error de 3 km... despreciable!

###################################################################
###################################################################
########################### FIN PARTE 3 ###########################
###################################################################
###################################################################
