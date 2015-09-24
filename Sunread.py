##########################
#Importar paquetes a usar#
##########################
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as cte # se agrega esto
from astropy import units as unit
from scipy import integrate as integrate #se agrega esto también
np.seterr(all='ignore')

############################
                           #
kronos=get_ipython().magic #
                           #
############################
#Agregado para la parte 4  #
############################

##############################################################################
#Arreglo de los datos, nótese que se pasa la longitud de onda a micrones y el#
#flujo a cgs, que según wikipedia en radiación es (erg·cm−2·s−1·um-1)        #
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



Integral_1=0 #Me demoré demasiado en notar que se pone afuera del 'for' u.u

for i in range(len(wavelength)-1):
    trapecitos= (wavelength[i+1]-wavelength[i])*(flujo[i]+flujo[i+1])/2.0
    Integral_1+=trapecitos
print 'Integral_1='+str(Integral_1)

#la luminosidad total corresponde a la superficie de la esfera x su flujo

luminosidad_total= 4*np.pi*Integral_1*(cte.au).to('cm')**2
print 'luminosidad_total= '+str(luminosidad_total)
estabien=luminosidad_total.to('W')-cte.L_sun #nótese que la diferencia es de 3 órdenes de magnitud menor...

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
#|#|#|#|#|#| ERROR DE #$&%€ ACA!!!! ### solucionado a la mala en linea9

#Integramos usando el método del trapecio
Integral_2=0

for i in range(len(planx)-1):
    x2=planx[i+1]
    x1=planx[i]
    y2=blackbody_function[i+1]
    y1=blackbody_function[i]
    kkk=((x2-x1)/2)*(y2+y1)
    Integral_2+=kkk

###################VERIFICANDO SI ESTÁ BIEN HASTA AHORA########################
#print 'Integral_2= '+str(print Integral_2)                                   #
#estabien2= np.pi**4/15.0 -Integral_2 #el error es de 7 magnitudes menor!!    #
##################################LO ESTÁ #####################################

#ahora la integral falta multiplicarla por las constantes omitidas y convertirla
# a unidades adecuadas

megaconstante=2*np.pi*((cte.k_B*5778*unit.K)**4)/((cte.c**2)*(cte.h**3))
bbrcgs=(Integral_2*megaconstante).to('erg / cm2 s') #radiacion en cgs
print 'radiacion en cgs= '+str(bbrcgs)

#cálculo del radio solar
radio_solar=((Integral_1/bbrcgs)*(cte.au).to('cm')**2)**0.5
print 'radio_solar en cm= '+ str(radio_solar)

#comprobando si está bien o no
estabien3=radio_solar-cte.R_sun #error de 3 km... despreciable!

###################################################################
###################################################################
########################### FIN PARTE 3 ###########################
###################################################################
###################################################################

#Mi integral 1 se llama energia_total (ahora Integral_1) y la segunda
#se llama blackbody_sum (ahora se llama Integral_2)
#I_1 es de parametros wavelength, flujo, len(wavelength)-1
#I_2 es de parametros planx, blackbody_function, len(planx)-1

# esto sólo tiene relevancia para mí, es para no perderme
'''
______________
ejemplo quad
Calculate ∫x^2dx from 0 to 4 and compare with an analytic result
>>>
>>> from scipy import integrate
>>> x2 = lambda x: x**2
>>> integrate.quad(x2, 0, 4)
(21.333333333333332, 2.3684757858670003e-13)
>>> print(4**3 / 3.)  # analytical result
21.3333333333
'''
I1_trapz=integrate.trapz(flujo,wavelength)
print 'I1_trapz= '+str(I1_trapz)
#print 'estabien4= '+str(Integral_1-I1_trapz) #error de 10^-9

#la Integral_1 no la hice con funcion... no puedo aplicar quad predeterminado

I2_trapz=integrate.trapz(blackbody_function,planx)*unit.erg/(unit.cm**2*unit.s)
print 'I2_trapz= '+str(I2_trapz)
#print 'estabien5= '+str(Integral_2-I2_trapz)

#quad calculation, crear función para el quad
funcion=lambda x: (x**3)/(np.exp(x) - 1)

I2_quad= integrate.quad(funcion,0.0,np.inf)[0]*unit.erg/(unit.cm**2*unit.s)
print 'I2_quad= '+str(I2_quad) # Me da un vector, lo saqué a la mala, ya que el primer número era el correcto
#print 'estabien41= '+str(Integral_1-I1_quad) #error de 10^-9
print 'Comparar 1'
print 'algoritmo scipy tarda en I1_trapz= ', kronos('%timeit I1_trapz')
print 'algoritmo mio tarda en Integral_1= ', kronos('%timeit Integral_1')
print 'Comparar 2'
print 'algoritmo scipy tarda en I2_trapz= ', kronos('%timeit I2_trapz')
print 'algoritmo scipy tarda en I2_quad= ', kronos('%timeit I2_quad')
print 'algoritmo mio tarda en bbrcgs= ', kronos('%timeit bbrcgs')

'''
  ___     ___________.__             ___________           .___  ___
 / _ \_/\ \__    ___/|  |__   ____   \_   _____/ ____    __| _/ / _ \_/\
 \/ \___/   |    |   |  |  \_/ __ \   |    __)_ /    \  / __ |  \/ \___/
            |    |   |   Y  \  ___/   |        \   |  \/ /_/ |
            |____|   |___|  /\___  > /_______  /___|  /\____ |
                          \/     \/          \/     \/      \/


                        o
                         _'
                        {_}
                        |=|
              .  '      | |
         o  .   o   o   |@|
         . o  _o_._'_  /___\
        o_.__'|~~~~~/  |=2 |
       |~~~~~/ '-.-'   |=0 |
        '-.-'    |     |=1 |
          |     _|_    |=5_|
         _|_   `"""`   |_._|
        `"""`          `"""`
'''
