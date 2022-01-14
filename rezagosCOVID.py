#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Determinación de días en que se rezaga la gráfica de hospitalizaciones o muertes de la de casos..

Modelo:
    t     = día
    C(t) = casos día t
    H(t) = Hospitalizaciones o muertes t
    Modelo: H(t) ~ h C(t+d)
            donde h factor de escala
                  d días de retraso
    S(d,h) = \suma_{t} (h C(t+d) - H(t))^2
    hm(d) =  valor de h que minimiza S(d,h) para d fija
    E(d) = S(d,hm(d))

Algoritmo:
    por cada retraso d
        calcular hm(d)
        calcular E(d)
    reportar d que minimiza E(d)

Cómputo de hm(d):
    
    d S(d,h)                   \sum_{t} H(t) C(t+d)
    -------- = 0   =>  hm(d) = --------------------
       dh                       \sum_{t} C^2(t+d)



Created on Thu Dec 30 08:48:32 2021
Posibilidad de seleccionar varias opciones de hospitalizaciones o muertes el 14/ene/2022

@author: José O. Sotero Esteva
@email: jose.sotero@upr.edu

Advertencia: Este código ni los resultados que produce tienen
el propósito de ser usados para tomar decisiones médicas o 
epidemiológicas de ningún tipo.

Fuentes de datos: https://github.com/rafalab/pr-covid/tree/master/dashboard

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# seleccionar:
REPORTE      = "prom. sem. muertes"  # seleccione un índice del diccionario REPORTES
FECHA_CORTE  = "2021-12-9"
FECHA_INICIO = "2021-04-01"
ARCH_CASOS   = "casos-2022-01-14_08 00 56.csv"
ARCH_HOSP    = "hosp-mort-2022-01-14_08 00 56.csv"

REPORTES = {"hospitalizaciones adultos": ["fecha","hospitalizaciones"],
           "hospitalizaciones (adultos+ped)": ["fecha","hospitalizaciones", "hospitalizaciones_ped"],
           "muertes diarias": ["fecha","muertes"],
           "prom. sem. muertes": ["fecha","muertes_week_avg"]}

# encabezado
print("==================================================")
print("===== {:^38s} =====".format("Reporte "+REPORTE))
print("==================================================")

# lectura
casos = pd.read_csv(ARCH_CASOS)
hosp_mort = pd.read_csv(ARCH_HOSP)

# cambia formato de fechas
casos["fecha"] = pd.to_datetime(casos["date"])
hosp_mort["fecha"] = pd.to_datetime(hosp_mort["dates"])


# seleccionar columnas y filas y unir tablas
casos = casos[casos["testType"] == "Molecular+Antigens"][["fecha","cases_week_avg"]].set_index("fecha")
casos.rename(columns={"cases_week_avg":"C"}, inplace=True)

# Selecciona datos degún tipo de reporte
#hosp_mort['H'] = hosp_mort.hospitalizaciones + hosp_mort.hospitalizaciones_ped
#hosp = hosp_mort[["fecha","H"]].dropna().set_index("fecha")
hosp_mort = hosp_mort[REPORTES[REPORTE]].set_index("fecha").dropna()
hosp_mort['H'] = hosp_mort.sum(axis=1)
hosp = hosp_mort[["H"]]

#hosp = hosp_mort[["fecha","hospitalizaciones"]].dropna().set_index("fecha")
#hosp.rename(columns={"hospitalizaciones":"H"}, inplace=True)

casos_hosp = casos.join(hosp)
casos_hosp = casos_hosp[casos_hosp.index >= pd.to_datetime(FECHA_INICIO)]
casos_hosp = casos_hosp[casos_hosp.index <= pd.to_datetime(FECHA_CORTE)]
casos_hosp.fillna(0, inplace=True)

# acopio de escalas y errores
hs   = []
errs = []
for d in range(90): # retrasos hasta tres meses
    SHC = (casos_hosp.C.shift(d) * casos_hosp.H).sum()
    SCC = (casos_hosp.C.shift(d)**2).sum()
    h   = SHC / SCC
    err = ((h * casos_hosp.C.shift(d) - casos_hosp.H)**2).sum()
    hs.append(h)
    errs.append(err)

escalas = pd.DataFrame(data={'h':hs, 'E':errs})
escalas.index.name = 'd'

min_err = escalas.E.min()
dmin, hmin, emin = escalas[escalas.E == min_err].reset_index().values[0]
dmin = int(dmin)
print("Rezago (d): {} días\nescala (h): {}\nerror: {}".format(dmin, hmin, emin))

# comparación entre datos reales y modelo
comp = pd.DataFrame(casos_hosp.H)
comp['Modelo'] = hmin * casos.C.shift(int(dmin))

# proyección
futuro = casos[casos.index > pd.to_datetime(FECHA_CORTE)].copy()
new_dates = [casos.index[-1] + pd.DateOffset(d) for d in range(1,dmin+1)]
futuro = futuro.append(pd.DataFrame(data={'C':dmin*[np.nan]}, index=new_dates))
futuro['H'] = (hmin * futuro.C.shift(int(dmin)))
futuro.drop(columns=['C'],inplace=True)
futuro.dropna(inplace=True)
#futuro.index = futuro.index + pd.DateOffset(dmin)
max_hosp = futuro.max()
p_dia, p_hosp = futuro[futuro == max_hosp].dropna().reset_index().values[0]
#p_dia, p_hosp = futuro.reset_index().iloc[-1].values
print("Proyección: máximo de {} {} el {}".format(int(p_hosp), REPORTE, str(p_dia)))

futuro['H-real'] = hosp.H

# ajuste de proyección
#T_INFECCION_REL = futuro['H-real'].sum() / futuro.dropna().H.sum()
T_INFECCION_REL = 0.16
newcolname = '{:0.2f}H'.format(T_INFECCION_REL)
futuro[newcolname] = T_INFECCION_REL * futuro.H
max_hosp = futuro[newcolname].max()
p_dia, p_hosp = futuro[futuro[newcolname] == max_hosp][newcolname].reset_index().values[0]
#p_dia, p_hosp = futuro[newcolname].reset_index().iloc[-1].values
print("Proyección ajustada con factor {:0.2f}: máximo de {} {} el {}".format(T_INFECCION_REL, int(p_hosp), REPORTE, str(p_dia)))

# gráficas
casos.plot()
plt.title("Casos COVID en PR (antígeno + molecular)")
plt.ylabel('cantidad')
plt.show()

hosp.plot()
plt.title("{} COVID en PR (totales)".format(REPORTE.title()))
plt.ylabel('cantidad')
plt.show()

escalas.E.plot()
plt.title("Error del rezago de {}".format(REPORTE.title()))
plt.xlabel('dias')
plt.ylabel('error S(d)')
plt.show()

futuro.rename(columns={"H":"Est-delta","H-real":"real", newcolname:"ajustado ("+newcolname+")"}, inplace=True)
futuro.plot(color=['b','r','g'])
plt.title("Proyección {}".format(REPORTE.title()))
plt.xlabel('dias')
plt.ylabel('cantidad (prom. sem.)')
plt.show()


comp.plot()
plt.title("Comparación entre datos reales y modelo")
plt.ylabel('{} (prom. sem.)'.format(REPORTE))
plt.show()

