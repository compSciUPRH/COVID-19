#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Determinación de días en que se rezaga la gráfica de hospitalizaciones de la de casos..

Modelo:
    t     = día
    C(t) = casos día t
    H(t) = Hospitalizaciones t
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

@author: José O. Sotero Esteva
@email: jose.sotero@upr.edu

Advertencia: Este código ni los resultados que produce tienen
el propósito de ser usados para tomar decisiones médicas o 
epidemiológicas de ningún tipo.

Fuentes de datos: https://github.com/rafalab/pr-covid/tree/master/dashboard

"""

import pandas as pd
import matplotlib.pyplot as plt

FECHA_CORTE = "2021-12-07"

casos = pd.read_csv("casos-2021-12-30_07 45 59.csv")
hosp_mort = pd.read_csv("hosp-mort-2021-12-30_07 45 59.csv")

# cambia formato de fechas
casos["fecha"] = pd.to_datetime(casos["date"])
hosp_mort["fecha"] = pd.to_datetime(hosp_mort["dates"])


# seleccionar columnas y filas y unir tablas
casos = casos[casos["testType"] == "Molecular+Antigens"][["fecha","cases_week_avg"]].set_index("fecha")
casos.rename(columns={"cases_week_avg":"C"}, inplace=True)
hosp = hosp_mort[["fecha","hospitalizaciones_week_avg"]].dropna().set_index("fecha")
hosp.rename(columns={"hospitalizaciones_week_avg":"H"}, inplace=True)

casos_hosp = casos.join(hosp)
casos_hosp = casos_hosp[casos_hosp.index >= pd.to_datetime("2020-04-22")]
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
print("Rezago (d): {} días\nescala (h): {}\nerror: {}".format(dmin, hmin, emin))

# proyección
futuro = casos[casos.index > pd.to_datetime(FECHA_CORTE)].copy()
futuro['H'] = (hmin * futuro.C.shift(int(dmin)))
futuro.dropna(inplace=True)
futuro.index = futuro.index + pd.DateOffset(dmin)
p_dia, x, p_hosp = futuro.reset_index().iloc[-1].values
print("Proyección: {} hospitalizaciones el {}".format(int(p_hosp), str(p_dia)))

# gráficas
casos.plot()
plt.title("Casos COVID en PR (antígeno + molecular)")
plt.ylabel('cantidad')
plt.show()

hosp.plot()
plt.title("Hospitalizaciones COVID en PR")
plt.ylabel('cantidad')
plt.show()

escalas.E.plot()
plt.title("Error del modelo de rezago")
plt.xlabel('dias')
plt.ylabel('error S(d)')
plt.show()

futuro.H.plot()
plt.title("Proyección hospitalizaciones")
plt.xlabel('dias')
plt.ylabel('cantidad')

plt.show()
