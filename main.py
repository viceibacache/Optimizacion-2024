from gurobipy import GRB, Model
import gurobipy as gp
import pandas as pd

#####################
# Carga de Archivos #
#####################
concentración_notas_estudante =
ranking_establecimiento =
capacidad_total_establecimiento =
estudiante_con_hermanos_en_establecimiento =
estudiante_hijo_de_funcionario_en_establecimiento =
estudiante_es_PESD =
estudante_postula_a_establecimiento =
estudiante_reside_en_la_misma_comuna_que_el_establecimiento =
estudiante_situacion_economica_inferior_al_promedio_en_comuna_de_establecimiento =
estudiante_vuelve_a_establecimiento =
ponderadores =

#############
# Conjuntos #
#############
I = 100 #Reemplazar por número de estudiantes
J = 20 #Reemplazar por número de colegios

##############
# Parametros #
##############
a = concentración_notas_estudante
b = ranking_establecimiento
c = capacidad_total_establecimiento
d = estudiante_con_hermanos_en_establecimiento
e = estudiante_hijo_de_funcionario_en_establecimiento
f = estudiante_es_PESD
h = estudante_postula_a_establecimiento
ip = estudiante_reside_en_la_misma_comuna_que_el_establecimiento
j = estudiante_situacion_economica_inferior_al_promedio_en_comuna_de_establecimiento
k = estudiante_vuelve_a_establecimiento

################
# Ponderadores #
################
alpha = ponderadores[0] #ponderador_notas
beta = ponderadores[1] #ponderador_estudiante_hijo_funcionario
gamma = ponderadores[2] #poderador_hermanos_en_establecimiento
delta = ponderadores[3] #poderador_situacion_economica_menor
epsilon = ponderadores[4] #ponderador_misma_comuna
zeta = ponderadores[5] #ponderador_estudiante_vuelve
eta = ponderadores[6] #ponderador_estudiante_es_PESD

###################################
# Creación del Modelo Y Variables #
###################################
model = Model('Admisión Sistema Escolar')

P = model.addVars(I, J, name = 'P', vtype = GRB.CONTINUOUS) #Parametro auxiliar
X = model.addVars(I, J, name = 'X', vtype = GRB.BINARY) #Se entrega vacante a I en J
Y = model.addVars(J, name = 'Y', vtype = GRB.INTEGER) #Cantidad vacantes ocupadas en J
Z = model.addVars(J, name = 'Z', vtype = GRB.INTEGER) #Cantidad estudiantes prioritarios en J
Omega = model.addVars(I, J, name = 'Omega', vtype = GRB.BINARY) #P es mayor o igual a b
Phi = model.addVars(I, J, name = 'Phi', vtype = GRB.CONTINUOUS) #Variable auxiliar de diferencia Absoluta
Psi = model.addVars(I, J, name = 'Psi', vtype = GRB.CONTINUOUS) #Variable auxiliar linealidad

#################
# Restricciones #
#################
#Cantidad de alumnos admitidos
for j in range(J):
    model.addConstr(Y[j] == gp.quicksum(X[i,j] for i in range(I)))

#Limitar cantidad de almumnos por vacantes disponibles por establecimiento
for j in range(J):
    model.addConstr(Y[j] <= c[j])

#Cada estudiante que postula se queda en exactamente un colegio
for i in range(I):
    model.addConstr(gp.quicksum(X[i,j] for j in J) == 1)

#Cantidad de alumnos con prioridad por ser PESD
for j in range(J):
    model.addConstr(Z[j] == gp.quicksum(X[i,j]*k[i,j] for i in I))

#Limitar cantidad minima de PESD
for j in range(J):
    model.addConstr(Z[j] >= c[j] * 0.15)

#Definición de parametro auxiliar
for i in range(I):
    for j in range(J):
        model.addConstr(P[i,j] == (alpha*a[i] + (b[j] - a[i]) *
                                    (beta*e[i][j] + gamma*d[i][j] + delta*[i][j] + 
                                     epsilon*ip[i][j] + zeta*k[i][j] + eta*f[i][j])) * h[i][j])

#Asignar valor de variable auxiliar Phi
for i in range(I):
    for j in range(J):
        model.addConstr(Psi[i,j] == 
                        Omega[i,j]*(P[i,j] - b[j]) + (1 - Omega[i,j])*(b[j] - P[i,j]))

#Psi no puede exceder Phi
for i in range(I):
    for j in range(J):
        model.addConstr(Psi[i,j] <= Omega[i,j])

#Psi no puede exceder el valor M = 5000 * Xij
for i in range(I):
    for j in range(J):
        model.addConstr(Psi[i,j] <= 5000 * X[i,j])

#Psi al menos igual a la expresion Phi - (1 - M) * Xij
for i in range(I):
    for j in range(J):
        model.addConstr(Psi[i,j] <= Phi[i,j] - (1 - 5000) * X[i,j])

####################
# Función Objetivo #
####################
objetivo = gp.quicksum(gp.quicksum(Psi[i,j] for i in I) for j in J)
model.setObjective(objetivo, GRB.MINIMIZE)

#########################
# Recopilación de Datos #
#########################
model.optimize()
if model.Status == GRB.OPTIMAL:
    print("El modelo es optimo")
    tabla = pd.DataFrame(0, index = [f"Alumno {i}" for i in range(I)], columns = ["Establecimiento"])
    for i in range(I):
        for j in range(J):
            if X[i,j].X == 1:
                tabla.at[f"Alumno {i}", "Establecimiento"] = j
    tabla.to_excel('tabla_colegios_de_estudiantes.xlsx')

else:
    print("El modelo no es optimo")
    print("Guz, hice algo mal D:")