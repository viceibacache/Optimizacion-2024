from gurobipy import GRB, Model
import gurobipy as gp
import pandas as pd
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

ruta_completa_a = os.path.join(script_dir, 'A.csv')
ruta_completa_b = os.path.join(script_dir, 'B.csv')
ruta_completa_c = os.path.join(script_dir, 'C.csv')
ruta_completa_d = os.path.join(script_dir, 'D.csv')
ruta_completa_e = os.path.join(script_dir, 'E.csv')
ruta_completa_f = os.path.join(script_dir, 'F.csv')
ruta_completa_h = os.path.join(script_dir, 'H.csv')
ruta_completa_i = os.path.join(script_dir, 'I.csv')
ruta_completa_j = os.path.join(script_dir, 'J.csv')
ruta_completa_k = os.path.join(script_dir, 'K.csv')

#####################
# Carga de Archivos #
#####################

concentración_notas_estudante = pd.read_csv(ruta_completa_a, header = None).squeeze("columns").to_numpy() # A_i
ranking_establecimiento = pd.read_csv(ruta_completa_b, header = None).squeeze("columns").to_numpy() # B_j
capacidad_total_establecimiento = pd.read_csv(ruta_completa_c, header = None).squeeze("columns").to_numpy() # C_j
estudiante_con_hermanos_en_establecimiento = pd.read_csv(ruta_completa_d, header = None).to_numpy() # D_i_j
estudiante_hijo_de_funcionario_en_establecimiento = pd.read_csv(ruta_completa_e, header = None).to_numpy() # E_i_j
estudiante_es_PESD = pd.read_csv(ruta_completa_f, header = None).squeeze("columns").to_numpy() # F_i
estudante_postula_a_establecimiento = pd.read_csv(ruta_completa_h, header = None).to_numpy() # H_i_j
estudiante_reside_en_la_misma_comuna_que_el_establecimiento = pd.read_csv(ruta_completa_i, header = None).to_numpy() # I_i_j
estudiante_situacion_economica_inferior_al_promedio_en_comuna_de_establecimiento = pd.read_csv(ruta_completa_j, header = None).to_numpy() # J_i_j
estudiante_vuelve_a_establecimiento = pd.read_csv(ruta_completa_k, header = None).to_numpy() # K_i_j

# Ponderadores
ponderadores = [1, 0.143, 0.143, 0.143, 0.143, 0.143, 0.143]

#############
# Conjuntos #
#############
I = 100 #Reemplazar por número de estudiantes
J = 50 #Reemplazar por número de colegios

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
jp = estudiante_situacion_economica_inferior_al_promedio_en_comuna_de_establecimiento
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

Punt = [[0 for j in range(J)] for i in range(I)]

#Definición de parametro auxiliar
for i in range(I):
    for j in range(J):
        Punt[i][j] = (alpha*a[i] + (b[j] - a[i]) * (beta*e[i][j] + gamma*d[i][j] + delta*jp[i][j] + epsilon*ip[i][j] + zeta*k[i][j] + eta*f[i])) * h[i][j]

###################################
# Creación del Modelo Y Variables #
###################################

model = Model('Admisión Sistema Escolar')

X = model.addVars(I, J, name = 'X', vtype = GRB.BINARY) #Se entrega vacante a I en J
Y = model.addVars(J, name = 'Y', vtype = GRB.INTEGER) #Cantidad vacantes ocupadas en J
Z = model.addVars(J, name = 'Z', vtype = GRB.INTEGER) #Cantidad estudiantes prioritarios en J
Omega = model.addVars(I, J, name = 'Omega', vtype = GRB.BINARY) #P es mayor o igual a b
Phi = model.addVars(I, J, name = 'Phi', vtype = GRB.CONTINUOUS) #Variable auxiliar de diferencia Absoluta
Psi = model.addVars(I, J, name = 'Psi', vtype = GRB.CONTINUOUS) #Variable auxiliar linealidad

#################
# Restricciones #
#################

# Prueba

for i in range(I):
    for j in range(J):
        if Punt[i][j] - b[j] >= 0:
            model.addConstr(Omega[i, j] == 1, "Omega_constraint")
        else:
            model.addConstr(Omega[i, j] == 0, "Omega_constraint")

#Cantidad de alumnos admitidos
for j in range(J):
    model.addConstr(Y[j] == gp.quicksum(X[i,j] for i in range(I)), name=f"R cantidad alumnos admitidos {j}")

#Limitar cantidad de almumnos por vacantes disponibles por establecimiento
for j in range(J):
    model.addConstr(Y[j] <= c[j], name=f"R limitar vacantes {j}")

#Cada estudiante que postula se queda en exactamente un colegio
for i in range(I):
    model.addConstr(gp.quicksum(X[i,j] for j in range(J)) == 1, name=f"R 1 colegio por alumno admitidos {i}")


# Esto lo saque porque puede que nunca se cumpla si hay pocos pesd hay que revisarlo

# #Cantidad de alumnos con prioridad por ser PESD
# for j in range(J):
#     model.addConstr(Z[j] == gp.quicksum(X[i,j]*k[i,j] for i in range(I)))

# #Limitar cantidad minima de PESD
# for j in range(J):
#     model.addConstr(Z[j] >= c[j] * 0.15)

#Asignar valor de variable auxiliar Phi
for i in range(I):
    for j in range(J):
        model.addConstr(Psi[i,j] == 
                        Omega[i,j]*(Punt[i][j] - b[j]) + (1 - Omega[i,j])*(b[j] - Punt[i][j]), name=f"Phi {i} {j}")

#Psi no puede exceder Phi
for i in range(I):
    for j in range(J):
        model.addConstr(Psi[i,j] <= Omega[i,j], name=f"Psi < Phi {i} {j}")

#Psi no puede exceder el valor M = 5000 * Xij
for i in range(I):
    for j in range(J):
        model.addConstr(Psi[i,j] <= 5000 * X[i,j], name=f"Psi < M*X {i} {j}")

#Psi al menos igual a la expresion Phi - (1 - M) * Xij
for i in range(I):
    for j in range(J):
        model.addConstr(Psi[i,j] >= Phi[i,j] - (1 - 5000) * X[i,j], name=f"Psi > Phi - M*X {i} {j}")

####################
# Función Objetivo #
####################
objetivo = gp.quicksum(gp.quicksum(Psi[i,j] for i in range(I)) for j in range(J))
model.setObjective(objetivo, GRB.MINIMIZE)

model.update()

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

elif model.status == GRB.INFEASIBLE:
    # Compute the IIS
    model.computeIIS()
    model.write("model.ilp")

    # Print the IIS
    for c in model.getConstrs():
        if c.IISConstr:
            print(f"Infeasible constraint: {c.constrName}")

    for v in model.getVars():
        if v.IISLB:
            print(f"Infeasible lower bound: {v.varName}")
        if v.IISUB:
            print(f"Infeasible upper bound: {v.varName}")


else:
    print("El modelo no es optimo")
    print("Guz, hice algo mal D:")

