import numpy as np
import pandas as pd

data = pd.read_excel("monacemebi.ods", engine="odf")

data["X1^2"] = data["X1"]**2
data["X2^2"] = data["X2"]**2
data["Y^2"] = data["Y"]**2
data["X1*X2"] = data["X1"] * data["X2"]
data["X1*Y"] = data["X1"] * data["Y"]
data["X2*Y"] = data["X2"] * data["Y"]

X = np.array([np.ones(len(data["X1"])), data["X1"], data["X2"]]).T
Y = np.array(data["Y"]).T

B_hat = (np.linalg.pinv(X.T @ X) @ X.T @ Y).round(4)

sum_of_data = pd.DataFrame(data.sum(axis = 0))
average_of_data = pd.DataFrame(sum_of_data/len(data))

var_x1 = average_of_data.loc["X1^2"] - average_of_data.loc["X1"]**2
var_x2 = average_of_data.loc["X2^2"] - average_of_data.loc["X2"]**2
cov_x1_x2 = average_of_data.loc["X1*X2"] - average_of_data.loc["X1"]*average_of_data.loc["X2"]
cov_x1_y = average_of_data.loc["X1*Y"] - average_of_data.loc["X1"]*average_of_data.loc["Y"]
cov_x2_y = average_of_data.loc["X2*Y"] - average_of_data.loc["X2"]*average_of_data.loc["Y"]

var_x1 = round(var_x1, 4)
var_x2 = round(var_x2, 4)
cov_x1_x2 = round(cov_x1_x2, 4)
cov_x1_y = round(cov_x1_y, 4)
cov_x2_y = round(cov_x2_y, 4)

b1 = (cov_x1_y * var_x2 - cov_x2_y * cov_x1_x2)/(var_x1 * var_x2 - cov_x1_x2**2)
b2 = (cov_x2_y * var_x1 - cov_x1_y * cov_x1_x2)/(var_x1 * var_x2 - cov_x1_x2**2)
b0 = (average_of_data.loc["Y"] - b1*average_of_data.loc["X1"] - b2*average_of_data.loc["X2"])[0]

b1 = round(b1, 4)
b2 = round(b2, 4)
b0 = round(b0, 4)

var_y = (average_of_data.loc["Y^2"] - average_of_data.loc["Y"]**2)[0]
var_y = round(var_y, 4)

S_x1 = round(var_x1[0]**(1/2), 4)
S_x2 = round(var_x2[0]**(1/2), 4)
S_y = round(var_y**(1/2), 4)

b1_talgovani = round(b1*S_x1/S_y, 4)
b2_talgovani = round(b2*S_x2/S_y, 4)

qsi_1 = b1*average_of_data.loc["X1"][0]/average_of_data.loc["Y"][0]
qsi_2 = b2*average_of_data.loc["X2"][0]/average_of_data.loc["Y"][0]

qsi_1 = round(qsi_1, 4)
qsi_2 = round(qsi_2, 4)

data["Y_hat"] = (b0 + b1*data["X1"] + b2*data["X2"]).round(4)

data["e"] = (data["Y"] - data["Y_hat"]).round(4)
data["e^2"] = (data["e"]**2).round(4)

sum_of_data = (pd.DataFrame(data.sum(axis = 0))).round(4)

average_of_data = pd.DataFrame(sum_of_data/len(data))

n=len(data)
m=2

Se_kvadrati = (1/(n-m-1) * sum_of_data.loc["e^2"])[0].round(4)

variaciebi = np.linalg.pinv(X.T @ X).round(4)

var_b0 = variaciebi[0][0]
var_b1 = variaciebi[1][1]
var_b2 = variaciebi[2][2]

S_b0 = (Se_kvadrati*var_b0)**(1/2)
S_b1 = (Se_kvadrati*var_b1)**(1/2)
S_b2 = (Se_kvadrati*var_b2)**(1/2)

S_b0 = round(S_b0, 4)
S_b1 = round(S_b1, 4)
S_b2 = round(S_b2, 4)

std_shecdomebi = [S_b0, S_b1, S_b2]

t_b0 = round(b0/S_b0, 4)
t_b1 = round(b1/S_b1, 4)
t_b2 = round(b2/S_b2, 4)

t_stat = [t_b0, t_b1, t_b2]

t_krit= 2.36

c = 1

pi = np.array([[0],[1],[1]])
pi_T = pi.T

Se = round(Se_kvadrati**(1/2), 4)

zeda = pi_T@B_hat - c
qveda = Se*(pi_T@np.linalg.pinv(X.T@X)@pi)**(1/2)

zeda = zeda.round(4)[0]
qveda = qveda.round(4)[0][0]

t = abs(zeda/qveda)

tb_talg = (B_hat[2]+B_hat[1] - 1)/B_hat[2]
tb_talg = abs(round(tb_talg, 4))

r_yx1 = cov_x1_y/(var_x1*var_y)**(1/2)
r_yx2 = cov_x2_y/(var_x2*var_y)**(1/2)
r_x1x2 = cov_x1_x2/(var_x1*var_x2)**(1/2)

r_yx1 = round(r_yx1, 4)[0]
r_yx2 = round(r_yx2, 4)[0]
r_x1x2 = round(r_x1x2, 4)[0]

Ryx_matrix = np.array([[1, r_yx1, r_yx2],
                      [r_yx1, 1, r_x1x2],
                      [r_yx2, r_x1x2, 1]])
 
Ryx = ((r_yx1**2 + r_yx2**2 - 2*r_yx1*r_yx2*r_x1x2)/(1-r_x1x2**2))**(1/2)
Ryx = round(Ryx, 4)

Q = round(var_y, 4)
Q_E = sum_of_data.loc["e^2"][0]
Q_R = Q - Q_E

R_2 = round(Q_R/Q, 4)

R_2_koreq = R_2 - (m/(n-m-1))*(1 - R_2)
R_2_koreq = round(R_2_koreq, 4)

F_krit = 4.74

F = (n-m-1)*R_2/(m*(1-R_2))
F = round(F, 4)

# prognozi

x1_fifqi = 70
x2_fifqi = 4

y_qudiani_fifqi = round(b0 + b1*x1_fifqi + b2*x2_fifqi, 4)

print(y_qudiani_fifqi)

















