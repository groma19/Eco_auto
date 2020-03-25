import numpy as np
import pandas as pd
import os
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.title("Econometrics Auto")
root.geometry("105x50")


t_table_ = """12.706
4.303
3.182
2.776
2.571
2.447
2.365
2.306
2.262
2.228
2.201
2.179
2.16
2.145
2.131
2.12
2.11
2.101
2.093
2.086
2.08
2.074
2.069
2.064
2.06
2.056
2.052
2.048
2.045
2.042""".split("\n")
t_table = list(map(float, t_table_))

f_table_ = """199.5
19
9.55
6.94
5.79
5.14
4.74
4.46
4.26
4.1
3.98
3.89
3.81
3.74
3.68
3.63
3.59
3.55
3.52
3.49
3.47
3.44
3.42
3.4
3.39
3.37
3.35
3.34
3.33
3.32""".split("\n")
f_table = list(map(float, f_table_))


path = os.getcwd()
filename = ""
def readFile():
    global filename
    filename = filedialog.askopenfilename(initialdir="path", title="Select File",
                                          filetypes=(("LibreOfficeCalc","*.ods"), ("all files", ".")))


def mainApp():
    if filename[-3:] == "ods":
        data = pd.read_excel(filename, engine="odf")
    else:
        data = pd.read_excel(filename)

    n = len(data)
    m = 2

    data["X1^2"] = data["X1"]**2
    data["X2^2"] = data["X2"]**2
    data["Y^2"] = data["Y"]**2
    data["X1*X2"] = data["X1"] * data["X2"]
    data["X1*Y"] = data["X1"] * data["Y"]
    data["X2*Y"] = data["X2"] * data["Y"]

    X = np.array([np.ones(len(data["X1"])), data["X1"], data["X2"]]).T
    Y = np.array(data["Y"]).T

    B_hat = (np.linalg.pinv(X.T @ X) @ X.T @ Y).round(4)

    b0 = B_hat[0]
    b1 = B_hat[1]
    b2 = B_hat[2]

    data["Y_hat"] = (b0 + b1*data["X1"] + b2*data["X2"]).round(4)
    data["e"] = (data["Y"] - data["Y_hat"]).round(4)
    data["e^2"] = (data["e"]**2).round(4)

    sum_of_data = pd.DataFrame(data.sum(axis = 0)).round(4)
    average_of_data = pd.DataFrame(sum_of_data/len(data)).round(4)

    sum_x1 = round(sum_of_data.loc["X1"][0], 4)
    sum_x2 = round(sum_of_data.loc["X2"][0], 4)
    sum_y = round(sum_of_data.loc["Y"][0], 4)
    sum_x1_kv = round(sum_of_data.loc["X1^2"][0], 4)
    sum_x2_kv = round(sum_of_data.loc["X2^2"][0], 4)
    sum_x1_x2 = round(sum_of_data.loc["X1*X2"][0], 4)
    sum_x1_y = round(sum_of_data.loc["X1*Y"][0], 4)
    sum_x2_y = round(sum_of_data.loc["X2*Y"][0], 4)

    var_x1 = (average_of_data.loc["X1^2"] - average_of_data.loc["X1"]**2)[0]
    var_x2 = (average_of_data.loc["X2^2"] - average_of_data.loc["X2"]**2)[0]
    cov_x1_x2 = (average_of_data.loc["X1*X2"] - average_of_data.loc["X1"]*average_of_data.loc["X2"])[0]
    cov_x1_y = (average_of_data.loc["X1*Y"] - average_of_data.loc["X1"]*average_of_data.loc["Y"])[0]
    cov_x2_y = (average_of_data.loc["X2*Y"] - average_of_data.loc["X2"]*average_of_data.loc["Y"])[0]
    var_y = (average_of_data.loc["Y^2"] - average_of_data.loc["Y"]**2)[0]

    var_x1 = round(var_x1, 4)
    var_x2 = round(var_x2, 4)
    cov_x1_x2 = round(cov_x1_x2, 4)
    cov_x1_y = round(cov_x1_y, 4)
    cov_x2_y = round(cov_x2_y, 4)
    var_y = round(var_y, 4)

    avg_x1 = average_of_data.loc["X1"][0]
    avg_x2 = average_of_data.loc["X2"][0]
    avg_y = average_of_data.loc["Y"][0]

    S_x1 = round(var_x1**(1/2), 4)
    S_x2 = round(var_x2**(1/2), 4)
    S_y = round(var_y**(1/2), 4)

    b1_talg = round(b1*S_x1/S_y, 4)
    b2_talg = round(b2*S_x2/S_y, 4)

    qsi_1 = b1*average_of_data.loc["X1"][0]/average_of_data.loc["Y"][0]
    qsi_2 = b2*average_of_data.loc["X2"][0]/average_of_data.loc["Y"][0]

    qsi_1 = round(qsi_1, 4)
    qsi_2 = round(qsi_2, 4)

    Se_kv = (1/(n-m-1) * sum_of_data.loc["e^2"])[0].round(4)
    Se = round(Se_kv**(1/2), 4)

    variaciebi = np.linalg.pinv(X.T @ X).round(4)

    var_b0 = variaciebi[0][0]
    var_b1 = variaciebi[1][1]
    var_b2 = variaciebi[2][2]

    S_b0 = (Se_kv*var_b0)**(1/2)
    S_b1 = (Se_kv*var_b1)**(1/2)
    S_b2 = (Se_kv*var_b2)**(1/2)

    S_b0 = round(S_b0, 4)
    S_b1 = round(S_b1, 4)
    S_b2 = round(S_b2, 4)

    t_krit = t_table[n-m-1]

    t_b0 = round(b0/S_b0, 4)
    t_b1 = round(b1/S_b1, 4)
    t_b2 = round(b2/S_b2, 4)

    std_shecdomebi = [S_b0, S_b1, S_b2]
    t_stat = [t_b0, t_b1, t_b2]

    c = 1

    PI = np.array([[0],[1],[1]])
    PI_T = PI.T

    zeda = PI_T@B_hat - c
    qveda = Se*(PI_T@np.linalg.pinv(X.T@X)@PI)**(1/2)

    zeda = zeda.round(4)[0]
    qveda = qveda.round(4)[0][0]

    t_abs = abs(zeda/qveda)

    tb_talg = (B_hat[2]+B_hat[1] - 1)/B_hat[2]
    tb_talg = round(tb_talg, 4)

    r_yx1 = cov_x1_y/(var_x1*var_y)**(1/2)
    r_yx2 = cov_x2_y/(var_x2*var_y)**(1/2)
    r_x1x2 = cov_x1_x2/(var_x1*var_x2)**(1/2)

    r_yx1 = round(r_yx1, 4)
    r_yx2 = round(r_yx2, 4)
    r_x1x2 = round(r_x1x2, 4)

    Ryx = ((r_yx1**2 + r_yx2**2 - 2*r_yx1*r_yx2*r_x1x2)/(1-r_x1x2**2))**(1/2)
    Ryx = round(Ryx, 4)

    Q = round(var_y, 4)
    Q_E = sum_of_data.loc["e^2"][0]
    Q_R = Q - Q_E

    R_kv = round(Q_R/Q, 4)

    R_kv_koreq = R_kv - (m/(n-m-1))*(1 - R_kv)
    R_kv_koreq = round(R_kv_koreq, 4)

    F_krit = f_table[n-m-2]

    F = (n-m-1)*R_kv/(m*(1-R_kv))
    F = round(F, 4)

    print("1) ნორმალურ განტოლებათა სისტემა ზოგადად:")
    print("""
     / n*b0 + b1*sum(x_i1) + b2*sum(x_i2) = sum(y_i)
     |
    <  b0*sum(x_i1) + b1*sum(x_i1^2) + b2*sum(x_i1*x_i2) = sum(y_i*x_i1)
     |
     \ b0*sum(x_i2) + b1*sum(x_i1*x_i2) + b2*sum(x_i2^2) = sum(y_i*x_i2) 
    """)

    print("ნორმალურ განტოლებათა სისტემა ჩვენს შემთხვევაში:")
    print(f"{n}*b0 + {sum_x1}*b1 + {sum_x2}*b2 = {sum_y}")
    print(f"{sum_x1}*b0 + {sum_x1_kv}*b1 + {sum_x1_x2}*b2 = {sum_x1_y}")
    print(f"{sum_x2}*b0 + {sum_x1_x2}*b1 + {sum_x2_kv}*b2 = {sum_x2_y}")
    print()

    print("კოეფიციენტების გამოთვლა ჩვეულებრივი ხერხით: ")
    print(f"var_x1 = (1/n)*sum(x1^2) - avg_x1^2 = {var_x1}")
    print(f"var_x2 = (1/n)*sum(x2^2) - avg_x2^2 = {var_x2}")
    print(f"cov_x1_x2 = (1/n)*x1*x2 - avg_x1*avg_x2 = {cov_x1_x2}")
    print(f"cov_x1_y = (1/n)*x1*y - avg_x1*avg_y = {cov_x1_y}")
    print(f"cov_x2_y = (1/n)*x2*y - avg_x2*avg_y = {cov_x2_y}")
    print(f"var_y = (1/n)*sum(y^2) - avg_y^2 = {var_y}")
    print()
    print(f"""b1 = (cov_x1_y*var_x2 - cov_x2_y*cov_x1_x2)/
                    (var_x1*var_x2 - (cov_x1_x2)^2) = {b1}""")

    print(f"""b2 = (cov_x2_y*var_x1 - cov_x1_y*cov_x1_x2)/
                    (var_x1*var_x2 - (cov_x1_x2)^2) = {b2}""")

    print(f"b0 = avg_y - b1*avg_x1 - b2*avg_x2 = {b0}")
    print()
    print("კოეფიციენტების გამოთვლა მატრიცული ხერხით:")
    print(f"B_hat = (X_t*X)^-1 * (X_t*Y) = {B_hat}.T")
    print()
    print("უმცირეს კვადრატთა მეთოდით განტოლების შეფასება:")
    print("y_hat = {0:.4f} + {1:.4f}*x1 + {2:.4f}*x2".format(b0, b1, b2))
    print()
    print("2) მიღებული მოდელისა დაკოფიციენტების ეკონომიკური ინტერპრეტაცია:")
    print("---------------------------------------------------------------")
    print()
    print("3) სტანდარტიზებული კოეფიციენტების გამოთვლა: ")
    print(f"""S_x1 = sqrt(var_x1) = {S_x1}
    S_x2 = sqrt(var_x2) = {S_x2}
    S_y = sqrt(var_y) = {S_y}""")
    print()
    print(f"""b1_talg = b1*S_x1/S_y = {b1_talg}
    b2_talg = b2*S_x2/S_y = {b2_talg}""")
    print()
    print(f"""qsi_1 = b1*avg_x1/avg_y = {qsi_1}
    qsi_2 = b2*avg_x2/avg_y = {qsi_2}""")
    print()
    print("სტანდარტიზებული განტოლება:")
    print(f"w = {b1_talg}*z1 + {b2_talg}*z2")
    print()
    print("4) კოეფიციენტების სტანდარტული შეცდომების გამოთვლა:")
    print(f"Se_kv = (1/(n-m-1))*sum(e^2) = {Se_kv}")
    print(f"""var_b0 = {var_b0}
    var_b1 = {var_b1}
    var_b2 = {var_b2}""")

    print(f"""S_b0 = sqrt(Se_kv*var_b0) = {S_b0}
    S_b1 = sqrt(Se_kv*var_b1) = {S_b1}
    S_b2 = sqrt(Se_kv*var_b2) = {S_b2}""")

    print()
    print("მნიშვნელოვნების შეფასება სტიუდენტის 5% მნიშვნელოვნების დონისთვის:")
    print(f"tკრ = {t_krit}")
    print("H0: beta == 0")
    print("H1: beta != 0")
    print(f"t_b0 = b0/S_b0 = {t_b0}")
    print(f"t_b1 = b1/S_b1 = {t_b1}")
    print(f"t_b2 = b2/S_b2 = {t_b2}")
    print()

    punct = ","
    for i in range(len(t_stat)):
        if i == 2:
            punct = "."
        if t_stat[i] > t_krit:
            print(f"b{i} არ არის მნიშვნელოვანი კოეფიციენტი{punct}")
        else:
            print(f"b{i} მნიშვნელოვანი კოეფიციენტია{punct}")
            
    print()
    print("95% ნდობის ინტერვალი:")

    for i in range(len(B_hat)): 
        print(f"b{i} - t_krit * S_b{i} < beta{i} < b{i} + t_krit * S_b{i}")
        
    print()
    for i in range(len(B_hat)):
        print("{0:.4f} < beta{2} < {1:.4f}".format(B_hat[i] - t_krit*std_shecdomebi[i],
                                                    B_hat[i] + t_krit*std_shecdomebi[i],
                                                    i))
                                                   
    print()
    print("5) ჰიპოთეზის შემოწმება წრფივი შეზღუდვის შესახებ:")
    print(f"5%-იანი მნიშვნელოვნების დონისთვის tკრ={t_krit}")
    print(f"H0: beta1 + beta2 = {c}")
    print()
    print("t = (PI*B_hat - c)/(Se*sqrt(PI_t * (X_t*X)^-1 * PI))")
    print(f"c={c}; PI_T={PI_T[0]}; B_hat = {B_hat}.T")
    print()
    print(f"|t|={t_abs}")

    if t_abs>=t_krit:
        print("|t|>tკრ")
    else:
        print("|t|<tკრ")

    if t_abs>=t_krit:
        print("H0 არ არის სამართლიანი.")
    else:
        print("H0 სამართლიანია.")
        
    print()
    print("მეორე გზა:")
    print("""
    y = beta0 + beta1*x1 + b2*x2 + u
    y = beta0 + beta1*(x1-x2) + (beta2+beta1)*x2 + u
    """)
    print("y_hat = {0} + {1}*(x1-x2) + {2}*x2 + u".format(B_hat[0], B_hat[1], B_hat[2]+B_hat[1]))
    print()
    print(f"tb_talg = ({B_hat[2]+B_hat[1]} - {c})/{qveda} = {tb_talg}")

    if t_abs>=t_krit:
        print("|t_talg|>tკრ")
    else:
        print("|t_talg|<tკრ")

    if t_abs>=t_krit:
        print("H0 ჰიპოთეზა არ არის სამართლიანი.")
    else:
        print("H0 ჰიპოთეზას ვერ უარვყოფთ.")
        
    print()
    print("6) მრავლობითი კორელაციის კოეფიციენტის გამოთვლა:")
    print(f"""
          /  1    r_yx1   r_yx2\\
    Ryx = |r_x1y    1    r_x1x2|
          \\r_x2y  r_x2x1    1  /
    """)
    print(f"""
    r_yx1 = sqrt(cov_x1_y/(var_x1*var_y)) = {r_yx1}
    r_yx2 = sqrt(cov_x2_y/(var_x2*var_y)) = {r_yx2}
    r_x1x2 = sqrt(cov_x1_x2/(var_x1*var_x2)) = {r_x1x2}
    """)

    print(f"Ryx = sqrt((r_yx1**2 + r_yx2**2 - 2*r_yx1*r_yx2*r_x1x2)/(1-r_x1x2**2)) = {Ryx}")
    print()
    print("მრავლობითი დეტერმინაციის კოეფიციენტის გამოთვლა:")
    print("Q = Q_R + Q_E")
    print(f"Q = var_y = {Q}")
    print(f"Q_E = sum(e^2) = {Q_E}")
    print(f"Q_R = Q - Q_E = {Q_R}")
    print(f"R_kv = Q_R/Q = 1 - Q_E/Q = {R_kv}")
    print()
    print("დეტერმინაციის გადაუადგილებელი (კორექტირებული) მნიშვნელობის გამოთვლა:")
    print(f"R_kv_koreq = R_kv_koreq = R_kv - (m/(n-m-1))*(1 - R_kv) = {R_kv_koreq}")
    print()
    print("განტოლების მნიშვნელოვნების შეფასება:")
    print(f"5% მნიშვნელოვნების დონისთვის Fკრ={F_krit}")

    print("H0: beta1 = beta2 = 0")

    print(f"F = (n-m-1)*R_kv/(m*(1-R_kv)) = {F}")
    print()
    if F < F_krit:
        print("F < Fკრ, ამიტომ ნულოვან ჰიპოთეზას ვერ უარვყოფთ -->")
        print("რეგრესიის განტოლების მოდელი არ არის მნიშვნელოვანი.")
        
    elif F > F_krit:
        print("F > Fკრ, შეგვიძლია უარვყოთ ნულოვანი ჰიპოთეზა.")
    else:
        print("F=Fკრ, საჭიროა დამატებითი კვლევის ჩატარება.")
        
    print()
    print("7) რეგრესიის მოდელით საპროგნოზო მნიშვნელობის განსაზღვრა:")



    x1_fifqi = int(input("შეიყვანეთ მნიშვნელობა x1 ცვლადისთვის: "))
    x2_fifqi = int(input("შეიყვანეთ მნიშვნელობა x2 ცვლადისთვის: "))

    y_qudiani_fifqi = round(b0 + b1*x1_fifqi + b2*x2_fifqi, 4)

    print(f"პროგნოზირებული მნიშვნელობაა {y_qudiani_fifqi}")


openFile = tk.Button(root, text="Open File", padx=10, pady=5,fg="white", bg="black", command=readFile, anchor="center")
openFile.pack()

runApp = tk.Button(root, text="Run App", padx=10, pady=5,fg="white", bg="black", command=mainApp, anchor="center")
runApp.pack()



root.mainloop()






