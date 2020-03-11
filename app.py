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

print(B_hat[0])
