import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

np.random.seed(7)

T = 250

A = np.array([[0.92, 0.10],
              [0.00, 0.95]])
B = np.array([[0.25],
              [0.10]])
C = np.array([[1.0, 0.0]])

D_true = np.array([[0.18, -0.10],
                   [0.05,  0.12]])

F_true = np.array([[0.995, 0.0],
                   [0.0,   0.990]])

Qx = np.diag([0.03**2, 0.03**2])
Qz = np.diag([0.01**2, 0.01**2])
R = np.array([[0.06**2]])

x = np.zeros((2, T + 1))
z = np.zeros((2, T + 1))
y = np.zeros(T + 1)
u = np.zeros((1, T))

for t in range(T):
    u[0, t] = 0.6 * np.sin(2 * np.pi * t / 35) + 0.3 * np.cos(2 * np.pi * t / 57)

z[:, 0] = np.array([0.8, -0.5])

for t in range(T):
    wx = np.random.multivariate_normal(np.zeros(2), Qx)
    wz = np.random.multivariate_normal(np.zeros(2), Qz)
    v = np.random.normal(0, np.sqrt(R[0, 0]))
    x[:, t + 1] = A @ x[:, t] + B @ u[:, t] + D_true @ z[:, t] + wx
    z[:, t + 1] = F_true @ z[:, t] + wz
    y[t] = (C @ x[:, t])[0] + v

y[T] = (C @ x[:, T])[0] + np.random.normal(0, np.sqrt(R[0, 0]))

A1 = A.copy()
B1 = B.copy()
C1 = C.copy()
Q1 = Qx.copy()
R1 = R.copy()

xhat1 = np.zeros((2, T + 1))
P1 = np.eye(2) * 0.5

y_pred1 = np.zeros(T + 1)
innov1 = np.zeros(T + 1)

for t in range(T):
    xpred = A1 @ xhat1[:, t] + B1 @ u[:, t]
    Ppred = A1 @ P1 @ A1.T + Q1

    y_pred1[t + 1] = (C1 @ xpred)[0]
    S = C1 @ Ppred @ C1.T + R1
    K = Ppred @ C1.T @ np.linalg.inv(S)

    innov1[t + 1] = y[t + 1] - y_pred1[t + 1]
    xhat1[:, t + 1] = xpred + (K.flatten() * innov1[t + 1])
    P1 = (np.eye(2) - K @ C1) @ Ppred

A2 = np.block([
    [A, D_true],
    [np.zeros((2, 2)), F_true]
])
B2 = np.vstack([B, np.zeros((2, 1))])
C2 = np.hstack([C, np.zeros((1, 2))])
Q2 = np.block([
    [Qx, np.zeros((2, 2))],
    [np.zeros((2, 2)), Qz]
])
R2 = R.copy()

xhat2 = np.zeros((4, T + 1))
P2 = np.eye(4) * 0.5

y_pred2 = np.zeros(T + 1)
innov2 = np.zeros(T + 1)

for t in range(T):
    xpred = A2 @ xhat2[:, t] + B2 @ u[:, t]
    Ppred = A2 @ P2 @ A2.T + Q2

    y_pred2[t + 1] = (C2 @ xpred)[0]
    S = C2 @ Ppred @ C2.T + R2
    K = Ppred @ C2.T @ np.linalg.inv(S)

    innov2[t + 1] = y[t + 1] - y_pred2[t + 1]
    xhat2[:, t + 1] = xpred + (K.flatten() * innov2[t + 1])
    P2 = (np.eye(4) - K @ C2) @ Ppred

zhat2 = xhat2[2:, :]

def autocorr(x, max_lag=15):
    x = np.asarray(x)
    x = x - x.mean()
    denom = np.dot(x, x)
    vals = []
    for lag in range(1, max_lag + 1):
        vals.append(np.dot(x[:-lag], x[lag:]) / denom)
    return np.array(vals)

mse1 = np.mean((y[1:] - y_pred1[1:]) ** 2)
mse2 = np.mean((y[1:] - y_pred2[1:]) ** 2)

acf1 = autocorr(innov1[1:], 15)
acf2 = autocorr(innov2[1:], 15)

summary = pd.DataFrame({
    "Metric": [
        "One-step prediction MSE",
        "Innovation variance",
        "Mean absolute innovation autocorr (lags 1-15)",
        "Lag-1 innovation autocorr",
        "Lag-5 innovation autocorr"
    ],
    "Model 1: no hidden context": [
        mse1,
        np.var(innov1[1:]),
        np.mean(np.abs(acf1)),
        acf1[0],
        acf1[4]
    ],
    "Model 2: hidden context KF": [
        mse2,
        np.var(innov2[1:]),
        np.mean(np.abs(acf2)),
        acf2[0],
        acf2[4]
    ]
})

print(summary)

plt.figure(figsize=(10, 4))
plt.plot(y, label="Observed y_t")
plt.plot(y_pred1, label="Predicted: Model 1")
plt.plot(y_pred2, label="Predicted: Model 2")
plt.legend()
plt.title("Observed series vs one-step predictions")
plt.xlabel("t")
plt.ylabel("y")
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(innov1[1:], label="Innovation: Model 1")
plt.plot(innov2[1:], label="Innovation: Model 2")
plt.legend()
plt.title("Innovation sequences")
plt.xlabel("t")
plt.ylabel("innovation")
plt.tight_layout()
plt.show()

lags = np.arange(1, 16)
plt.figure(figsize=(10, 4))
plt.plot(lags, acf1, marker='o', label="Model 1")
plt.plot(lags, acf2, marker='o', label="Model 2")
plt.axhline(0, linewidth=1)
plt.legend()
plt.title("Innovation autocorrelation")
plt.xlabel("lag")
plt.ylabel("ACF")
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(z[0, :], label="True z1")
plt.plot(zhat2[0, :], label="Estimated z1")
plt.legend()
plt.title("Hidden context component 1: true vs estimated")
plt.xlabel("t")
plt.ylabel("value")
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(z[1, :], label="True z2")
plt.plot(zhat2[1, :], label="Estimated z2")
plt.legend()
plt.title("Hidden context component 2: true vs estimated")
plt.xlabel("t")
plt.ylabel("value")
plt.tight_layout()
plt.show()