#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 13:00:27 2020

@author: will
"""
import numpy as np
from statsmodels.tsa.api import VAR
from numba import jit, float64, vectorize, njit, int64
from BigVAR.BigVARAlgorithms import pen_loop

"""
Internal Function to create VAR lag matrix, not to be used by end-user
"""


def Zmat(Y, p, k, intercept=True, oos=False, contemp=False, offset=0):
    T = Y.shape[0]
    if oos and contemp:
        T += 1
    Y2 = np.array(Y)
    Y2 = np.fliplr(Y2)
    if contemp:
        p += 1
    Y2a = Y2[offset:(p+offset), :]
    Y2a = np.transpose(Y2a)
    Y2a = Y2a.reshape((-1, 1), order="F")
    M = T-p-offset
    if contemp:
        M += 1
    Z = np.zeros((k*p, M))
    Z[:, 0] = np.flip(Y2a).ravel()
    for i in range(1, M):
        Y1M = np.transpose(Y2[(i+offset):(p+i+offset), 0:k])
        Y1M = (np.flip(Y1M))
        Z[:, i] = Y1M.reshape((-1, 1), order="F").ravel()
    ones = np.ones((1, M))
    if intercept:
        ZF = np.concatenate((ones, Z), axis=0)
    else:
        ZF = Z
    return ZF


def VARXCons(Y, X, k, p, m=0, s=0, oos=False, contemp=False):
    if k == 1:
        Y = np.reshape(Y, (Y.shape[0], 1))
    if m == 1:
        X = np.reshape(X, (X.shape[0], 1))
    if s == 0:
        Z1 = Zmat(Y, p, k, True, oos, contemp)
        return Z1
    elif p == 0:
        Z1 = Zmat(X, s, m, True, oos, contemp)
        return Z1
    else:
        offsetX = 0
        offsetE = 0
        if p > s:
            offsetX = p-s
        else:
            offsetE = s-p
        Z1 = Zmat(Y, p, k, True, oos, False, offset=offsetE)
        Z2 = Zmat(X, s, m, False, oos, contemp, offsetX)
        ZZ = np.concatenate((Z1, Z2), axis=0)
        return ZZ


def LGSearch(lamstart, Y, Z, B, k, p, MN, alpha, C, intercept, tol, k1, m, s):
    lambdah = lamstart
    lambdal = 0
    while abs(lambdah-lambdal) > 10*tol:
        lam = (lambdah+lambdal)/2
        B = pen_loop(Z, Y, B, lam, tol, alpha)
        B2 = np.array(B[:, 1:B.shape[1], :], copy=True)
        if np.max(abs(B2)) < tol:
            lambdah = lam
        else:
            lambdal = lam
    return lambdah


def LambdaGridCons(gran1, gran2, Y, Z, p, k, MN, C, alpha, intercept, tol, separate_lambdas, linear, k1, m, s):
    B = np.zeros((k, k*p+1, 1))

    if not separate_lambdas:
        lamstart = np.max(np.matmul(Y.transpose(), Z.transpose()))
        lamstart = LGSearch(lamstart, Y, Z, B, k, p, MN,
                            alpha, C, intercept, tol, k1, m, s)
    else:
        lamstart = np.zeros(gran2)
        for i in range(1, k):
            lamstart[i] = np.max(np.matmul(Y[:, i].transpose(), Z.transpose()))
            lamstart[i] = LGSearch(
                lamstart[i], Y, Z, B, k, p, MN, alpha, C, intercept, tol, k1, m, s)
    if(linear):
        lam_grid = np.linspace(lamstart, lamstart/gran1, gran2)
    else:
        lam_grid = np.exp(np.linspace(np.log(lamstart),
                                      np.log(lamstart/gran1), gran2))
    return lam_grid.flatten()


@jit(nopython=True)
def huber(delta, x):
    if np.abs(x) < delta:
        l = .5*np.abs(x)**2
    else:
        l = delta*(np.abs(x)-0.5*delta)
    return l


@jit(nopython=True)
def calc_loss(x, loss="L2", delta=2.5):
    if(loss == "L1"):
        l = np.sum(np.abs(x))
    elif loss == "Huber":
        l = np.sum(huber(delta, x))
    else:
        l = np.linalg.norm(x)**2
    return l


def eval_mean(Y, T1, T2, loss, p):
    MSFE = []
#    k=Y.shape
    for u in range(T1, T2):
        trainY = Y[p:(u+p), :]
        yhat = trainY.mean(0)
        MSFE_temp = calc_loss(Y[u+p, :]-yhat, 'loss')
        MSFE.append(MSFE_temp)
    MSFE = np.array(MSFE)
    return np.mean(MSFE)


def eval_ar(Y, T1, T2, ic, p, loss):
    MSFE = []
    for u in range(T1, T2):
        trainY = Y[p:(u+p), :]
        var_mod = VAR(trainY)
        mod = var_mod.fit(maxlags=p, ic=ic)
        lag_order = mod.k_ar
        yhat = mod.forecast(trainY[-lag_order:], 1)
        MSFE_temp = calc_loss(Y[u+p, :]-yhat, loss)
        MSFE.append(MSFE_temp)
    MSFE = np.array(MSFE)
    return(np.mean(MSFE))


def BigVAR_eval(Y, Z, T1, T2, alpha, opt_lambda, h, p, tol, B, loss, k1):
    counter = 0
    MSFE = []
    for v in range(T1, T2):
        trainY = Y[(p):(v+p), 0:k1]
        trainZ = Z[1:Z.shape[0], 0:(v)]
        B = pen_loop(trainZ, trainY, B, opt_lambda, tol, alpha)
        eZ = Z[:, [v]]
        preds = np.matmul(B[:, :, 0], eZ)
        MSFE.append(calc_loss(Y[v+p, :]-preds[:, 0], loss))
        counter += 1
    MSFE_mean = np.mean(np.array(MSFE))
    return(MSFE)

# convert coefficient matrix to VAR(1)


def CreateCoefMat(B, p, k):
    A = np.zeros((k*p, k*p))
    A[0:k, :] = B
    A[(k):A.shape[0], 0:(A.shape[1]-k)] = np.identity(k*p-k)
    np.linalg.eigvals(A)
    if np.max(np.absolute(np.linalg.eigvals(A))) > 1:
        raise ValueError("Generator Matrix is not stationary")
    return(A)


def MultVARSim(A, p, k, Sigma, T):
    if not np.allclose(Sigma.flatten(), Sigma.T.flatten()):
        raise ValueError("Sigma must be square and symmetric")
    Y = np.zeros((T+500+p, k))
    YY = Y.flatten()
    for i in range(k*p, Y.shape[0]*k, k):
        u = np.concatenate((np.random.multivariate_normal(
            np.zeros(k), Sigma, 1).flatten(), np.zeros(k*p-k))).flatten()
        YY[(i-k*p+k):(i+k)] = A@np.flipud(YY[(i-k*p):(i)])+u
        YY[(i-k*p+k):(i+k)] = np.flipud(YY[(i-k*p+k):(i+k)])
    YY = YY[YY != 0]
    N = int(YY.shape[0]/k)
    YY = YY.reshape((N, k))
    YY = np.flipud(YY)
    YY = YY[500:YY.shape[0], :]
    return(YY)
