#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 10:37:11 2020

@author: will
"""
import os
import numpy as np
import pdb
from numba import jit, njit, b1
from numba.typed import List


"""
Soft-thresholding operator - computes one-dimensional solution to the lasso problem
"""


@jit(nopython=True)
def ST(x, lam):
    if x > lam:
        return(x-lam)
    if x < -lam:
        return(x+lam)
    else:
        return float(0)


"""
inner kernel of lasso function (coordinate descent) with 3 for loops
"""


def lasso(Z, Y, lam, maxiter, tol, B, znorm2):
    pk = Z.shape[0]
    k = Y.shape[0]
    if B.shape[0] == 0:
        B = np.zeros((k, pk))
    for i in range(maxiter):
        BOLD = np.array(B, copy=True)
        for j in range(k):
            for m in range(pk):
                BS = B[j, :]
                YS = Y[j, :]
                r = YS-np.matmul(np.delete(BS, m), np.delete(Z, m, 0))
                test = np.asscalar(np.sum(np.squeeze(np.asarray(r))*Z[m, :]))
                B[j, m] = ST(test, lam)/znorm2[m]
        if np.max(abs(BOLD-B)) < tol:
            break

    return B


"""
elastic net- same thing but with alpha penalty
"""


@njit
def elastic_net(Z, Y, lam, maxiter, tol, B, znorm2, alpha):
    pk = Z.shape[0]
    k = Y.shape[0]
    znorm2 = List(znorm2)
    if B.shape[0] == 0:
        B = np.zeros((k, pk))
    for i in range(maxiter):
        BOLD = np.zeros((B.shape[0], B.shape[1]))
        for ii in range(B.shape[0]):
            for jj in range(B.shape[1]):
                BOLD[ii, jj] = B[ii, jj]
        for j in range(k):
            for m in range(pk):
                BS = B[j, :]
                YS = Y[j, :]
                ind = np.ones(BS.shape[0], dtype=b1)
                ind[m] = 0
                ind2 = np.ones(Z.shape[0], dtype=b1)
                ind2[m] = 0
                r1 = YS-np.dot(BS[ind], Z[ind2, :])
                # to convert to scalar# np.sum(np.squeeze(np.asarray(r1))*Z[m,:])
                test = np.dot(r1, Z[m, :]).item()
                B[j, m] = ST(test, lam*alpha)/(1+lam*(1-alpha))/znorm2[m]
        chg = B-BOLD
        if np.linalg.norm(chg, np.inf) < tol:
            break
    return B


"""
full estimation function
"""


def pen_loop(Z, Y, B, lambda_grid, tol, alpha=1):
    Y = np.transpose(Y)
    # subtract mean (added later as an intercept)
    Ymeans = np.reshape(np.mean(Y, axis=1), (Y.shape[0], 1))
    Zmeans = np.reshape(np.mean(Z, axis=1), (Z.shape[0], 1))
    Y = Y-np.matmul(Ymeans, np.ones((1, Y.shape[1])))
    Z = Z-np.matmul(Zmeans, np.ones((1, Z.shape[1])))

    znorm2 = List(np.sum(pow(Z, 2), axis=1))
    if lambda_grid.ndim == 0:
        lambda_grid = lambda_grid[None]
        nlambdas = 1
    else:
        nlambdas = lambda_grid.shape[0]
    for i in range(nlambdas):
        if alpha == 1:
            BOLD = np.matrix(B[:, 1:B.shape[1], i], copy=True)
            B[:, 1:B.shape[1], i] = lasso(Z, Y, np.asscalar(
                lambda_grid[i]), 10000, tol, BOLD, znorm2)
        else:
            B[:, 1:B.shape[1], i] = elastic_net(Z, Y, np.asscalar(
                lambda_grid[i]), 10000, tol, B[:, 1:B.shape[1], i], znorm2, alpha)
        nu = np.subtract(Ymeans, np.matmul(B[:, 1:B.shape[1], i], Zmeans))
        B[:, 0, i] = np.reshape(nu, (B.shape[0]))
    return B
