#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: will
"""
# VARX={'k':5,'s':20,'contemp':True}
import numpy as np
from BigVAR.BigVARSupportFunctions import VARXCons, LambdaGridCons, calc_loss, BigVAR_eval, eval_mean, eval_ar
from BigVAR.BigVARAlgorithms import pen_loop


# Right now, just going to support Lasso/Elastic Net
supported_structs = ['Basic', 'Basic_EN']
supported_cv = ['Rolling']


class BigVAR():

    ''' 
    Y - T x k multivariate time series
    p - lagmax, maximal lag order for modeled series
    struct - Penalty Structure
    gran - granularity of penalty grid
    T1 - index of time series to start CV
    T2 index of time series to start Forecast evaluation
    RVAR - indicator for relaxed VAR
    h - desired forecast horizon
    cv - validation procedure
    MN -Minnesota Prior Indicator
    verbose - indicator for verbose output
    IC - indicator for including AIC and BIC benchmarks
    VARX - VARX model specifications
    ONESE - 'One Standard Error' heuristic
    own_lambdas - user-supplied lambdas
    tf - transfer function
    alpha - grid of candidate alpha values (applies only to sparse VARX-L models)
    recursive - whether recursive multi-step forecasts are used (applies only to muliple horizon VAR models)
    C - vector of coefficients to shrink toward random walk
    dates - optional vector of dates corresponding to Y
    '''

    def __init__(self, Y, p, struct, gran, T1, T2, RVAR=False, h=1, cv='Rolling', MN=False, verbose=True, IC=True,
                 VARX={}, ONESE=False, own_lambdas=False, alpha=np.array([]), recursive=False, C=np.array([]), intercept=True, separate_lambdas=False, loss="L2"):

        if Y.shape[1] > Y.shape[0]:
            raise ValueError('k > T!')
        if p < 0:
            raise ValueError('p must be >= 0')
        #if p == 0 and struct !='Basic': raise ValueError('Only Basic VARX-L supports a transfer function')
        if struct not in supported_structs:
            raise ValueError(
                'penalty structure must be one of {}'.format(supported_structs))
        if h < 1:
            raise ValueError('h must be greater than 1!')
        if cv not in supported_cv:
            raise ValueError(
                'Cross-Validation must be one of {}'.format(supported_cv))
        if len(gran) != 2 and not own_lambdas:
            raise ValueError('granularity must have two parameters')
        if gran[0] <= 0 or gran[1] <= 0:
            raise ValueError('granularity parameters must be positive')

        if len(VARX) != 0:
            k = VARX['k']

            if k > Y.shape[1]:
                raise ValueError(
                    'k is greater than the number of columns in Y')
        else:
            k = Y.shape[1]
        self.m = Y.shape[1] - k
        self.n_series = Y.shape[1] - (self.m if self.m < Y.shape[1] else 0)
        self.tf = (p == 0)
      #  if self.n_series == 1 and struct not in ['Basic', 'Lag', 'HVARC']:
      #      raise ValueError('Univariate support is only available for Lasso, Lag Group, and Componentwise HVAR')
      #  if len(VARX) == 0 and struct=='EFX': raise ValueError('EFX is only supported in the VARX framework')
        # TODO check for contemporaneous dependence
      #  structs = ['HVARC', 'HVAROO', 'HVARELEM']
      #  if len(VARX) != 0 and struct in structs: raise ValueError('EFZ is the only nested model supported in the VARX framework')
        if T1 > Y.shape[0] or T2 > Y.shape[0] or T2 < T1:
            raise ValueError('Training dates exceed series length')

        # TODO verify VARX specifications entered correctly

        #if len(alpha) > 0 and any(a < 0 for a in alpha) and any(a > 1 for a in alpha): raise ValueError('alpha must be [0,1]')
        if len(C) != 0:
            if len(C) != k:
                raise ValueError('C must have length k')
            if not all(c == 0 or c == 1 for c in C):
                raise ValueError('Values of C must be either 0 or 1')
        else:
            self.C = [1]*k
        # TODO add logic for dates

        self.Y = Y
        self.p = p
        self.struct = struct
        self.gran = gran
        self.T1 = T1
        self.T2 = T2
        self.RVAR = RVAR
        self.h = h
        self.cv = cv
        self.MN = MN
        self.verbose = verbose
        self.IC = IC
        self.VARX = VARX
        self.ONESE = ONESE
        self.own_lambdas = own_lambdas
        self.alpha = alpha
        self.recursive = recursive
        self.intercept = intercept
        self.separate_lambdas = separate_lambdas
        self.loss = loss


def rolling_validate(self):

   # if len(self.alpha) == 0:
    #    if len(self.VARX) > 0:
    #      alpha = 1/(self.VARX/(self.k) + 1)
    #   else:
    #     alpha = 1/(self.k + 1)

   # dual = len(self.alpha > 1) and self.struct in ['SparseLag', 'SparseOO']
  #  jj = 0
    Y = self.Y
    alpha = self.alpha
    T1 = self.T1 if self.cv == 'Rolling' else self.p + 2
    T2 = self.T2
    MN = self.MN
    C = self.C
    struct = self.struct
    VARX_list = self.VARX
    s = self.VARX['s'] if len(self.VARX) != 0 else 0
    if self.own_lambdas:
        gamm = self.gran
        gran2 = len(gamm)
    ONESE = self.ONESE
    intercept = self.intercept
    separate_lambdas = self.separate_lambdas
    p = self.p
    loss = self.loss
    tol = 1e-4
    h = 1
    cv = self.cv
    verbose = self.verbose
    IC = self.IC
    RVAR = self.RVAR

    if (self.cv == 'Rolling'):
        T1 = T1 - np.max([self.p, s])
        T2 = T2 - np.max([self.p, s])
    if not self.own_lambdas:
        gran = self.gran
        gran2 = self.gran[1]
        gran1 = self.gran[0]

    MSFE = np.zeros((len(np.arange(T1-h, T2)), gran2))
    # constructing a lag matrix in VARX setting
    if len(self.VARX) != 0:
        VARX = True
        k1 = self.VARX['k']
        s = self.VARX['s']
        contemp = self.VARX.get('contemp')
        if contemp == None:
            contemp = False

        s1 = 0

        m = self.Y.shape[1] - k1  # k - k1
        Z1 = VARXCons(Y[:, 0:k1], Y[:, k1:], k1, self.p, m, s, False, contemp)
        trainZ = Z1[1:Z1.shape[0], 0:T2]
        trainY = np.matrix(self.Y[p:T2+p, 0:k1])
        lambdas = LambdaGridCons(gran1, gran2, trainY, trainZ, p,
                                 Y.shape[1], MN, C, alpha, intercept, tol, separate_lambdas, False, k1, m, s)
        B = np.zeros((k1, k1*p+m*s+1, gran2))

    else:  # VAR estimation
        contemp = False
        k = Y.shape[1]
        k1 = k
        Z1 = VARXCons(Y, np.matrix(
            np.zeros(shape=Y.shape)), Y.shape[1], p, 0, 0)
        trainZ = Z1[1:Z1.shape[0], 0:T2]
        trainY = np.matrix(Y[p:Y.shape[0], :], copy=True)
        trainY = np.matrix(trainY[0:T2, :], copy=True)
        lambdas = LambdaGridCons(gran1, gran2, trainY, trainZ, p, k,
                                 MN, C, alpha, intercept, tol, separate_lambdas, False, 0, 0, 0)
        B = np.zeros((k, k*p+1, gran2))

    counter = 0
    for v in range(T1-h, T2):
        trainY = Y[(p):(v+p), 0:k1]
        trainZ = Z1[1:Z1.shape[0], 0:(v)]
        B = pen_loop(trainZ, trainY, B, lambdas, tol, alpha)
        eZ = Z1[:, [v]]
        for j in range(lambdas.shape[0]):
            preds = np.matmul(B[:, :, j], eZ)
            MSFE[counter, j] = calc_loss(Y[v+p, :]-preds[:, 0], loss)
        counter += 1
    MSFE_mean = np.mean(MSFE, axis=0)
    index = np.argmin(MSFE_mean)
    opt_lambda = lambdas[index]
    oos_msfe = BigVAR_eval(
        Y, Z1, T2, Z1.shape[1], alpha, opt_lambda, h, p, tol, B[:, :, [index]], loss, k1)

    oos_aic = eval_ar(Y, T2, Z1.shape[1], 'aic', p, loss)

    oos_bic = eval_ar(Y, T2, Z1.shape[1], 'bic', p, loss)

    oos_mean = eval_mean(Y, T2, Z1.shape[1], loss, p)

    result = BigVAR_results(Y, p, struct, gran, T1, T2, RVAR, h, cv, MN, verbose, IC,
                            VARX_list, ONESE, self.own_lambdas, alpha, self.recursive, C, intercept, separate_lambdas, loss, oos_msfe, index, B, lambdas, opt_lambda, MSFE_mean, oos_aic, oos_bic, oos_mean)

    print(result)
    return(result)


class BigVAR_results(BigVAR):
    def __init__(self, Y, p, struct, gran, T1, T2, RVAR, h, cv, MN, verbose, IC,
                 VARX, ONESE, own_lambdas, alpha, recursive, C, intercept, separate_lambdas, loss, oos_msfe, index, B, lambdas, opt_lambda, is_msfe, oos_aic, oos_bic, oos_mean):
        super().__init__(Y, p, struct, gran, T1, T2, RVAR, h, cv, MN, verbose, IC,
                         VARX, ONESE, own_lambdas, alpha, recursive, C, intercept, separate_lambdas, loss)
        self.oos_msfe = oos_msfe
        self.index = index,
        self.B = B
        self.lambdas = lambdas
        self.opt_lambda = opt_lambda
        self.is_msfe = is_msfe
        self.oos_aic = oos_aic
        self.oos_bic = oos_bic
        self.oos_mean = oos_mean

    def __str__(self):
        str_res = ""
        str_res = str_res+"BigVAR.Results\n"
        str_res = str_res+'loss: ' + self.loss+"\n"
        str_res = str_res + \
            'In-sample loss: {}'.format(self.is_msfe[self.index])+"\n"

        str_res = str_res+'index: ' + str(self.index)+"\n"
        str_res = str_res+'optimal lambda: '+str(self.opt_lambda)+"\n"
        str_res = str_res + \
            'Out-of-sample loss: {}'.format(np.mean(self.oos_msfe))+"\n"
        str_res = str_res + \
            'Out-of-sample AIC benchmark: {}'.format(
                np.mean(self.oos_aic))+"\n"
        str_res = str_res + \
            'Out-of-sample BIC benchmark: {}'.format(
                np.mean(self.oos_bic))+"\n"
        str_res = str_res + \
            'Out-of-sample mean benchmark: {}'.format(
                np.mean(self.oos_mean))+"\n"
        return(str_res)
