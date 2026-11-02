#!/usr/bin/env python3
"""Exact translation-accessibility calculations for a thermal two-coordinate model.

Units used internally: nanometres, seconds, and zeptojoules (1 zJ = 1e-21 J).
A stiffness of 1 microN/m is 0.001 zJ/nm^2; a friction of
1 microN*s/m is 0.001 zJ*s/nm^2. No experimental parameter is fitted.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.special import ndtr

@dataclass(frozen=True)
class Parameters:
    k: float = 0.00147
    kb: float = 0.0018
    gamma: float = 0.00034
    gammab: float = 0.015
    kBT: float = 1.380649e-23 * 298 * 1e21
    def __post_init__(self):
        if min(self.k,self.kb,self.gamma,self.gammab,self.kBT)<=0:
            raise ValueError('All parameters must be strictly positive.')
    @property
    def H(self): return np.array([[self.k+self.kb,-self.kb],[-self.kb,self.kb]])
    @property
    def A(self): return np.diag([1/self.gamma,1/self.gammab]) @ self.H
    @property
    def C(self): return self.kBT*np.linalg.inv(self.H)

def riccati(h: float | np.ndarray, p: Parameters=Parameters()):
    h=np.asarray(h,dtype=float)
    if np.any(h<0): raise ValueError('The horizon must be nonnegative.')
    a=p.kb/p.gammab; omega=np.sqrt(a*a+p.kb*a/p.gamma)
    pp=2*p.gamma*(omega-a); pm=-2*p.gamma*(omega+a)
    ratio=(p.kb-pp)/(p.kb-pm)*np.exp(-2*omega*h)
    return (pp-ratio*pm)/(1-ratio)

def recovery(h,p=Parameters()): return 1-riccati(h,p)/p.kb

def accessible(mean, horizon, p=Parameters()):
    x,y=np.asarray(mean,dtype=float)
    return 0.5*p.k*x*x+0.5*(p.kb-riccati(horizon,p))*(x-y)**2

def constant_speed_work(mean,horizon,speed,p=Parameters()):
    B=np.array([[p.k/p.gamma,p.kb/p.gamma],
                [p.k/p.gamma,p.kb*(1/p.gamma+1/p.gammab)]])
    E=expm(-B*horizon); I=np.eye(2); e=np.array([1.,0.])
    state=np.array([mean[0],mean[0]-mean[1]])
    L=e@np.linalg.solve(B,(I-E)@state)
    D=e@(horizon*np.linalg.solve(B,e)-np.linalg.solve(B,np.linalg.solve(B,(I-E)@e)))
    if speed is None:
        speed=L/(2*D) if D>0 else 0.
    return float(p.k*(speed*L-speed*speed*D)),float(speed)

def binary_condition(delay:float,i:int,j:int,threshold=70.,p=Parameters()):
    """Stationary Gaussian conditioning by two binary coordinate observations.

    Returns exact quadrature probabilities, conditional means (x,y), and
    second moment of the observed coordinate. The conditional state itself
    is not assumed Gaussian after thresholding.
    """
    if delay<=0 or i not in (0,1) or j not in (0,1): raise ValueError('Invalid measurement specification.')
    C=p.C; E=expm(-p.A*delay); sx=np.sqrt(C[0,0]); r=(E@C)[0,0]/C[0,0]
    r=float(np.clip(r,-1+1e-14,1-1e-14)); sd=np.sqrt(1-r*r); h=threshold/sx
    lo,hi=(-np.inf,h) if i==0 else (h,np.inf); sign=1 if j==1 else -1
    phi=lambda u:np.exp(-u*u/2)/np.sqrt(2*np.pi)
    def integ(u,order):
        z=(h-r*u)/sd; prob=ndtr(-z) if j else ndtr(z); f=phi(z)
        if order==0: val=prob
        elif order==1: val=u*prob
        elif order==2: val=r*u*prob+sign*sd*f
        else: val=((r*u)**2+sd**2)*prob+sign*sd*(r*u+h)*f
        return phi(u)*val
    vals=[quad(lambda u:integ(u,o),lo,hi,epsabs=2e-12,epsrel=2e-11,limit=250)[0] for o in range(4)]
    P=vals[0]
    if P<=0: raise ArithmeticError('Conditional branch below numerical resolution.')
    mobs=sx*np.array(vals[1:3])/P
    Cobs=C[0,0]*np.array([[1.,r],[r,1.]])
    cross=np.column_stack(((E@C)[:,0],C[:,0]))
    mean=cross@np.linalg.solve(Cobs,mobs)
    return float(P),mean,float(sx*sx*vals[3]/P)

def optimal_path(mean,horizon,p=Parameters(),n=1001):
    x0,y0=map(float,mean); a=p.kb/p.gammab
    def rhs(t,s):
        x,d=s; v=-float(riccati(horizon-t,p))/(2*p.gamma)*d
        return [v,v-a*d]
    ts=np.linspace(0,horizon,n)
    sol=solve_ivp(rhs,(0,horizon),(x0,x0-y0),t_eval=ts,rtol=2e-11,atol=2e-11)
    if not sol.success:raise RuntimeError(sol.message)
    x,d=sol.y; v=-riccati(horizon-ts,p)*d/(2*p.gamma)
    lam=x+(p.gamma*v+p.kb*d)/p.k
    return ts,x,x-d,lam,v

def discrete_optimum(mean,horizon,N,p=Parameters()):
    """Independent exact optimization over N hold intervals and N+1 trap jumps.

    Each jump has a linear expected work contribution plus a quadratic cost.
    Coordinates relax by the exact matrix exponential between jumps.
    The final trap centre is free; no finite-difference dynamics is used.
    """
    E=expm(-p.A*horizon/N); one=np.ones(2); gain=(np.eye(2)-E)@one
    m=np.asarray(mean,float); sensitivity=np.zeros((2,N+1)); Q=np.zeros((N+1,N+1)); b=np.zeros(N+1)
    for j in range(N+1):
        # Absolute trap-centre variables lambda_j (after jump j).
        row=sensitivity[0].copy()
        if j:row[j-1]-=1
        b[j]+=p.k*m[0]
        if j:b[j-1]-=p.k*m[0]
        # W_j=k*m_x*(lambda_j-lambda_prev)-k/2*(lambda_j^2-lambda_prev^2).
        a=np.zeros(N+1);a[j]=1
        if j:a[j-1]=-1
        # sensitivity here is m_x, not r; use exact energy-jump expression.
        sx=sensitivity[0]
        Q-=p.k*(np.outer(sx,a)+np.outer(a,sx))
        Q[j,j]+=p.k
        if j:Q[j-1,j-1]-=p.k
        if j<N:
            m=E@m; sensitivity=E@sensitivity; sensitivity[:,j]+=gain
    eig=np.linalg.eigvalsh(Q)
    if eig.min()<=0:raise ArithmeticError('Discrete work-cost matrix is not positive definite.')
    lam=np.linalg.solve(Q,b)
    return float(0.5*b@lam),lam,float(eig.min())
