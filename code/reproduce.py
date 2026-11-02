#!/usr/bin/env python3
"""Reproduce all new numerical results and read the publisher's source workbooks.

Run from any working directory: python code/reproduce.py
The xlsx files are read directly as OOXML; no workbook is modified.
Outputs are written to results/ and figures/. No network connection is used.
"""
from __future__ import annotations
import csv, hashlib, json, platform, sys, zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from dataclasses import asdict, replace
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.integrate import solve_ivp, simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from accessibility import Parameters, riccati, recovery, accessible, constant_speed_work, binary_condition, optimal_path, discrete_optimum

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; FIG=ROOT/'figures'
OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
RNS='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'

def read_xlsx(path:Path):
    """Read cached numeric/text values, preserving blank columns and sheet names."""
    with zipfile.ZipFile(path) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            ss=ET.fromstring(z.read('xl/sharedStrings.xml'))
            strings=[''.join(e.itertext()) for e in ss]
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        targets={r.attrib['Id']:r.attrib['Target'] for r in rels}
        book=ET.fromstring(z.read('xl/workbook.xml'))
        sheets={}
        for sh in book.find('m:sheets',NS):
            target=targets[sh.attrib[RNS+'id']]
            name=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            root=ET.fromstring(z.read(name)); rows=[]
            for row in root.findall('m:sheetData/m:row',NS):
                vals={}
                for cell in row:
                    ref=cell.attrib.get('r','A1'); letters=''.join(c for c in ref if c.isalpha()); idx=0
                    for c in letters:idx=26*idx+ord(c)-64
                    v=cell.find('m:v',NS); typ=cell.attrib.get('t')
                    if typ=='inlineStr': value=''.join(cell.find('m:is',NS).itertext())
                    elif v is None: value=None
                    elif typ=='s':value=strings[int(v.text)]
                    elif typ in ('str','e'):value=v.text
                    else:
                        try:value=float(v.text)
                        except (TypeError,ValueError):value=v.text
                    vals[idx-1]=value
                rows.append([vals.get(i) for i in range(max(vals,default=-1)+1)])
            sheets[sh.attrib['name']]=rows
        return sheets

def save_csv(name,rows):
    with (OUT/name).open('w',newline='') as f:csv.writer(f).writerows(rows)

checks=[]
def check(name,value,tolerance):
    value=float(value)
    ok=bool(np.isfinite(value) and value<=tolerance)
    checks.append({'name':name,'residual':value,'tolerance':tolerance,'passed':ok})
    if not ok:raise AssertionError(checks[-1])

p=Parameters(); tau=3.
# Independent scalar ODE verification of the Riccati closed form.
for ratio in [.1,1.,10.,p.gammab/p.gamma,100.]:
    pr=replace(p,gammab=ratio*p.gamma); a=pr.kb/pr.gammab
    t=np.linspace(0,10,201)
    sol=solve_ivp(lambda t,y:[2*pr.kb*a-2*a*y[0]-y[0]**2/(2*pr.gamma)],(0,10),[pr.kb],t_eval=t,rtol=1e-11,atol=1e-14)
    check(f'Riccati_ODE_ratio_{ratio:g}',np.max(abs(sol.y[0]-riccati(t,pr)))/pr.kb,2e-9)
    eta=recovery(t,pr)
    check(f'eta_monotonic_ratio_{ratio:g}',max(0.,-np.min(np.diff(eta))),1e-12)
    check(f'eta_range_ratio_{ratio:g}',max(0.,-np.min(eta),np.max(eta)-1),1e-12)
# Direct energy and work integrals including both endpoint trap jumps.
for mean in [(20.,80.),(-20.,5.),(30.,30.),(0.,0.)]:
    for tauf in [.1,1.,3.,10.]:
        t,x,y,ell,v=optimal_path(mean,tauf,p,n=6001)
        d=x-y; a=p.kb/p.gammab; ph=riccati(tauf-t,p)
        # ell=x+kb*d/(2*k) along the optimum only at final horizon; use derivative fully.
        dp_h=2*p.kb*a-2*a*ph-ph**2/(2*p.gamma)
        dd=v-a*d
        dv=dp_h*d/(2*p.gamma)-ph*dd/(2*p.gamma)
        ell_dot=v+(p.gamma*dv+p.kb*dd)/p.k
        w0=p.k*mean[0]*ell[0]-.5*p.k*ell[0]**2
        wf=.5*p.k*(x[-1]-ell[-1])**2
        w=w0+simpson(p.k*(x-ell)*ell_dot,x=t)+wf
        exact=float(accessible(mean,tauf,p))
        check(f'optimal_work_integral_{mean}_{tauf:g}',abs(w-exact)/p.kBT,2e-8)
        e0=.5*p.k*mean[0]**2+.5*p.kb*(mean[0]-mean[1])**2
        residual=.5*p.kb*d[-1]**2
        loss=simpson(p.gamma*v*v+p.gammab*(a*d)**2,x=t)
        check(f'energy_balance_{mean}_{tauf:g}',abs(e0-residual-loss-exact)/p.kBT,2e-8)
# Exact interval propagation and independent quadratic optimization.
conv=[['N_intervals','W_zJ','analytic_zJ','absolute_error_zJ','relative_error']]
prev=np.inf; mean=(20.,80.); exact=float(accessible(mean,tau,p))
for n in [5,10,20,40,80,160]:
    val,_,mineig=discrete_optimum(mean,tau,n,p)
    err=exact-val
    check(f'discrete_lower_bound_N{n}',max(0.,-err),1e-9)
    check(f'discrete_convergence_N{n}',max(0.,err-prev),1e-9)
    prev=err
    conv.append([n,val,exact,err,err/exact])
check('discrete_relative_error_N160',prev/exact,5e-4)
save_csv('control_convergence.csv',conv)
# Gaussian equilibrium, FDT, and two-threshold conditioning checks.
check('equilibrium_Lyapunov',np.max(abs(p.A@p.C+p.C@p.A.T-2*p.kBT*np.diag([1/p.gamma,1/p.gammab]))),1e-9)
branch=[['state','probability','mean_x_nm','mean_y_nm','extension_nm','B_kBT','hidden_increment_kBT','W_opt_kBT','W_source_speed_kBT','V_initial_kBT','source_speed_nm_s']]
means=[]; probs=[]; second=[]
for i,j in [(0,0),(0,1),(1,0),(1,1)]:
    prob,mu,x2=binary_condition(1.,i,j,p=p)
    speed=25. if (i,j)==(1,0) else 45.
    b=.5*p.k*mu[0]**2; w=float(accessible(mu,tau,p)); wc,_=constant_speed_work(mu,tau,speed,p)
    branch.append([f'{i}{j}',prob,*mu,mu[0]-mu[1],b/p.kBT,(w-b)/p.kBT,w/p.kBT,wc/p.kBT,.5*p.k*x2/p.kBT,speed])
    means.append(mu);probs.append(prob);second.append(x2)
    check(f'constant_speed_below_optimum_{i}{j}',max(0.,wc-w),1e-10)
check('conditional_probabilities_sum',abs(sum(probs)-1),1e-10)
check('conditional_means_recombine',np.max(abs(np.array(probs)@np.array(means))),1e-8)
check('conditional_second_moment_recombines',abs(np.dot(probs,second)-p.C[0,0])/p.C[0,0],1e-10)
save_csv('binary_model_results.csv',branch)
# Exact measurement refinement identity, comparing two binary records to final bit only.
fine=np.array(means); probs=np.array(probs); finework=sum(probs[k]*accessible(fine[k],tau,p) for k in range(4))
coarsework=0.; decomposition=0.
for ids in ([0,2],[1,3]):
    w=probs[ids]; prob=w.sum(); mu=np.sum(w[:,None]*fine[ids],axis=0)/prob
    coarsework+=prob*accessible(mu,tau,p)
    decomposition+=.5*p.k*np.sum(w*(fine[ids,0]-mu[0])**2)
    d=fine[ids,0]-fine[ids,1]; dm=mu[0]-mu[1]
    decomposition+=.5*p.kb*recovery(tau,p)*np.sum(w*(d-dm)**2)
check('measurement_refinement_identity',abs(finework-coarsework-decomposition)/p.kBT,1e-11)
# Analytic recovery fractions and physical scales.
eta_inf=(np.sqrt(1+p.gammab/p.gamma)-1)/(np.sqrt(1+p.gammab/p.gamma)+1)
meta={'parameters_nm_s_zJ':asdict(p),'eta_3s':float(recovery(3,p)),'eta_infinity':float(eta_inf),
      'eigen_relaxation_times_s':sorted((1/np.linalg.eigvals(p.A)).tolist()),'bath_bare_time_s':p.gammab/p.kb,
      'sigma_x_nm':float(np.sqrt(p.C[0,0])),'threshold_over_sigma':float(70/np.sqrt(p.C[0,0])),
      'average_fine_work_kBT':float(finework/p.kBT),'average_coarse_work_kBT':float(coarsework/p.kBT),
      'refinement_gain_kBT':float(decomposition/p.kBT)}
(OUT/'model_results.json').write_text(json.dumps(meta,indent=2))
# Publisher-source inventory, preserving original column identities and limiting Fig.4 to 3 s.
bookmeta=[]; peakrows=[['workbook','sheet','column','header_as_stored','peak_time_s','peak_value_as_stored','adjacent_error_as_stored','n_points_0_to_3s']]
for no in range(3,7):
    path=ROOT/'data'/f'41567_2026_3442_MOESM{no}_ESM.xlsx'
    if not path.exists():raise FileNotFoundError(path)
    sheets=read_xlsx(path)
    bookmeta.append({'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'sheets':[{ 'name':name,'data_rows':len(rows)-1,'header':rows[0]} for name,rows in sheets.items()]})
    if no==6:
        for name in ['Figure_4a','Figure_4b','Inset_4b']:
            rows=sheets[name]; header=rows[0]
            for col in [1,3,5,7]:
                points=[r for r in rows[1:] if len(r)>col+1 and isinstance(r[0],(int,float)) and 0<=r[0]<=3 and isinstance(r[col],(int,float))]
                r=max(points,key=lambda r:r[col])
                peakrows.append([path.name,name,chr(65+col),header[col],r[0],r[col],r[col+1],len(points)])
    if no==5:
        rr=sheets['Figure_3d']; r=max(rr[1:],key=lambda r:r[1])
        meta['source_Phi_peak']={'delay_s':r[0],'value_s':r[1],'row_as_stored':r,'header':rr[0]}
save_csv('publisher_source_peaks.csv',peakrows)
(OUT/'source_inventory.json').write_text(json.dumps(bookmeta,indent=2,ensure_ascii=False))
(OUT/'model_results.json').write_text(json.dumps(meta,indent=2))
# Figures: one plot per figure, default cycle, no manually selected colors.
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.labelsize':11,'legend.fontsize':9})
fig,ax=plt.subplots(figsize=(6.2,3.6)); tt=np.geomspace(.003,100,350)
for ratio in [1.,10.,p.gammab/p.gamma,100.]:
    pr=replace(p,gammab=ratio*p.gamma)
    ax.plot(tt,recovery(tt,pr),label=rf'$\gamma_b/\gamma={ratio:.2f}$')
ax.set_xscale('log');ax.set_xlabel(r'Extraction duration $\tau$ (s)');ax.set_ylabel(r'Hidden-mean recovery fraction $\eta(\tau)$');ax.set_ylim(0,1);ax.legend(loc='upper left');ax.grid(alpha=.2);fig.tight_layout();fig.savefig(FIG/'recovery_fraction.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(6.2,3.6)); ds=np.geomspace(.02,100,80);scan=[['delay_s','state','probability','W_opt_kBT','B_kBT','hidden_increment_kBT']]
for i,j in [(0,0),(0,1),(1,0),(1,1)]:
    vals=[]
    for dt in ds:
        prob,mu,_=binary_condition(float(dt),i,j,p=p)
        b=.5*p.k*mu[0]**2;w=float(accessible(mu,3,p)); vals.append((w-b)/p.kBT)
        scan.append([dt,f'{i}{j}',prob,w/p.kBT,b/p.kBT,(w-b)/p.kBT])
    ax.plot(ds,vals,label=f'$m={i}{j}$')
ax.set_xscale('log');ax.set_xlabel(r'Measurement separation $\delta$ (s)');ax.set_ylabel('Optimal hidden-mean increment ($k_B T$)');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(FIG/'branch_hidden_gain.pdf');plt.close(fig);save_csv('delay_scan.csv',scan)
fig,ax=plt.subplots(figsize=(6.2,3.6));t=np.geomspace(.001,100,500);c=4.
q=np.array([1+(c-1)*p.k/p.kb*expm(-p.A*s)[0,1]**2 for s in t]);kl=.5*(q-1-np.log(q));h=.5*(c-1-np.log(c))
ax.plot(t,kl,label=r'$D(\rho_t^X\Vert q^X)$');ax.axhline(h,linestyle='--',label=r'$\mathcal{H}/(k_BT)$');ax.axhline(0,linestyle=':',label=r'$W_{\rm acc}/(k_BT)=0$');ax.set_xscale('log');ax.set_xlabel('Relaxation time $t$ (s)');ax.set_ylabel('Relative entropy or normalized work');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(FIG/'covariance_counterexample.pdf');plt.close(fig)
check('covariance_KL_contraction',max(0.,float(np.max(kl))-h),1e-12)
report={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'count':len(checks),'passed':sum(c['passed'] for c in checks),'checks':checks}
(OUT/'verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps({'checks':report['count'],'passed':report['passed'],'model':meta},indent=2))
