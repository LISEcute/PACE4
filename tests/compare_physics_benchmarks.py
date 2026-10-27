"""Compare reproducible PACE4 and portable FORTRAN suites without assuming agreement."""
import argparse
from collections import Counter
import csv
import html
import json
import math
from pathlib import Path
import numpy as np
from scipy.stats import chisquare, chi2_contingency
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ELEMENTS={25:'Mn',26:'Fe',27:'Co',28:'Ni'}
LABELS={'ni60_fixed_j10':'Fixed J = 10','ni60_equal_j0_10_20':'Equal J = 0, 10, 20','ni60_generated_j0_20':'Weights 2J+1, J = 0 to 20'}
def nuclide(k):
 if k=='fission': return 'Fission'
 z,n=map(int,k.split(',')); return str(z+n)+ELEMENTS.get(z,'Z'+str(z))
def aggregate(suite,name,field):
 runs=[r for r in suite['runs'] if r['case']['name']==name]; counts=Counter()
 for r in runs: counts.update(r[field])
 return counts,sum(r['events'] for r in runs),sum(r['sigma_input_mb'] for r in runs)/len(runs)
def wilson(k,n):
 z=1.959963984540054; p=k/n; den=1+z*z/n
 mid=(p+z*z/(2*n))/den; rad=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return max(0,mid-rad),min(1,mid+rad)
def ptext(v): return '<1e-300' if v<1e-300 else f'{v:.4g}'
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('root',type=Path); ap.add_argument('--baseline',type=Path)
 ap.add_argument('--label',default='PACE4 4.37 corrected'); args=ap.parse_args()
 root=args.root.resolve(); out=root/'report'; out.mkdir(exist_ok=True)
 suites={k:json.loads((root/v/'results.json').read_text()) for k,v in {'pace':'cpp_msvc','spins':'cpp_mingw','fortran':'fortran_independent'}.items()}
 if args.baseline: suites['baseline']=json.loads(args.baseline.read_text())
 for key in ['spins','fortran']: assert suites[key]['fyrst']==suites['pace']['fyrst'], 'FYRST mismatch'
 paired={(r['case']['name'],r['seed']):r for r in suites['spins']['runs']}
 build_match=all(r['residual_counts']==paired[(r['case']['name'],r['seed'])]['residual_counts'] for r in suites['pace']['runs'])
 cases=list(dict.fromkeys(r['case']['name'] for r in suites['pace']['runs']))
 results={'label':args.label,'fyrst':suites['pace']['fyrst'],'msvc_mingw_counts_identical':build_match,
 'reference':'Documented GNU port of PACE3.for; original source unchanged; not a historical VAX executable','cases':[]}
 rows=[]; sections=[]
 for name in cases:
  a,n,_=aggregate(suites['spins'],name,'initial_spin_counts'); f,nf,_=aggregate(suites['fortran'],name,'initial_spin_counts')
  assert sum(a.values())==n and sum(f.values())==nf
  keys=sorted(a.keys()|f.keys(),key=int); spin={}
  if len(keys)>1:
   weights=[1/3]*3 if name=='ni60_equal_j0_10_20' else [(2*int(k)+1)/441 for k in keys]
   for tag,c,nn in [('pace',a,n),('fortran',f,nf)]:
    stat,p=chisquare([c[k] for k in keys],[nn*w for w in weights]); spin[tag]={'chi_square':float(stat),'df':len(keys)-1,'p':float(p)}
   stat,p,df,_=chi2_contingency([[a[k] for k in keys],[f[k] for k in keys]],correction=False)
   spin['comparison']={'chi_square':float(stat),'df':int(df),'p':float(p)}
  ca,n,sigma=aggregate(suites['pace'],name,'residual_counts'); cf,nf,sigmaf=aggregate(suites['fortran'],name,'residual_counts')
  keys=sorted(ca.keys()|cf.keys(),key=lambda k:-(ca[k]+cf[k]))
  common=[k for k in keys if min(n,nf)*(ca[k]+cf[k])/(n+nf)>=5]
  rare=[k for k in keys if k not in common]
  left=[ca[k] for k in common]; right=[cf[k] for k in common]
  if rare and min(n,nf)*sum(ca[k]+cf[k] for k in rare)/(n+nf)>=5:
   left.append(sum(ca[k] for k in rare)); right.append(sum(cf[k] for k in rare)); excluded=[]
  else: excluded=rare
  stat,p,df,_=chi2_contingency([left,right],correction=False)
  result={'name':name,'events_per_program':n,'sigma_input_mb':sigma,'spin_tests':spin,
   'residue_test':{'chi_square':float(stat),'df':int(df),'p':float(p),'excluded_rare_residues':excluded},'residues':[]}
  body=[]
  for k in keys:
   lo,hi=wilson(ca[k],n); lf,hf=wilson(cf[k],nf)
   row={'case':name,'residue':nuclide(k),'pace_events':ca[k],'fortran_events':cf[k],
    'pace_mb':sigma*ca[k]/n,'pace_95_low':sigma*lo,'pace_95_high':sigma*hi,
    'fortran_mb':sigmaf*cf[k]/nf,'fortran_95_low':sigmaf*lf,'fortran_95_high':sigmaf*hf}
   if 'baseline' in suites:
    cb,nb,sb=aggregate(suites['baseline'],name,'residual_counts'); row['baseline_mb']=sb*cb[k]/nb
   result['residues'].append(row); rows.append(row)
   body.append(f"<tr><td>{row['residue']}</td><td>{ca[k]}</td><td>{cf[k]}</td><td>{row['pace_mb']:.4f} [{sigma*lo:.4f}, {sigma*hi:.4f}]</td><td>{row['fortran_mb']:.4f} [{sigmaf*lf:.4f}, {sigmaf*hf:.4f}]</td></tr>")
  results['cases'].append(result)
  spintext=f"Initial-spin comparison p = {ptext(spin['comparison']['p'])}." if spin else f'All {n:,} initial spins are J=10 in each program.'
  sections.append(f"<h2>{LABELS[name]}</h2><p>{html.escape(spintext)} Residue comparison p = {html.escape(ptext(p))}; chi-square {stat:.3f}, df {df}. {n:,} events per program; normalization {sigma:g} mb.</p><table><tr><th>Residue</th><th>PACE4 count</th><th>FORTRAN count</th><th>PACE4 mb and 95% interval</th><th>FORTRAN mb and 95% interval</th></tr>{''.join(body)}</table>")
 (out/'summary.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
 with (out/'cross_sections.csv').open('w',newline='') as handle:
  writer=csv.DictWriter(handle,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
 fig,axes=plt.subplots(1,3,figsize=(15,4.7),constrained_layout=True)
 for i,name in enumerate(cases[1:]):
  a,n,_=aggregate(suites['spins'],name,'initial_spin_counts'); f,nf,_=aggregate(suites['fortran'],name,'initial_spin_counts')
  js=[0,10,20] if i==0 else list(range(21)); x=np.arange(len(js))
  if i==0:
   axes[i].bar(x-.18,[100*a[str(j)]/n for j in js],.36,label=args.label,color='#2166ac')
   axes[i].bar(x+.18,[100*f[str(j)]/nf for j in js],.36,label='FORTRAN port',color='#d6604d')
   axes[i].axhline(100/3,ls='--',color='#333333',lw=1); axes[i].set_xticks(x,js)
  else:
   axes[i].plot(js,[100*(2*j+1)/441 for j in js],color='#333333',lw=1,label='Expected 2J+1')
   axes[i].scatter(js,[100*a[str(j)]/n for j in js],s=22,label=args.label,color='#2166ac')
   axes[i].scatter(js,[100*f[str(j)]/nf for j in js],s=22,marker='x',label='FORTRAN port',color='#d6604d')
  axes[i].set(xlabel='Initial spin J',ylabel='Events (%)',title=LABELS[name]); axes[i].legend(fontsize=7)
 keys=['28,31','27,32','27,31','26,30']; x=np.arange(4)
 series=[('pace',args.label,'#2166ac'),('fortran','FORTRAN port','#d6604d')]
 if 'baseline' in suites: series.insert(0,('baseline','Before corrections','#999999'))
 width=.8/len(series)
 for i,(key,label,color) in enumerate(series):
  c,n,sigma=aggregate(suites[key],cases[0],'residual_counts'); vals=np.array([sigma*c[k]/n for k in keys]); bounds=np.array([wilson(c[k],n) for k in keys])*sigma
  axes[2].bar(x+(i-(len(series)-1)/2)*width,vals,width,label=label,color=color,yerr=np.vstack([np.maximum(0,vals-bounds[:,0]),np.maximum(0,bounds[:,1]-vals)]),capsize=2)
 axes[2].set(xticks=x,xticklabels=[nuclide(k) for k in keys],ylabel='Cross section (mb)',title='Fixed J=10 residual cross sections'); axes[2].legend(fontsize=7)
 for ax in axes: ax.spines[['top','right']].set_visible(False)
 fig.suptitle(f"60Ni at 30.125 MeV | FYRST={results['fyrst']:g} | 25,000 events per case/program",fontsize=12)
 fig.savefig(out/'comparison.png',dpi=180); plt.close(fig)
 failures=[LABELS[c['name']] for c in results['cases'] if c['residue_test']['p']<.05/len(cases)]
 conclusion=('Residual distributions still disagree beyond sampling uncertainty in: '+', '.join(failures)+'.') if failures else 'No residual comparison rejects agreement at the familywise 5% threshold in these cases. This is not general physics validation.'
 report=f'''<!doctype html><meta charset="utf-8"><title>Corrected PACE4 physics comparison</title><style>body{{font:16px/1.5 system-ui;max-width:1200px;margin:32px auto;padding:0 24px}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{border:1px solid #ddd;padding:7px;text-align:right}}td:first-child{{text-align:left}}img{{width:100%}}</style>
 <h1>{html.escape(args.label)} versus the portable FORTRAN reference</h1><p>26 October 2026 | ASTRA physics corrections | FYRST={results['fyrst']:g}</p><p><b>{html.escape(conclusion)}</b></p>
 <p>MSVC and instrumented MinGW residue counts identical across all matched runs: {build_match}. Five independent 5,000-event replicas per case. Initial-spin histograms come from the instrumented MinGW executable.</p><img src="comparison.png" alt="Initial-spin distributions and fixed-spin residue comparisons">
 <p>Ni-60 at excitation 30.125 MeV, zero recoil, BARFAC=-1, ARATIO=1, FACLA=10, NOSHL=1, IDIST=0. Cross sections use supplied compound-nucleus normalization, not calculated entrance-channel fusion cross sections. Intervals are individual 95% Wilson intervals and exclude systematic/model uncertainty.</p>
 <p>FORTRAN uses documented GNU portability adaptations: static local storage, portable masked Park-Miller RNG with REAL*4 output, retained 30 warm-up draws, and a guarded zero-weight boundary access. Original source is unchanged. Independent reference seeds are offset by 1,000,003. See the reference manifest for compiler flags and source hashes. This is not a pristine historical VAX execution.</p>
 <p>Residual chi-square tests pool rare categories where expected counts permit; any insufficient pooled remainder is excluded and listed in summary.json. Conclusions use 0.05 divided by three cases. Full unpooled counts and intervals remain in the CSV.</p>
 {''.join(sections)}<h2>Evidence</h2><p><a href="summary.json">Summary JSON</a> | <a href="cross_sections.csv">Cross sections CSV</a> | <a href="../cpp_msvc/results.json">MSVC runs</a> | <a href="../cpp_mingw/results.json">MinGW runs</a> | <a href="../fortran_independent/results.json">FORTRAN runs and provenance</a></p>'''
 (out/'report.html').write_text(report,encoding='utf-8')
 print(json.dumps({k:v for k,v in results.items() if k!='cases'}))
 for case in results['cases']: print(case['name'],case['residue_test'])
if __name__=='__main__': main()
