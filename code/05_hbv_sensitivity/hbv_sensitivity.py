"""Idealized rainfall-only HBV soil, groundwater response and triangular routing.

Independent, uncalibrated implementation of the equations described by Seibert
(2005), HBV-light User's Manual pp. 4-5. This is not an official HBV-light run.
All forcings are rain, PET is zero; snow processes are inactive. Daily flux
ordering and routing discretisation are specified below and in methods.txt.
FAI/FCI are fitted in the separate fit_indices.py postprocessing script.
No threshold-based event-boundary selection is used for these idealized runs.
"""
from pathlib import Path
import csv,json,hashlib
import numpy as np

OUT=Path(__file__).resolve().parent / "outputs" / "analysis"
N=730
START=5
PARAMS=dict(FC_mm=200.,BETA=2.,LP=0.7,PERC_mm_day=1.,UZL_mm=10.,
            K0_day_inv=0.3,K1_day_inv=0.1,K2_day_inv=0.02,PET_mm_day=0.,
            initial_SUZ_mm=0.,initial_SLZ_mm=0.,time_step_days=1.)
BASE=dict(rainfall_total_mm=120.,rainfall_duration_days=6,SM0_fraction=0.5,MAXBAS_days=3)
TRIBUTARY_LAGS=(0,2,3)
SOURCE_URL='https://www.geo.uzh.ch/dam/jcr:c8afa73c-ac90-478e-a8c7-929eed7b1b62/hbv_manual_2005.pdf'

def csv_write(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as h:
        writer=csv.DictWriter(h,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def triangle_weights(maxbas):
    """Integrate a unit triangular travel-time density over one-day lag bins.

    Bin j=[j,j+1] receives weight CDF(j+1)-CDF(j), and is applied to
    discrete lag j. For MAXBAS=1 the kernel is [1] (no daily-bin lag).
    This is a daily discrete routing approximation, not continuous-time exactness.
    """
    edges=np.arange(int(np.ceil(maxbas))+1,dtype=float)
    x=np.clip(edges,0.,maxbas)
    cdf=np.where(x<=maxbas/2,2*(x/maxbas)**2,1-2*((maxbas-x)/maxbas)**2)
    weights=np.diff(cdf)
    assert np.all(weights>=0) and abs(weights.sum()-1)<1e-13
    return weights

def soil_rainfall_update(sm,p):
    """Integrate dSM/dP=1-(SM/FC)^BETA over rain depth, for BETA=2.

    The exact flow map avoids artificial dependence of total recharge on the
    subdivision of a fixed storm depth into daily cells. This implementation
    choice is documented; it is not claimed to reproduce official software.
    """
    assert PARAMS['BETA']==2. and p>=0
    if p==0:return sm,0.
    fc=PARAMS['FC_mm'];r=sm/fc;u=float(np.tanh(p/fc))
    new=fc*(r+u)/(1+r*u)
    recharge=p-(new-sm)
    assert 0<=new<=fc+1e-12 and recharge>=-1e-12
    return new,recharge

def hbv(config,start=START):
    p=np.zeros(N)
    duration=int(config['rainfall_duration_days'])
    j=np.arange(duration)
    rain_weights=np.minimum(j+1,duration-j).astype(float)
    p[start:start+duration]=config['rainfall_total_mm']*rain_weights/rain_weights.sum()
    sm=config['SM0_fraction']*PARAMS['FC_mm'];uz=PARAMS['initial_SUZ_mm'];lz=PARAMS['initial_SLZ_mm']
    initial=sm+uz+lz
    records=[];raw=np.zeros(N)
    for day in range(N):
        old_sm,old_uz,old_lz=sm,uz,lz
        sm,recharge=soil_rainfall_update(old_sm,p[day])
        overflow=0. # Integrated storage approaches FC without exceeding it.
        aet=min(sm,PARAMS['PET_mm_day']*min(sm/(PARAMS['FC_mm']*PARAMS['LP']),1.))
        sm-=aet
        uz=old_uz+recharge
        perc=min(PARAMS['PERC_mm_day'],uz)
        uz-=perc;lz=old_lz+perc
        # Both upper-store fluxes use the same state, after recharge/percolation.
        q0=PARAMS['K0_day_inv']*max(uz-PARAMS['UZL_mm'],0.)
        q1=PARAMS['K1_day_inv']*uz
        q2=PARAMS['K2_day_inv']*lz
        uz-=q0+q1;lz-=q2
        raw[day]=q0+q1+q2
        residual=p[day]-aet-raw[day]-(sm+uz+lz-old_sm-old_uz-old_lz)
        records.append(dict(day_index=day,P_mm_day=float(p[day]),recharge_mm_day=float(recharge),
                            overflow_mm_day=float(overflow),AET_mm_day=float(aet),percolation_mm_day=float(perc),
                            Q0_mm_day=float(q0),Q1_mm_day=float(q1),Q2_mm_day=float(q2),
                            Q_generated_mm_day=float(raw[day]),SM_end_mm=float(sm),SUZ_end_mm=float(uz),SLZ_end_mm=float(lz),
                            daily_balance_error_mm=float(residual)))
        assert -1e-12<=sm<=PARAMS['FC_mm']+1e-12 and uz>=-1e-12 and lz>=-1e-12
    kernel=triangle_weights(config['MAXBAS_days'])
    routed=np.convolve(raw,kernel,mode='full')
    final=sm+uz+lz
    summary=dict(initial_total_storage_mm=float(initial),final_total_storage_mm=float(final),
                 final_soil_storage_mm=float(sm),final_upper_storage_mm=float(uz),final_lower_storage_mm=float(lz),
                 rainfall_total_check_mm=float(p.sum()),generated_runoff_volume_mm=float(raw.sum()),
                 routed_runoff_volume_mm=float(routed.sum()),routing_volume_error_mm=float(routed.sum()-raw.sum()),
                 whole_balance_error_mm=float(initial+p.sum()-routed.sum()-sum(r['AET_mm_day'] for r in records)-final),
                 maximum_daily_balance_error_mm=float(max(abs(r['daily_balance_error_mm']) for r in records)),
                 minimum_Q_mm_day=float(routed.min()),minimum_state_mm=float(min(min(r['SM_end_mm'],r['SUZ_end_mm'],r['SLZ_end_mm']) for r in records)),
                 routing_weights=kernel.tolist(),mean_discrete_routing_lag_days=float(np.dot(kernel,np.arange(len(kernel)))))
    assert abs(summary['whole_balance_error_mm'])<1e-9
    assert abs(summary['routing_volume_error_mm'])<1e-9
    assert summary['maximum_daily_balance_error_mm']<1e-9
    return routed,p,records,summary

def scenarios():
    result=[]
    def add(panel,code,label,changes=None,baseline=False,lag=None):
        config={**BASE,**(changes or {})}
        result.append(dict(panel=panel,scenario_id=code,label=label,baseline=baseline,tributary_lag_days=lag,**config))
    add('a','rain60','P = 60 mm, D = 6 days',{'rainfall_total_mm':60.})
    add('a','baseline_a','P = 120 mm, D = 6 days',baseline=True)
    add('a','rain180','P = 180 mm, D = 6 days',{'rainfall_total_mm':180.})
    add('a','duration3','P = 120 mm, D = 3 days',{'rainfall_duration_days':3})
    add('a','duration12','P = 120 mm, D = 12 days',{'rainfall_duration_days':12})
    for f in [.2,.5,.8]:add('b',f'soil{int(f*100)}',f'SM0 / FC = {f:.1f}',{'SM0_fraction':f},baseline=f==.5)
    for m in [1,3,7]:add('c',f'routing{m}',f'MAXBAS = {m} days',{'MAXBAS_days':m},baseline=m==3)
    for lag in TRIBUTARY_LAGS:add('d',f'lag{lag}',f'Tributary lag = {lag} days',baseline=lag==0,lag=lag)
    return result

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    all_series=[];all_states=[];metrics=[];details={};curves={};kernels=[]
    for s in scenarios():
        q,p,states,balance=hbv(s)
        individual=[(1,q,p,states,balance)]
        if s['panel']=='d':
            q2,p2,states2,bal2=hbv(s,START+s['tributary_lag_days'])
            q=.5*(q+q2);p=.5*(p+p2)
            individual.append((2,q2,p2,states2,bal2))
            balance={k: .5*(balance[k]+bal2[k]) if isinstance(balance[k],(int,float)) else balance[k] for k in balance}
        curves[s['scenario_id']]=q
        tp=int(np.argmax(q));t=np.arange(len(q))+.5-START
        row={**s,'peak_Q_mm_day':float(q[tp]),'peak_time_from_first_rainfall_days':float(t[tp]),
             'runoff_volume_mm':float(q.sum()),'final_upper_storage_mm':balance['final_upper_storage_mm'],
             'final_lower_storage_mm':balance['final_lower_storage_mm'],
             'outflow_plus_remaining_groundwater_mm':float(q.sum()+balance['final_upper_storage_mm']+balance['final_lower_storage_mm']),
             'whole_balance_error_mm':balance['whole_balance_error_mm'],
             'routing_volume_error_mm':balance['routing_volume_error_mm'],
             'maximum_daily_balance_error_mm':balance['maximum_daily_balance_error_mm'],
             'volume_fraction_displayed_first_20_days':float(q[(t>=0)&(t<=20)].sum()/q.sum())}
        metrics.append(row);details[s['scenario_id']]=dict(settings=s,balance=balance,metrics=row)
        for j in range(len(q)):
            all_series.append(dict(panel=s['panel'],scenario_id=s['scenario_id'],label=s['label'],time_days=float(t[j]),
                                   Q_mm_day=float(q[j]),P_mm_day=float(p[j]) if j<len(p) else 0.,baseline=s['baseline']))
        for tributary,_,_,state_rows,_ in individual:
            for record in state_rows:all_states.append(dict(scenario_id=s['scenario_id'],tributary=tributary,**record))
    # OFAT baseline is identical across all four panels; the synchronised pair
    # is area weighted, not summed to double the catchment-normalised discharge.
    base=curves['baseline_a']
    for key in ['soil50','routing3','lag0']:assert np.array_equal(base,curves[key])
    # Routing changes no generated runoff volume: full convolution retains tails.
    routed_volumes=[next(m['runoff_volume_mm'] for m in metrics if m['scenario_id']==f'routing{n}') for n in [1,3,7]]
    assert max(routed_volumes)-min(routed_volumes)<1e-9
    lag_volumes=[next(m['runoff_volume_mm'] for m in metrics if m['scenario_id']==f'lag{n}') for n in TRIBUTARY_LAGS]
    assert max(lag_volumes)-min(lag_volumes)<1e-5 # finite-record residual lower-store tails
    lag_including_storage=[next(m['outflow_plus_remaining_groundwater_mm'] for m in metrics if m['scenario_id']==f'lag{n}') for n in TRIBUTARY_LAGS]
    assert max(lag_including_storage)-min(lag_including_storage)<1e-9
    # Single-peak check uses all original daily values, before any interpolation.
    unimodality={}
    for lag in TRIBUTARY_LAGS:
        q=curves[f'lag{lag}']; peak=int(np.argmax(q)); tol=1e-10
        pre=np.diff(q[:peak+1]); post=np.diff(q[peak:])
        pre_dips=int(np.sum(pre < -tol)); post_rises=int(np.sum(post > tol))
        assert pre_dips==0 and post_rises==0, f'Non-unimodal lag={lag}'
        unimodality[f'lag{lag}']=dict(peak_Q_mm_day=float(q[peak]),
            peak_time_from_first_rainfall_days=float(peak+.5-START),
            pre_peak_decrease_count=pre_dips,post_peak_increase_count=post_rises,
            tolerance_mm_day=tol,full_daily_series_checked=True,single_peaked=True)
    duration_soil_final=[details[key]['balance']['final_soil_storage_mm'] for key in ['duration3','baseline_a','duration12']]
    assert max(duration_soil_final)-min(duration_soil_final)<1e-10
    duration_releasable=[next(m['outflow_plus_remaining_groundwater_mm'] for m in metrics if m['scenario_id']==key) for key in ['duration3','baseline_a','duration12']]
    assert max(duration_releasable)-min(duration_releasable)<1e-9
    # Independent fourth-order Runge-Kutta integration in rain-depth space.
    soil_exact,_=soil_rainfall_update(100.,120.)
    soil_rk=100.;h=0.01
    f=lambda s:1-(s/PARAMS['FC_mm'])**2
    for _ in range(12000):
        k1=f(soil_rk);k2=f(soil_rk+.5*h*k1);k3=f(soil_rk+.5*h*k2);k4=f(soil_rk+h*k3)
        soil_rk+=h*(k1+2*k2+2*k3+k4)/6
    assert abs(soil_exact-soil_rk)<1e-9
    expected_kernels={1:np.array([1.]),3:np.array([2.,5.,2.])/9,7:np.array([2.,6.,10.,13.,10.,6.,2.])/49}
    # An independently computed continuous triangle integral checks the CDF bins.
    for maxbas in [1,3,7]:
        w=triangle_weights(maxbas)
        assert np.allclose(w,expected_kernels[maxbas],rtol=0,atol=1e-14)
        x=(np.arange(maxbas*10000)+.5)/10000
        density=np.maximum(0.,(2/maxbas)*(1-np.abs((x-maxbas/2)/(maxbas/2))))
        approx=np.bincount(np.floor(x).astype(int),weights=density/10000,minlength=maxbas)
        assert np.max(np.abs(w-approx))<1e-7
        for lag,value in enumerate(w):kernels.append(dict(MAXBAS_days=maxbas,lag_days=lag,weight=float(value)))
    checks=dict(model='Independent rainfall-only HBV soil-groundwater-triangular-routing implementation',
                official_software_run=False,calibrated=False,source=SOURCE_URL,
                parameters=PARAMS,baseline=BASE,simulation_days=N,first_rain_day_index=START,
                plotted_hydrographs=14,unique_single_catchment_settings=9,
                four_panel_baselines_identical=True,full_convolution=True,
                validation_scope='Model simulation; index fitting is separately checked in index_fit_validation.json',
                FAI_FCI_fitting_in_this_model_script=False,
                tributary_lags_days=list(TRIBUTARY_LAGS),panel_d_single_peak_checks=unimodality,
                routing_volume_spread_mm=float(max(routed_volumes)-min(routed_volumes)),
                tributary_lag_volume_spread_mm=float(max(lag_volumes)-min(lag_volumes)),
                tributary_lag_outflow_plus_residual_storage_spread_mm=float(max(lag_including_storage)-min(lag_including_storage)),
                duration_soil_final_spread_mm=float(max(duration_soil_final)-min(duration_soil_final)),
                duration_outflow_plus_groundwater_spread_mm=float(max(duration_releasable)-min(duration_releasable)),
                analytical_soil_vs_independent_RK4_error_mm=float(abs(soil_exact-soil_rk)),
                max_daily_balance_error_mm=float(max(m['maximum_daily_balance_error_mm'] for m in metrics)),
                max_whole_balance_error_mm=float(max(abs(m['whole_balance_error_mm']) for m in metrics)),
                max_remaining_groundwater_mm=float(max(m['final_upper_storage_mm']+m['final_lower_storage_mm'] for m in metrics)),
                minimum_Q_mm_day=float(min(r['Q_mm_day'] for r in all_series)),
                scenarios=details,all_checks_passed=True)
    csv_write('hydrographs.csv',all_series);csv_write('states_and_fluxes.csv',all_states)
    csv_write('scenario_metrics.csv',metrics);csv_write('routing_weights.csv',kernels)
    (OUT/'validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    methods='''HBV sensitivity experiments

We implemented the soil, two-store groundwater-response and triangular-routing routines described in the HBV-light manual (Seibert, 2005, pp. 4–5). This is an independent, uncalibrated rainfall-only implementation, not a run of the official HBV-light software. It is an idealized process sensitivity experiment and is not calibrated to a particular river or flood. All precipitation is rain; snow processes are inactive and potential evapotranspiration is zero. No simulated discharge is peak-normalised or amplitude-rescaled.

Time step and initial conditions. Each simulation has 730 one-day cells, with the first rainfall starting at cell index 5. Soil capacity FC is 200 mm and BETA is 2. The default initial soil storage is 100 mm (SM0/FC=0.5); initial upper and lower groundwater stores are zero. LP is 0.7; with PET=0 this parameter has no effect. State values are imposed for these controlled event experiments, without a calibration or warm-up period.

Within-day sequence. The HBV soil-recharge relation is integrated over rainfall depth: dSM/dP=1-(SM/FC)^BETA. For BETA=2, let r=SMstart/FC and u=tanh(P/FC); the stable exact update is SMrain=FC*(r+u)/(1+r*u). Recharge is P-(SMrain-SMstart). This prevents artificial changes in total recharge when a fixed storm depth is split among different daily bins. Soil storage approaches FC without exceeding it; an initially saturated store sends all rainfall to recharge. This analytical rain-depth integration is our numerical implementation choice, not a claim about the official HBV-light program. An independent Runge-Kutta integration and equal-total-rainfall comparisons verify the update. Actual evaporation is min(available soil water, PET*min(SM/(LP*FC),1)). Recharge enters the upper groundwater store. Percolation min(PERC,SUZ), with PERC=1 mm/day, is then transferred to the lower store. From the same post-percolation upper-store state, Q0=K0*max(SUZ-UZL,0) and Q1=K1*SUZ are computed simultaneously; both are subtracted together. Q2=K2*SLZ is computed after receipt of percolation. K0=0.3/day, K1=0.1/day, K2=0.02/day and UZL=10 mm. Generated runoff is Q0+Q1+Q2. Storage and fluxes are expressed as depths over catchment area. These are numerical settings of this experiment, not recommended or calibrated parameter values.

Routing. Generated daily runoff is convolved with non-negative, unit-sum triangular weights. The triangle is supported on [0,MAXBAS], with its apex at MAXBAS/2. Weight j is the integral of this density from j to j+1 and is applied to discrete lag j, using the full convolution. Thus MAXBAS is the triangular support width, not mean lag, and MAXBAS=1 day has daily kernel [1]. No continuous-time exactness is claimed for this daily-bin routing approximation. Full convolution retains all routed mass, including the end tail; residual groundwater at day 730 is explicitly included in mass-balance checks.

Controlled comparisons. The common baseline is triangular rainfall with total 120 mm over 6 days, SM0/FC=0.5 and MAXBAS=3 days. Discrete storm weights are proportional to min(j+1,D-j) for j=0,...,D-1 and sum to the prescribed total. Panel a varies rainfall total (60,120,180 mm at D=6 days) or duration (3,6,12 days at P=120 mm); the baseline is shown once. Panel b varies SM0/FC over 0.2,0.5,0.8. Panel c varies MAXBAS over 1,3,7 days. In panel d, two identical equal-area catchments receive the baseline storm with onset separation 0,2,3 days. These lags were selected to retain single-peaked outlet responses: the full original daily series are checked to be non-decreasing up to their maximum and non-increasing thereafter, with tolerance 1e-10 mm/day. Each catchment is simulated independently and their routed discharge depths are averaged with weights 0.5 and 0.5; they are not added as unweighted depths. The time origin is the first storm onset. Equal storm totals are retained in each tributary, so catchment-average rainfall totals are constant. Any tiny difference in finite-record outlet volume is reported with the corresponding remaining groundwater storage.

Display and validation. Figures use a common time and discharge scale, with time measured from the first rainfall onset rather than each curve's peak. Shape-preserving interpolation is a display guide through the daily means at bin centres; calculations use the original daily series. No threshold-based event extraction is performed. All prescribed cases are retained. Per-day and whole-record water balances, non-negative states, normalized routing weights, conservative routing and identical common baselines are checked. The full source data contain all simulation cells, even though only days 0-20 after rainfall onset are plotted.

FAI and FCI. The separate fit_indices.py script uses the same integer durations T=1,...,20 days for every hydrograph. For each T, W(T) is the maximum sum of T consecutive original daily discharge depths, searched over the entire available simulation including complete routing tails. With Qref=1 mm/day and Tref=1 day, unweighted ordinary least squares fits log10[W(T)/(Qref*Tref)] = FAI + (1-FCI)*log10(T/Tref). Thus FAI is the fitted intercept and FCI is one minus the fitted volume exponent. The intercept is referenced to the stated units, rather than assumed to equal log10 of the actual maximum daily discharge. These indices describe the specified 1-20-day fitting range and are not full-event estimates. The calendar interval displayed in the figure does not restrict the maximum-window search. No smoothed samples, new event-end threshold or fit-quality exclusion is used. fitted_indices.csv contains the indices and R-squared for all 14 plotted curves; maximum_window_volumes.csv contains the 20 maximizing windows and fitting residuals per curve. index_fit_validation.json checks the calculation separately from the model water-balance checks.

Source: Seibert, J. (2005). HBV light version 2 User's Manual, pp. 4–5.
'''+SOURCE_URL+'\n'
    (OUT/'methods.txt').write_text(methods,encoding='utf-8')
    readme='''Reproduce: run run_all.py from the HBV code directory (Python, NumPy, ReportLab, pypdfium2 and Pillow). See package/docs/HBV_README.md.
The CSV files contain all daily values, states/fluxes, scenario settings/metrics and routing weights.
validation.json documents full-record mass conservation and finite-record residual storage.
index_fit_validation.json documents the separate 1-20-day FAI/FCI fitting; fitted_indices.csv contains all indices and R-squared values.
maximum_window_volumes.csv records the maximizing windows over the full source series, independently of the displayed 0-20-day interval.
methods.txt specifies equations, flux ordering, parameters, assumptions and sources.
No observations, external datasets, random samples or calibration were used.
The resulting hydrographs are simulated experiments, not observed flood events.
'''
    (OUT/'README.txt').write_text(readme,encoding='utf-8')
    print(json.dumps({k:v for k,v in checks.items() if k not in ['scenarios','parameters']},indent=2))

if __name__=='__main__':main()
