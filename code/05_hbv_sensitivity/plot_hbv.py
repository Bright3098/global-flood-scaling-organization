"""Plot one four-panel HBV sensitivity figure from supplied daily model output.

No model runs, rescaling, normalization, or peak alignment are performed here.
PCHIP is a shape-preserving display guide; the source CSV remains unchanged.
Run: python plot_hbv.py, after the model and fitting scripts.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pypdfium2 as pdfium
from reportlab.graphics import renderPDF, renderSVG
from reportlab.graphics.shapes import Drawing, Group, Line, PolyLine, Rect, String
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'outputs' / 'analysis' / 'hydrographs.csv'
INDICES = ROOT / 'outputs' / 'analysis' / 'fitted_indices.csv'
OUT = ROOT / 'outputs' / 'figures'
X_MIN, X_MAX = 0.0, 20.0
TOL = 1e-10
TITLES = {
    'a': 'Precipitation amount and duration',
    'b': 'Antecedent soil moisture',
    'c': 'Routing duration',
    'd': 'Tributary asynchrony',
}
BASELINE_COLOR = '#3F3F3F'
BLUE, ORANGE, GREEN, PURPLE = '#0072B2', '#D55E00', '#009E73', '#CC79A7'


def pchip_tangents(x, y):
    h = np.diff(x)
    delta = np.diff(y) / h
    m = np.zeros_like(y)
    for i in range(1, len(y) - 1):
        if delta[i - 1] * delta[i] > 0:
            w1, w2 = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
            m[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i])

    def endpoint(h0, h1, d0, d1):
        value = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
        if np.sign(value) != np.sign(d0):
            return 0.0
        if np.sign(d0) != np.sign(d1) and abs(value) > 3 * abs(d0):
            return 3 * d0
        return value

    m[0] = endpoint(h[0], h[1], delta[0], delta[1])
    m[-1] = endpoint(h[-1], h[-2], delta[-1], delta[-2])
    return m


def coefficients(x, y, m):
    h = np.diff(x)
    return (2*y[:-1] - 2*y[1:] + h*(m[:-1] + m[1:]),
            -3*y[:-1] + 3*y[1:] - h*(2*m[:-1] + m[1:]),
            h*m[:-1], y[:-1])


def interpolate(x, y, m, query):
    if np.min(query) < x[0] or np.max(query) > x[-1]:
        raise ValueError('Extrapolation is not allowed')
    i = np.clip(np.searchsorted(x, query, side='right') - 1, 0, len(x)-2)
    u = (query - x[i]) / (x[i+1] - x[i])
    a, b, c, d = coefficients(x, y, m)
    return ((a[i]*u + b[i])*u + c[i])*u + d[i]


def interpolation_audit(x, y, m):
    a, b, c, d = coefficients(x, y, m)
    overshoot = 0.0
    monotonicity = 0.0
    for i in range(len(x)-1):
        candidates = [0.0, 1.0]
        poly = np.trim_zeros(np.asarray([3*a[i], 2*b[i], c[i]]), 'f')
        if len(poly) > 1:
            for root in np.roots(poly):
                if abs(np.imag(root)) < 1e-12 and 0 < np.real(root) < 1:
                    candidates.append(float(np.real(root)))
        values = [((a[i]*u + b[i])*u + c[i])*u + d[i] for u in candidates]
        overshoot = max(overshoot, max(values)-max(y[i], y[i+1]),
                        min(y[i], y[i+1])-min(values))
        derivative_points = [0.0, 1.0]
        if abs(a[i]) > 1e-30:
            vertex = -b[i] / (3*a[i])
            if 0 < vertex < 1:
                derivative_points.append(float(vertex))
        direction = np.sign(y[i+1] - y[i])
        monotonicity = max(monotonicity, *[
            -direction*(3*a[i]*u*u + 2*b[i]*u + c[i]) for u in derivative_points])
    knot_error = float(np.max(np.abs(interpolate(x, y, m, x)-y)))
    derivative_jump = float(np.max(np.abs(
        ((3*a + 2*b + c)/np.diff(x))[:-1] - (c/np.diff(x))[1:])))
    assert max(overshoot, monotonicity, knot_error, derivative_jump) < TOL
    return dict(maximum_interval_overshoot_mm_day=float(overshoot),
                maximum_monotonicity_violation=float(monotonicity),
                maximum_knot_error_mm_day=knot_error,
                maximum_derivative_jump=derivative_jump,
                all_checks_passed=True)


def load_data():
    groups = defaultdict(list)
    with SOURCE.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        required = {'panel','scenario_id','label','time_days','Q_mm_day','P_mm_day','baseline'}
        assert required.issubset(reader.fieldnames), reader.fieldnames
        for row in reader:
            groups[row['panel'].strip().lower(), row['scenario_id']].append(row)
    panels = defaultdict(list)
    for (panel, sid), rows in groups.items():
        assert panel in TITLES, panel
        rows.sort(key=lambda r: float(r['time_days']))
        x = np.asarray([float(r['time_days']) for r in rows])
        y = np.asarray([float(r['Q_mm_day']) for r in rows])
        rain = np.asarray([float(r['P_mm_day']) for r in rows])
        assert len(x) >= 3 and np.all(np.diff(x) > 0)
        assert np.all(np.isfinite(x)) and np.all(np.isfinite(y))
        assert np.min(y) >= 0
        assert len({r['label'] for r in rows}) == 1
        baseline_flags = {r['baseline'].strip().lower() for r in rows}
        assert len(baseline_flags) == 1
        baseline = baseline_flags.pop() in ('1','true','yes')
        tangent = pchip_tangents(x, y)
        checks = interpolation_audit(x, y, tangent)
        lo, hi = max(X_MIN, x[0]), min(X_MAX, x[-1])
        assert hi > lo
        query = np.unique(np.concatenate((
            np.linspace(lo, hi, int(math.ceil((hi-lo)*80))+1), x[(x>=lo)&(x<=hi)])))
        smooth_y = interpolate(x, y, tangent, query)
        assert np.min(smooth_y) >= -TOL
        assert np.max(smooth_y) <= np.max(y) + TOL
        panels[panel].append(dict(
            scenario_id=sid, label=rows[0]['label'], baseline=baseline,
            x=x, y=y, query=query, smooth_y=smooth_y,
            rain_positive_day_count=int(np.count_nonzero(rain > 1e-12)),
            checks=checks))
    for panel, expected in zip('abcd', [5,3,3,3]):
        assert len(panels[panel]) == expected, (panel,len(panels[panel]))
        assert sum(r['baseline'] for r in panels[panel]) == 1
    return panels


def style_curves(panels):
    for panel, rows in panels.items():
        baseline = next(r for r in rows if r['baseline'])
        baseline_duration = baseline['rain_positive_day_count']
        regular = [r for r in rows if not r['baseline'] and
                   not (panel=='a' and r['rain_positive_day_count'] != baseline_duration)]
        regular_ids = {r['scenario_id'] for r in regular}
        additional = [r for r in rows if not r['baseline'] and r['scenario_id'] not in regular_ids]
        for row, color in zip(regular, [BLUE, ORANGE]):
            row.update(color=color, dash=None, linewidth=1.2)
        for row, color in zip(additional, [GREEN, PURPLE]):
            row.update(color=color, dash=[4,2], linewidth=1.2)
        baseline.update(color=BASELINE_COLOR, dash=None, linewidth=1.5)
        for row in rows:
            assert 'color' in row
            display_label = row['label'].replace('SM0', 'SM₀').replace('= 1 days', '= 1 day')
            display_label = display_label.replace(' (baseline)', '')
            if panel == 'a':
                match = re.fullmatch(r'P\s*=\s*([\d.]+)\s*mm,\s*D\s*=\s*([\d.]+)\s*days?', display_label)
                assert match, display_label
                display_label = f'{match[1]} mm / {match[2]} days'
            elif panel == 'd':
                display_label = display_label.replace('Tributary lag', 'Lag')
            display_label = display_label.replace('SM₀ / FC', 'SM₀/FC')
            row['legend_label'] = display_label + (' (base)' if row['baseline'] else '')


def attach_indices(panels):
    """Read approved values only; this renderer never fits the hydrographs."""
    with INDICES.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        assert {'scenario_id','FAI','FCI'}.issubset(reader.fieldnames), reader.fieldnames
        rows = list(reader)
    by_id = {r['scenario_id']: r for r in rows}
    assert len(by_id) == len(rows), 'Duplicate scenario_id in fitted_indices.csv'
    source_ids = {r['scenario_id'] for group in panels.values() for r in group}
    assert source_ids == set(by_id), (source_ids - set(by_id), set(by_id) - source_ids)
    for group in panels.values():
        for record in group:
            values = {key: float(by_id[record['scenario_id']][key]) for key in ('FAI','FCI')}
            assert all(math.isfinite(v) for v in values.values()), record['scenario_id']
            record.update(values)
    baselines = [next(r for r in panels[p] if r['baseline']) for p in 'abcd']
    assert all(abs(r[k] - baselines[0][k]) < TOL
               for r in baselines for k in ('FAI','FCI')), 'Baseline indices differ across panels'


def y_axis(panels):
    maximum = max(float(r['y'].max()) for rows in panels.values() for r in rows)
    raw_step = maximum * 1.06 / 4.5
    power = 10 ** math.floor(math.log10(raw_step))
    step = next(v*power for v in [1,2,2.5,5,10] if v*power >= raw_step)
    top = math.ceil(maximum * 1.06 / step) * step
    return top, np.arange(0, top + step/2, step)


def render(panels):
    OUT.mkdir(parents=True, exist_ok=True)
    width, height = 180*72/25.4, 181*72/25.4
    drawing = Drawing(width,height)
    drawing.add(Rect(0,0,width,height,fillColor=colors.white,strokeColor=None))
    y_max, ticks = y_axis(panels)

    def text(x,y,value,size=8.6,anchor='start',bold=False,rotate=0):
        obj = String(0 if rotate else x,0 if rotate else y,value,
                     fontName='Helvetica-Bold' if bold else 'Helvetica',fontSize=size,
                     textAnchor=anchor,fillColor=colors.HexColor('#202020'))
        if rotate:
            group=Group(obj); group.translate(x,y); group.rotate(rotate); drawing.add(group)
        else:
            drawing.add(obj)

    def line(x1,y1,x2,y2,color='#303030',weight=.6,dash=None):
        drawing.add(Line(x1,y1,x2,y2,strokeColor=colors.HexColor(color),strokeWidth=weight,strokeDashArray=dash))

    def label_width(value,size):
        # Render the zero as a positioned subscript rather than relying on a
        # Unicode glyph that is absent from some versions of the Helvetica font.
        if '₀' not in value:
            return pdfmetrics.stringWidth(value,'Helvetica',size)
        before,after=value.split('₀',1)
        return (pdfmetrics.stringWidth(before+after,'Helvetica',size)
                +pdfmetrics.stringWidth('0','Helvetica',size*.72))

    def legend_text(x,y,value,size):
        if '₀' not in value:
            text(x,y,value,size)
            return
        before,after=value.split('₀',1)
        text(x,y,before,size)
        x+=pdfmetrics.stringWidth(before,'Helvetica',size)
        text(x,y-size*.23,'0',size*.72)
        x+=pdfmetrics.stringWidth('0','Helvetica',size*.72)
        text(x,y,after,size)

    text(12,height-17,'HBV rainfall–runoff sensitivity',12,bold=True)
    lefts=[48.0,299.0]
    panel_width,panel_height=195.0,110.0
    bottoms=[height-249,41.0]
    for index,panel in enumerate('abcd'):
        row_index,col=index//2,index%2
        x0,y0=lefts[col],bottoms[row_index]
        title_y=height-46 if row_index==0 else height-294
        header_y=height-63 if row_index==0 else height-311
        legend_y=header_y-14
        text(x0-16,title_y,panel,10.5,bold=True)
        text(x0,title_y,TITLES[panel],9.3,bold=True)
        text(x0+21,header_y,'P / D' if panel == 'a' else 'Scenario',8.2,bold=True)
        text(x0+154,header_y,'FAI',8.2,anchor='middle',bold=True)
        text(x0+183,header_y,'FCI',8.2,anchor='middle',bold=True)
        line(x0,header_y-4,x0+panel_width,header_y-4,color='#B8B8B8',weight=.4)
        for j,record in enumerate(panels[panel]):
            ly=legend_y-j*12
            line(x0,ly+2,x0+15,ly+2,record['color'],record['linewidth'],record['dash'])
            font_size=8.2
            while label_width(record['legend_label'],font_size)>117:
                font_size=round(font_size-.1,1)
            assert font_size >= 7.5, record['legend_label']
            legend_text(x0+21,ly,record['legend_label'],font_size)
            text(x0+154,ly,f"{record['FAI']:.3f}",8.2,anchor='middle')
            text(x0+183,ly,f"{record['FCI']:.3f}",8.2,anchor='middle')
        line(x0,y0,x0+panel_width,y0)
        line(x0,y0,x0,y0+panel_height)
        for tick in [0,5,10,15,20]:
            px=x0+(tick-X_MIN)/(X_MAX-X_MIN)*panel_width
            line(px,y0,px,y0-3)
            text(px,y0-13,f'{tick:g}',8.2,anchor='middle')
        for tick in ticks:
            py=y0+tick/y_max*panel_height
            line(x0-3,py,x0,py)
            text(x0-6,py-2.8,f'{tick:g}',8.2,anchor='end')
        text(x0+panel_width/2,y0-27,'Time from rainfall onset (days)',8.5,anchor='middle')
        if col==0:
            text(x0-32,y0+panel_height/2,'Discharge (mm/day)',8.7,anchor='middle',rotate=90)
        for record in sorted(panels[panel],key=lambda r:r['baseline']):
            px=x0+(record['query']-X_MIN)/(X_MAX-X_MIN)*panel_width
            py=y0+record['smooth_y']/y_max*panel_height
            pts=[float(v) for pair in zip(px,py) for v in pair]
            drawing.add(PolyLine(pts,strokeColor=colors.HexColor(record['color']),
                                strokeWidth=record['linewidth'],fillColor=None,
                                strokeDashArray=record['dash']))
    stem=OUT/'Extended_Data_Figure_5_HBV_hydrographs'
    renderPDF.drawToFile(drawing,str(stem.with_suffix('.pdf')),title='HBV rainfall–runoff sensitivity',author='')
    renderSVG.drawToFile(drawing,str(stem.with_suffix('.svg')))
    pdf=pdfium.PdfDocument(str(stem.with_suffix('.pdf')))
    page=pdf[0]
    page.render(scale=600/72).to_pil().save(stem.with_suffix('.png'),dpi=(600,600))
    page.render(scale=180/72).to_pil().save(OUT/(stem.name+'_preview.png'))
    page.close();pdf.close()
    audit=dict(
        source=SOURCE.name,source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        fitted_indices_source=INDICES.name,
        fitted_indices_source_sha256=hashlib.sha256(INDICES.read_bytes()).hexdigest(),
        indices_fitted_by_plotter=False,indices_display_decimal_places=3,
        panel_count=4,curve_counts={p:len(rows) for p,rows in panels.items()},
        shared_time_range_days=[X_MIN,X_MAX],shared_discharge_range_mm_day=[0,float(y_max)],
        shared_time_ticks_days=[0,5,10,15,20],
        shared_discharge_ticks=ticks.tolist(),time_origin='Rainfall onset; earliest tributary onset in panel d',
        peak_alignment=False,discharge_rescaling=False,normalization=False,
        interpolation='Shape-preserving cubic Hermite interpolation through supplied daily values; display only',
        extrapolation=False,negative_value_clipping=False,numerical_tolerance=TOL,
        all_interpolation_checks_passed=True,
        curves=[dict(panel=p,scenario_id=r['scenario_id'],label=r['legend_label'],baseline=r['baseline'],
                     FAI=r['FAI'],FCI=r['FCI'],
                     source_daily_count=len(r['x']),source_time_range=[float(r['x'][0]),float(r['x'][-1])],
                     display_point_count=len(r['query']),daily_peak_mm_day=float(r['y'].max()),
                     color=r['color'],line_dash=r['dash'],**r['checks'])
                for p,rows in panels.items() for r in rows])
    (OUT/'figure_data_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps({k:audit[k] for k in ['panel_count','curve_counts','shared_time_range_days',
                     'shared_discharge_range_mm_day','all_interpolation_checks_passed']},indent=2))


def main():
    panels=load_data()
    baselines=[next(r for r in panels[p] if r['baseline']) for p in 'abcd']
    assert all(np.array_equal(r['x'],baselines[0]['x']) for r in baselines)
    assert all(np.allclose(r['y'],baselines[0]['y'],rtol=0,atol=TOL) for r in baselines)
    style_curves(panels)
    attach_indices(panels)
    render(panels)


if __name__=='__main__':
    main()
