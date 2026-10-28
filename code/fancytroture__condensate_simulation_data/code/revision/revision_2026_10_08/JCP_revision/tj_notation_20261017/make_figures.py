"""Redraw only the two TJ figures with unchanged numerical inputs and styles."""
from pathlib import Path
import importlib.util
import json

W = Path(__file__).resolve().parent
P = W.parent
OUT = W / 'figures'
OUT.mkdir(exist_ok=True)
source = P / 'response_revision_20261015_v3/plot_figures.py'
spec = importlib.util.spec_from_file_location('ratio_figures', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.OUT = OUT
original_axes_setup = module.axes_setup
def tj_axes_setup(ax, title, **kwargs):
    kwargs['ylabel'] = r'$\tau^{\mathrm{TJ}}(T_{\mathrm{H}}^{\mathrm{TJ}})/\tau^{\mathrm{TJ}}(T_{\mathrm{L}}^{\mathrm{TJ}})$'
    return original_axes_setup(ax, title, **kwargs)
module.axes_setup = tj_axes_setup
module.controls()
old = json.loads((P / 'response_revision_20261015_v3/plotted_values.json').read_text())
old = [row for row in old if row['figure'] == 'R1']
assert module.VALUES == old, 'TJ plotted values or uncertainty intervals changed'
(W / 'plotted_values.json').write_text(json.dumps(module.VALUES, indent=2) + '\n')
(W / 'layout_bounds.json').write_text(json.dumps(module.LAYOUT, indent=2) + '\n')

source = P / 'density_labels_20261016/make_density.py'
code = source.read_text()
assert (source.parent / 'density_summary.npz').exists()
assert code.count('label=style.temp_label(role)') == 1
code = code.replace('label=style.temp_label(role)', "label=(r'$T_{\\mathrm{L}}^{\\mathrm{TJ}}$' if role == 'Tl' else r'$T_{\\mathrm{H}}^{\\mathrm{TJ}}$')")
code = code.replace("W/'figures'/", "TJ_OUT/")
code = code.replace("W/'density_report.json'", "TJ_WORK/'density_report.json'")
exec(compile(code, str(source), 'exec'), {'__file__': str(source), '__name__': '__main__', 'TJ_OUT': OUT, 'TJ_WORK': W})
assert json.loads((W / 'density_report.json').read_text()) == json.loads((source.parent / 'density_report.json').read_text())
print('Both figures regenerated; plotted values and density summaries match their prior versions.')
