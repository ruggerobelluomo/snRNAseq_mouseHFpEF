"""Execute a notebook top-to-bottom in one in-process IPython shell (no Jupyter/ZMQ kernel), cells in order, stopping at
the first error, and write the outputs back into the notebook. Usage (from notebooks/): python run_notebook.py notebook.ipynb"""
import sys, time, traceback
import nbformat as nbf
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output

path = sys.argv[1]
nb = nbf.read(path, as_version=4)
from IPython.core.displayhook import DisplayHook
class QuietDisplayHook(DisplayHook):
    # keep result capture (fill_exec_result) but do not print Out[n] to stdout; the runner formats it
    def write_output_prompt(self): pass
    def write_format_data(self, format_dict, md_dict=None): pass
    def finish_displayhook(self): pass
class HeadlessShell(InteractiveShell):
    displayhook_class = QuietDisplayHook
    def enable_gui(self, gui=None):
        pass
import matplotlib; matplotlib.use('module://matplotlib_inline.backend_inline')
shell = HeadlessShell.instance(colors='NoColor')
shell.display_formatter.active_types = ['text/plain', 'text/html', 'image/png', 'image/svg+xml']
count = 0
t_all = time.time()
for i, cell in enumerate(nb.cells):
    if cell.cell_type != 'code':
        continue
    count += 1
    cell.outputs = []
    t0 = time.time()
    with capture_output(stdout=True, stderr=True, display=True) as cap:
        res = shell.run_cell(cell.source, store_history=True, silent=False)
    outs = []
    if cap.stdout:
        outs.append(nbf.v4.new_output('stream', name='stdout', text=cap.stdout))
    if cap.stderr:
        outs.append(nbf.v4.new_output('stream', name='stderr', text=cap.stderr))
    for o in cap.outputs:
        outs.append(nbf.v4.new_output('display_data', data=o.data, metadata=o.metadata or {}))
    if res.result is not None:
        data, md = shell.display_formatter.format(res.result)
        outs.append(nbf.v4.new_output('execute_result', data=data, metadata=md or {}, execution_count=count))
    err = res.error_before_exec or res.error_in_exec
    if err is not None:
        tb = traceback.format_exception(type(err), err, err.__traceback__)
        outs.append(nbf.v4.new_output('error', ename=type(err).__name__, evalue=str(err), traceback=tb))
    cell.outputs = outs
    cell.execution_count = count
    nbf.write(nb, path)
    print(f'[cell {i} / exec {count}] {time.time() - t0:.0f}s  total {(time.time() - t_all) / 60:.1f} min', flush=True)
    if err is not None:
        print(f'ERROR in cell {i}: {type(err).__name__}: {err}', flush=True)
        print(''.join(tb[-15:]), flush=True)
        sys.exit(1)
print('DONE', flush=True)
