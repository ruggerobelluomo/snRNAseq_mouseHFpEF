# 01b_integration_benchmark.py
# Code-only export of the analysis pipeline (code cells of the Jupyter notebook). Data are not included; see README.md.

# %%
import os, time, logging, warnings, subprocess
import numpy as np
import pandas as pd
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import anndata as ad
from anndata.io import read_elem
import scanpy as sc
import harmonypy
import scanorama
import torch
import scvi
import scib_metrics
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from IPython.display import SVG, display
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
logging.getLogger('harmonypy').setLevel(logging.WARNING)
logging.getLogger('scvi').setLevel(logging.WARNING)
logging.getLogger('lightning.pytorch').setLevel(logging.WARNING)
%matplotlib inline
sc.settings.n_jobs = 8
sc.settings.verbosity = 0
sc.set_figure_params(vector_friendly=True, dpi_save=300)

mpl.rcParams.update({
    'svg.fonttype': 'none', 'savefig.format': 'svg', 'pdf.fonttype': 42, 'savefig.dpi': 300,
    'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.grid': False, 'grid.alpha': 0.15, 'grid.linewidth': 0.4,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': '0.6', 'axes.linewidth': 0.6, 'xtick.color': '0.6', 'ytick.color': '0.6'})
from importlib.metadata import version
packages = ['scanpy', 'anndata', 'numpy', 'pandas', 'scipy', 'scikit-learn', 'matplotlib', 'harmonypy', 'scanorama',
            'scvi-tools', 'torch', 'lightning', 'scib-metrics', 'jax', 'jaxlib', 'pynndescent', 'leidenalg', 'igraph', 'plottable', 'h5py']
versions = pd.Series({pkg: version(pkg) for pkg in packages}, name='version')
print(versions.to_string())
print('torch CUDA', torch.version.cuda, '| GPU:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')

# %%
PROJECT_ROOT = os.path.normpath(os.environ.get('PROJECT_ROOT', os.path.join(os.getcwd(), '..')))
path_atlas  = os.path.join(PROJECT_ROOT, 'data', 'objects', 'atlas_annotated.h5ad')
path_out    = os.path.join(PROJECT_ROOT, 'results', '01b_integration_benchmark')
path_tables = os.path.join(path_out, 'tables')
path_figs   = os.path.join(path_out, 'figures')
for folder_path in [path_tables, path_figs]:
    os.makedirs(folder_path, exist_ok=True)
print('project root:', PROJECT_ROOT)

N_PCS = 50
N_NEIGHBORS = 15
UMAP_MIN_DIST = 0.2
HARMONY_MAX_ITER = 20
MIN_TOP_SCORE = 0.3
MIN_MARGIN = 0.1
sample_list = ['B2FYH','B3MYH','B3MOC','B5MYC','B5MOH','B6FYC','B6FOH','B7MOC',
               'B4FYH','B5MYH','B5MOC','B6FYC_2','B6FOC','B7MYC','B7MOH','B8FOH']

SEED = 0
SCVI_N_LATENT = 30
N_JOBS = 8

CELL_TYPE_ORDER = ['Cardiomyocyte', 'Endothelial', 'Fibroblast', 'Pericyte', 'VSMC', 'Macrophage',
                   'Lymphocyte', 'Lymphatic_EC', 'Glia', 'Proliferating', 'B_cell', 'Unassigned']
CELL_TYPE_PALETTE = dict(zip(CELL_TYPE_ORDER, ['#D55E00', '#0072B2', '#009E73', '#CC79A7', '#E69F00', '#56B4E9',
                                               '#F0E442', '#8C564B', '#17BECF', '#BCBD22', '#000000', '#7F7F7F']))
CELL_TYPE_PALETTE['Adipocyte'] = '#B8BDC4'
METHODS = ['Unintegrated PCA', 'Harmony (sample)', 'Harmony (sample + seq)', 'ComBat', 'Scanorama', 'scVI', 'scANVI']
METHOD_PALETTE = dict(zip(METHODS, ['#7F7F7F', '#0072B2', '#56B4E9', '#E69F00', '#009E73', '#D55E00', '#CC79A7']))

def save_fig(fig, name):
    fig.savefig(os.path.join(path_figs, name + '.svg'), bbox_inches='tight')
    fig.savefig(os.path.join(path_figs, name + '.png'), bbox_inches='tight', dpi=200)
    display(SVG(filename=os.path.join(path_figs, name + '.svg')))
    plt.close(fig)

runtimes = {}

# %%
t0 = time.time()
with h5py.File(path_atlas, 'r') as f:
    obs = read_elem(f['obs'])
    var = read_elem(f['var'])
    hvg = np.flatnonzero(var['highly_variable'].to_numpy())
    lognorm = read_elem(f['X'])[:, hvg]
    counts = read_elem(f['layers/counts'])[:, hvg]
    obsm = {k: np.asarray(read_elem(f['obsm'][k])) for k in ['X_pca', 'PC_harmony', 'UMAP_PC_harmony']}
adata = ad.AnnData(X=lognorm, obs=obs, var=var.iloc[hvg].copy(), layers={'counts': counts}, obsm=obsm)
adata.obs['cell_type'] = adata.obs['cell_type'].cat.reorder_categories(
    [c for c in CELL_TYPE_ORDER if c in adata.obs['cell_type'].cat.categories])
adata.obs['sample'] = adata.obs['sample'].cat.reorder_categories(sorted(sample_list))
del lognorm, counts
print(f'loaded in {time.time() - t0:.0f}s:', adata)
summary = pd.DataFrame({'nuclei': [adata.n_obs], 'HVGs': [adata.n_vars], 'samples': [adata.obs['sample'].nunique()],
                        'sequencing runs': [adata.obs['seq'].nunique()], 'cell types': [adata.obs['cell_type'].nunique()]})
print(summary.to_string(index=False))
print(pd.crosstab(adata.obs['seq'], adata.obs['sample']).to_string())

# %%
def run_harmony(X_pca, obs, batch_keys):

    harmony_out = harmonypy.run_harmony(X_pca, obs, batch_keys, max_iter_harmony=HARMONY_MAX_ITER, verbose=False)
    harmony_z = harmony_out.Z_corr
    if harmony_z.shape[0] != X_pca.shape[0]:
        harmony_z = harmony_z.T
    return np.asarray(harmony_z, dtype=np.float32)

def agreement(new, stored):
    r = np.array([np.corrcoef(new[:, j], stored[:, j])[0, 1] for j in range(stored.shape[1])])
    aligned = new * np.sign(r)
    return pd.Series({'min |r| across components': np.abs(r).min(), 'median |r|': np.median(np.abs(r)),
                      'max |difference|': np.abs(aligned - stored).max(), 'stored value range': np.ptp(stored)})

t0 = time.time()
scaled = ad.AnnData(adata.X.copy())
sc.pp.scale(scaled, max_value=10)
sc.tl.pca(scaled, n_comps=N_PCS, svd_solver='arpack')
pca_new = scaled.obsm['X_pca']
runtimes['Unintegrated PCA'] = time.time() - t0
del scaled

t0 = time.time()
harmony_new = run_harmony(adata.obsm['X_pca'], adata.obs, ['sample'])
runtimes['Harmony (sample)'] = time.time() - t0

reproduction = pd.DataFrame({'PCA (X_pca)': agreement(pca_new, adata.obsm['X_pca']),
                             'Harmony (PC_harmony)': agreement(harmony_new, adata.obsm['PC_harmony'])}).T
reproduction.to_csv(os.path.join(path_tables, 'reproduction_of_stored_embeddings.csv'))
print(reproduction.round(6).to_string())

# %%
t0 = time.time()
adata.obsm['X_harmony_sample_seq'] = run_harmony(adata.obsm['X_pca'], adata.obs, ['sample', 'seq'])
runtimes['Harmony (sample + seq)'] = time.time() - t0

t0 = time.time()
combat = ad.AnnData(adata.X.copy(), obs=adata.obs[['sample']].copy())
combat.X = np.array(sc.pp.combat(combat, key='sample', inplace=False), dtype=np.float32)
sc.pp.scale(combat, max_value=10)
sc.tl.pca(combat, n_comps=N_PCS, svd_solver='arpack')
adata.obsm['X_combat'] = combat.obsm['X_pca']
runtimes['ComBat'] = time.time() - t0
del combat

t0 = time.time()
per_sample = [adata[adata.obs['sample'] == s].copy() for s in sample_list]
scanorama.integrate_scanpy(per_sample, dimred=N_PCS, seed=SEED, verbose=False)
scanorama_emb = pd.concat([pd.DataFrame(a.obsm['X_scanorama'], index=a.obs_names) for a in per_sample])
adata.obsm['X_scanorama'] = scanorama_emb.loc[adata.obs_names].to_numpy(dtype=np.float32)
runtimes['Scanorama'] = time.time() - t0
del per_sample, scanorama_emb
print({k: round(v) for k, v in runtimes.items()})

# %%
scvi.settings.seed = SEED
scvi_data = ad.AnnData(X=adata.layers['counts'].copy(), obs=adata.obs[['sample', 'seq', 'cell_type']].copy(), var=adata.var[[]].copy())
scvi.model.SCVI.setup_anndata(scvi_data, batch_key='sample', categorical_covariate_keys=['seq'])
t0 = time.time()
scvi_model = scvi.model.SCVI(scvi_data, n_latent=SCVI_N_LATENT)
scvi_model.train(accelerator='gpu', devices=1, enable_progress_bar=False)
adata.obsm['X_scVI'] = scvi_model.get_latent_representation().astype(np.float32)
runtimes['scVI'] = time.time() - t0

t0 = time.time()
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category='Unassigned', labels_key='cell_type')
scanvi_model.train(accelerator='gpu', devices=1, enable_progress_bar=False)
adata.obsm['X_scANVI'] = scanvi_model.get_latent_representation().astype(np.float32)
runtimes['scANVI'] = time.time() - t0

training = pd.DataFrame({
    'model': ['scVI', 'scANVI'],
    'epochs': [len(scvi_model.history['elbo_train']), len(scanvi_model.history['elbo_train'])],
    'final train ELBO': [scvi_model.history['elbo_train'].iloc[-1, 0], scanvi_model.history['elbo_train'].iloc[-1, 0]],
    'runtime (s)': [runtimes['scVI'], runtimes['scANVI']]})
print(training.round(1).to_string(index=False))

# %%
EMBEDDINGS = {'Unintegrated PCA': 'X_pca', 'Harmony (sample)': 'PC_harmony', 'Harmony (sample + seq)': 'X_harmony_sample_seq',
              'ComBat': 'X_combat', 'Scanorama': 'X_scanorama', 'scVI': 'X_scVI', 'scANVI': 'X_scANVI'}
device = {m: 'GPU' if m in ('scVI', 'scANVI') else 'CPU (8 cores)' for m in METHODS}
runtime_table = pd.DataFrame({'method': METHODS,
                              'embedding': [EMBEDDINGS[m] for m in METHODS],
                              'dimensions': [adata.obsm[EMBEDDINGS[m]].shape[1] for m in METHODS],
                              'device': [device[m] for m in METHODS],
                              'integration runtime (s)': [runtimes[m] for m in METHODS]})
runtime_table['note'] = ['PCA recomputed (Step 2)', 'harmonypy re-run (Step 2); stored embedding benchmarked', '', 'ComBat + PCA',
                         '', 'training + latent', 'training from scVI + latent']
runtime_table.to_csv(os.path.join(path_tables, 'integration_runtime.csv'), index=False)
emb_out = ad.AnnData(obs=adata.obs[['sample', 'seq', 'batch', 'cell_type']].copy(),
                     obsm={EMBEDDINGS[m]: adata.obsm[EMBEDDINGS[m]] for m in METHODS})
emb_out.write(os.path.join(path_out, 'embeddings_by_method.h5ad'), compression='gzip')
print(runtime_table.round(0).to_string(index=False))

# %%
bench = ad.AnnData(obs=adata.obs[['sample', 'seq', 'cell_type']].copy(),
                   obsm={m: adata.obsm[EMBEDDINGS[m]] for m in METHODS})
bio_metrics = BioConservation(isolated_labels=True, nmi_ari_cluster_labels_leiden={'n_jobs': N_JOBS},
                              nmi_ari_cluster_labels_kmeans=True, silhouette_label=True, clisi_knn=True)
batch_metrics = BatchCorrection(bras=True, ilisi_knn=True, kbet_per_label=True, graph_connectivity=True, pcr_comparison=True)
t0 = time.time()
bm_sample = Benchmarker(bench, batch_key='sample', label_key='cell_type', embedding_obsm_keys=METHODS,
                        bio_conservation_metrics=bio_metrics, batch_correction_metrics=batch_metrics,
                        pre_integrated_embedding_obsm_key='Unintegrated PCA', n_jobs=N_JOBS, progress_bar=False)
bm_sample.benchmark()
metric_runtime = {'sample': time.time() - t0}
print(f'benchmark (batch_key = sample): {metric_runtime["sample"] / 60:.0f} min')

# %%
def results_tables(bm, tag):
    tables = {}
    for scaled in (False, True):
        res = bm.get_results(min_max_scale=scaled)
        metric_type = res.loc['Metric Type']
        res = res.drop('Metric Type').astype(float).sort_values('Total', ascending=False)
        res.index.name = 'method'
        name = f'scib_results_{tag}_{"minmax" if scaled else "unscaled"}'
        out = res.copy(); out.loc['Metric Type'] = metric_type
        out.to_csv(os.path.join(path_tables, name + '.csv'))
        tables[scaled] = res
    ranked = pd.concat({'unscaled': tables[False][['Bio conservation', 'Batch correction', 'Total']],
                        'min-max scaled': tables[True][['Bio conservation', 'Batch correction', 'Total']]}, axis=1)
    ranked[('unscaled', 'rank')] = ranked[('unscaled', 'Total')].rank(ascending=False).astype(int)
    ranked[('min-max scaled', 'rank')] = ranked[('min-max scaled', 'Total')].rank(ascending=False).astype(int)
    ranked = ranked.sort_values(('unscaled', 'Total'), ascending=False)
    ranked.to_csv(os.path.join(path_tables, f'scib_ranked_{tag}.csv'))
    return tables, metric_type, ranked

def scib_table_plot(bm, name):
    for scaled in (False, True):
        tab = bm.plot_results_table(min_max_scale=scaled, show=False)
        suffix = 'minmax' if scaled else 'unscaled'
        tab.ax.figure.savefig(os.path.join(path_figs, f'{name}_{suffix}.svg'), bbox_inches='tight')
        tab.ax.figure.savefig(os.path.join(path_figs, f'{name}_{suffix}.png'), bbox_inches='tight', dpi=200)
        plt.close(tab.ax.figure)
    display(SVG(filename=os.path.join(path_figs, f'{name}_unscaled.svg')))

res_sample, metric_type, ranked_sample = results_tables(bm_sample, 'batch_sample')
scib_table_plot(bm_sample, 'scib_table_batch_sample')
print(res_sample[False].round(3).to_string())
print()
print(ranked_sample.round(3).to_string())

# %%
labels = adata.obs['cell_type'].to_numpy()
silhouette_batch = pd.Series({m: scib_metrics.silhouette_batch(adata.obsm[EMBEDDINGS[m]], labels, adata.obs['sample'].to_numpy())
                              for m in METHODS}, name='silhouette batch (sample)')
silhouette_batch.to_csv(os.path.join(path_tables, 'silhouette_batch_sample.csv'))
print(pd.concat([silhouette_batch, res_sample[False]['BRAS'].rename('BRAS (sample)')], axis=1).loc[METHODS].round(3).to_string())

# %%
t0 = time.time()
bm_seq = Benchmarker(bench, batch_key='seq', label_key='cell_type', embedding_obsm_keys=METHODS,
                     bio_conservation_metrics=BioConservation(isolated_labels=True, nmi_ari_cluster_labels_leiden=False,
                                                              nmi_ari_cluster_labels_kmeans=False, silhouette_label=False, clisi_knn=True),
                     batch_correction_metrics=batch_metrics,
                     pre_integrated_embedding_obsm_key='Unintegrated PCA', n_jobs=N_JOBS, progress_bar=False)
bm_seq.benchmark()
metric_runtime['seq'] = time.time() - t0

carried = bm_sample._results.loc[[i for i in bm_sample._results.index if i not in bm_seq._results.index]]
bm_seq._results = pd.concat([bm_seq._results, carried]).loc[bm_sample._results.index]
bm_seq._bio_conservation_metrics = bio_metrics
res_seq, _, ranked_seq = results_tables(bm_seq, 'batch_seq')
scib_table_plot(bm_seq, 'scib_table_batch_seq')
print(f'benchmark (batch_key = seq): {metric_runtime["seq"] / 60:.0f} min')
print(res_seq[False].round(3).to_string())
print()
print(ranked_seq.round(3).to_string())

# %%
def clean_table(res, metric_type, title, name):
    bio_cols = [c for c in res.columns if metric_type.get(c) == 'Bio conservation']
    batch_cols = [c for c in res.columns if metric_type.get(c) == 'Batch correction']
    agg_cols = ['Bio conservation', 'Batch correction', 'Total']
    cols = bio_cols + batch_cols + agg_cols
    data = res[cols]
    fig, ax = plt.subplots(figsize=(0.62 * len(cols) + 1.8, 0.32 * len(data) + 1.1))
    x_pos, x = [], 0.0
    for j, c in enumerate(cols):
        if j in (len(bio_cols), len(bio_cols) + len(batch_cols)):
            x += 0.35
        x_pos.append(x); x += 1
    for j, c in enumerate(cols):
        values = data[c].to_numpy()
        span = np.ptp(values) or 1.0
        cmap = plt.get_cmap('YlGnBu' if c in agg_cols else 'Greens')
        for i, v in enumerate(values):
            shade = 0.15 + 0.6 * (v - values.min()) / span
            ax.add_patch(plt.Rectangle((x_pos[j] - 0.46, i - 0.42), 0.92, 0.84, color=cmap(shade), lw=0))
            ax.text(x_pos[j], i, f'{v:.2f}', ha='center', va='center', fontsize=7, color='white' if shade > 0.6 else '0.15')
    ax.set_xlim(-0.6, x_pos[-1] + 0.6); ax.set_ylim(len(data) - 0.5, -1.6)
    ax.set_yticks(range(len(data))); ax.set_yticklabels(data.index, fontsize=8)
    for lab in ax.get_yticklabels():
        lab.set_color(METHOD_PALETTE[lab.get_text()])
    ax.set_xticks(x_pos); ax.set_xticklabels([c.replace('Bio conservation', 'Bio').replace('Batch correction', 'Batch') for c in cols],
                                             rotation=45, ha='left', fontsize=7)
    ax.xaxis.tick_top(); ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    for group, members in [('Bio conservation', bio_cols), ('Batch correction', batch_cols), ('Aggregate', agg_cols)]:
        xs = [x_pos[cols.index(c)] for c in members]
        ax.plot([xs[0] - 0.4, xs[-1] + 0.4], [-0.75, -0.75], color='0.6', lw=0.6)
        ax.text(np.mean(xs), -1.0, group, ha='center', va='bottom', fontsize=7, color='0.3')
    ax.set_title(title, fontsize=8, pad=56, loc='left')
    save_fig(fig, name)

for tag, res in [('sample', res_sample), ('seq', res_seq)]:
    for scaled in (False, True):
        scaling = 'min-max scaled' if scaled else 'unscaled'
        clean_table(res[scaled], metric_type, f'scib-metrics, batch_key = {tag}, {scaling} (total = 0.6 bio + 0.4 batch)',
                    f'scib_table_clean_batch_{tag}_{"minmax" if scaled else "unscaled"}')

# %%
fig, axes = plt.subplots(1, 4, figsize=(13, 3.2))
for ax, (tag, scaled) in zip(axes, [('sample', False), ('sample', True), ('seq', False), ('seq', True)]):
    res = (res_sample if tag == 'sample' else res_seq)[scaled]
    for m in METHODS:
        ax.scatter(res.loc[m, 'Batch correction'], res.loc[m, 'Bio conservation'], s=36, color=METHOD_PALETTE[m],
                   edgecolor='white', linewidth=0.5, zorder=3)
    ax.set_xlabel('Batch correction score'); ax.set_ylabel('Bio conservation score')
    ax.set_title(f'batch_key = {tag}, {"min-max scaled" if scaled else "unscaled"}')
    pad_x = 0.15 * np.ptp(res['Batch correction']); pad_y = 0.15 * np.ptp(res['Bio conservation'])
    ax.set_xlim(res['Batch correction'].min() - pad_x, res['Batch correction'].max() + pad_x)
    ax.set_ylim(res['Bio conservation'].min() - pad_y, res['Bio conservation'].max() + pad_y)
axes[-1].legend(handles=[Line2D([], [], marker='o', ls='', color=METHOD_PALETTE[m], label=m) for m in METHODS],
                frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5))
fig.tight_layout()
save_fig(fig, 'batch_vs_bio_scatter')

# %%
score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
panel_names = np.array([c.replace('score_', '') for c in score_cols])
scores = adata.obs[score_cols].to_numpy()
order = np.argsort(scores, axis=1)
panel_label = panel_names[order[:, -1]]
top = scores[np.arange(len(scores)), order[:, -1]]
margin = top - scores[np.arange(len(scores)), order[:, -2]]
confident = (top >= MIN_TOP_SCORE) & (margin >= MIN_MARGIN)
adata.obs['panel_label'] = panel_label

labelled = adata.obs['cell_type'] != 'Unassigned'
agree = (adata.obs['panel_label'] == adata.obs['cell_type'].astype(str))
panel_summary = pd.DataFrame({
    'nuclei': adata.obs['panel_label'].value_counts(),
    'confident nuclei': adata.obs.loc[confident, 'panel_label'].value_counts()}).fillna(0).astype(int)
panel_summary['confident (%)'] = 100 * panel_summary['confident nuclei'] / panel_summary['nuclei']
panel_summary.to_csv(os.path.join(path_tables, 'panel_label_summary.csv'))
crosstab = pd.crosstab(adata.obs['cell_type'], adata.obs['panel_label'])
crosstab.to_csv(os.path.join(path_tables, 'celltype_vs_panel_label_counts.csv'))
print(f'confident nuclei: {confident.sum():,} of {len(confident):,} ({100 * confident.mean():.1f} %)')
print(f'panel label = cell_type: {100 * agree[labelled].mean():.1f} % of labelled nuclei, '
      f'{100 * agree[labelled & confident].mean():.1f} % of confident labelled nuclei')
print(panel_summary.sort_values('nuclei', ascending=False).round(1).to_string())

# %%
def nearest_without_self(indices, k=N_NEIGHBORS):

    first = indices[:, :k + 1]
    keep = first != np.arange(len(first))[:, None]
    keep[keep.all(axis=1), k] = False
    return first[keep].reshape(len(first), k)

def knn_purity(indices, lab, subset):
    neighbours = nearest_without_self(indices)
    per_nucleus = (lab[neighbours] == lab[:, None]).mean(axis=1)
    per_label = pd.Series(per_nucleus[subset]).groupby(lab[subset]).mean()
    return per_label.mean(), per_label

label_free, purity_by_panel = {}, {}
for m in METHODS:
    neighbours = bm_sample._emb_adatas[m].uns['15_neighbor_res']
    wide = bm_sample._emb_adatas[m].uns['50_neighbor_res'].indices
    purity_all, _ = knn_purity(wide, panel_label, np.ones(len(panel_label), bool))
    purity_conf, per_panel = knn_purity(wide, panel_label, confident)
    nmi_ari = scib_metrics.nmi_ari_cluster_labels_leiden(neighbours, panel_label, n_jobs=N_JOBS)
    label_free[m] = {'kNN panel purity, all nuclei': purity_all, 'kNN panel purity, confident nuclei': purity_conf,
                     'Leiden NMI vs panel label': nmi_ari['nmi'], 'Leiden ARI vs panel label': nmi_ari['ari']}
    purity_by_panel[m] = per_panel
label_free = pd.DataFrame(label_free).T
label_free['iLISI (sample)'] = res_sample[False]['iLISI']
label_free['PCR comparison (sample)'] = res_sample[False]['PCR comparison']
bio_free = ['kNN panel purity, all nuclei', 'kNN panel purity, confident nuclei', 'Leiden NMI vs panel label', 'Leiden ARI vs panel label']
label_free['label-free bio'] = label_free[bio_free].mean(axis=1)
label_free['label-free batch'] = label_free[['iLISI (sample)', 'PCR comparison (sample)']].mean(axis=1)
label_free['label-free total'] = 0.6 * label_free['label-free bio'] + 0.4 * label_free['label-free batch']
label_free['label-based bio (scib)'] = res_sample[False]['Bio conservation']
label_free['rank label-free bio'] = label_free['label-free bio'].rank(ascending=False).astype(int)
label_free['rank label-based bio'] = label_free['label-based bio (scib)'].rank(ascending=False).astype(int)
label_free = label_free.sort_values('label-free total', ascending=False)
label_free.to_csv(os.path.join(path_tables, 'label_free_check.csv'))
purity_by_panel = pd.DataFrame(purity_by_panel).T
purity_by_panel.to_csv(os.path.join(path_tables, 'knn_panel_purity_by_panel_confident.csv'))
print(label_free.round(3).to_string())

# %%
panel_order = panel_summary.sort_values('confident nuclei', ascending=False).index.tolist()
data = purity_by_panel.loc[label_free.index, [p for p in panel_order if p in purity_by_panel.columns]]
fig, ax = plt.subplots(figsize=(0.55 * data.shape[1] + 2, 0.32 * data.shape[0] + 1.3))
im = ax.imshow(data.to_numpy(), cmap='Greens', vmin=0, vmax=1, aspect='auto')
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        v = data.iat[i, j]
        ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=6, color='white' if v > 0.65 else '0.2')
ax.set_xticks(range(data.shape[1])); ax.set_xticklabels([f'{p} ({panel_summary.loc[p, "confident nuclei"]:,})' for p in data.columns],
                                                        rotation=45, ha='right')
ax.set_yticks(range(data.shape[0])); ax.set_yticklabels(data.index)
ax.tick_params(length=0)
for side in ax.spines.values():
    side.set_visible(False)
cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cbar.set_label('kNN panel purity (proportion)'); cbar.outline.set_visible(False)
ax.set_title('Proportion of the 15 nearest neighbours with the same top marker panel (confident nuclei; n per panel)')
save_fig(fig, 'knn_panel_purity_by_panel')

# %%
umaps = {'Harmony (sample)': adata.obsm['UMAP_PC_harmony']}
t0 = time.time()
for m in METHODS:
    if m in umaps:
        continue
    graph = ad.AnnData(obs=adata.obs[[]].copy(), obsm={'emb': adata.obsm[EMBEDDINGS[m]]})
    sc.pp.neighbors(graph, n_neighbors=N_NEIGHBORS, use_rep='emb', metric='cosine', random_state=SEED)
    sc.tl.umap(graph, min_dist=UMAP_MIN_DIST, random_state=SEED)
    umaps[m] = graph.obsm['X_umap']
print(f'UMAPs: {(time.time() - t0) / 60:.0f} min')

rng = np.random.RandomState(SEED)
draw_order = rng.permutation(adata.n_obs)
sample_cats = adata.obs['sample'].cat.categories
sample_palette = {s: mpl.colors.to_hex(plt.get_cmap('tab20')(k % 20)) for k, s in enumerate(sample_cats)}
colour_rows = [('cell_type', adata.obs['cell_type'].astype(str).map(CELL_TYPE_PALETTE).to_numpy()),
               ('sample', adata.obs['sample'].astype(str).map(sample_palette).to_numpy())]
fig, axes = plt.subplots(2, len(METHODS), figsize=(2.3 * len(METHODS) + 2.2, 5.0))
for j, m in enumerate(METHODS):
    xy = umaps[m][draw_order]
    for i, (key, colours) in enumerate(colour_rows):
        ax = axes[i, j]
        ax.scatter(xy[:, 0], xy[:, 1], s=0.3, c=colours[draw_order], linewidths=0, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        for side in ax.spines.values():
            side.set_visible(False)
        if i == 0:
            ax.set_title(m, color=METHOD_PALETTE[m])
        if j == 0:
            ax.set_ylabel(f'coloured by {key}')
axes[0, -1].legend(handles=[Line2D([], [], marker='o', ls='', color=CELL_TYPE_PALETTE[c], label=c) for c in adata.obs['cell_type'].cat.categories],
                   frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5), markerscale=0.9)
axes[1, -1].legend(handles=[Line2D([], [], marker='o', ls='', color=sample_palette[s], label=s) for s in sample_cats],
                   frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=2, markerscale=0.9)
fig.suptitle(f'UMAP per integration method ({adata.n_obs:,} nuclei, {len(sample_cats)} samples)', y=1.0, fontsize=9)
fig.tight_layout()
save_fig(fig, 'umap_grid_by_method')

# %%
runtime_table['scib benchmark, batch_key = sample (s, all methods)'] = metric_runtime['sample']
runtime_table['scib benchmark, batch_key = seq (s, all methods)'] = metric_runtime['seq']
runtime_table.to_csv(os.path.join(path_tables, 'integration_runtime.csv'), index=False)
versions.to_csv(os.path.join(path_tables, 'software_versions.csv'))
print(runtime_table[['method', 'device', 'integration runtime (s)']].round(0).to_string(index=False))
print(f"scib benchmark: {metric_runtime['sample'] / 60:.0f} min (sample), {metric_runtime['seq'] / 60:.0f} min (seq)")
