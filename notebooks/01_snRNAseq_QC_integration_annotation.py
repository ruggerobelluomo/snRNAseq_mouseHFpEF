# 01_snRNAseq_QC_integration_annotation.py
# Code-only export of the analysis pipeline (code cells of the Jupyter notebook). Data are not included; see README.md.

# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import os
import logging
import scanpy as sc
import anndata
import harmonypy
import scrublet as scr
import scipy.sparse as sp
from scipy.stats import entropy, mannwhitneyu
import warnings
from IPython.display import SVG, display
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
logging.getLogger('harmonypy').setLevel(logging.WARNING)
%matplotlib inline
sc.settings.n_jobs = 8
sc.settings.file_format_figs = 'svg'
sc.set_figure_params(vector_friendly=True, dpi_save=300)

mpl.rcParams.update({
    'svg.fonttype': 'none', 'savefig.format': 'svg', 'pdf.fonttype': 42, 'savefig.dpi': 300,
    'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.grid': False, 'grid.alpha': 0.15, 'grid.linewidth': 0.4,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': '0.6', 'axes.linewidth': 0.6, 'xtick.color': '0.6', 'ytick.color': '0.6'})
from importlib.metadata import version
for pkg in ['scanpy', 'anndata', 'numpy', 'pandas', 'scipy', 'matplotlib', 'harmonypy', 'scrublet', 'leidenalg', 'igraph', 'h5py']:
    print(pkg, version(pkg))

# %%
PROJECT_ROOT = os.path.normpath(os.environ.get('PROJECT_ROOT', os.path.join(os.getcwd(), '..')))
path_samples = os.path.join(PROJECT_ROOT, 'data', 'raw')
path_out     = PROJECT_ROOT
path_objects = os.path.join(path_out, 'data', 'objects')
path_tables  = os.path.join(path_out, 'results', '01_qc_atlas', 'tables')
path_figs    = os.path.join(path_out, 'results', '01_qc_atlas', 'figures')
qc_dir       = os.path.join(path_figs, 'qc_supplementary')
scrublet_dir = os.path.join(path_out, 'results', '01_qc_atlas', 'scrublet')
for folder_path in [path_objects, path_tables, path_figs, qc_dir, scrublet_dir]:
    os.makedirs(folder_path, exist_ok=True)
print('project root:', PROJECT_ROOT)

sample_list = ['B2FYH','B3MYH','B3MOC','B5MYC','B5MOH','B6FYC','B6FOH','B7MOC',
               'B4FYH','B5MYH','B5MOC','B6FYC_2','B6FOC','B7MYC','B7MOH','B8FOH']
individual_map = {sample_id: sample_id for sample_id in sample_list}

MIN_GENES = 150
MIN_UMI = 200
MAX_PCT_MITO = 5.0
IQR_MULT = 1.5
OVERCLUSTER_RES = 2.0
MAP_RES = 0.8
N_PCS = 50
N_HVG = 2000
N_NEIGHBORS = 15
DOUBLET_FALLBACK = 0.30
TECH_EXON_PROP_PCTL = 95
FIXED_LIMIT_PCTL = 99

SCRUBLET_MIN_CALLED = 0.005
SCRUBLET_MIN_DETECTABLE = 0.10
SCRUBLET_THRESHOLDS = {}

LOW_YIELD_LIBRARIES = {'B5MOH': 'low_yield_library'}

CELLBENDER_DROPLETS = 25000
EMPTY_PLATEAU_UMI = 200

MIN_TOP_SCORE = 0.3
MIN_MARGIN = 0.1

LOWC_REL = 0.5
INCOH_REL = 0.25

COLUMN_UNITS = {
    'med_mito': ('med_mito_pct', '%'), 'med_mito_cellbender': ('med_mito_cellbender_pct', '%'),
    'median_mito': ('median_mito_pct', '%'), 'median_pct_mito': ('median_mito_pct', '%'),
    'med_exon_prop': ('med_exon_prop', 'proportion'), 'frac_pred_dbl': ('pred_doublet_prop', 'proportion'),
    'frac_pred_doublet': ('pred_doublet_prop', 'proportion'), 'largest_sample_fraction': ('largest_sample_prop', 'proportion'),
    'frac_removed': ('removed_prop', 'proportion'), 'frac_mito_positive': ('mito_positive_prop', 'proportion'),
    'predicted_fraction': ('predicted_prop', 'proportion'), 'detectable_fraction': ('detectable_prop', 'proportion'),
    'estimated_overall_fraction': ('estimated_overall_prop', 'proportion'),
    'cb_frac_droplets_nonempty': ('cb_droplets_nonempty_prop', 'proportion'), 'cb_frac_counts_removed': ('cb_counts_removed_prop', 'proportion')}
EXPORT_NAMES = {column: exported for column, (exported, unit) in COLUMN_UNITS.items()}
DISPLAY_NAMES = {column: f'{column} ({unit})' for column, (exported, unit) in COLUMN_UNITS.items()}

GROUP_ORDER = ['YC', 'OC', 'YH', 'OH']
GROUP_FILL = {'YC': '#9ecae1', 'OC': '#1f78b4', 'YH': '#fdbf6f', 'OH': '#e31a1c'}
CELL_TYPE_ORDER = ['Cardiomyocyte', 'Endothelial', 'Fibroblast', 'Pericyte', 'VSMC', 'Macrophage', 'Dendritic_cell',
            'Lymphocyte', 'Lymphatic_EC', 'Glia', 'Proliferating', 'B_cell']
CELL_TYPE_PALETTE = dict(zip(CELL_TYPE_ORDER, ['#D55E00', '#0072B2', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#882255',
                                 '#F0E442', '#8C564B', '#17BECF', '#BCBD22', '#000000']))
CELL_TYPE_PALETTE.update({'Doublet': '#B8BDC4', 'Unassigned': '#7F7F7F'})
GREY = '#B8BDC4'; KEEP = '#2E6FB0'; REMOVE = '#C0392B'

def save_svg(fig, name, folder=None):
    path = os.path.join(folder or qc_dir, name + '.svg')
    fig.savefig(path, bbox_inches='tight')
    display(SVG(filename=path))

# %%
sex_map = {'B2FYH':'F','B3MYH':'M','B3MOC':'M','B5MYC':'M','B5MOH':'M','B6FYC':'F',
           'B6FOH':'F','B7MOC':'M','B4FYH':'F','B5MYH':'M','B5MOC':'M','B6FYC_2':'F',
           'B6FOC':'F','B7MYC':'M','B7MOH':'M','B8FOH':'F'}

disease_map = {'B2FYH':'YH','B3MYH':'YH','B3MOC':'OC','B5MYC':'YC','B5MOH':'OH','B6FYC':'YC',
               'B6FOH':'OH','B7MOC':'OC','B4FYH':'YH','B5MYH':'YH','B5MOC':'OC','B6FYC_2':'YC',
               'B6FOC':'OC','B7MYC':'YC','B7MOH':'OH','B8FOH':'OH'}
batch_map = {'B2FYH':'B2','B3MYH':'B3','B3MOC':'B3','B5MYC':'B5','B5MOH':'B5','B6FYC':'B6',
             'B6FOH':'B6','B7MOC':'B7','B4FYH':'B4','B5MYH':'B5','B5MOC':'B5','B6FYC_2':'B6',
             'B6FOC':'B6','B7MYC':'B7','B7MOH':'B7','B8FOH':'B8'}
seq_map = {'B2FYH':'S2','B3MYH':'S2','B3MOC':'S2','B5MYC':'S2','B5MOH':'S2','B6FYC':'S2',
           'B6FOH':'S2','B7MOC':'S2','B4FYH':'S3','B5MYH':'S3','B5MOC':'S3','B6FYC_2':'S3',
           'B6FOC':'S3','B7MYC':'S3','B7MOH':'S3','B8FOH':'S3'}

age_map = {k:('young' if v[0]=='Y' else 'old') for k,v in disease_map.items()}
cond_map = {k:('HFpEF' if v[1]=='H' else 'control') for k,v in disease_map.items()}

phenos = pd.DataFrame({'sample':sample_list})
phenos['sex'] = phenos['sample'].map(sex_map)
phenos['disease'] = phenos['sample'].map(disease_map)
phenos['age'] = phenos['sample'].map(age_map)
phenos['condition'] = phenos['sample'].map(cond_map)
phenos['batch'] = phenos['sample'].map(batch_map)
phenos['seq'] = phenos['sample'].map(seq_map)
phenos

# %%
def counts_entropy(counts):

    counts = sp.csr_matrix(counts, dtype=np.float64)
    total_counts = np.asarray(counts.sum(1)).ravel()
    proportions = sp.diags(1 / np.clip(total_counts, 1, None)) @ counts
    p_log_p = proportions.copy(); p_log_p.data = proportions.data * np.log2(proportions.data)
    return -np.asarray(p_log_p.sum(1)).ravel()

def import_sample(sample, path_samples, phenos):

    adata_cellranger = sc.read_10x_h5(os.path.join(path_samples, sample, 'raw_feature_bc_matrix.h5'))
    adata_cellranger.var_names_make_unique()

    adata_cellbender = sc.read_10x_h5(os.path.join(path_samples, sample, 'cellbender_filtered.h5'))
    adata_cellbender.var_names_make_unique()

    cellranger_gene_ids = pd.Index(adata_cellranger.var['gene_ids'])
    cellbender_gene_ids = pd.Index(adata_cellbender.var['gene_ids'])
    assert cellranger_gene_ids.is_unique and cellbender_gene_ids.is_unique, 'Duplicate gene IDs in input'
    gene_order = cellranger_gene_ids.get_indexer(cellbender_gene_ids)
    assert (gene_order >= 0).all(), 'CellBender genes missing from CellRanger input'
    adata_cellranger = adata_cellranger[adata_cellbender.obs_names, gene_order]
    adata_cellbender.layers['cellranger_raw'] = adata_cellranger.X.copy()

    adata_cellbender.layers['counts'] = adata_cellbender.X.copy()

    adata_cellbender.obs['sample'] = sample
    adata_cellbender.obs['individual'] = individual_map[sample]
    for phenotype in [c for c in phenos.columns if c != 'sample']:
        adata_cellbender.obs[phenotype] = phenos.loc[phenos['sample'] == sample, phenotype].iloc[0]

    adata_cellbender.obs['cellranger_ncount'] = np.asarray(adata_cellbender.layers['cellranger_raw'].sum(1)).ravel()
    adata_cellbender.obs['cellbender_ncount'] = np.asarray(adata_cellbender.X.sum(1)).ravel()
    adata_cellbender.obs['cellranger_ngenes'] = np.asarray((adata_cellbender.layers['cellranger_raw'] > 0).sum(1)).ravel()
    adata_cellbender.obs['cellbender_ngenes'] = np.asarray((adata_cellbender.X > 0).sum(1)).ravel()

    mito_genes = adata_cellbender.var_names.str.startswith('mt-')
    adata_cellbender.obs['percent_mito'] = 100 * np.asarray(adata_cellbender[:, mito_genes].layers['cellranger_raw'].sum(1)).ravel() / \
                                   np.clip(adata_cellbender.obs['cellranger_ncount'].values, 1, None)
    adata_cellbender.obs['percent_mito_cellbender'] = 100 * np.asarray(adata_cellbender[:, mito_genes].X.sum(1)).ravel() / \
                                              np.clip(adata_cellbender.obs['cellbender_ncount'].values, 1, None)

    adata_cellbender = adata_cellbender[adata_cellbender.obs['cellbender_ncount'] > 0].copy()
    adata_cellbender.obs['cellbender_entropy'] = counts_entropy(adata_cellbender.X)
    adata_cellbender.obs['cellranger_entropy'] = counts_entropy(adata_cellbender.layers['cellranger_raw'])
    adata_cellbender.obs['log_genes_x_entropy'] = np.log(adata_cellbender.obs['cellbender_ngenes']) * adata_cellbender.obs['cellbender_entropy']

    intron_metrics = pd.read_csv(os.path.join(path_samples, sample, 'intron_metrics.csv'), index_col='barcode')
    adata_cellbender.obs['exon_prop'] = intron_metrics['exon_prop'].reindex(adata_cellbender.obs_names).values
    assert adata_cellbender.obs['exon_prop'].notna().all(), 'nuclei missing from intron_metrics.csv'
    return adata_cellbender

# %%
for sample_id in sample_list:
    print(sample_id)
    adata = import_sample(sample_id, path_samples, phenos)
    adata.write(os.path.join(path_objects, sample_id + '.preqc.h5ad'))

# %%
scrublet_runs = {}
for sample_id in sample_list:
    print(sample_id)
    adata = sc.read_h5ad(os.path.join(path_objects, sample_id + '.preqc.h5ad'))
    for source, layer in [('cellbender', 'counts'), ('cellranger', 'cellranger_raw')]:
        scrublet_obj = scr.Scrublet(adata.layers[layer], random_state=0)
        scores, predicted = scrublet_obj.scrub_doublets(verbose=False)
        sim = scrublet_obj.doublet_scores_sim_
        pd.DataFrame({'observed': pd.Series(scores), 'simulated': pd.Series(sim)}).to_csv(
            os.path.join(path_tables, sample_id + '_' + source + '_scrublet_scores.csv.gz'), index=False)
        scrublet_runs[(sample_id, source)] = {'scores': scores, 'sim': sim, 'n': adata.n_obs,
                                      'threshold_auto': getattr(scrublet_obj, 'threshold_', np.nan)}

def bimodal(run):

    auto_threshold = run['threshold_auto']
    return np.isfinite(auto_threshold) and np.mean(run['sim'] > auto_threshold) >= SCRUBLET_MIN_DETECTABLE

def threshold_rule(run, key, fallback):
    auto_threshold = run['threshold_auto']
    if key in SCRUBLET_THRESHOLDS:
        return SCRUBLET_THRESHOLDS[key], 'manual override'
    if bimodal(run) and np.mean(run['scores'] > auto_threshold) >= SCRUBLET_MIN_CALLED:
        return auto_threshold, 'automatic'
    if bimodal(run):
        return auto_threshold, 'simulated-histogram minimum'
    return fallback, 'fixed fallback (median automatic threshold)'

auto_ok = {k: bimodal(scrub_run) and np.mean(scrub_run['scores'] > scrub_run['threshold_auto']) >= SCRUBLET_MIN_CALLED for k, scrub_run in scrublet_runs.items()}
fallback_cutoff = {src: float(np.median([scrublet_runs[k]['threshold_auto'] for k in scrublet_runs if k[1] == src and auto_ok[k]]))
                   for src in ['cellbender', 'cellranger']}
print('fixed fallback cutoff per source:', {k: round(v, 3) for k, v in fallback_cutoff.items()})

scrublet_summary = []
for sample_id in sample_list:
    adata = sc.read_h5ad(os.path.join(path_objects, sample_id + '.preqc.h5ad'))
    manual = False
    for source in ['cellbender', 'cellranger']:
        run = scrublet_runs[(sample_id, source)]
        threshold, method = threshold_rule(run, (sample_id, source), fallback_cutoff[source])
        predicted = run['scores'] > threshold
        detectable = float(np.mean(run['sim'] > threshold))
        adata.obs[source + '_doublet_scores'] = run['scores']
        adata.obs[source + '_predicted_doublets'] = predicted
        manual |= (source == 'cellranger') and (method != 'automatic')
        scrublet_summary.append({'sample': sample_id, 'source': source, 'n': run['n'],
            'threshold_auto': run['threshold_auto'], 'threshold_used': threshold, 'method': method,
            'predicted_fraction': float(np.mean(predicted)), 'detectable_fraction': detectable,
            'estimated_overall_fraction': np.mean(predicted) / detectable if detectable > 0 else np.nan})
        fig, axes = plt.subplots(1, 2, figsize=(8, 3))
        for ax, values, title in zip(axes, [run['scores'], run['sim']], ['Observed nuclei', 'Simulated doublets']):
            ax.hist(values, bins=50, color=KEEP)
            ax.axvline(threshold, color=REMOVE, ls='--', label=f'used {threshold:.2f} ({method})')
            if np.isfinite(run['threshold_auto']):
                ax.axvline(run['threshold_auto'], color='k', ls=':', label=f"automatic {run['threshold_auto']:.2f}")
            ax.set(title=title, xlabel='Scrublet score', ylabel='Nuclei')
        axes[0].legend(frameon=False)
        fig.suptitle(sample_id + ' - ' + source)
        fig.tight_layout()
        fig.savefig(os.path.join(scrublet_dir, sample_id + '_' + source + '.svg'), bbox_inches='tight')
        plt.close(fig)
    adata.obs['qc_flag'] = LOW_YIELD_LIBRARIES.get(sample_id, 'manual_doublet_threshold' if manual else 'pass')
    adata.write(os.path.join(path_objects, sample_id + '.preqc.scrub.h5ad'))

scrublet_summary = pd.DataFrame(scrublet_summary)
scrublet_summary['qc_flag'] = scrublet_summary['sample'].map(
    {sample_id: LOW_YIELD_LIBRARIES.get(sample_id, 'manual_doublet_threshold' if needs_manual else 'pass')
     for sample_id, needs_manual in scrublet_summary[scrublet_summary['source'] == 'cellranger'].set_index('sample')['method'].ne('automatic').items()})
scrublet_summary.to_csv(os.path.join(path_tables, 'scrublet_summary.csv'), index=False)
print(scrublet_summary[scrublet_summary['source'] == 'cellranger'].rename(columns=DISPLAY_NAMES).round(3).to_string())

# %%
toaggr = [sc.read_h5ad(os.path.join(path_objects, sample_id + '.preqc.scrub.h5ad')) for sample_id in sample_list]
adata = anndata.concat(toaggr, join='outer', merge='first', index_unique='-', label='batch_concat')
adata.obs_names_make_unique()
adata.write(os.path.join(path_objects, 'data.PreQC.Raw.h5ad'))
print(adata.shape)
print(adata.obs.groupby('sample', observed=True)['qc_flag'].first().to_string())

# %%
def run_harmony(X_pca, obs, batch_key):
    harmony_out = harmonypy.run_harmony(X_pca, obs, [batch_key], max_iter_harmony=20, verbose=False)
    harmony_z = harmony_out.Z_corr
    if harmony_z.shape[0] != X_pca.shape[0]:
        harmony_z = harmony_z.T
    assert harmony_z.shape == X_pca.shape, (harmony_z.shape, X_pca.shape)
    return np.asarray(harmony_z)

# %%
adata = sc.read_h5ad(os.path.join(path_objects, 'data.PreQC.Raw.h5ad'))

sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor='seurat_v3', batch_key='sample')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver='arpack')
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

# %%
adata.obs['batch_harmony'] = adata.obs['sample'].astype('category')
adata.obsm['PC_harmony'] = run_harmony(adata.obsm['X_pca'], adata.obs, 'batch_harmony')

sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep='PC_harmony',
                metric='cosine', key_added='PC_harmony')
sc.tl.umap(adata, neighbors_key='PC_harmony', min_dist=0.2)
adata.obsm['UMAP_PC_harmony'] = adata.obsm['X_umap'].copy()
sc.tl.leiden(adata, resolution=OVERCLUSTER_RES, key_added='leiden_overcluster',
             neighbors_key='PC_harmony')
adata.write(os.path.join(path_objects, 'PreQC.Leiden.h5ad'))
print(adata.obs['leiden_overcluster'].value_counts())

# %%
adata = sc.read_h5ad(os.path.join(path_objects, 'PreQC.Leiden.h5ad'))
adata.obsm['X_umap'] = adata.obsm['UMAP_PC_harmony'].copy()
fig, ax = plt.subplots(1, 2, figsize=(16, 7))
sc.pl.umap(adata, color='sample', ax=ax[0], show=False, frameon=False, title='sample (mixing check)')
sc.pl.umap(adata, color='leiden_overcluster', ax=ax[1], show=False, frameon=False, legend_loc='on data',
           legend_fontsize=6, title='leiden res %.1f' % OVERCLUSTER_RES)
fig.savefig(os.path.join(path_figs, 'step4_preqc_overcluster.svg'), bbox_inches='tight')

# %%
sc.pl.dotplot(adata,
              ['Ttn','Ryr2','Myh6',
               'Vwf','Pecam1','Emcn',
               'Dcn','Gsn','Gucy1a2',
               'Cd163','F13a1',
               'Pdgfrb','Rgs5',
               'Myh11',
               'Skap1','Il7r',
               'Cd79a','Ms4a1',
               'Mmrn1','Reln',
               'Plp1','Cdh19',
               'Adipoq'],
              groupby='leiden_overcluster', standard_scale='var')

# %%
adata = sc.read_h5ad(os.path.join(path_objects, 'PreQC.Leiden.h5ad'))
cluster_key = 'leiden_overcluster'

cluster_profiles = adata.obs.groupby(cluster_key, observed=True).agg(
    n=(cluster_key, 'size'),
    med_dbl=('cellranger_doublet_scores', 'median'),
    frac_pred_dbl=('cellranger_predicted_doublets', 'mean'),
    med_mito=('percent_mito', 'median'),
    med_mito_cellbender=('percent_mito_cellbender', 'median'),
    med_entropy=('cellbender_entropy', 'median'),
    med_genes=('cellbender_ngenes', 'median'),
    med_exon_prop=('exon_prop', 'median'),
)

global_entropy_floor = np.percentile(adata.obs['cellbender_entropy'], 5)
global_exon_prop_ceiling = np.percentile(adata.obs['exon_prop'], TECH_EXON_PROP_PCTL)
cluster_profiles['is_technical'] = ((cluster_profiles['frac_pred_dbl'] > 0.5) |
                        (cluster_profiles['med_mito'] > MAX_PCT_MITO) |
                        (cluster_profiles['med_genes'] < MIN_GENES) |
                        (cluster_profiles['med_entropy'] < global_entropy_floor) |
                        (cluster_profiles['med_exon_prop'] > global_exon_prop_ceiling))

FIXED_LIMITS = {metric: float(np.percentile(adata.obs[metric], FIXED_LIMIT_PCTL)) for metric in ['exon_prop', 'log_genes_x_entropy']}
print('exon_prop ceiling of technical clusters (proportion): %.3f; fixed limits (exon_prop proportion, log_genes_x_entropy score):' % global_exon_prop_ceiling,
      {k: round(v, 3) for k, v in FIXED_LIMITS.items()})
bad_clusters = cluster_profiles[cluster_profiles['is_technical']].index.tolist()
print(cluster_profiles.round(3).sort_values('frac_pred_dbl', ascending=False).rename(columns=DISPLAY_NAMES).to_string())
print('\nTechnical clusters flagged:', bad_clusters)
cluster_profiles.to_csv(os.path.join(path_tables, 'preqc_cluster_profiles.csv'))

# %%
def iqr_out(vector, high, fallback):

    q1, q3 = np.percentile(vector, [25, 75])
    use_iqr = len(vector) > 30 and q3 > q1
    cutoff = fallback
    if use_iqr:
        cutoff = q3 + IQR_MULT * (q3-q1) if high else q1 - IQR_MULT * (q3-q1)
    mask = vector > cutoff if high else vector < cutoff
    return mask, cutoff, use_iqr

# %%
def per_group_qc(group_adata):
    nucleus_qc = group_adata.obs.copy()
    nucleus_qc['iqr_outlier'] = False
    nucleus_qc['fixed_doublet_outlier'] = False
    nucleus_qc['fixed_exon_complexity_outlier'] = False
    limits = []

    metrics = [('cellbender_ngenes', False, -np.inf, 'flag_cellbender_ngenes'),
               ('cellbender_ngenes', True, np.inf, 'flag_cellbender_ngenes_high'),
               ('cellbender_ncount', False, -np.inf, 'flag_cellbender_ncount'),
               ('cellbender_ncount', True, np.inf, 'flag_cellbender_ncount_high'),
               ('cellbender_entropy', False, -np.inf, 'flag_cellbender_entropy'),
               ('percent_mito', True, MAX_PCT_MITO, 'flag_percent_mito'),
               ('cellranger_doublet_scores', True, DOUBLET_FALLBACK, 'flag_cellranger_doublet_scores'),
               ('exon_prop', True, FIXED_LIMITS['exon_prop'], 'flag_exon_prop'),
               ('log_genes_x_entropy', True, FIXED_LIMITS['log_genes_x_entropy'], 'flag_log_genes_x_entropy')]
    for metric, high, fallback, flag in metrics:
        mask, cutoff, use_iqr = iqr_out(nucleus_qc[metric].values, high, fallback)
        nucleus_qc[flag] = mask
        if use_iqr:
            nucleus_qc['iqr_outlier'] |= mask
        elif metric == 'cellranger_doublet_scores':
            nucleus_qc['fixed_doublet_outlier'] = mask
        elif metric in FIXED_LIMITS:
            nucleus_qc['fixed_exon_complexity_outlier'] |= mask
        limits.append({'sample': nucleus_qc['sample'].iloc[0], 'cluster': nucleus_qc[cluster_key].iloc[0],
                       'n': len(nucleus_qc), 'metric': metric, 'side': 'high' if high else 'low', 'cutoff': cutoff,
                       'rule': 'IQR' if use_iqr else 'fixed/no low-side fence'})
    nucleus_qc['hard_low_genes'] = nucleus_qc['cellbender_ngenes'] < MIN_GENES
    nucleus_qc['hard_low_umi'] = nucleus_qc['cellbender_ncount'] < MIN_UMI
    nucleus_qc['hard_high_mito'] = nucleus_qc['percent_mito'] > MAX_PCT_MITO
    nucleus_qc['hard_doublet'] = nucleus_qc['cellranger_predicted_doublets'].astype(bool)
    nucleus_qc['outlier'] = (nucleus_qc['iqr_outlier'] | nucleus_qc['fixed_doublet_outlier'] | nucleus_qc['fixed_exon_complexity_outlier'] | nucleus_qc['hard_low_genes'] |
                    nucleus_qc['hard_low_umi'] | nucleus_qc['hard_high_mito'] | nucleus_qc['hard_doublet'])
    return nucleus_qc, pd.DataFrame(limits)

# %%
group_qc_list, thresholds = [], []
for sample_id in adata.obs['sample'].unique():
    for cluster in adata.obs.loc[adata.obs['sample'] == sample_id, cluster_key].unique():
        group_adata = adata[(adata.obs['sample'] == sample_id) & (adata.obs[cluster_key] == cluster)]
        group_qc, group_thresholds = per_group_qc(group_adata)
        group_qc_list.append(group_qc)
        thresholds.append(group_thresholds)
qc_results = pd.concat(group_qc_list)
pd.concat(thresholds).to_csv(os.path.join(path_tables, 'sample_cluster_qc_thresholds.csv'), index=False)

qc_results['technical_cluster'] = qc_results[cluster_key].isin(bad_clusters)
qc_results['outlier'] = qc_results['outlier'] | qc_results['technical_cluster']
qc_results.to_csv(os.path.join(path_tables, 'cellQC_results.txt.gz'), sep='\t', compression='gzip')
print(qc_results['outlier'].value_counts())
print('retained: %d / %d (%.1f%%)' % ((~qc_results['outlier']).sum(), qc_results.shape[0],
                                       100*(~qc_results['outlier']).mean()))
print('\nScrublet-predicted doublets removed per library:')
print(qc_results.groupby('sample', observed=True)['hard_doublet'].agg(n_doublets='sum', proportion='mean').round(3).to_string())

# %%
for keys, filename in [(['disease'], 'qc_retention_by_condition.csv'),
                        (['sample'], 'qc_retention_by_sample.csv'),
                        (['sample', cluster_key], 'qc_retention_by_sample_cluster.csv')]:
    retention_tab = qc_results.groupby(keys, observed=True)['outlier'].agg(n_preqc='size', n_removed='sum')
    retention_tab['n_retained'] = retention_tab['n_preqc'] - retention_tab['n_removed']
    retention_tab['frac_removed'] = retention_tab['n_removed'] / retention_tab['n_preqc']
    retention_tab.rename(columns=EXPORT_NAMES).to_csv(os.path.join(path_tables, filename))
    if keys != ['sample', cluster_key]:
        print(retention_tab.rename(columns=DISPLAY_NAMES).round(3).to_string(), '\n')

# %%
qc_fig_dir = os.path.join(path_figs, 'qc_per_sample')
os.makedirs(qc_fig_dir, exist_ok=True)

umap_coords = pd.DataFrame(adata.obsm['UMAP_PC_harmony'], index=adata.obs_names, columns=['UMAP1', 'UMAP2'])
qc_results[['UMAP1', 'UMAP2']] = umap_coords.reindex(qc_results.index).values
metrics = ['cellbender_ncount', 'cellbender_ngenes', 'percent_mito', 'cellbender_entropy',
           'cellranger_doublet_scores', 'exon_prop', 'log_genes_x_entropy']

def plot_sample_qc(qc_res, sample, savefig):
    sample_qc = qc_res[qc_res['sample'] == sample]
    order = sorted(sample_qc[cluster_key].unique(), key=int)
    cluster_codes = pd.Categorical(sample_qc[cluster_key], categories=order).codes
    fig = plt.figure(figsize=(25, 8))
    grid = gridspec.GridSpec(2, 5, wspace=0.4, hspace=0.45)
    for pos, metric in zip([(0, 0), (0, 1), (1, 0), (1, 1), (0, 2), (0, 3), (1, 3)], metrics):
        ax = fig.add_subplot(grid[pos])
        ax.boxplot([sample_qc.loc[cluster_codes == k, metric].values for k in range(len(order))], positions=range(len(order)),
                   showfliers=False, widths=0.7, medianprops={'color': 'k'})
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=6)
        ax.set_title(metric); ax.set_xlabel('pre-QC cluster')
    ax = fig.add_subplot(grid[1, 2])
    frac_removed = sample_qc.groupby(cluster_codes)['outlier'].mean()
    ax.bar(frac_removed.index, frac_removed.values, color=REMOVE)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_ylim(0, 1); ax.set_ylabel('fraction removed (proportion)'); ax.set_xlabel('pre-QC cluster'); ax.set_title('removal per cluster')
    ax = fig.add_subplot(grid[0, 4])
    ax.scatter(sample_qc['UMAP1'], sample_qc['UMAP2'], s=2, c=np.where(sample_qc['outlier'], REMOVE, GREY), linewidths=0, rasterized=True)
    ax.set_title('removed nuclei (red) on the pre-QC UMAP'); ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    ax = fig.add_subplot(grid[1, 4]); ax.axis('off')
    n_removed, n_retained = int(sample_qc['outlier'].sum()), int((~sample_qc['outlier']).sum())
    ax.text(0, .7, 'Removed = %d' % n_removed, fontsize=16, color=REMOVE)
    ax.text(0, .55, 'Retained = %d' % n_retained, fontsize=16)
    ax.text(0, .4, '%.1f%% retained' % (100*n_retained/(n_retained+n_removed)), fontsize=13)
    fig.suptitle('SAMPLE: %s' % sample, fontsize=16)
    fig.savefig(savefig, bbox_inches='tight', dpi=150); plt.close(fig)

for sample_id in qc_results['sample'].unique():
    plot_sample_qc(qc_results, sample_id, os.path.join(qc_fig_dir, str(sample_id) + '_qc.png'))
print('wrote %d per-sample QC figures to %s' % (qc_results['sample'].nunique(), qc_fig_dir))

# %%
adata = sc.read_h5ad(os.path.join(path_objects, 'data.PreQC.Raw.h5ad'))
qc_table = pd.read_csv(os.path.join(path_tables, 'cellQC_results.txt.gz'), sep='\t', index_col=0)
keep = qc_table[~qc_table['outlier']].index
adata = adata[adata.obs.index.isin(keep)].copy()
print('clean nuclei:', adata.shape)

sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor='seurat_v3', batch_key='sample')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver='arpack')
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

# %%
adata.obs['batch_harmony'] = adata.obs['sample'].astype('category')
adata.obsm['PC_harmony'] = run_harmony(adata.obsm['X_pca'], adata.obs, 'batch_harmony')
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep='PC_harmony',
                metric='cosine', key_added='PC_harmony')
sc.tl.umap(adata, neighbors_key='PC_harmony', min_dist=0.2)
adata.obsm['UMAP_PC_harmony'] = adata.obsm['X_umap'].copy()
sc.tl.leiden(adata, resolution=MAP_RES, key_added='leiden', neighbors_key='PC_harmony')
adata.write(os.path.join(path_objects, 'PostQC.Map.h5ad'))
print(adata.obs['leiden'].value_counts())

# %%
adata = sc.read_h5ad(os.path.join(path_objects, 'PostQC.Map.h5ad'))
adata.obsm['X_umap'] = adata.obsm['UMAP_PC_harmony'].copy()

panels = {
    'Cardiomyocyte':      ['Ttn','Ryr2','Myh6','Tnnt2','Mlip'],
    'Endothelial':        ['Vwf','Pecam1','Emcn','Cdh5','Fabp4'],
    'Fibroblast':         ['Dcn','Gsn','Gucy1a2','Pdgfra','Col1a1'],
    'Pericyte':           ['Pdgfrb','Rgs5','Kcnj8','Abcc9'],
    'VSMC':               ['Myh11','Tagln','Acta2','Cnn1'],
    'Macrophage':         ['Cd163','F13a1','Mrc1','Csf1r'],
    'Lymphocyte':         ['Skap1','Il7r','Themis','Cd3e'],
    'B_cell':             ['Cd79a','Cd79b','Ms4a1','Cd19','Pax5','Ebf1'],
    'Lymphatic_EC':       ['Mmrn1','Reln','Prox1','Flt4'],
    'Glia':               ['Plp1','Cdh19','Scn7a','Kcna1','Sox10','Nrxn1'],
    'Adipocyte':          ['Adipoq','Gpam','Plin1'],
    'Proliferating':      ['Top2a','Mki67'],
}
panel_description = {'Glia': 'peripheral glia: non-myelinating Schwann-cell markers', 'Proliferating': 'cell-cycle state, not a lineage'}
panels = {k: [g for g in v if g in adata.raw.var_names] for k, v in panels.items()}
pd.DataFrame([{'cell_type': ct, 'genes': ';'.join(genes), 'n_genes': len(genes),
               'description': panel_description.get(ct, 'lineage')} for ct, genes in panels.items()]).to_csv(
    os.path.join(path_tables, 'annotation_panel_coverage.csv'), index=False)
for ct, genes in panels.items():
    sc.tl.score_genes(adata, genes, score_name='score_'+ct, use_raw=True)

# %%
score_cols = ['score_'+ct for ct in panels]
clmean = adata.obs.groupby('leiden', observed=True)[score_cols].mean()
ranked = np.sort(clmean.values, axis=1)
top_score = pd.Series(ranked[:, -1], index=clmean.index)
top_two_margin = pd.Series(ranked[:, -1] - ranked[:, -2], index=clmean.index)
cluster_label = clmean.idxmax(axis=1).str.replace('score_', '', regex=False)
cluster_label[(top_score < MIN_TOP_SCORE) | (top_two_margin < MIN_MARGIN)] = 'Unassigned'
dbl_frac = adata.obs.groupby('leiden', observed=True)['cellranger_predicted_doublets'].mean()
cluster_label[dbl_frac > 0.5] = 'Doublet'
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_label).astype('category')

cluster_scores = adata.obs.groupby('leiden', observed=True).agg(
    n=('sample', 'size'), n_samples=('sample', 'nunique'),
    n_individuals=('individual', 'nunique'),
    median_genes=('cellbender_ngenes', 'median'), median_umi=('cellbender_ncount', 'median'),
    median_mito=('percent_mito', 'median'), median_entropy=('cellbender_entropy', 'median'))
support = pd.crosstab(adata.obs['leiden'], adata.obs['sample'])
support.to_csv(os.path.join(path_tables, 'postqc_cluster_sample_counts.csv'))
cluster_scores['largest_sample'] = support.idxmax(axis=1)
cluster_scores['largest_sample_fraction'] = support.max(axis=1) / support.sum(axis=1)
cluster_scores['frac_pred_doublet'] = dbl_frac
cluster_scores['cell_type'] = cluster_label
cluster_scores['top_score'] = top_score
cluster_scores['top_two_margin'] = top_two_margin
cluster_scores['candidate_noncycling_lineage'] = clmean.drop(
    columns='score_Proliferating').idxmax(axis=1).str.replace('score_', '', regex=False)
cluster_scores = cluster_scores.join(clmean)
cluster_scores.to_csv(os.path.join(path_tables, 'cluster_annotation_before_exclusions.csv'))
print(cluster_scores.drop(columns=score_cols).rename(columns=DISPLAY_NAMES).round(3).to_string())
adata.uns['cell_type_colors'] = [CELL_TYPE_PALETTE[c] for c in adata.obs['cell_type'].cat.categories]
adata.write(os.path.join(path_objects, 'atlas_before_final_exclusions.h5ad'))
sc.pl.umap(adata, color=['leiden', 'cell_type', 'sample', 'cellbender_ngenes'], frameon=False,
           legend_loc='on data', ncols=2)
sc.pl.dotplot(adata, panels, groupby='leiden', standard_scale='var', use_raw=True)

# %%
lineage_scores = clmean.drop(columns='score_Proliferating')
second_lineage = lineage_scores.apply(lambda row: row.drop(row.idxmax()).max(), axis=1)
label_median_genes = adata.obs.groupby('cell_type', observed=True)['cellbender_ngenes'].median()
exclusion_rules = pd.DataFrame({
    'n': cluster_scores['n'], 'cell_type': cluster_label, 'median_genes': cluster_scores['median_genes'],
    'label_median_genes': cluster_label.map(label_median_genes).astype(float),
    'median_mito': cluster_scores['median_mito'], 'largest_sample': cluster_scores['largest_sample'],
    'largest_sample_fraction': cluster_scores['largest_sample_fraction'],
    'doublet_dominated': dbl_frac > 0.5,
    'high_mito': cluster_scores['median_mito'] > MAX_PCT_MITO,
    'low_complexity': cluster_scores['median_genes'] < LOWC_REL * cluster_label.map(label_median_genes).astype(float),
    'mixed_lineage': second_lineage > INCOH_REL * lineage_scores.max(axis=1)})
exclusion_rules['exclude'] = exclusion_rules['doublet_dominated'] | exclusion_rules['high_mito'] | (exclusion_rules['low_complexity'] & exclusion_rules['mixed_lineage'])
exclusion_rules['reason'] = np.select([exclusion_rules['doublet_dominated'], exclusion_rules['high_mito'], exclusion_rules['exclude']],
                            ['doublet-dominated cluster', 'high-mito cluster', 'low-complexity mixed-lineage cluster'], 'retained')
exclusion_rules.rename(columns=EXPORT_NAMES).to_csv(os.path.join(path_tables, 'cluster_exclusion_rules.csv'))
print(exclusion_rules.rename(columns=DISPLAY_NAMES).round(3).to_string())

reason = adata.obs['leiden'].map(exclusion_rules['reason']).astype(str)

inclusion = pd.read_csv(os.path.join(path_tables, 'cellQC_results.txt.gz'), sep='\t', index_col=0)
inclusion['initial_qc_pass'] = ~inclusion['outlier']
inclusion['final_exclusion_reason'] = reason.reindex(inclusion.index).fillna('failed initial QC')
inclusion['cell_type'] = adata.obs['cell_type'].reindex(inclusion.index).astype(object).fillna('')
inclusion['final_atlas'] = inclusion['final_exclusion_reason'] == 'retained'
inclusion.to_csv(os.path.join(path_tables, 'atlas_inclusion_by_barcode.csv.gz'))

print('\nFinal-stage exclusions:', int((reason != 'retained').sum()), 'nuclei in clusters',
      exclusion_rules.index[exclusion_rules['exclude']].tolist())
adata = adata[(reason == 'retained').values].copy()
adata.obs['cell_type'] = adata.obs['cell_type'].cat.remove_unused_categories()
adata.obs['cell_type'] = adata.obs['cell_type'].cat.reorder_categories(
    [c for c in CELL_TYPE_ORDER + ['Adipocyte', 'Unassigned'] if c in adata.obs['cell_type'].cat.categories])
adata.uns['cell_type_colors'] = [CELL_TYPE_PALETTE[c] for c in adata.obs['cell_type'].cat.categories]
cluster_scores.loc[cluster_scores.index.isin(adata.obs['leiden'])].rename(columns=EXPORT_NAMES).to_csv(os.path.join(path_tables, 'cluster_annotation.csv'))
adata.write(os.path.join(path_objects, 'atlas_annotated.h5ad'))
print(adata.obs['cell_type'].value_counts().to_string())

# %%
sc.tl.rank_genes_groups(adata, 'cell_type', method='wilcoxon', tie_correct=False, n_genes=adata.raw.shape[1], use_raw=True)
marker_rows = []
names = adata.uns['rank_genes_groups']['names'].dtype.names
counts = adata.obs['cell_type'].value_counts()
for group_pos, cl in enumerate(names):
    group_markers = {k: [field_values[group_pos] for field_values in adata.uns['rank_genes_groups'][k].tolist()]
         for k in adata.uns['rank_genes_groups'] if k != 'params'}
    marker_df = pd.DataFrame(group_markers); marker_df['cell_type'] = cl
    n_in_group = counts[cl]; n_out_group = adata.shape[0] - n_in_group
    marker_df['U'] = marker_df['scores'] * np.sqrt(n_in_group*n_out_group*(n_in_group+n_out_group+1)/12) + n_in_group*n_out_group/2
    marker_df['AUC'] = marker_df['U'] / (n_in_group*n_out_group)
    marker_rows.append(marker_df)
markers = pd.concat(marker_rows)
markers.to_csv(os.path.join(path_tables, 'broad_cluster_markers.txt.gz'), sep='\t',
               index=False, compression='gzip')
markers[markers['AUC'] > 0.6].groupby('cell_type').head(5)[['cell_type','names','AUC','logfoldchanges']]

# %%
ct_counts = adata.obs['cell_type'].value_counts()
fig, ax = plt.subplots(1, 1, figsize=(9, 7))
sc.pl.umap(adata, color='cell_type', ax=ax, show=False, frameon=False, palette=CELL_TYPE_PALETTE, legend_loc='none',
           title='Broad cell types (n=%d nuclei)' % adata.shape[0])
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=CELL_TYPE_PALETTE[c], label=f'{c} (n={ct_counts[c]:,})')
                   for c in adata.obs['cell_type'].cat.categories],
          frameon=False, loc='center left', bbox_to_anchor=(1.01, 0.5))
fig.savefig(os.path.join(path_figs, 'step7_annotated_atlas.svg'), bbox_inches='tight')
fig.savefig(os.path.join(path_figs, 'step7_annotated_atlas.png'), bbox_inches='tight', dpi=300)

top_genes = [markers[markers['cell_type']==ct].sort_values('AUC', ascending=False).head(3)['names'].tolist()
        for ct in adata.obs['cell_type'].cat.categories]
top_genes = [g for gene_list in top_genes for g in gene_list]
sc.pl.dotplot(adata, top_genes, groupby='cell_type', standard_scale='var')

# %%
adata = sc.read_h5ad(os.path.join(path_objects, 'atlas_annotated.h5ad'))

fig = plt.figure(figsize=(11, 14))
grid = fig.add_gridspec(2, 1, height_ratios=[1.4, 1], hspace=0.3)
ax_umap = fig.add_subplot(grid[0])
ax_dot  = fig.add_subplot(grid[1])

sc.pl.umap(adata, color='cell_type', ax=ax_umap, show=False, frameon=False, legend_loc='on data', palette=CELL_TYPE_PALETTE,
           legend_fontsize=7, title='Broad cell types (n=%d nuclei)' % adata.shape[0])
sc.pl.dotplot(adata, top_genes, groupby='cell_type', standard_scale='var',
              ax=ax_dot, show=False)

fig.savefig(os.path.join(path_figs, 'step7_annotated_atlas_combined.svg'), bbox_inches='tight')

# %%
dotplot_fig = sc.pl.dotplot(adata, top_genes, groupby='cell_type', standard_scale='var', return_fig=True)
dotplot_fig.savefig(os.path.join(path_figs, 'step7_annotated_atlas_refined_markers.svg'), bbox_inches='tight')

# %%
sample_counts = pd.crosstab(adata.obs['sample'], adata.obs['cell_type'])
sample_counts.to_csv(os.path.join(path_tables, 'final_celltype_counts_by_sample.csv'))
animal_counts = pd.crosstab(adata.obs['individual'], adata.obs['cell_type'])
prop_individual = animal_counts.div(animal_counts.sum(axis=1), axis=0)
animal_meta = adata.obs[['individual', 'disease', 'qc_flag']].drop_duplicates().set_index('individual')
assert animal_meta.index.is_unique, 'An individual has conflicting condition labels'
prop_individual = prop_individual.join(animal_meta)
prop_individual['n_nuclei'] = animal_counts.sum(axis=1)
prop_individual.rename(columns={ct: ct + '_prop' for ct in animal_counts.columns}).to_csv(os.path.join(path_tables, 'celltype_proportions_by_individual.csv'))
mean_prop_by_condition = prop_individual.groupby('disease', observed=True)[animal_counts.columns].mean()
mean_prop_by_condition.add_suffix('_prop').to_csv(os.path.join(path_tables, 'celltype_proportions_by_condition.csv'))
pooled = pd.crosstab(adata.obs['disease'], adata.obs['cell_type'], normalize='index')
pooled.add_suffix('_prop').to_csv(os.path.join(path_tables, 'celltype_proportions_pooled_nuclei.csv'))
adata.obs[['sample','individual','disease','sex','age','condition','batch','seq','qc_flag',
           'cellbender_ncount','cellbender_ngenes','percent_mito','percent_mito_cellbender','exon_prop',
           'cellranger_doublet_scores','leiden','cell_type']].to_csv(os.path.join(path_tables, 'cell_metadata.csv'))
print(prop_individual[['disease', 'qc_flag', 'n_nuclei']].to_string())
mean_prop_by_condition.rename_axis(columns='cell type (proportion, mean of animals)')

# %%
atlas_adata = sc.read_h5ad(os.path.join(path_objects, 'atlas_annotated.h5ad'))
atlas_adata.obsm['X_umap'] = atlas_adata.obsm['UMAP_PC_harmony']
order = list(atlas_adata.obs['cell_type'].cat.categories)
n_samples = atlas_adata.obs['sample'].nunique()

rng = np.random.RandomState(0)
lisi_idx = rng.choice(atlas_adata.n_obs, min(30000, atlas_adata.n_obs), replace=False)
lisi_meta = atlas_adata.obs.iloc[lisi_idx][['sample']].copy()
lisi_pre = harmonypy.compute_lisi(atlas_adata.obsm['X_pca'][lisi_idx, :], lisi_meta, ['sample'])[:, 0]
lisi_post = harmonypy.compute_lisi(atlas_adata.obsm['PC_harmony'][lisi_idx, :], lisi_meta, ['sample'])[:, 0]
lisi = atlas_adata.obs.iloc[lisi_idx][['sample', 'cell_type']].copy()
lisi['pre_harmony'] = lisi_pre
lisi['post_harmony'] = lisi_post
lisi.to_csv(os.path.join(path_tables, 'ilisi_per_nucleus.csv.gz'))
lisi.groupby('cell_type', observed=True)[['pre_harmony', 'post_harmony']].median().to_csv(
    os.path.join(path_tables, 'ilisi_by_celltype.csv'))
lisi_by_sample = lisi.groupby('sample', observed=True)[['pre_harmony', 'post_harmony']].median()
lisi_by_sample.to_csv(os.path.join(path_tables, 'ilisi_by_sample.csv'))

draw_order = rng.permutation(atlas_adata.n_obs)
fig = plt.figure(figsize=(16, 4.6)); grid = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.9], wspace=0.3)
for j, key in enumerate(['sample', 'batch', 'seq']):
    ax = fig.add_subplot(grid[j]); categories = atlas_adata.obs[key].astype('category')
    cmap = plt.get_cmap('tab20' if categories.cat.categories.size > 10 else 'tab10')
    colour_map = {c: mpl.colors.to_hex(cmap(k % cmap.N)) for k, c in enumerate(categories.cat.categories)}
    ax.scatter(atlas_adata.obsm['X_umap'][draw_order, 0], atlas_adata.obsm['X_umap'][draw_order, 1], s=1, c=categories.map(colour_map).astype(str).values[draw_order], linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f'UMAP coloured by {key}')
    ax.spines['left'].set_visible(False); ax.spines['bottom'].set_visible(False)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=colour_map[c], label=c) for c in categories.cat.categories],
              frameon=False, fontsize=7, ncol=2 if len(colour_map) > 8 else 1, loc='upper left', bbox_to_anchor=(0, 0), markerscale=1.2)
ax = fig.add_subplot(grid[3])
boxes = ax.boxplot([lisi_pre, lisi_post], tick_labels=['pre-Harmony\n(PCA)', 'post-Harmony'], showfliers=False, patch_artist=True, widths=0.6)
for patch, c in zip(boxes['boxes'], [GREY, KEEP]): patch.set_facecolor(c)
for median_line in boxes['medians']: median_line.set_color('k')
ax.axhline(n_samples, ls=':', color='grey', lw=0.8); ax.text(2.42, n_samples, f' {n_samples} samples', fontsize=7, va='center')
ax.set_ylabel('iLISI (effective number of samples per neighbourhood)'); ax.set_title('Sample mixing before / after')
ax.text(1, np.median(lisi_pre), f" {np.median(lisi_pre):.1f}", fontsize=7, va='bottom')
ax.text(2, np.median(lisi_post), f" {np.median(lisi_post):.1f}", fontsize=7, va='bottom')
fig.suptitle(f'Harmony integration QC: {atlas_adata.n_obs:,} nuclei, {n_samples} samples - mixing across sample / batch / run', y=1.03, fontsize=9)
save_svg(fig, 'FigS6_integration_mixing'); plt.close(fig)
print(lisi_by_sample.round(2).to_string())

# %%
panel = {'Cardiomyocyte': ['Mhrt', 'Atp2a2', 'Nppa'], 'Endothelial': ['Flt1', 'Egfl7', 'Npr3'],
         'Fibroblast': ['Abca8a', 'Col3a1', 'Lama2'], 'Pericyte': ['Notch3', 'Vtn', 'Higd1b'],
         'VSMC': ['Flna', 'Lmod1', 'Mylk'], 'Macrophage': ['Ptprc', 'Lyz2', 'Adgre1'],
         'Lymphocyte': ['Cd247', 'Bcl11b', 'Prkcq'], 'B_cell': ['Ighm', 'Bank1', 'Fcmr'],
         'Lymphatic_EC': ['Ccl21a', 'Lyve1', 'Pdpn'], 'Glia': ['Sox10', 'Adam23', 'Gfra3'],
         'Proliferating': ['Hjurp', 'Knl1', 'Cenpf']}
panel = {k: [g for g in v if g in atlas_adata.raw.var_names] for k, v in panel.items() if k in order}
fig = sc.pl.dotplot(atlas_adata, panel, groupby='cell_type', standard_scale='var', use_raw=True, return_fig=True)
fig.savefig(os.path.join(qc_dir, 'FigS7_marker_dotplot.svg'), bbox_inches='tight')

# %%
composition = pd.crosstab(atlas_adata.obs['individual'], atlas_adata.obs['cell_type'], normalize='index')
animal_meta = atlas_adata.obs[['individual', 'disease', 'qc_flag']].drop_duplicates().set_index('individual')
composition = composition.join(animal_meta); cell_types = [c for c in order if c in composition.columns]
fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={'width_ratios': [1, 2.2]})
group_mean = composition.groupby('disease', observed=True)[cell_types].mean().reindex(GROUP_ORDER); bottom = np.zeros(len(group_mean))
for ct in cell_types:
    axs[0].bar(range(len(group_mean)), group_mean[ct], bottom=bottom, color=CELL_TYPE_PALETTE[ct], label=ct, width=0.8); bottom += group_mean[ct].values
axs[0].set_xticks(range(len(group_mean))); axs[0].set_xticklabels(group_mean.index); axs[0].set_ylabel('mean proportion of nuclei per animal'); axs[0].set_ylim(0, 1)
axs[0].set_xlabel('group'); axs[0].set_title('Composition by group (mean of per-animal proportions)')
axs[0].legend(frameon=False, ncol=1, bbox_to_anchor=(-0.45, 1), loc='upper left')
x = np.arange(len(cell_types)); bar_width = 0.19
for group_idx, group in enumerate(GROUP_ORDER):
    group_comp = composition[composition['disease'] == group]
    for type_idx, ct in enumerate(cell_types):
        x_pos = x[type_idx] + (group_idx - 1.5) * bar_width
        axs[1].scatter(np.repeat(x_pos, len(group_comp)) + rng.uniform(-0.05, 0.05, len(group_comp)), group_comp[ct], s=14, c=GROUP_FILL[group],
                       edgecolor=np.where(group_comp['qc_flag'] == 'pass', 'none', 'k'), linewidths=0.6, zorder=3)
        axs[1].plot([x_pos - bar_width * 0.4, x_pos + bar_width * 0.4], [group_comp[ct].median()] * 2, color='k', lw=1, zorder=4)
axs[1].set_yscale('log'); axs[1].set_ylim(composition[cell_types].values[composition[cell_types].values > 0].min() * 0.7, 1)
axs[1].set_xticks(x); axs[1].set_xticklabels(cell_types, rotation=40, ha='right')
axs[1].set_ylabel('proportion per animal (log scale)'); axs[1].set_title('Per-animal cell-type proportions by group (black outline = qc_flag != pass)')
axs[1].legend(handles=[Line2D([], [], marker='o', ls='', color=GROUP_FILL[group], label=group) for group in GROUP_ORDER], frameon=False, ncol=4, loc='lower left')
fig.suptitle('Cell-type composition across age x HFpEF (descriptive; each animal weighted equally)', y=1.02, fontsize=9); fig.tight_layout()
save_svg(fig, 'FigS8_composition_by_condition'); plt.close(fig)

# %%
qc_metric_specs = [('cellbender_ncount', 'log10 UMIs per nucleus (CellBender)', True), ('cellbender_ngenes', 'log10 genes per nucleus (CellBender)', True),
       ('percent_mito', '% mitochondrial UMIs (raw counts)', False), ('cellranger_doublet_scores', 'Scrublet doublet score (raw counts)', False),
       ('exon_prop', 'exonic read proportion (cytoplasmic carry-over)', False)]
fig, axs = plt.subplots(1, 5, figsize=(21, 4.2))
for ax, (metric, title_text, log_scale) in zip(axs, qc_metric_specs):
    v = atlas_adata.obs[metric].values.astype(float)
    if log_scale: v = np.log10(np.clip(v, 1, None))
    draw_order = np.argsort(v)
    scatter = ax.scatter(atlas_adata.obsm['X_umap'][draw_order, 0], atlas_adata.obsm['X_umap'][draw_order, 1], c=v[draw_order], s=1.2, cmap='viridis', linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title_text)
    ax.spines['left'].set_visible(False); ax.spines['bottom'].set_visible(False)
    colorbar = fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.02); colorbar.ax.tick_params(labelsize=7); colorbar.set_label(title_text.split(' (')[0], fontsize=7)
fig.suptitle('Final atlas: QC metrics on the UMAP', y=1.03, fontsize=9)
fig.tight_layout(); save_svg(fig, 'FigS9_qc_on_umap'); plt.close(fig)

# %%
genes_amb = [g for g in ['Ttn', 'Myh6', 'Malat1', 'Dcn', 'Pecam1', 'Cd163'] if g in atlas_adata.var_names]
gene_idx = [atlas_adata.var_names.get_loc(g) for g in genes_amb]
cellranger_counts = sp.csc_matrix(atlas_adata.layers['cellranger_raw']); cellbender_counts = sp.csc_matrix(atlas_adata.layers['counts'])
detection_frac = lambda mat, j: float((mat[:, j] > 0).sum()) / mat.shape[0]
ambient_detection = pd.DataFrame({'gene': genes_amb, 'CellRanger': [100 * detection_frac(cellranger_counts, j) for j in gene_idx], 'CellBender': [100 * detection_frac(cellbender_counts, j) for j in gene_idx]})
fig, ax = plt.subplots(figsize=(6.5, 3.6)); x = np.arange(len(ambient_detection)); bar_width = 0.38
ax.bar(x - bar_width / 2, ambient_detection['CellRanger'], bar_width, color=GREY, label='Raw (STARsolo)')
ax.bar(x + bar_width / 2, ambient_detection['CellBender'], bar_width, color=KEEP, label='CellBender (adjusted counts)')
ax.set_xticks(x); ax.set_xticklabels(ambient_detection['gene'], style='italic'); ax.set_xlabel('gene')
ax.set_ylabel('% nuclei with detected counts'); ax.set_title('Gene detection before / after CellBender'); ax.legend(frameon=False)
fig.tight_layout(); save_svg(fig, 'FigS10_cellbender_ambient'); plt.close(fig)

# %%
y_genes = [g for g in ['Ddx3y', 'Uty', 'Eif2s3y', 'Kdm5d'] if g in atlas_adata.raw.var_names]
raw_adata = atlas_adata.raw.to_adata()
xist_expr = np.asarray(raw_adata[:, 'Xist'].X.todense()).ravel()
y_expr = np.asarray(raw_adata[:, y_genes].X.todense()).mean(1).ravel()
sex_check = pd.DataFrame({'sample': atlas_adata.obs['sample'].values, 'sex': atlas_adata.obs['sex'].values, 'Xist': xist_expr, 'Y': y_expr})
sex_check_mean = sex_check.groupby('sample', observed=True).agg(Xist=('Xist', 'mean'), Y=('Y', 'mean'), sex=('sex', 'first')).reset_index()
sex_check_mean.to_csv(os.path.join(path_tables, 'sex_check_by_sample.csv'), index=False)
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, (sex_code, c, title_text) in zip(axs, [('M', KEEP, 'recorded Male'), ('F', REMOVE, 'recorded Female')]):
    sex_subset = sex_check_mean[sex_check_mean['sex'] == sex_code].sort_values('Y'); y_pos = np.arange(len(sex_subset))
    ax.scatter(sex_subset['Y'], y_pos, c=c, s=45, zorder=3, label='Y-gene mean'); ax.scatter(sex_subset['Xist'], y_pos, c='k', marker='x', s=30, zorder=3, label='Xist')
    ax.set_yticks(y_pos); ax.set_yticklabels(sex_subset['sample']); ax.set_xlabel('mean log-normalised expression per nucleus'); ax.set_title(title_text, color=c)
    ax.legend(frameon=False, loc='lower right')
fig.suptitle('Sex check: Y-gene (' + '/'.join(y_genes) + ') vs Xist per sample', y=1.02, fontsize=9); fig.tight_layout()
save_svg(fig, 'FigS11_sex_check'); plt.close(fig)

# %%
thousands_fmt = mpl.ticker.FuncFormatter(lambda x, _: f'{x/1000:g}k' if x >= 1000 else f'{x:g}')
qc_table = pd.read_csv(os.path.join(path_tables, 'cellQC_results.txt.gz'), sep='\t', index_col=0)
stages = {'CellBender nuclei\n(pre-QC)': len(qc_table), 'Passed initial QC': int((~qc_table['outlier']).sum()), 'Final atlas': atlas_adata.n_obs}
fig, ax = plt.subplots(figsize=(4.4, 3.6)); x_pos = range(len(stages)); stage_counts = list(stages.values())
ax.bar(x_pos, stage_counts, color=[GREY, '#80B1D3', KEEP], width=0.6)
for i, v in enumerate(stage_counts): ax.text(i, v + stage_counts[0] * 0.02, f'{v:,}\n({100*v/stage_counts[0]:.0f}%)', ha='center', fontsize=7)
ax.set_xticks(list(x_pos)); ax.set_xticklabels(stages.keys()); ax.set_ylabel('nuclei'); ax.set_ylim(0, stage_counts[0] * 1.15); ax.yaxis.set_major_formatter(thousands_fmt)
ax.set_title('Nuclei recovered through QC'); fig.tight_layout(); save_svg(fig, 'FigS12_cell_funnel'); plt.close(fig)
print('wrote atlas supplementary SVGs to', qc_dir)

# %%
import h5py
qc_table = pd.read_csv(os.path.join(path_tables, 'cellQC_results.txt.gz'), sep='\t', index_col=0)
qc_table['barcode'] = [i.rsplit('-', 1)[0] for i in qc_table.index]
retention_by_sample = qc_table.groupby('sample').agg(n_cellbender=('outlier', 'size'), n_removed=('outlier', 'sum'), qc_flag=('qc_flag', 'first'))
retention_by_sample['n_retained'] = retention_by_sample['n_cellbender'] - retention_by_sample['n_removed']
retention_by_sample['frac_removed'] = retention_by_sample['n_removed'] / retention_by_sample['n_cellbender']

FLAGGED = retention_by_sample[retention_by_sample['frac_removed'] > 0.5].index.tolist()
print('flagged (>50% removed):', FLAGGED)
sample_colour = {sample_id: (REMOVE if sample_id in FLAGGED else KEEP) for sample_id in sample_list}
def color_ticklabels(ax):
    for t in ax.get_yticklabels() + ax.get_xticklabels():
        if t.get_text() in sample_colour: t.set_color(sample_colour[t.get_text()])

library_rows = []
for sample_id in sample_list:
    cb_metrics = pd.read_csv(os.path.join(path_samples, sample_id, 'cellbender_metrics.csv'), header=None, index_col=0)[1]
    with h5py.File(os.path.join(path_samples, sample_id, 'filtered_feature_bc_matrix.h5')) as h5:
        n_cellranger_cells = h5['matrix/barcodes'].shape[0]
    library_rows.append({'sample': sample_id, 'n_cellranger': n_cellranger_cells, 'cb_found_cells': int(cb_metrics['found_cells']),
                 'cb_frac_droplets_nonempty': cb_metrics['fraction_of_analyzed_droplets_that_are_nonempty'],
                 'cb_frac_counts_removed': cb_metrics['fraction_counts_removed'],
                 'cb_ratio_found_expected': cb_metrics['ratio_of_found_cells_to_expected_cells'],
                 'total_raw_counts': cb_metrics['total_raw_counts']})
library_metrics = pd.DataFrame(library_rows).set_index('sample').join(retention_by_sample)
library_metrics['cellbender_over_cellranger'] = library_metrics['n_cellbender'] / library_metrics['n_cellranger']
library_metrics['retained_over_cellranger'] = library_metrics['n_retained'] / library_metrics['n_cellranger']
library_metrics.rename(columns=EXPORT_NAMES).to_csv(os.path.join(path_tables, 'qc_validation_library_metrics.csv'))
print(library_metrics.round(3).sort_values('frac_removed', ascending=False).rename(columns=DISPLAY_NAMES).to_string())

sample_order = library_metrics.sort_values('frac_removed').index; y = np.arange(len(sample_order)); bar_height = 0.27
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.8), gridspec_kw={'width_ratios': [2.1, 1.2, 0.9]})
panel_ax = axs[0]
panel_ax.barh(y + bar_height, library_metrics.loc[sample_order, 'n_cellbender'], bar_height, color=GREY, label='CellBender non-empty droplets (QC input)')
panel_ax.barh(y,     library_metrics.loc[sample_order, 'n_cellranger'], bar_height, color='#7A7F87', label='STARsolo cell calls')
panel_ax.barh(y - bar_height, library_metrics.loc[sample_order, 'n_retained'],   bar_height, color=KEEP, label='retained after QC')
for yi, sample_id in zip(y, sample_order):
    panel_ax.text(library_metrics.loc[sample_id, 'n_cellbender'] + 250, yi + bar_height, f"{library_metrics.loc[sample_id,'cellbender_over_cellranger']:.1f}x STARsolo", va='center', fontsize=6, color='#555')
    panel_ax.text(library_metrics.loc[sample_id, 'n_retained'] + 250, yi - bar_height, f"{library_metrics.loc[sample_id,'retained_over_cellranger']:.2f}x CR", va='center', fontsize=6, color=KEEP)
panel_ax.set_yticks(y); panel_ax.set_yticklabels(sample_order); color_ticklabels(panel_ax)
panel_ax.set_xlabel('barcodes'); panel_ax.xaxis.set_major_formatter(thousands_fmt); panel_ax.legend(frameon=False, loc='lower right')
panel_ax.set_title('Barcodes per sample: CellBender input vs STARsolo call vs QC-retained')
panel_ax = axs[1]
panel_ax.scatter(library_metrics.loc[sample_order, 'cb_frac_droplets_nonempty'], y, c=[sample_colour[sample_id] for sample_id in sample_order], s=30, zorder=3,
          label='fraction of analysed droplets\ncalled non-empty by CellBender')
panel_ax.scatter(library_metrics.loc[sample_order, 'cb_frac_counts_removed'], y, marker='x', c='k', s=22, zorder=3,
          label='fraction of counts removed\nas ambient by CellBender')
for yi in y: panel_ax.axhline(yi, color=GREY, lw=0.4, zorder=1)
panel_ax.set_yticks(y); panel_ax.set_yticklabels([]); panel_ax.set_xlim(0, 1.02); panel_ax.set_xlabel('proportion')
panel_ax.set_title('CellBender run metrics'); panel_ax.legend(frameon=False, loc='center right')
panel_ax = axs[2]
panel_ax.barh(y, library_metrics.loc[sample_order, 'total_raw_counts'] / 1e6, color=[sample_colour[sample_id] for sample_id in sample_order])
panel_ax.set_yticks(y); panel_ax.set_yticklabels([]); panel_ax.set_xlabel('raw UMIs in called droplets (millions)'); panel_ax.set_title('Library size')
fig.suptitle('Droplet calling per sample (red labels = samples with >50% of CellBender droplets removed by QC)', y=1.02, fontsize=9)
fig.tight_layout(); save_svg(fig, 'FigS13_droplet_calling_per_sample'); plt.close(fig)

# %%
def barcode_umis(path):
    with h5py.File(path) as h5:
        indptr = h5['matrix/indptr'][:]; umi_data = h5['matrix/data'][:]; barcodes = h5['matrix/barcodes'][:].astype(str)
    barcode_totals = np.zeros(len(barcodes)); nonzero = np.diff(indptr) > 0
    barcode_totals[nonzero] = np.add.reduceat(umi_data, indptr[:-1][nonzero])
    return pd.Series(barcode_totals, index=barcodes)

rng = np.random.RandomState(0); umi_at_window = {}
fig, axs = plt.subplots(4, 4, figsize=(13, 12), sharex=True, sharey=True); axs = axs.ravel()
for ax, sample_id in zip(axs, sample_order):
    barcode_totals = barcode_umis(os.path.join(path_samples, sample_id, 'raw_feature_bc_matrix.h5'))
    barcode_totals = barcode_totals[barcode_totals > 0].sort_values(ascending=False); barcode_rank = np.arange(1, len(barcode_totals) + 1)
    umi_at_window[sample_id] = int(barcode_totals.iloc[CELLBENDER_DROPLETS - 1])
    called_barcodes = set(pd.read_csv(os.path.join(path_samples, sample_id, 'cellbender_cell_barcodes.csv'), header=None)[0])
    retained_barcodes = set(qc_table.loc[(qc_table['sample'] == sample_id) & ~qc_table['outlier'], 'barcode'])
    barcode_fate = np.where(barcode_totals.index.isin(retained_barcodes), 'retained', np.where(barcode_totals.index.isin(called_barcodes), 'called, removed by QC', 'not called'))
    for fate_label, colour, point_size in [('not called', GREY, 1.5), ('called, removed by QC', REMOVE, 2.5), ('retained', KEEP, 2.5)]:
        fate_idx = np.where(barcode_fate == fate_label)[0]; n = len(fate_idx)
        if fate_label == 'not called' and n > 20000: fate_idx = np.sort(rng.choice(fate_idx, 20000, replace=False))
        ax.scatter(barcode_rank[fate_idx], barcode_totals.values[fate_idx], s=point_size, c=colour, linewidths=0, rasterized=True, label=f'{fate_label} (n={n:,})')
    ax.axvline(library_metrics.loc[sample_id, 'n_cellranger'], color='k', ls='--', lw=0.7)
    ax.text(library_metrics.loc[sample_id, 'n_cellranger'] * 1.15, barcode_totals.max() * 0.4, f"STARsolo calls\nn={library_metrics.loc[sample_id,'n_cellranger']:,}", fontsize=6)
    ax.axhline(MIN_UMI, color='k', ls=':', lw=0.6); ax.axvline(CELLBENDER_DROPLETS, color='k', ls='-.', lw=0.6)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.grid(True, which='major')
    ax.set_title(f"{sample_id}   ({100*(1-library_metrics.loc[sample_id,'frac_removed']):.0f}% retained)", color=sample_colour[sample_id])
    ax.legend(frameon=False, fontsize=6, markerscale=3, loc='lower left', handletextpad=0.2)
for ax in axs[12:]: ax.set_xlabel('barcode rank')
for ax in axs[::4]: ax.set_ylabel('raw UMIs per barcode')
fig.suptitle('Barcode-rank curves: fate of every droplet (dashed = STARsolo cell count, dash-dot = end of the CellBender window, dotted = UMI floor)', y=1.005, fontsize=9)
fig.tight_layout(); save_svg(fig, 'FigS14_barcode_rank_curves'); plt.close(fig)

# %%
cellranger_metrics = pd.DataFrame({sample_id: pd.read_csv(os.path.join(path_samples, sample_id, 'cellranger', 'metrics_summary.csv'), thousands=',').iloc[0]
                                   for sample_id in sample_list}).T
percent = lambda column: cellranger_metrics[column].str.rstrip('%').astype(float)
sex_by_sample = sex_check_mean.set_index('sample')
sample_qc = pd.DataFrame({
    'estimated_cells': cellranger_metrics['Estimated Number of Cells'].astype(int),
    'reads': cellranger_metrics['Number of Reads'].astype(int),
    'median_genes_per_cell': cellranger_metrics['Median Genes per Cell'].astype(int),
    'reads_in_cells_pct': percent('Fraction Reads in Cells'),
    'intronic_reads_pct': percent('Reads Mapped Confidently to Intronic Regions'),
    'umi_last_cellbender_droplet': pd.Series(umi_at_window),
    'recorded_sex': sex_by_sample['sex'],
    'inferred_sex': sex_by_sample['Xist'].gt(sex_by_sample['Y']).map({True: 'F', False: 'M'})}).reindex(sample_list)
sample_qc['empty_droplet_flag'] = sample_qc['umi_last_cellbender_droplet'] >= EMPTY_PLATEAU_UMI
sample_qc['sex_flag'] = sample_qc['inferred_sex'] != sample_qc['recorded_sex']
sample_qc['qc_flag'] = library_metrics['qc_flag']
sample_qc.rename_axis('sample').to_csv(os.path.join(path_tables, 'sample_qc_flags.csv'))
print(sample_qc.rename(columns={'reads_in_cells_pct': 'reads in cells (%)', 'intronic_reads_pct': 'intronic reads (%)'}).to_string())

# %%
removed = qc_table['outlier'].values
primary_reason = np.full(len(qc_table), '', dtype=object)
for mask, reason_name in [((qc_table['cellbender_ngenes'] < MIN_GENES).values, 'low genes'),
                   ((qc_table['cellbender_ncount'] < MIN_UMI).values, 'low UMIs'),
                   ((qc_table['percent_mito'] > MAX_PCT_MITO).values, 'high %mito'),
                   (qc_table['technical_cluster'].values, 'technical cluster'),
                   (qc_table['hard_doublet'].values, 'predicted doublet'),
                   ((qc_table['flag_cellranger_doublet_scores'] & (qc_table['iqr_outlier'] | qc_table['fixed_doublet_outlier'])).values, 'doublet score fence'),
                   (((qc_table['flag_exon_prop'] | qc_table['flag_log_genes_x_entropy']) & (qc_table['iqr_outlier'] | qc_table['fixed_exon_complexity_outlier'])).values, 'exon_prop / complexity fence'),
                   (qc_table['iqr_outlier'].values, 'other adaptive (IQR)')]:
    assign = removed & mask & (primary_reason == ''); primary_reason[assign] = reason_name
qc_table['primary_reason'] = primary_reason
reasons = ['low genes', 'low UMIs', 'high %mito', 'technical cluster', 'predicted doublet', 'doublet score fence',
           'exon_prop / complexity fence', 'other adaptive (IQR)']
reason_palette = dict(zip(reasons, ['#0072B2', '#56B4E9', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#F0E442', '#999999']))
reason_table = (pd.crosstab(qc_table['sample'], qc_table['primary_reason']).reindex(columns=reasons, fill_value=0)
        .div(retention_by_sample['n_cellbender'], axis=0)).loc[sample_order]
reason_table.add_suffix('_prop').to_csv(os.path.join(path_tables, 'qc_validation_primary_reason_by_sample.csv'))

dbl_any = qc_table.groupby('sample').apply(lambda values: ((values['flag_cellranger_doublet_scores'] | values['hard_doublet']) & values['outlier']).sum() / values['outlier'].sum()).loc[sample_order]

cluster_profiles = pd.read_csv(os.path.join(path_tables, 'preqc_cluster_profiles.csv'), index_col=0)
lowc_clusters = cluster_profiles[cluster_profiles['med_genes'] < 50].sort_values('med_genes').index
comp_lowc = pd.crosstab(qc_table['leiden_overcluster'], qc_table['sample']).loc[lowc_clusters]
comp_lowc = comp_lowc[FLAGGED].assign(**{f'other {len(sample_list) - len(FLAGGED)} samples': comp_lowc.drop(columns=FLAGGED).sum(1)}).div(comp_lowc.sum(1), axis=0)

fig = plt.figure(figsize=(14, 8.6)); grid = gridspec.GridSpec(2, 2, width_ratios=[1.5, 1], hspace=0.5, wspace=0.25)
panel_ax = fig.add_subplot(grid[0, 0])
data = [np.log10(np.clip(qc_table.loc[qc_table['sample'] == sample_id, 'cellbender_ngenes'].values, 1, None)) for sample_id in sample_order]
parts = panel_ax.violinplot(data, positions=np.arange(len(sample_order)), showextrema=False, widths=0.85)
for pc, sample_id in zip(parts['bodies'], sample_order): pc.set_facecolor(sample_colour[sample_id]); pc.set_alpha(0.85); pc.set_edgecolor('none')
for i, values in enumerate(data): panel_ax.plot(i, np.median(values), 'o', color='w', mec='k', ms=3, zorder=4)
panel_ax.axhline(np.log10(MIN_GENES), color='k', ls='--', lw=0.8); panel_ax.text(len(sample_order) - 0.5, np.log10(MIN_GENES) + 0.04, f'floor {MIN_GENES}', fontsize=6, ha='right')
panel_ax.set_xticks(np.arange(len(sample_order))); panel_ax.set_xticklabels(sample_order, rotation=90); color_ticklabels(panel_ax)
panel_ax.set_yticks([0, 1, 2, 3, 4]); panel_ax.set_yticklabels(['$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$']); panel_ax.set_ylabel('genes / nucleus (pre-QC)')
panel_ax.set_title('(a) Pre-QC genes per nucleus, all CellBender droplets')
panel_ax = fig.add_subplot(grid[0, 1]); left = np.zeros(len(sample_order))
for reason in reasons:
    panel_ax.barh(np.arange(len(sample_order)), reason_table[reason].values, left=left, color=reason_palette[reason], label=reason); left += reason_table[reason].values
panel_ax.set_yticks(np.arange(len(sample_order))); panel_ax.set_yticklabels(sample_order); color_ticklabels(panel_ax)
panel_ax.set_xlabel('fraction of CellBender droplets removed (proportion)'); panel_ax.set_xlim(0, 0.7); panel_ax.legend(frameon=False, fontsize=6, loc='lower right')
panel_ax.set_title('(b) Primary removal reason (hard thresholds take priority)')
panel_ax = fig.add_subplot(grid[1, 0])
panel_ax.bar(np.arange(len(sample_order)), 100 * dbl_any.values, color=[sample_colour[sample_id] for sample_id in sample_order])
panel_ax.set_xticks(np.arange(len(sample_order))); panel_ax.set_xticklabels(sample_order, rotation=90); color_ticklabels(panel_ax)
panel_ax.set_ylabel('% of removed nuclei carrying\nany doublet flag'); panel_ax.set_title('(c) Contribution of Scrublet-dependent criteria to removal')
panel_ax = fig.add_subplot(grid[1, 1]); left = np.zeros(len(lowc_clusters))
for sample_col, c in zip(comp_lowc.columns, [REMOVE, '#E07B6C', '#F2B8AE', '#F7D6D0', GREY][-len(comp_lowc.columns):]):
    panel_ax.barh(np.arange(len(lowc_clusters)), comp_lowc[sample_col].values, left=left, color=c, label=sample_col); left += comp_lowc[sample_col].values
panel_ax.set_yticks(np.arange(len(lowc_clusters))); panel_ax.set_yticklabels([f"cl {c}  (n={int(cluster_profiles.loc[c,'n']):,}, med genes {int(cluster_profiles.loc[c,'med_genes'])})" for c in lowc_clusters], fontsize=6)
panel_ax.set_xlabel('fraction of cluster (proportion)'); panel_ax.set_xlim(0, 1); panel_ax.legend(frameon=False, fontsize=6, loc='lower right')
share_flagged = qc_table['sample'].isin(FLAGGED).mean()
panel_ax.set_title(f'(d) Sample of origin of low-complexity pre-QC clusters (median genes < 50)\nthe {len(FLAGGED)} flagged samples are {100*share_flagged:.0f}% of all pre-QC nuclei')
fig.suptitle('Why nuclei were removed, by library (flagged libraries in red)', y=0.96, fontsize=10)
save_svg(fig, 'FigS15_removal_anatomy_flagged_samples'); plt.close(fig)
print(reason_table.rename_axis(columns='primary reason (proportion of CellBender droplets)').round(3).to_string()); print('\nany doublet flag among removed (proportion):\n', dbl_any.round(3).to_string())

# %%
scrublet_tab = pd.read_csv(os.path.join(path_tables, 'scrublet_summary.csv'))
fixed_cutoff = float(np.median(scrublet_tab.loc[(scrublet_tab['source'] == 'cellranger') & (scrublet_tab['method'] == 'automatic'), 'threshold_auto']))
separation_rows = []
for sample_id in sample_list:
    for src in ['cellranger', 'cellbender']:
        t = pd.read_csv(os.path.join(path_tables, f'{sample_id}_{src}_scrublet_scores.csv.gz'))
        obs, sim = t['observed'].dropna().values, t['simulated'].dropna().values
        auc = mannwhitneyu(sim, obs, alternative='greater').statistic / (len(sim) * len(obs))
        separation_rows.append({'sample': sample_id, 'source': src, 'auroc_sim_gt_obs': auc,
                         'obs_median': np.median(obs), 'sim_median': np.median(sim), 'sim_p90': np.percentile(sim, 90)})
separation = pd.DataFrame(separation_rows).merge(scrublet_tab[['sample', 'source', 'threshold_auto', 'threshold_used', 'method', 'predicted_fraction', 'detectable_fraction']], on=['sample', 'source'])
separation['flagged'] = separation['sample'].isin(FLAGGED)
separation.rename(columns=EXPORT_NAMES).to_csv(os.path.join(path_tables, 'qc_validation_scrublet_separation.csv'), index=False)
print(separation[separation['source'] == 'cellranger'].round(3).sort_values('auroc_sim_gt_obs').rename(columns=DISPLAY_NAMES).to_string())

fig = plt.figure(figsize=(14, 12)); grid = gridspec.GridSpec(5, 4, height_ratios=[1.25, 1, 1, 1, 1], hspace=0.6, wspace=0.3, top=0.93)

for j, (y_col, y_label) in enumerate([('threshold_auto', 'automatic Scrublet threshold'), ('detectable_fraction', 'detectable doublet proportion (threshold used)')]):
    panel_ax = fig.add_subplot(grid[0, 2*j:2*j+2])
    for src, marker_shape in [('cellranger', 'o'), ('cellbender', 's')]:
        source_sep = separation[separation['source'] == src]
        panel_ax.scatter(source_sep['auroc_sim_gt_obs'], source_sep[y_col], marker=marker_shape, s=28, c=[sample_colour[x] for x in source_sep['sample']], edgecolor='k', linewidths=0.3, label=f'{src} counts')
    for _, sep_row in separation[separation['flagged'] & (separation['source'] == 'cellranger')].iterrows():
        panel_ax.annotate(sep_row['sample'], (sep_row['auroc_sim_gt_obs'], sep_row[y_col]), xytext=(5, 3), textcoords='offset points', fontsize=6, color=REMOVE)
    if y_col == 'threshold_auto':
        panel_ax.axhline(fixed_cutoff, color='k', ls=':', lw=0.7); panel_ax.text(0.51, fixed_cutoff + 0.01, f"fixed fallback cutoff {fixed_cutoff:.2f} (CellRanger)", fontsize=6)
    panel_ax.axvline(0.5, color=GREY, lw=0.8); panel_ax.set_xlim(0.45, 1.0); panel_ax.set_xlabel('AUROC (simulated doublets score > observed nuclei)'); panel_ax.set_ylabel(y_label)
    panel_ax.set_title(('(a) ' if j == 0 else '(b) ') + f'Score separation vs {y_label}')
    if j == 0: panel_ax.legend(frameon=False, loc='upper right')

for k, sample_id in enumerate(sample_order):
    panel_ax = fig.add_subplot(grid[1 + k // 4, k % 4])
    t = pd.read_csv(os.path.join(path_tables, f'{sample_id}_cellranger_scrublet_scores.csv.gz'))
    obs, sim = t['observed'].dropna().values, t['simulated'].dropna().values
    bins = np.linspace(0, 1, 61)
    panel_ax.hist(obs, bins=bins, density=True, color=KEEP, alpha=0.75, label='observed nuclei')
    panel_ax.hist(sim, bins=bins, density=True, histtype='step', color='k', lw=0.9, label='simulated doublets')
    sep_row = separation[(separation['sample'] == sample_id) & (separation['source'] == 'cellranger')].iloc[0]
    panel_ax.axvline(sep_row['threshold_used'], color=REMOVE, ls='--', lw=0.9, label='threshold used')
    panel_ax.set_title(f"{sample_id}   AUROC={sep_row['auroc_sim_gt_obs']:.2f}  thr={sep_row['threshold_used']:.2f} ({sep_row['method'].split(' (')[0]})  called {100*sep_row['predicted_fraction']:.1f}%", color=sample_colour[sample_id], fontsize=6.5)
    panel_ax.set_xlim(0, 1); panel_ax.set_yticks([]); panel_ax.set_ylabel('density')
    if k == 0: panel_ax.legend(frameon=False, fontsize=6)
    if k >= 12: panel_ax.set_xlabel('Scrublet score')
fig.suptitle('Scrublet diagnostics: automatic threshold vs simulated/observed score separation, and the threshold used per library', y=0.985, fontsize=9)
save_svg(fig, 'FigS16_scrublet_separation_all_samples'); plt.close(fig)

# %%
good_libraries = library_metrics.drop(index=FLAGGED)['frac_removed']
REF_SAMPLE = good_libraries.index[np.argsort(good_libraries.values)[len(good_libraries) // 2]]
print('reference sample:', REF_SAMPLE)

def scrublet_run(counts, seed=0):
    scrublet_obj = scr.Scrublet(counts, random_state=seed); scores, predicted = scrublet_obj.scrub_doublets(verbose=False)
    return {'n': counts.shape[0], 'threshold': getattr(scrublet_obj, 'threshold_', np.nan),
            'predicted_fraction': np.nan if predicted is None else float(np.mean(predicted)),
            'detectable_fraction': scrublet_obj.detectable_doublet_fraction_,
            'auroc_sim_gt_obs': mannwhitneyu(scrublet_obj.doublet_scores_sim_, scores, alternative='greater').statistic / (len(scores) * len(scrublet_obj.doublet_scores_sim_)),
            'scores': scores, 'sim': scrublet_obj.doublet_scores_sim_}
def thin_counts(counts, target_median, rng):
    counts = sp.csr_matrix(counts); thin_frac = min(1.0, target_median / np.median(np.asarray(counts.sum(1)).ravel()))
    thinned = counts.copy(); thinned.data = rng.binomial(thinned.data.astype(int), thin_frac).astype(counts.dtype); thinned.eliminate_zeros()
    return thinned, thin_frac

rng = np.random.default_rng(0); controls = []
reference_adata = sc.read_h5ad(os.path.join(path_objects, REF_SAMPLE + '.preqc.scrub.h5ad')); reference_counts = reference_adata.layers['cellranger_raw']
median_umis = {sample_id: float(np.median(qc_table.loc[qc_table['sample'] == sample_id, 'cellranger_ncount'])) for sample_id in FLAGGED + [REF_SAMPLE]}
control_row = scrublet_run(reference_counts); control_row.update(control=f'{REF_SAMPLE} native (median UMI {median_umis[REF_SAMPLE]:.0f})', kind='reference'); controls.append(control_row)
for sample_id in FLAGGED:
    if median_umis[sample_id] >= median_umis[REF_SAMPLE]: continue
    thinned, thin_frac = thin_counts(reference_counts, median_umis[sample_id], rng)
    control_row = scrublet_run(thinned); control_row.update(control=f'{REF_SAMPLE} thinned to {sample_id} depth (median UMI {median_umis[sample_id]:.0f}, p={thin_frac:.2f})', kind='depth'); controls.append(control_row)
N_DRAWS = 5
small_sizes = sorted({int(retention_by_sample.loc[sample_id, 'n_cellbender']) for sample_id in FLAGGED if retention_by_sample.loc[sample_id, 'n_cellbender'] < 5000} |
               {int(retention_by_sample.loc[sample_id, 'n_retained']) for sample_id in FLAGGED if retention_by_sample.loc[sample_id, 'n_retained'] < 5000})
for n in small_sizes:
    draws = [scrublet_run(reference_counts[rng.choice(reference_counts.shape[0], n, replace=False)], seed=draw) for draw in range(N_DRAWS)]
    n_detectable = sum(draw['detectable_fraction'] > 0.2 for draw in draws)
    control_row = sorted(draws, key=lambda draw: draw['detectable_fraction'])[N_DRAWS // 2]
    control_row.update(control=f'{REF_SAMPLE} subsampled to n={n:,} (threshold found in {n_detectable}/{N_DRAWS} draws)', kind='size'); controls.append(control_row)
for sample_id in FLAGGED:
    sample_adata = sc.read_h5ad(os.path.join(path_objects, sample_id + '.preqc.scrub.h5ad'))
    with h5py.File(os.path.join(path_samples, sample_id, 'filtered_feature_bc_matrix.h5')) as h5: cellranger_barcodes = set(h5['matrix/barcodes'][:].astype(str))
    control_row = scrublet_run(sample_adata.layers['cellranger_raw'][sample_adata.obs_names.isin(cellranger_barcodes)]); control_row.update(control=f'{sample_id} on CellRanger cell calls only', kind='input'); controls.append(control_row)
controls_tab = pd.DataFrame(controls).drop(columns=['scores', 'sim'])
controls_tab.rename(columns=EXPORT_NAMES).to_csv(os.path.join(path_tables, 'qc_validation_scrublet_controls.csv'), index=False)
print(controls_tab.rename(columns=DISPLAY_NAMES).round(3).to_string())

kind_colour = {'reference': KEEP, 'depth': '#4C9F70', 'size': '#E0A03C', 'input': REMOVE}
n_cols = int(np.ceil(len(controls) / 2))
fig = plt.figure(figsize=(6 + 2.2 * n_cols, 6.2)); grid = gridspec.GridSpec(2, 1 + n_cols, width_ratios=[1.8] + [1] * n_cols, wspace=0.35, hspace=0.6)
panel_ax = fig.add_subplot(grid[:, 0]); y_pos = np.arange(len(controls))[::-1]
panel_ax.barh(y_pos, controls_tab['auroc_sim_gt_obs'], color=[kind_colour[k] for k in controls_tab['kind']])
for yi, (_, control_row) in zip(y_pos, controls_tab.iterrows()):
    panel_ax.text(control_row['auroc_sim_gt_obs'] + 0.01, yi, f"thr {control_row['threshold']:.2f} · called {100*control_row['predicted_fraction']:.1f}%", va='center', fontsize=6)
panel_ax.set_yticks(y_pos); panel_ax.set_yticklabels(controls_tab['control'], fontsize=6.5); panel_ax.set_xlim(0.4, 1.15); panel_ax.axvline(0.5, color=GREY, lw=0.8)
panel_ax.set_xlabel('AUROC (simulated > observed)'); panel_ax.set_title('(a) Same automatic Scrublet call on controlled inputs')
panel_ax.legend(handles=[mpl.patches.Patch(color=c, label=k) for k, c in kind_colour.items()], frameon=False, ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.1))
bins = np.linspace(0, 1, 61)
for k, control_row in enumerate(controls):
    pos = (k // n_cols, 1 + k % n_cols)
    panel_ax = fig.add_subplot(grid[pos[0], pos[1]])
    panel_ax.hist(control_row['scores'], bins=bins, density=True, color=kind_colour[control_row['kind']], alpha=0.75); panel_ax.hist(control_row['sim'], bins=bins, density=True, histtype='step', color='k', lw=0.9)
    if not np.isnan(control_row['threshold']): panel_ax.axvline(control_row['threshold'], color=REMOVE, ls='--', lw=0.9)
    panel_ax.set_title(control_row['control'].replace(' (', '\n('), fontsize=6.5); panel_ax.set_yticks([]); panel_ax.set_xlim(0, 1)
    if pos[0] == 1: panel_ax.set_xlabel('Scrublet score')
fig.suptitle(f'Scrublet degradation controls: {REF_SAMPLE} thinned / subsampled with the same code, and the flagged libraries restricted to STARsolo cell calls', y=1.0, fontsize=9)
save_svg(fig, 'FigS17_scrublet_degradation_controls'); plt.close(fig)

# %%
atlas_meta = pd.read_csv(os.path.join(path_tables, 'cell_metadata.csv'), index_col=0)
retained_qc = qc_table[~qc_table['outlier']]
retained_metric_specs = [('cellbender_ngenes', 'genes / nucleus', True), ('cellbender_ncount', 'UMIs / nucleus', True),
        ('percent_mito', '% mitochondrial (raw counts)', False), ('cellbender_entropy', 'transcriptome entropy (bits)', False),
        ('cellranger_doublet_scores', 'Scrublet score', False)]
fig = plt.figure(figsize=(14, 8.5)); grid = gridspec.GridSpec(2, 5, hspace=0.75, wspace=0.35)
for j, (metric, metric_label, log_scale) in enumerate(retained_metric_specs):
    panel_ax = fig.add_subplot(grid[0, j])
    data = [retained_qc.loc[retained_qc['sample'] == sample_id, metric].values for sample_id in sample_order]
    ref_med = np.median(retained_qc.loc[~retained_qc['sample'].isin(FLAGGED), metric].values)
    if log_scale: data = [np.log10(np.clip(values, 1, None)) for values in data]; ref_med = np.log10(max(ref_med, 1))
    parts = panel_ax.violinplot(data, positions=np.arange(len(sample_order)), showextrema=False, widths=0.85)
    for pc, sample_id in zip(parts['bodies'], sample_order): pc.set_facecolor(sample_colour[sample_id]); pc.set_alpha(0.85); pc.set_edgecolor('none')
    for i, values in enumerate(data): panel_ax.plot(i, np.median(values), 'o', color='w', mec='k', ms=2.5, zorder=4)
    panel_ax.axhline(ref_med, color='k', ls=':', lw=0.8)
    panel_ax.set_xticks(np.arange(len(sample_order))); panel_ax.set_xticklabels(sample_order, rotation=90, fontsize=6); color_ticklabels(panel_ax)
    all_values = np.concatenate(data); panel_ax.set_ylim(np.percentile(all_values, 0.1), np.percentile(all_values, 99.9) * 1.05)
    if log_scale:
        y_lo, y_hi = int(np.floor(panel_ax.get_ylim()[0])), int(np.ceil(panel_ax.get_ylim()[1]))
        panel_ax.set_yticks(range(y_lo, y_hi + 1)); panel_ax.set_yticklabels([f'$10^{{{k}}}$' for k in range(y_lo, y_hi + 1)])
    panel_ax.set_title(metric_label); panel_ax.set_xlabel('library')
fig.text(0.5, 0.94, '(a) QC metrics of retained nuclei per library (dotted = median of the non-flagged libraries)', ha='center', fontsize=8.5)
panel_ax = fig.add_subplot(grid[1, :3])
composition = pd.crosstab(atlas_meta['sample'], atlas_meta['cell_type'], normalize='index').reindex(index=sample_order, columns=order).fillna(0)
left = np.zeros(len(sample_order))
for c in order:
    panel_ax.barh(np.arange(len(sample_order)), composition[c].values, left=left, color=CELL_TYPE_PALETTE[c], label=c); left += composition[c].values
panel_ax.set_yticks(np.arange(len(sample_order))); panel_ax.set_yticklabels(sample_order); color_ticklabels(panel_ax); panel_ax.set_xlim(0, 1)
panel_ax.set_xlabel('proportion of nuclei in the final atlas'); panel_ax.set_title('(b) Cell-type composition per library (final atlas)')
panel_ax.legend(frameon=False, fontsize=6, ncol=6, loc='upper center', bbox_to_anchor=(0.5, -0.22))
panel_ax = fig.add_subplot(grid[1, 3:])
n_final = atlas_meta.groupby('sample').size().reindex(sample_order)
panel_ax.barh(np.arange(len(sample_order)), n_final.values, color=[sample_colour[sample_id] for sample_id in sample_order])
for yi, v in enumerate(n_final.values): panel_ax.text(v + 100, yi, f'{v:,}', va='center', fontsize=6)
panel_ax.set_yticks(np.arange(len(sample_order))); panel_ax.set_yticklabels([]); panel_ax.set_xlabel('nuclei in final atlas'); panel_ax.xaxis.set_major_formatter(thousands_fmt)
panel_ax.set_title('(c) Nuclei contributed to the atlas')
fig.suptitle('QC metrics and composition of retained nuclei per library; libraries with >50% QC attrition highlighted in red', y=1.0, fontsize=9)
save_svg(fig, 'FigS18_retained_nuclei_flagged_samples'); plt.close(fig)

retained_metrics = retained_qc.groupby('sample').agg(
    qc_flag=('qc_flag', 'first'), n_retained_initial_qc=('outlier', 'size'),
    median_genes=('cellbender_ngenes', 'median'), median_umi=('cellbender_ncount', 'median'),
    median_pct_mito=('percent_mito', 'median'), frac_mito_positive=('percent_mito', lambda v: float((v > 0).mean())),
    median_entropy=('cellbender_entropy', 'median'), median_scrublet=('cellranger_doublet_scores', 'median'))
retained_metrics['n_final_atlas'] = n_final
retained_metrics = retained_metrics.join(composition.round(3).add_suffix('_prop')).loc[sample_order]
retained_metrics.rename(columns=EXPORT_NAMES).to_csv(os.path.join(path_tables, 'qc_validation_retained_metrics_by_sample.csv'))
print(retained_metrics.drop(columns=composition.add_suffix('_prop').columns).rename(columns=DISPLAY_NAMES).round(3).to_string())
print('\ncomposition of the flagged libraries:\n', composition.rename_axis(columns='cell type (proportion)').round(3).loc[FLAGGED].to_string())

# %%
IFN_PANEL = ['Gbp6', 'Gbp10', 'Gbp4', 'Iigp1', 'Igtp', 'Irgm1', 'Gm4841', 'Stat1', 'Cd274', 'Ifit3']
cluster_annotation = pd.read_csv(os.path.join(path_tables, 'cluster_annotation_before_exclusions.csv'), index_col=0)
single_animal_tab = cluster_annotation[(cluster_annotation['cell_type'] == 'Cardiomyocyte') & (cluster_annotation['largest_sample_fraction'] > 0.5)]
SINGLE_ANIMAL = [str(c) for c in single_animal_tab.index]
FOCUS = single_animal_tab['largest_sample'].mode().iloc[0]
print('single-animal cardiomyocyte clusters:', SINGLE_ANIMAL, '| animal:', FOCUS)
cm_ec_adata = atlas_adata[atlas_adata.obs['cell_type'].isin(['Cardiomyocyte', 'Endothelial'])].copy()
sc.tl.score_genes(cm_ec_adata, [group_label for group_label in IFN_PANEL if group_label in cm_ec_adata.var_names], score_name='score_IFN', use_raw=False, random_state=0)
cm_obs = cm_ec_adata.obs[cm_ec_adata.obs['cell_type'] == 'Cardiomyocyte']
cm_clusters = cm_obs['leiden'].value_counts().loc[lambda n: n >= 100].index
cm_clusters = sorted(cm_clusters, key=int)
focus_group = atlas_adata.obs.loc[atlas_adata.obs['individual'] == FOCUS, 'disease'].iloc[0]
highlight_colour = GROUP_FILL[focus_group]
cluster_colour = {c: (highlight_colour if c in SINGLE_ANIMAL else GREY) for c in cm_clusters}
xlab = [f"{c}\n{(cm_obs['leiden'] == c).sum():,}" for c in cm_clusters]

fig = plt.figure(figsize=(15, 8.6)); grid = gridspec.GridSpec(2, 4, height_ratios=[1, 1], hspace=0.6, wspace=0.4)

ax = fig.add_subplot(grid[0, :2])
share = pd.crosstab(cm_obs['leiden'], cm_obs['individual'], normalize='index').loc[cm_clusters]
focus_share = share[FOCUS]; other_share = 1 - focus_share
ax.bar(range(len(cm_clusters)), focus_share, color=highlight_colour, width=0.8, label=f'{FOCUS} ({focus_group})')
ax.bar(range(len(cm_clusters)), other_share, bottom=focus_share, color=GREY, width=0.8, label=f'other {len(sample_list) - 1} animals')
ax.axhline(1 / len(sample_list), color='0.4', ls=':', lw=0.8); ax.text(len(cm_clusters) - 0.4, 1 / len(sample_list), f' 1/{len(sample_list)}', fontsize=7, va='center')
ax.set_xticks(range(len(cm_clusters))); ax.set_xticklabels(xlab); ax.set_xlabel('post-QC cluster (cardiomyocyte label) and n nuclei')
ax.set_ylabel(f'fraction of nuclei from {FOCUS} (proportion)'); ax.set_ylim(0, 1)
ax.set_title(f"a  Clusters {', '.join(SINGLE_ANIMAL)}: > 50 % of nuclei from one animal", loc='left', pad=20)
ax.legend(frameon=False, loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2, borderaxespad=0.2)

lowc_threshold = LOWC_REL * atlas_adata.obs.loc[atlas_adata.obs['cell_type'] == 'Cardiomyocyte', 'cellbender_ngenes'].median()
for j, (metric, metric_label) in enumerate([('cellbender_ngenes', 'genes per nucleus'), ('cellranger_doublet_scores', 'Scrublet doublet score')]):
    ax = fig.add_subplot(grid[0, 2 + j])
    data = [cm_obs.loc[cm_obs['leiden'] == c, metric].values for c in cm_clusters]
    boxes = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.7, medianprops={'color': 'k', 'lw': 0.8},
                    whiskerprops={'color': '0.6', 'lw': 0.6}, capprops={'color': '0.6', 'lw': 0.6}, boxprops={'lw': 0})
    for body, c in zip(boxes['boxes'], cm_clusters): body.set_facecolor(cluster_colour[c])
    ax.set_xticks(range(1, len(cm_clusters) + 1)); ax.set_xticklabels(cm_clusters, fontsize=7); ax.set_ylabel(metric_label); ax.set_xlabel('cluster')
    ax.set_title(f"{'b' if j == 0 else ' '}  {metric_label}", loc='left')
ax_genes = fig.axes[1]; ax_genes.axhline(lowc_threshold, color='0.4', ls=':', lw=0.8)
ax_genes.text(0.6, lowc_threshold, f'low-complexity threshold ({lowc_threshold:.0f})', fontsize=6, va='bottom', color='0.35')

ax = fig.add_subplot(grid[1, :2])
animal_ifn = (cm_ec_adata.obs.groupby(['individual', 'cell_type'], observed=True)['score_IFN'].mean().unstack()
        .join(cm_ec_adata.obs[['individual', 'disease']].drop_duplicates().set_index('individual')))
animal_ifn = animal_ifn.sort_values('Cardiomyocyte')
x = np.arange(len(animal_ifn))
for k, (ct, marker_shape) in enumerate([('Cardiomyocyte', 'o'), ('Endothelial', 's')]):
    ax.scatter(x + (k - 0.5) * 0.25, animal_ifn[ct], marker=marker_shape, s=28, c=[GROUP_FILL[group_label] for group_label in animal_ifn['disease']], edgecolor='0.3', linewidths=0.4, label=ct, zorder=3)
ax.set_xticks(x); ax.set_xticklabels(animal_ifn.index, rotation=90, fontsize=7)
for t in ax.get_xticklabels(): t.set_fontweight('bold' if t.get_text() == FOCUS else 'normal')
ax.set_ylabel('mean interferon score per animal'); ax.set_title(f'c  Interferon score per animal in cardiomyocytes and endothelium ({FOCUS} in bold)', loc='left', pad=20)
ax.legend(handles=[Line2D([], [], marker='o', ls='', mfc='w', mec='0.3', label='cardiomyocytes'), Line2D([], [], marker='s', ls='', mfc='w', mec='0.3', label='endothelium')]
          + [Line2D([], [], marker='o', ls='', color=GROUP_FILL[group_label], label=group_label) for group_label in GROUP_ORDER], frameon=False, ncol=6, loc='lower left', bbox_to_anchor=(0, 1.0), borderaxespad=0.2)

ax = fig.add_subplot(grid[1, 2:])
focus_cm = cm_obs[cm_obs['individual'] == FOCUS]
nucleus_group = np.where(focus_cm['leiden'].astype(str).isin(SINGLE_ANIMAL), 'cluster ' + focus_cm['leiden'].astype(str), 'other CM clusters')
order_d = sorted(set(nucleus_group), key=lambda group_label: np.median(focus_cm['score_IFN'].values[nucleus_group == group_label]))
data = [focus_cm['score_IFN'].values[nucleus_group == group_label] for group_label in order_d]
parts = ax.violinplot(data, showextrema=False, widths=0.8)
ifn_median_others = cm_obs.loc[cm_obs['individual'] != FOCUS, 'score_IFN'].median()
ax.axhline(ifn_median_others, color='0.4', ls=':', lw=0.8); ax.text(len(order_d) + 0.45, ifn_median_others, f' median,\n other {len(sample_list) - 1}\n animals', fontsize=6, va='center')
for body, group_label in zip(parts['bodies'], order_d): body.set_facecolor(highlight_colour if group_label != 'other CM clusters' else GREY); body.set_alpha(0.85); body.set_edgecolor('none')
for i, values in enumerate(data): ax.plot(i + 1, np.median(values), 'o', color='w', mec='k', ms=3, zorder=3)
ax.set_xticks(range(1, len(order_d) + 1)); ax.set_xticklabels([f'{group_label}\n(n = {len(values):,})' for group_label, values in zip(order_d, data)]); ax.set_xlim(0.4, len(order_d) + 0.9)
ax.set_ylabel('interferon score per nucleus'); ax.set_title(f'd  Within {FOCUS}: interferon score by cluster', loc='left')
fig.suptitle('Single-animal cardiomyocyte clusters: Step 7 exclusion metrics', y=1.01, fontsize=9)
save_svg(fig, 'FigS19_single_animal_clusters'); plt.close(fig)

single_animal_check = pd.DataFrame({
    'n': cm_obs['leiden'].value_counts().reindex(cm_clusters),
    f'{FOCUS}_prop': focus_share.round(3),
    'median_genes': cm_obs.groupby('leiden', observed=True)['cellbender_ngenes'].median().reindex(cm_clusters),
    'median_mito_pct': cm_obs.groupby('leiden', observed=True)['percent_mito'].median().reindex(cm_clusters).round(3),
    'median_scrublet': cm_obs.groupby('leiden', observed=True)['cellranger_doublet_scores'].median().reindex(cm_clusters).round(3),
    'mean_CM_score': cm_obs.groupby('leiden', observed=True)['score_Cardiomyocyte'].mean().reindex(cm_clusters).round(2),
    'mean_IFN_score': cm_obs.groupby('leiden', observed=True)['score_IFN'].mean().reindex(cm_clusters).round(3)})
single_animal_check.to_csv(os.path.join(path_tables, 'single_animal_cluster_check.csv'))
print(single_animal_check.to_string())

# %%
from scipy import io, sparse
import gzip, shutil
adata = sc.read_h5ad(os.path.join(path_objects, 'atlas_annotated.h5ad'))
path_export = os.path.join(path_objects, 'for_seurat')
os.makedirs(path_export, exist_ok=True)

counts = sparse.csr_matrix(adata.layers['counts'])
assert np.isfinite(counts.data).all() and (counts.data >= 0).all(), 'Invalid counts'
assert np.equal(counts.data, np.floor(counts.data)).all(), 'Export requires integer counts'
io.mmwrite(os.path.join(path_export, 'counts.mtx'), counts.T.astype(np.int64))
with open(os.path.join(path_export, 'counts.mtx'), 'rb') as file_in, gzip.open(os.path.join(path_export, 'counts.mtx.gz'), 'wb') as file_out:
    shutil.copyfileobj(file_in, file_out)
os.remove(os.path.join(path_export, 'counts.mtx'))

pd.DataFrame({'gene': adata.var_names}).to_csv(os.path.join(path_export, 'genes.csv.gz'), index=False)
pd.DataFrame({'barcode': adata.obs_names}).to_csv(os.path.join(path_export, 'barcodes.csv.gz'), index=False)
adata.obs.to_csv(os.path.join(path_export, 'metadata.csv.gz'))
adata.var[['gene_ids']].rename_axis('gene_symbol').to_csv(os.path.join(path_export, 'gene_annotations.csv.gz'))
pd.crosstab(adata.obs['sample'], adata.obs['cell_type']).to_csv(
    os.path.join(path_export, 'final_celltype_counts_by_sample.csv'))
print('Final export:', adata.n_obs, 'nuclei;', adata.n_vars, 'genes;', adata.obs.shape[1], 'metadata columns')
for reduction, key in [('pca', 'X_pca'), ('harmony', 'PC_harmony'), ('umap', 'UMAP_PC_harmony')]:
    pd.DataFrame(adata.obsm[key], index=adata.obs_names).to_csv(
        os.path.join(path_export, reduction + '.csv.gz'), header=False)
print('exported for Seurat ->', path_export)
print(sorted(os.listdir(path_export)))

# %%
from importlib.metadata import version
packages = ['scanpy', 'anndata', 'numpy', 'pandas', 'scipy', 'matplotlib',
            'harmonypy', 'scrublet', 'scikit-learn', 'scikit-misc',
            'leidenalg', 'igraph', 'h5py']
versions = pd.Series({pkg: version(pkg) for pkg in packages}, name='version')
versions.to_csv(os.path.join(path_tables, 'software_versions.csv'))
print(versions.to_string())
