"""Recover veloVI's own pre-min-max moment units for numerical evaluation."""

import hashlib
import warnings

import numpy as np
from scipy import sparse
from sklearn.preprocessing import MinMaxScaler


def dense(values):
    return values.toarray() if sparse.issparse(values) else np.asarray(values)


def restore_velovi_coordinates(adata, gene_indices, preprocessed=None):
    """Return own moment arrays, velocity multipliers and exact-forward evidence.

    The saved result does not retain its scaler. Recover the moment source from
    its saved abundance and binary/self-connected, row-normalized graph, as in
    scVelo 0.2.5. For capo_pre the caller supplies its actual preprocessed input,
    because its saved final graph is not the graph that produced its moments.
    Reject any source that cannot reproduce saved Ms/Mu bit for bit. Nothing in
    the input AnnData objects is modified.
    """
    gene_indices = np.asarray(gene_indices, dtype=int)
    genes = np.asarray(adata.var_names)[gene_indices].astype(str)
    original_graph = adata.obsp["connectivities"].tocsr()
    graph = (original_graph > 0).tocsr()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", sparse.SparseEfficiencyWarning)
        graph.setdiag(1)
    graph = graph.multiply(1.0 / graph.sum(1)).tocsr().astype(np.float32)
    graph_copy = original_graph.copy()
    graph_copy.sort_indices()
    digest = hashlib.sha256()
    for value in (np.asarray(graph_copy.shape), graph_copy.data,
                  graph_copy.indices, graph_copy.indptr):
        digest.update(value.tobytes())
    graph_hash = digest.hexdigest()
    source_kind = "own_saved_abundance_and_graph"
    if preprocessed is not None:
        if not preprocessed.obs_names.is_unique or not preprocessed.var_names.is_unique:
            raise ValueError("veloVI preprocessed input identifiers are not unique")
        rows = preprocessed.obs_names.astype(str).get_indexer(adata.obs_names.astype(str))
        columns = preprocessed.var_names.astype(str).get_indexer(genes)
        if np.any(rows < 0) or np.any(columns < 0):
            raise ValueError("veloVI preprocessed input does not cover saved cell/gene identifiers")
        source_kind = "actual_preprocessed_input_moments; saved_final_graph_not_moment_origin"

    moments, multipliers, evidence = [], [], []
    for axis, abundance in (("s", "spliced"), ("u", "unspliced")):
        saved_abundance = dense(adata.layers[abundance][:, gene_indices])
        restored = graph.dot(sparse.csr_matrix(saved_abundance)).astype(np.float32).toarray()
        native = dense(adata.layers["M" + axis][:, gene_indices])
        graph_error = float(np.max(np.abs(MinMaxScaler().fit_transform(restored) - native)))
        if preprocessed is not None:
            source_abundance = dense(preprocessed.layers[abundance][:, columns])[rows]
            if not np.array_equal(source_abundance, saved_abundance):
                raise ValueError("veloVI preprocessed input abundance differs: " + abundance)
            restored = dense(preprocessed.layers["M" + axis][:, columns])[rows]
        scaler = MinMaxScaler().fit(restored)
        forward = scaler.transform(restored)
        if not np.array_equal(forward, native) or np.any(scaler.data_range_ <= 0):
            raise ValueError("veloVI " + axis + " moments do not exactly reproduce saved min-max coordinates")
        inverse = 1.0 / scaler.scale_.astype(np.float64)
        moments.append(restored.astype(np.float64))
        multipliers.append(inverse)
        for j, gene in enumerate(genes):
            evidence.append({"gene": gene, "axis": axis,
                             "data_min": float(scaler.data_min_[j]),
                             "data_max": float(scaler.data_max_[j]),
                             "data_range": float(scaler.data_range_[j]),
                             "scale": float(scaler.scale_[j]), "offset": float(scaler.min_[j]),
                             "velocity_inverse_multiplier": float(inverse[j]),
                             "forward_exact": True, "forward_max_abs_error": 0.0,
                             "moment_source_kind": source_kind,
                             "source_graph_sha256": graph_hash,
                             "saved_graph_forward_error_before_source_resolution": graph_error})
    return moments[0], moments[1], multipliers[0], multipliers[1], evidence


def velocity_cosine(velocity, truth):
    """Historical simulation CosSim with the same >1e-6 norm mask."""
    nv, nt = np.linalg.norm(velocity, axis=1), np.linalg.norm(truth, axis=1)
    valid = (nv > 1e-6) & (nt > 1e-6)
    if not valid.any():
        return 0.0, 0
    values = np.sum(velocity[valid] * truth[valid], axis=1) / (nv[valid] * nt[valid])
    return float(np.clip(values, -1.0, 1.0).mean()), int(valid.sum())
