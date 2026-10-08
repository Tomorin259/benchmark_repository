"""Shared gene-evaluation inputs and scientific case definitions."""

METHOD_ORDER = [
    "cellDancer", "DeepVelo", "DeepVelo_scv", "Graphvelo_scv", "MultiVelo",
    "scVelo_dyn", "TIvelo", "UniTVelo", "veloVI",
]

METHOD_KEYS = {
    "cellDancer": "celldancer",
    "DeepVelo": "deepvelo",
    "DeepVelo_scv": "deepvelo2022",
    "Graphvelo_scv": "graphvelo_scv",
    "MultiVelo": "multivelo",
    "scVelo_dyn": "scvelo_dyn",
    "TIvelo": "tivelo",
    "UniTVelo": "unitvelo_with_vu",
    "veloVI": "veloVI",
}

SIM_METHOD_KEYS = {
    "cellDancer": "celldancer",
    "DeepVelo": "deepvelo",
    "DeepVelo_scv": "deepvelo2022",
    "Graphvelo_scv": "graphvelo",
    "MultiVelo": "multivelo",
    "scVelo_dyn": "scvelo",
    "TIvelo": "tivelo",
    "UniTVelo": "unitvelo",
    "veloVI": "velovi",
}

VELOCITY_KEYS = {
    "cellDancer": ("velocity", "velocity_u"),
    "DeepVelo": ("velocity", "velocity_unspliced"),
    "DeepVelo_scv": ("velocity_reconstructed", "velocity_u"),
    "Graphvelo_scv": ("velocity_gv", "velocity_u"),
    "MultiVelo": ("velo_s", "velo_u"),
    "scVelo_dyn": ("velocity", "velocity_u"),
    "TIvelo": ("velocity", "velocity_u"),
    "UniTVelo": ("velocity", "velocity_u"),
    "veloVI": ("velocity", "velocity_u"),
}

REAL_CONFIGS = {
    "retina": {
        "neighbor_asset": "retina.npz",
        "label_reference": "results/deepvelo/results_with_vu/deepvelo_retina.h5ad",
        "genes": ["Anks1b", "Cacna1b", "Calm2", "Celf4", "Hcn1", "Myt1l", "Rasgef1b", "Slc17a6", "Snap25", "Zfhx3"],
        "edges": [("Neuroblast", "PR"), ("Neuroblast", "AC/HC"), ("Neuroblast", "RGC")],
        "paths": {
            "cellDancer": "results/celldancer/celldancer_retina.h5ad",
            "DeepVelo": "results/deepvelo/results_with_vu/deepvelo_retina.h5ad",
            "DeepVelo_scv": "results/deepvelo2022/deepvelo_retina_full_recon.h5ad",
            "Graphvelo_scv": "results/graphvelo_scv/graphveloretina.h5ad",
            "MultiVelo": "results/multivelo/multivelo_retina.h5ad",
            "scVelo_dyn": "results/scvelo_dynamic/scv_dynamic_retina.h5ad",
            "TIvelo": "results/TIvelo/tivelo_retina.h5ad",
            "UniTVelo": "results/unitvelo/results_with_vu/unitvelo_retina_with_velocity_u_mode1.h5ad",
            "veloVI": "results/veloVI/veloVI_retina_full_all.h5ad",
        },
    },
    "DentateGyrus_lamanno": {
        "neighbor_asset": "DentateGyrus_lamanno.npz",
        "label_reference": "results/deepvelo/results_with_vu/deepvelo_DentateGyrus_lamanno.h5ad",
        "genes": ["Gpm6b", "Ncald", "Diaph3", "Dab1", "Dpp6", "Malat1", "Nnat", "Tubb3", "Plppr4", "Mmd2"],
        "edges": [
            ("nIPC", "GlialProg"), ("GlialProg", "ImmAstro"), ("GlialProg", "OPC"),
            ("Nbl1", "Nbl2"), ("Nbl2", "ImmGranule1"),
            ("ImmGranule1", "ImmGranule2"), ("ImmGranule2", "Granule"),
            ("Nbl2", "CA"), ("CA", "CA1-Sub"), ("CA", "CA2-3-4"),
        ],
        "paths": {
            "cellDancer": "results/celldancer/celldancer_DentateGyrus_lamanno.h5ad",
            "DeepVelo": "results/deepvelo/results_with_vu/deepvelo_DentateGyrus_lamanno.h5ad",
            "DeepVelo_scv": "results/deepvelo2022/deepvelo_DentateGyrus_lamanno_recon.h5ad",
            "Graphvelo_scv": "results/graphvelo_scv/graphveloDentateGyrus_lamanno.h5ad",
            "MultiVelo": "results/multivelo/multivelo_DentateGyrus_lamanno.h5ad",
            "scVelo_dyn": "results/scvelo_dynamic/scv_dynamic_DentateGyrus_lamanno.h5ad",
            "TIvelo": "results/TIvelo/tivelo_DentateGyrus_lamanno.h5ad",
            "UniTVelo": "results/unitvelo/results_with_vu/unitvelo_DentateGyrus_lamanno_with_velocity_u_mode1.h5ad",
            "veloVI": "results/veloVI/veloVI_DentateGyrus_lamanno_all.h5ad",
        },
    },
    "gastrulation": {
        "neighbor_asset": "gastrulation.npz",
        "label_reference": "results/deepvelo/results_with_vu/deepvelo_Gastrulation.h5ad",
        "genes": [
            "Abcg2", "Alas2", "Bpgm", "Cdc25b", "Dhrs11", "Fam210b", "Fech",
            "Gclm", "Gpx1", "Gypa", "Hba-x", "Hbb-bh1", "Hbb-y", "Hemgn",
            "Klf3", "Mllt3", "Ncf4", "Rhag",
        ],
        "edges": [
            ("Blood progenitors 1", "Blood progenitors 2"),
            ("Blood progenitors 2", "Erythroid1"),
            ("Erythroid1", "Erythroid2"),
            ("Erythroid2", "Erythroid3"),
        ],
        "paths": {
            "cellDancer": "results/celldancer/celldancer_Gastrulation.h5ad",
            "DeepVelo": "results/deepvelo/results_with_vu/deepvelo_Gastrulation.h5ad",
            "DeepVelo_scv": "results/deepvelo2022/deepvelo_Gastrulation_recon.h5ad",
            "Graphvelo_scv": "results/graphvelo_scv/graphveloGastrulation.h5ad",
            "MultiVelo": "results/multivelo/multivelo_Gastrulation.h5ad",
            "scVelo_dyn": "results/scvelo_dynamic/scv_dynamic_Gastrulation.h5ad",
            "TIvelo": "results/TIvelo/tivelo_Gastrulation.h5ad",
            "UniTVelo": "results/unitvelo/results_with_vu/unitvelo_Gastrulation_with_velocity_u_mode1.h5ad",
            "veloVI": "results/veloVI/veloVI_Gastrulation_all.h5ad",
        },
    },
    "capo_pre": {
        "neighbor_asset": "capofibro.npz",
        "label_reference": "results/deepvelo/results_with_vu/deepvelo_capo_pre.h5ad",
        "genes": ["ALPK2", "C4orf21", "C19orf48", "CLSPN", "DEK", "E2F7", "FAM111A"],
        "edges": [("G1", "S"), ("S", "G2M"), ("G2M", "G1")],
        "paths": {
            "cellDancer": "results/celldancer/celldancer_capo_pre.h5ad",
            "DeepVelo": "results/deepvelo/results_with_vu/deepvelo_capo_pre.h5ad",
            "DeepVelo_scv": "results/deepvelo2022/deepvelo_capo_pre_recon.h5ad",
            "Graphvelo_scv": "results/graphvelo_scv/graphvelocapo_pre.h5ad",
            "MultiVelo": "results/multivelo/multivelo_capo_pre.h5ad",
            "scVelo_dyn": "results/scvelo_dynamic/scv_dynamic_capo_pre.h5ad",
            "TIvelo": "results/TIvelo/tivelo_capo_pre.h5ad",
            "UniTVelo": "results/unitvelo/results_with_vu/unitvelo_capo_pre_with_velocity_u_mode1.h5ad",
            "veloVI": "results/veloVI/veloVI_capo_pre_all.h5ad",
        },
    },
    "bonemarrow": {
        "neighbor_asset": "bonemarrow.npz",
        "label_reference": "results/deepvelo/results_with_vu/deepvelo_bonemarrow.h5ad",
        "genes": ["FAM178B", "LMNA", "REXO2"],
        "edges": [("HSC_1", "Ery_1"), ("HSC_1", "HSC_2"), ("Ery_1", "Ery_2")],
        "paths": {
            "cellDancer": "results/celldancer/celldancer_bonemarrow.h5ad",
            "DeepVelo": "results/deepvelo/results_with_vu/deepvelo_bonemarrow.h5ad",
            "DeepVelo_scv": "results/deepvelo2022/deepvelo_bonemarrow_recon.h5ad",
            "Graphvelo_scv": "results/graphvelo_scv/graphvelobonemarrow.h5ad",
            "MultiVelo": "results/multivelo/multivelo_bonemarrow.h5ad",
            "scVelo_dyn": "results/scvelo_dynamic/scv_dynamic_bonemarrow.h5ad",
            "TIvelo": "results/TIvelo/tivelo_bonemarrow.h5ad",
            "UniTVelo": "results/unitvelo/results_with_vu/unitvelo_bonemarrow_with_velocity_u_mode1.h5ad",
            "veloVI": "results/veloVI/veloVI_bonemarrow_all.h5ad",
        },
    },
    "mousebrain": {
        "neighbor_asset": "mousebrain.npz",
        "label_reference": "results/deepvelo/results_with_vu/deepvelo_mousebrain.h5ad",
        "genes": ["Adcy1", "Cacna1e", "Flrt2", "Gria2"],
        "edges": [
            ("RG, Astro, OPC", "IPC"), ("IPC", "V-SVZ"),
            ("V-SVZ", "Upper Layer"), ("Upper Layer", "Deeper Layer"),
        ],
        "paths": {
            "cellDancer": "results/celldancer/celldancer_mousebrain.h5ad",
            "DeepVelo": "results/deepvelo/results_with_vu/deepvelo_mousebrain.h5ad",
            "DeepVelo_scv": "results/deepvelo2022/deepvelo_mousebrain_recon.h5ad",
            "Graphvelo_scv": "results/graphvelo_scv/graphvelomousebrain.h5ad",
            "MultiVelo": "results/multivelo/multivelo_mousebrain.h5ad",
            "scVelo_dyn": "results/scvelo_dynamic/scv_dynamic_mousebrain.h5ad",
            "TIvelo": "results/TIvelo/tivelo_mousebrain.h5ad",
            "UniTVelo": "results/unitvelo/results_with_vu/unitvelo_mousebrain_with_velocity_u_mode1.h5ad",
            "veloVI": "results/veloVI/veloVI_mousebrain_all.h5ad",
        },
    },
    "pancreas": {
        "neighbor_asset": "pancreas.npz",
        "label_reference": "results/deepvelo/results_with_vu/deepvelo_Pancreas.h5ad",
        "genes": ["Atp2a2", "Chst8", "Emb", "G6pc2", "Igf1r", "Nxph1", "Rbbp7"],
        "edges": [
            ("Ductal", "Ngn3 low EP"), ("Ngn3 low EP", "Ngn3 high EP"),
            ("Ngn3 high EP", "Pre-endocrine"), ("Pre-endocrine", "Beta"),
            ("Pre-endocrine", "Alpha"), ("Pre-endocrine", "Delta"),
            ("Pre-endocrine", "Epsilon"),
        ],
        "paths": {
            "cellDancer": "results/celldancer/celldancer_Pancreas.h5ad",
            "DeepVelo": "results/deepvelo/results_with_vu/deepvelo_Pancreas.h5ad",
            "DeepVelo_scv": "results/deepvelo2022/deepvelo_Pancreas_recon.h5ad",
            "Graphvelo_scv": "results/graphvelo_scv/graphveloPancreas.h5ad",
            "MultiVelo": "results/multivelo/multivelo_Pancreas.h5ad",
            "scVelo_dyn": "results/scvelo_dynamic/scv_dynamic_Pancreas.h5ad",
            "TIvelo": "results/TIvelo/tivelo_Pancreas.h5ad",
            "UniTVelo": "results/unitvelo/results_with_vu/unitvelo_Pancreas_with_velocity_u_mode1.h5ad",
            "veloVI": "results/veloVI/veloVI_Pancreas_all.h5ad",
        },
    },
}

SIM_CONFIGS = {
    "branching": {
        "reference_source": "indicators/gene/simulation/branching_analyze/dataset/s1_branch.h5ad",
        "genes": [f"MultiForward_{i}" for i in range(10)] + [f"BoostALL_{i}" for i in range(10)],
        "edges": [(f"BranchA_{i}", f"BranchA_{i + 1}") for i in range(9)] + [(f"BranchB_{i}", f"BranchB_{i + 1}") for i in range(9)],
    },
    "murk": {
        "reference_source": "indicators/gene/simulation/murk_analyze/data/s1_only_murk.h5ad",
        "genes": [f"OneSidedBoost_{i}" for i in range(10)],
        "edges": [(f"BranchB_{i}", f"BranchB_{i + 1}") for i in range(9)],
    },
    "all": {
        "reference_source": "indicators/gene/simulation/sim_all_analyze/data/sim_all.h5ad",
        "genes": (
            [f"MultiForward_{i}" for i in range(10)]
            + [f"BoostALL_{i}" for i in range(10)]
            + [f"OneSidedBoost_{i}" for i in range(10)]
            + [f"OneSidedBoost2_{i}" for i in range(10)]
        ),
        "edges": [(f"BranchA_{i}", f"BranchA_{i + 1}") for i in range(9)] + [(f"BranchB_{i}", f"BranchB_{i + 1}") for i in range(9)],
    },
}

SIM_PATHS = {
    "branching": {
        "cellDancer": "indicators/gene/simulation/results_for_analyze/celldancer/celldancer_branching.h5ad",
        "DeepVelo": "indicators/gene/simulation/results_for_analyze/deepvelo/deepvelo_branching.h5ad",
        "DeepVelo_scv": "indicators/gene/simulation/results_for_analyze/deepvelo2022/deepvelo2022_branching_recon.h5ad",
        "Graphvelo_scv": "indicators/gene/simulation/results_for_analyze/graphvelo_scv/graphvelo_branching.h5ad",
        "MultiVelo": "indicators/gene/simulation/results_for_analyze/multivelo/multivelo_branching.h5ad",
        "scVelo_dyn": "indicators/gene/simulation/results_for_analyze/scvelo/scv_dynamic_branching.h5ad",
        "TIvelo": "indicators/gene/simulation/results_for_analyze/tivelo/tivelo_branching.h5ad",
        "UniTVelo": "indicators/gene/simulation/results_for_analyze/unitvelo/resultstouse/unitvelo_branching.h5ad",
        "veloVI": "indicators/gene/simulation/results_for_analyze/velovi/veloVI_branching.h5ad",
    },
    "murk": {
        "cellDancer": "indicators/gene/simulation/results_for_analyze/celldancer/celldancer_murk.h5ad",
        "DeepVelo": "indicators/gene/simulation/results_for_analyze/deepvelo/deepvelo_murk.h5ad",
        "DeepVelo_scv": "indicators/gene/simulation/results_for_analyze/deepvelo2022/deepvelo2022_murk_recon.h5ad",
        "Graphvelo_scv": "indicators/gene/simulation/results_for_analyze/graphvelo_scv/graphvelo_murk.h5ad",
        "MultiVelo": "indicators/gene/simulation/results_for_analyze/multivelo/multivelo_murk.h5ad",
        "scVelo_dyn": "indicators/gene/simulation/results_for_analyze/scvelo/scv_dynamic_murk.h5ad",
        "TIvelo": "indicators/gene/simulation/results_for_analyze/tivelo/tivelo_murk.h5ad",
        "UniTVelo": "indicators/gene/simulation/results_for_analyze/unitvelo/resultstouse/unitvelo_murk.h5ad",
        "veloVI": "indicators/gene/simulation/results_for_analyze/velovi/veloVI_murk.h5ad",
    },
    "all": {
        "cellDancer": "indicators/gene/simulation/results_for_analyze/celldancer/celldancer_all.h5ad",
        "DeepVelo": "indicators/gene/simulation/results_for_analyze/deepvelo/deepvelo_all.h5ad",
        "DeepVelo_scv": "indicators/gene/simulation/results_for_analyze/deepvelo2022/deepvelo2022_all_recon.h5ad",
        "Graphvelo_scv": "indicators/gene/simulation/results_for_analyze/graphvelo_scv/graphvelo_all.h5ad",
        "MultiVelo": "indicators/gene/simulation/results_for_analyze/multivelo/multivelo_all.h5ad",
        "scVelo_dyn": "indicators/gene/simulation/results_for_analyze/scvelo/scv_dynamic_all.h5ad",
        "TIvelo": "indicators/gene/simulation/results_for_analyze/tivelo/tivelo_all.h5ad",
        "UniTVelo": "indicators/gene/simulation/results_for_analyze/unitvelo/resultstouse/unitvelo_all.h5ad",
        "veloVI": "indicators/gene/simulation/results_for_analyze/velovi/veloVI_all.h5ad",
    },
}
