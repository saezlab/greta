## Utils:

import yaml
from pathlib import Path
import ast
import argparse
from scdori import ppConfig, trainConfig # type: ignore
import numpy as np
import pandas as pd
import logging

def aggregate_grn_max_val(grn_act, grn_rep, rna_metacell):
    """
    Aggregate GRN activator and repressor matrices by taking the maximum absolute value
    across the two matrices for each TF-gene pair.
    Keeps the sign of the values.
    """
    stacked = np.stack([grn_rep, grn_act], axis=0)  # shape (2, num_topics, num_tfs, num_genes)

    # Collapse along first two dims
    stacked = stacked.reshape(-1, 300, 4000) 

    # Find indices of max abs value along 0th axis
    idx = np.abs(stacked).argmax(axis=0)

    # Extract actual values with sign
    grn = np.take_along_axis(stacked, idx[None, :, :], axis=0).squeeze(0)  # shape (num_tfs, num_genes)]
    tf_names = rna_metacell.var[rna_metacell.var.gene_type == "TF"].index.values
    grn_df = pd.DataFrame(
        grn, index=tf_names, columns=rna_metacell.var_names
    )
    grn_long = grn_df.reset_index(
        ).melt(
            id_vars="index", var_name="column", value_name="value"
        ).sort_values(
            'value', ascending=True
        ).rename(
            columns={'index': 'source', 'column': 'target', 'value': 'score'}
        )

    return grn_long

def load_config_from_yaml_or_cmdline(yaml_path: str, cmdline_args: dict = None):
    """
    Load configuration from a YAML file and override with command line arguments if provided.
    """

    with open(yaml_path, 'r') as f:
        cfg = yaml.safe_load(f)

    shared_cfg = cfg.get("scDoRIConfig", {})
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=shared_cfg['logging_level'])

    # Override yaml if command line arguments are provided
    if cmdline_args:
        for key, val in cmdline_args.items():
            if hasattr(ppConfig, key) or hasattr(trainConfig, key) or key in shared_cfg:
                logger.info(f"Overriding parameter {key} from commandline: {shared_cfg[key]} -> {val}")
                shared_cfg[key] = val
            else:
                raise ValueError(f"Key {key} not found in configuration.")
    
    data_dir = Path(shared_cfg.get("data_dir", "."))

    # Resolve paths
    for key in ["genome_dir", "motif_directory", "model_dir"]:
        if key in shared_cfg:
            shared_cfg[key] = Path(data_dir, shared_cfg[key])


    if "weight_dir" in shared_cfg:
        shared_cfg["weight_dir"] = Path(shared_cfg["weight_dir"])
        for k in ["weights_folder_scdori", "weights_folder_grn"]:
            shared_cfg[k] = Path(shared_cfg["weight_dir"], shared_cfg[k])

    if "model_dir" in shared_cfg:
        for k in ["best_scdori_model_path", "best_grn_model_path"]:
            shared_cfg[k] = Path(shared_cfg["model_dir"], shared_cfg[k])
    # Take care of naming mismatches
    shared_cfg['batch_col'] = shared_cfg.get('batch_key') # Different naming in preprocessing and training config
    shared_cfg['output_subdir'] = Path(shared_cfg.get('output_subdir_name'))

    # Override any default values in ppConfig and trainConfig with shared_cfg
    for key, val in shared_cfg.items():
        if hasattr(ppConfig, key):
            setattr(ppConfig, key, val)
        if hasattr(trainConfig, key):
            setattr(trainConfig, key, val)
        elif key == 'grn_file_out':
            setattr(trainConfig, key, val)
        else:
            # Adding custom preprocessing parameters to ppConfig
            if not hasattr(ppConfig, key):
                logger.info(f"Added additional parameter {key} to ppConfig: -> {val}")
                setattr(ppConfig, key, val)
        # if key == 'mudata_file_name':
        #     setattr(ppConfig, key, val)
        
    return ppConfig, trainConfig

def parse_cmdline_args():
    """Parse command line arguments and return the YAML file path and a dictionary of overrides.
    Command line arguments should be in the format key=value, separated by a whitespace.
    """
    cmdline_args = {}
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None, help="Path to the YAML configuration file.")
    known_args, overrides = parser.parse_known_args()

    yaml_file = vars(known_args).get("config")
    for override in overrides:
        if "=" not in override:
            raise ValueError(f"Invalid override: {override}. Use key=value format.")
        key, val = override.split("=", 1)
        # Parse value to correct type
        try:
            val = ast.literal_eval(val)
        except:
            pass
        cmdline_args[key] = val
    return yaml_file, cmdline_args
