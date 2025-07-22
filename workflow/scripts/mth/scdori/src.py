
from utils import load_config_from_yaml_or_cmdline, parse_cmdline_args
from pp_scDORI import wrapper_scdori_preprocessing
from train_scDORI import wrapper_scdori_training
from grns_scDORI import wrapper_scdori_grns
#from scdori import ppConfig, trainConfig # type: ignore

def main():
    # Load configuration from YAML file and command line arguments
    yaml_file, cmdline_args = parse_cmdline_args()
    if yaml_file:
        ppConfig, trainConfig = load_config_from_yaml_or_cmdline(yaml_file, cmdline_args)
    else:
        raise ValueError("Configuration file not provided. Use --config to specify the YAML file.")

    # Run preprocessing
    wrapper_scdori_preprocessing(ppConfig)

    # Run training
    wrapper_scdori_training(trainConfig)

    # Run GRN analysis
    wrapper_scdori_grns(trainConfig)

if __name__ == "__main__":
    main()