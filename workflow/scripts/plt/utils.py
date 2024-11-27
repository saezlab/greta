

def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config


def savefigs(lst_figs, path_fname):
    import matplotlib.backends.backend_pdf
    pdf = matplotlib.backends.backend_pdf.PdfPages(path_fname)
    for fig in lst_figs:
        pdf.savefig(fig, bbox_inches='tight')
    pdf.close()
