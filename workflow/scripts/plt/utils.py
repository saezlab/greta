

def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config


def savefigs(lst_figs, path_fname, index_pngs=[]):
    import matplotlib.backends.backend_pdf
    import io
    from PIL import Image
    import matplotlib.pyplot as plt
    pdf = matplotlib.backends.backend_pdf.PdfPages(path_fname)
    for i, fig in enumerate(lst_figs):
        if i not in index_pngs:
            pdf.savefig(fig, bbox_inches='tight')
        else:
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            buf.seek(0)
            image = Image.open(buf)
            new_fig, ax = plt.subplots(figsize=(image.width / 100, image.height / 100), dpi=300)
            ax.imshow(image)
            ax.axis("off")
            pdf.savefig(new_fig, bbox_inches='tight')
            plt.close(new_fig)
    pdf.close()
