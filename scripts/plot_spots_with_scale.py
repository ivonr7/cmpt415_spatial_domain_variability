from hest import iter_hest
import scanpy as sc
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.pyplot as plt
def plot_with_scale_bar(adata:sc.AnnData,pixel_size:float, key:str = "total_counts",show:bool = False):
    fig = sc.pl.spatial(
        adata,show=show,img_key="downscaled_fullres",
        color = [key],title = "Spots within Tissue",return_fig=True
    )
    scale = ScaleBar(pixel_size,'um',box_alpha=0.5)
    fig.legend()
    fig.gca().add_artist(scale)

    return fig
if __name__ == "__main__":
    for sample in iter_hest('../hest_data', id_list=['MISC1']):

        fig = plot_with_scale_bar(sample.adata,sample.pixel_size,show=False)
        fig.savefig("splots_with_scale.png")