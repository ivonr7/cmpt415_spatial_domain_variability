from hest import iter_hest



# Read 1 sample of spatialLIBD dataset
for sample in iter_hest('../hest_data', id_list=[f'MISC{i}' for i in range(1,9)]):
    print(sample.meta['id'],sample.meta['subseries'])
    # adata = sample.adata
    # print('\n* Scanpy adata:')
    # print(adata)
    # # WSI:
    # wsi = sample.wsi
    # print('\n* WSI:')
    # print(wsi)
    
    # # Shapes:
    # shapes = sample.shapes
    # print('\n* Shapes:')
    # print(shapes)
    
    # # Tissue segmentation
    # tissue_contours = sample.tissue_contours
    # print('\n* Tissue contours:')
    # print(tissue_contours)
    
    # # Conversion to SpatialData
    # sdata = sample.to_spatial_data()
    # print('\n* SpatialData conversion:')
    # print(sdata)

    # print(adata.obsm['spatial'])


    # # visualize the spots over a downscaled version of the full resolution image
    # save_dir = '.'
    # sample.save_spatial_plot(save_dir)

