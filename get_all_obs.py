import pandas as pd
import anndata as ad
from collections import Counter
import cellxgene_apply


def main():

    cellxgene_data_table_file = "./data_tables/cellxgene_all_data_latest.tsv"

    results = cellxgene_apply.apply_function_portal_h5ads(
        cellxgene_data_table_file, make_all_obs
    )

    big_obs = pd.concat(list(results.values()))


def make_all_obs(adata: ad.AnnData, max_categories: int = 1000):

    # Remove columns with high number of items
    columns_to_drop = []
    for column in adata.obs.columns:

        if isinstance(adata.obs[column].dtype, pd.api.types.CategoricalDtype):

            categories = adata.obs[column].cat.categories
            n_categories = len(categories)

            if n_categories > max_categories:
                columns_to_drop.append(column)

    adata.obs = adata.obs.drop(columns=columns_to_drop)

    return adata.obs.drop_duplicates().reset_index()


if __name__ == "__main__":
    main()