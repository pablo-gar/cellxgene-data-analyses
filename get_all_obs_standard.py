import pandas as pd
import anndata as ad
from collections import Counter
import cellxgene_apply


def main():

    cellxgene_data_table_file = "./data_tables/cellxgene_all_data_prod_09142022.json"

    results = cellxgene_apply.apply_function_portal_h5ads(
        cellxgene_data_table_file, make_all_obs
    )

    for key, values in results.items():
        results[key]["explorer_url"] = key

    big_obs = pd.concat(list(results.values()), ignore_index=True)

    big_obs.to_csv("./results/all_obs_standard.tsv.gz", sep="\t", compression="gzip")


def make_all_obs(adata: ad.AnnData, max_categories: int = 1000):

    columns = [
        "assay_ontology_term_id",
        "tissue_ontology_term_id",
        "cell_type_ontology_term_id",
        "development_stage_ontology_term_id",
        "donor_id",
        "is_primary_data",
        "organism_ontology_term_id",
        "self_reported_ethnicity_ontology_term_id",
        "ethnicity_ontology_term_id",
        "sex_ontology_term_id",
        "suspension_type",
        "assay",
        "cell_type",
        "development_stage",
        "disease",
        "organism",
        "self_reported_ethnicity",
        "ethnicity",
        "sex",
        "tissue",
    ]

    # Remove columns with high number of items
    columns_to_drop = []
    for column in adata.obs.columns:

        if column not in columns:
            columns_to_drop.append(column)

    obs = adata.obs.drop(columns=columns_to_drop)
    obs = count_unique_row_pandas(obs)
    return obs


def count_unique_row_pandas(df: pd.DataFrame):

    df_result = df.iloc[0:1, :].copy()
    df_result["count"] = 0

    for i in range(len(df)):
        row = df.iloc[i:i+1, :]
        which = (df_result.drop(columns=["count"]) == row.to_numpy()).all(1)
        if which.any():
            df_result.loc[which, "count"] += 1
        else:
            new_row = row.copy()
            new_row["count"] = 1
            df_result = pd.concat([df_result, new_row])

    return df_result





if __name__ == "__main__":
    main()