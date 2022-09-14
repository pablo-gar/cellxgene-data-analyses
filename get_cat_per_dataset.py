import pandas as pd
import anndata as ad
from collections import Counter
import cellxgene_apply


def main():

    cellxgene_data_table_file = "./data_tables/cellxgene_all_data_prod_09142022.json"
    results = cellxgene_apply.apply_function_portal_h5ads(
        cellxgene_data_table_file, count_categories_h5ad
    )

    final_results = {
        "exp_url": [],
        "n_col_one_category": [],
        "n_col_two_or_more_category": [],
        "largest_category_n": [],
        "largest_category_label": [],
    }

    n_cat_per_column = []

    for exp_url in results:

        final_results["exp_url"].append(exp_url)
        final_results["n_col_one_category"].append(
            results[exp_url]["n_col_one_category"]
        )
        final_results["n_col_two_or_more_category"].append(
            results[exp_url]["n_col_two_or_more_category"]
        )
        final_results["largest_category_n"].append(
            results[exp_url]["largest_category_n"]
        )
        final_results["largest_category_label"].append(
            results[exp_url]["largest_category_label"]
        )
        n_cat_per_column.extend(results[exp_url]["n_categories_per_column"])

    n_cat_per_column = Counter(n_cat_per_column)
    n_cat_per_column = pd.DataFrame(n_cat_per_column.most_common())
    n_cat_per_column.rename(
        columns={0: "n_categories", 1: "counts"}, inplace=True
    )

    n_cat_per_column.to_csv(
        "./results/category_counts.tsv", sep="\t", index=False
    )
    pd.DataFrame(final_results).to_csv(
        "./results/n_categories.tsv", sep="\t", index=False
    )


def count_categories_h5ad(adata: ad.AnnData):

    results = {
        "n_col_one_category": 0,
        "n_col_two_or_more_category": 0,
        "largest_category_n": 0,
        "largest_category_label": 0,
        "n_categories_per_column": [],
    }

    for column in adata.obs.columns:

        if "ontology_term_id" in column:
            continue

        if isinstance(adata.obs[column].dtype, pd.api.types.CategoricalDtype):

            categories = adata.obs[column].cat.categories
            n_categories = len(categories)

            results["n_categories_per_column"].append(n_categories)

            if n_categories > results["largest_category_n"]:
                results["largest_category_n"] = n_categories
                results["largest_category_label"] = column

            if len(categories) > 1:
                results["n_col_two_or_more_category"] += 1
            else:
                results["n_col_one_category"] += 1

    return results


if __name__ == "__main__":
    main()
