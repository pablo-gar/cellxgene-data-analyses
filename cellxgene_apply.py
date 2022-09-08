import boto3
import pandas as pd
import anndata as ad
import os
import logging
from typing import Dict

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
S3 = boto3.resource("s3")


def apply_function_portal_h5ads(
    data_table_file: str, fun: callable, sep: str = "\t"
) -> Dict:
    """
     Applies a function to all h5ads from CZ CELLxGENE Discover as defined
     in a data table file.

     :param str data_table_file: Path to tabular file containing CZ CELLxGENE
     Discover, it is obtained from
     here "https://github.com/chanzuckerberg/single-cell-data-portal/blob/
     810a64812a9942784f939d7fc76444145ded2946/scripts/cxg_admin.py#L289-L295"

     :param callable fun: function to be applied to h5ads, must only have
     1 positional argument the anndata.AnnData obj
     :param str sep="\t": separator for columns of `data_table_file`
    cellxgene_data_table_file
     :return Dict: Dictionary with keys as the s3 uris for each h5ad and
     items are the corresponding return of `fun`
    """

    data_info = get_all_data_info(data_table_file)

    results = {}
    counter = 1
    total = len(data_info)

    for dev_uri in data_info:

        logger.info(f"Working on data {counter} of {total}")

        try:
            download_from_s3(dev_uri, "temp.h5ad")
        except:
            logger.warning(f"Could not find {dev_uri}")
            continue

        logger.info(f"Applying function on {dev_uri}")
        adata = ad.read("temp.h5ad", "r")

        results[dev_uri] = fun(adata)

        os.remove("temp.h5ad")
        if counter > 5:
            break
        counter += 1

    return results


def get_all_data_info(
    data_table_file: str, sep: str = "\t"
) -> Dict[str, Dict]:
    """
    Parses information of tabular file containing CZ CELLxGENE Discover,
    it is obtained from
    here "https://github.com/chanzuckerberg/single-cell-data-portal/blob/
    810a64812a9942784f939d7fc76444145ded2946/scripts/cxg_admin.py#L289-L295"

    :return Dict[Dict]:  nested dictionary, keys for top dictionary are s3
    uris of each h5ad, keys for the nested
    dictionary are the following:
        - "dataset_uri_prod": str
        - "explorer_url": str
        - "collection_id": str
        - "collection_name": str
        - "cell_count": int
        - "organisms": List[str]
        - "tissues": List[str]
        - "assays": List[str]
    """

    def converter_fun(x): return x.strip("[]").replace("'", "").split(", ")

    all_data = pd.read_csv(
        data_table_file,
        sep=sep,
        index_col=0,
        converters={
            "S3 URIs": converter_fun,
            "Organisms": converter_fun,
            "Tissues": converter_fun,
            "Assays": converter_fun,
        },
    )

    data_dict = {}

    for row in all_data.iterrows():

        row = row[1]

        for dataset_uri in row["S3 URIs"]:

            if "h5ad" not in dataset_uri:
                continue

            data_dict[dataset_uri] = {
                "dataset_uri_prod": dataset_uri.replace(
                    "corpora-data-dev", "corpora-data-prod"
                ),
                "explorer_url": row["Explorer URL"],
                "collection_id": row["ID"],
                "collection_name": row["Name"],
                "cell_count": row["Cell Count"],
                "organisms": row["Organisms"],
                "tissues": row["Tissues"],
                "assays": row["Assays"],
            }

    return data_dict


def download_from_s3(s3_uri, filename):

    s3_uri = s3_uri.split("/")
    s3_bucket = s3_uri[2]
    s3_key = "/".join(s3_uri[3:])

    S3.meta.client.download_file(s3_bucket, s3_key, filename)
