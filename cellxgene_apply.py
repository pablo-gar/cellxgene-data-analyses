import boto3
import pandas as pd
import anndata as ad
import os
import logging
from typing import Dict
import json

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
S3 = boto3.resource("s3")


def apply_function_portal_h5ads(
    data_json_file: str, fun: callable, sep: str = "\t"
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

    data_info = get_all_data_info_from_json(data_json_file)

    results = {}
    counter = 1
    total = len(data_info)

    for uri in data_info:

        logger.info(f"Working on data {counter} of {total}")

        try:
            download_from_s3(uri, "temp.h5ad")
        except:
            logger.warning(f"Could not find {uri}")
            continue

        logger.info(f"Applying function on {uri}")
        explorer_url = data_info[uri]["explorer_url"]

        try:
            adata = ad.read("temp.h5ad", "r")
        except:
            logger.warning(f"Could not find read {explorer_url}")
            continue

        try:
            results[explorer_url] = fun(adata)
        except:
            logger.warning(f"Could not apply function to {explorer_url}")
            continue

        os.remove("temp.h5ad")
        #if counter > 15:
        #    break
        counter += 1

    return results


def get_all_data_info_from_json(
        data_json_file: str, sep: str = "\t"
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
    json_file_obj = open(data_json_file, "r")
    all_data = json.load(json_file_obj)
    json_file_obj.close()

    data_dict = {}

    for current_dict in all_data.values():

        dataset_uri = ""
        for uri in current_dict["s3_uris"]:
            if "h5ad" in uri:
                dataset_uri = uri

        data_dict[dataset_uri] = {
            "dataset_uri_prod": dataset_uri,
            "explorer_url": current_dict["explorer_url"],
            "cell_count": current_dict["cell_count"],
            "organisms": current_dict["organisms"],
            "tissues": current_dict["tissue"],
            "assays": current_dict["assay"],
        }

    return data_dict


def download_from_s3(s3_uri, filename):

    s3_uri = s3_uri.split("/")
    s3_bucket = s3_uri[2]
    s3_key = "/".join(s3_uri[3:])

    S3.meta.client.download_file(s3_bucket, s3_key, filename)
