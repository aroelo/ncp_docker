from flask import request, Blueprint
import pandas as pd

# Define the blueprint
get_json_data_blueprint = Blueprint('get_json_data', __name__)


@get_json_data_blueprint.route("/get_json_data")
def create_json_data():
    """Create json from pandas dataframe of primary and secondary alignments of reads"""
    df_pickle = request.args.get('df_pickle')
    taxid_list = request.args.get('taxid_list')
    taxid_list = taxid_list.split(',')

    df = pd.read_pickle(df_pickle)
    df_taxid = df[df['taxid0'].isin(taxid_list)]
    json_taxid = df_taxid.to_json(orient='split')
    return json_taxid
