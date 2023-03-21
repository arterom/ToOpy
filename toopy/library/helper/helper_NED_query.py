# IMPORTS
from astroquery.exceptions import RemoteServiceError
from astroquery.ned import Ned
import numpy as np


def query_NED_object(name, verbose=False, print_header=False):
    """
    Query SIMBAD by name.
    """
    try:
        q = Ned.query_object(name)
        #main_id  = str(q['Redshift'][0], encoding='utf-8')
        z        = q['Redshift']
        #print(z)
    except RemoteServiceError:
        print('This query results in a RemoteServiceError')
        z = np.nan
    return z
