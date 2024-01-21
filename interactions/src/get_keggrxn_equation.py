#!/usr/bin/env python

import re
import time

from functools import wraps
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen



def queue_pop_column_from_df(dataframe, column):
    return [rxn for rxn in dataframe[column]]


def retry(ExceptionToCheck, tries=4, delay=3, backoff=2, logger=None):
    """Retry calling the decorated function using an exponential backoff.

    http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

    :param ExceptionToCheck: the exception to check. may be a tuple of
        exceptions to check
    :type ExceptionToCheck: Exception or tuple
    :param tries: number of times to try (not retry) before giving up
    :type tries: int
    :param delay: initial delay between retries in seconds
    :type delay: int
    :param backoff: backoff multiplier e.g. value of 2 will double the delay
        each retry
    :type backoff: int
    :param logger: logger to use. If None, print
    :type logger: logging.Logger instance
    """
    def deco_retry(f):

        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except ExceptionToCheck as e:
                    msg = "%s, Retrying in %d seconds..." % (str(e), mdelay)
                    if logger:
                        logger.warning(msg)
                    else:
                        print(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry



@retry((TimeoutError, URLError, HTTPError), tries=4, delay=15, backoff=2)
def urlopen_with_retry(kegg_rxnid):
    rxn_req = Request('https://rest.kegg.jp/get/%s' %kegg_rxnid, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/99.0.4844.82 Safari/537.36'})
    return urlopen(rxn_req)


def get_keggrxn_equation(kegg_rxnid):

    temp = urlopen_with_retry(kegg_rxnid)
    rxn_lines = temp.readlines()

    results = {}

    for num in range(len(rxn_lines)):
        rxn_info = rxn_lines[num].decode('utf-8')
        sep_info = re.split('\s\s+', rxn_info.strip())
        table = dict(enumerate(sep_info))

        if 'ENTRY' == table.get(0):
            results['ENTRY'] = table.get(1)
        if 'DEFINITION' == table.get(0):
            results['DEFINITION'] = table.get(1)
        if 'EQUATION' == table.get(0):
            results['EQUATION'] = table.get(1)

    return results