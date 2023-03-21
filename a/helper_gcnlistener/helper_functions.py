# IMPORTS
from urllib.error import HTTPError
import traceback
from astropy.io import fits
from time import sleep
import urllib.request
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
def url_tester(url, attempts, sleeptime):
    while attempts > 0:
        request = urllib.request.Request(url)
        try:
            #code with possible error
            print('# '+str(attempts)+' attempts left')
            urllib.request.urlopen(request)
        except HTTPError:
            print('attempt failed')
            attempts -= 1
            print('Sleep for '+str(sleeptime)+' sec')
            sleep(sleeptime)
            continue
        except:
            print(traceback.format_exc())

        #the rest of the code
        break

    print(url)
    return url

def track_REV_tester(url, rev, attempts, sleeptime):
    while attempts > 0:
        request = urllib.request.Request(url)
        try:
            #code with possible error
            print('# '+str(attempts)+' attempts left')
            urllib.request.urlopen(request)
            data = urllib.request.urlopen(url).read() 
            data=data.decode()
            rev_str='REVISION:         '+str(rev)
            if rev_str not in data:
                raise ValueError
        except (HTTPError, ValueError):
            print('attempt failed')
            attempts -= 1
            print('Sleep for '+str(sleeptime)+' sec')
            sleep(sleeptime)
            continue
        except:
            print(traceback.format_exc())

        #the rest of the code
        break

    print(url)
    return url

def find_between_tags(lst, start_tag, end_tag): 
    start_index = lst.index(start_tag) 
    end_index = lst.index(end_tag, start_index) 
    return lst[start_index + 1: end_index] 
