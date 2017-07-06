import os
import ConfigParser
import subprocess


PATHS = dict()
config = ConfigParser.ConfigParser()
config.read(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '../..', 'config/paths.cfg'))
for key, value in config.items('paths'):
    PATHS[key] = value


def quiet_call(call_list, stdout=os.devnull):
    '''
    Call outside program while suppressing messages to stdout and stderr.

    stdout_file -- captures stdout (default os.devnull)
    '''
    with open(stdout, 'w') as OUT:
        return subprocess.call(
            call_list, stdout=OUT, stderr=open(os.devnull, 'w'))
