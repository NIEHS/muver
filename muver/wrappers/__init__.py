import os
import ConfigParser


PATHS = dict()
config = ConfigParser.ConfigParser()
config.read(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '../..', 'config/paths.cfg'))
for key, value in config.items('paths'):
    PATHS[key] = value
