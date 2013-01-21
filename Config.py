#!/usr/bin/env python

#/AM/home-0/shared/python/Python-2.7.1/python

import ConfigParser
import os
import re

class Config():
    def __init__(self, openedConfigFile, section):
        self._config = ConfigParser.ConfigParser()
        #self._config.readfp(open(os.path.normpath(configFilePath)))
        self._config.readfp(openedConfigFile)
        self._section = section

    def get(self, option):
        return self._config.get(self._section, option)


class Config2():
    def __init__(self, config, section):
        self._config = config._config
        self._section = section

    def get(self, option):
        return self._config.get(self._section, option)


def test():
    config = Config(open(os.path.normpath('D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//config01.cfg')), 'pPPS')
    pattern = config.get('scaffoldPattern')
    print pattern
    sn = re.findall(pattern, 'Scaffold_1_4082367.lucy.pga.C2');
    print sn

    print config.get('fastaLineMaxChar')
    print config.get('databaseFile')
    print config.get('taxonomicRanks')

    print config.get('rankIdAll')
    print config.get('rankIdCut')
    print config.get('rankIdCutMinBp')


if __name__ == "__main__":
  test()